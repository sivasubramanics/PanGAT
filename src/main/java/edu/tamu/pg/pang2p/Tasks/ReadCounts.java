package edu.tamu.pg.pang2p.Tasks;

import edu.tamu.pg.pang2p.Data.*;
import edu.tamu.pg.pang2p.Utils.*;
import picocli.CommandLine.*;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.concurrent.atomic.AtomicLong;

@Command(name = "count-reads", description = "Count reads to each anchors based on the kmers", mixinStandardHelpOptions = true)
public class ReadCounts implements Callable<Integer>, Runnable {
    @Option(names = {"-r", "--reference"}, description = "Input reference fasta file path", required = true)
    private Path refFastaFile;

    @Option(names = {"-k", "--kmer"}, description = "Kmer size", required = true)
    private int kmerSize;

    @Option(names = {"-i", "--input"}, description = "Input reads file path(s)", required = true, split = ",", arity = "1..*")
    private List<Path> inputFiles;

    @Option(names = {"-t", "--threads"}, description = "Number of threads to use", required = false)
    private int nthreads = 1;

    @Option(names = {"-o", "--output"}, description = "Output file path", required = true)
    private Path outputFile;

    @Option(names = {"--relax"}, description = "Skip remaining kmers in a read once a matching kmer is found", required = false)
    private boolean relax = false;

    @Option(names = {"--skip"}, description = "Step interval for sampling reference kmers (e.g., 5 for every 5th kmer: 1st, 6th, 11th)", required = false)
    private int skip = 1;

    private String CLASS_NAME = ReadCounts.class.getSimpleName();

    public ReadCounts() {}

    public ReadCounts(Path refFastaFile, List<Path> inputFiles, Path outputFile, int kmerSize, int nthreads, boolean relax, int skip) {
        this.refFastaFile = refFastaFile;
        this.kmerSize = kmerSize;
        this.inputFiles = inputFiles;
        this.nthreads = nthreads;
        this.outputFile = outputFile;
        this.relax = relax;
        this.skip = skip;
    }

    @Override
    public void run() {
        try {
            call();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Override
    public Integer call() throws Exception {
        countReads();
        return 0;
    }

    private void countReads() throws IOException {
        if (!FileUtils.isFastaFile(refFastaFile)) {
            Logger.error(CLASS_NAME, "Reference fasta file should have a valid fasta extension (.fasta, .fa, .fna)");
            return;
        }

        Path kmcPrefix = FileUtils.replaceExtension(refFastaFile, ".k" + kmerSize + ".kmc_pre");
        Path kmcSuffix = FileUtils.replaceExtension(refFastaFile, ".k" + kmerSize + ".kmc_suf");
        Path kmcDB = FileUtils.replaceExtension(refFastaFile, ".k" + kmerSize);

        if (!kmcPrefix.toFile().exists() || !kmcSuffix.toFile().exists()) {
            Logger.info(CLASS_NAME, "KMC database files not found for the reference fasta. It will be generated now.");
            KMCConfig kmcConfig = KMCConfig.builder()
                    .inputPath(refFastaFile)
                    .outputPrefix(kmcDB)
                    .kmerSize(kmerSize)
                    .minCount(0)
                    .threads(nthreads)
                    .build();

            Logger.info(CLASS_NAME, "Running KMC for reference fasta: " + refFastaFile.toString());
            KMCExecutor.run(kmcConfig);
        }

        Map<Path, String> validFilesMap = parseInputFiles(inputFiles);
        if (validFilesMap == null) {
            Logger.error(CLASS_NAME, "Input files are not valid. Please check the file paths and extensions.");
            return;
        }

        try (KMC kmc = new KMC(kmcDB, true)) {
            long totalKmers = kmc.getTotalKmers();
            int prefixLength = kmc.getPrefixLength();
            int kmerSize = kmc.getKmerLength();
            int chunkSize = 1_000_000;
            int numChunks = (int) Math.ceil((double) totalKmers / chunkSize);

            AtomicIntegerArray[] kmerCountChunks = new AtomicIntegerArray[numChunks];
            for (int i = 0; i < numChunks; i++) {
                int currentSize = (int) Math.min(chunkSize, totalKmers - ((long) i * chunkSize));
                kmerCountChunks[i] = new AtomicIntegerArray(currentSize);
            }

            AtomicLong totalReadsProcessed = new AtomicLong(0);
            AtomicLong matchedReads = new AtomicLong(0); // <--- Added counter for matched reads
            long reportingInterval = 1_000_000;

            ExecutorService executor = Executors.newFixedThreadPool(nthreads);

            for (Map.Entry<Path, String> entry : validFilesMap.entrySet()) {
                Path inputFile = entry.getKey();
                String format = entry.getValue();

                Logger.info(CLASS_NAME, "Processing input file sequentially: " + inputFile + " [Format: " + format + "]");

                List<Future<Void>> readFutures = new ArrayList<>();
                List<Read> batch = new ArrayList<>();
                int batchSize = 10_000;

                try (SeqReader reader = new SeqReader(inputFile, format)) {
                    for (Read read : reader) {
                        batch.add(read);

                        if (batch.size() >= batchSize) {
                            List<Read> batchToProcess = new ArrayList<>(batch);
                            readFutures.add(executor.submit(() -> {
                                processReadBatch(batchToProcess, kmc, kmerCountChunks, chunkSize, totalReadsProcessed, matchedReads, reportingInterval, relax);
                                return null;
                            }));
                            batch.clear();
                        }
                    }

                    if (!batch.isEmpty()) {
                        List<Read> batchToProcess = new ArrayList<>(batch);
                        readFutures.add(executor.submit(() -> {
                            processReadBatch(batchToProcess, kmc, kmerCountChunks, chunkSize, totalReadsProcessed, matchedReads, reportingInterval, relax);
                            return null;
                        }));
                    }

                    for (Future<Void> future : readFutures) {
                        future.get();
                    }

                } catch (Exception e) {
                    Logger.error(CLASS_NAME, "Error reading input file " + inputFile + ": " + e.getMessage());
                }
            }

            executor.shutdown();

            // --- PRINT READ STATISTICS SUMMARY ---
            long total = totalReadsProcessed.get();
            long matched = matchedReads.get();
            long unmatched = total - matched;
            double matchPct = total > 0 ? (double) matched / total * 100.0 : 0.0;
            double unmatchPct = total > 0 ? (double) unmatched / total * 100.0 : 0.0;

            Logger.info(CLASS_NAME, "==================================================");
            Logger.info(CLASS_NAME, "               READ COUNTING SUMMARY              ");
            Logger.info(CLASS_NAME, "==================================================");
            Logger.info(CLASS_NAME, "Total reads processed       : " + String.format("%,d", total));
            Logger.info(CLASS_NAME, "Reads with k-mer matching   : " + String.format("%,d (%.2f%%)", matched, matchPct));
            Logger.info(CLASS_NAME, "Unmatched reads             : " + String.format("%,d (%.2f%%)", unmatched, unmatchPct));
            Logger.info(CLASS_NAME, "==================================================");

            // 1. Set up a thread pool for processing the reference fasta records with correct sequence metrics
            ExecutorService writerExecutor = Executors.newFixedThreadPool(nthreads);
            List<Future<String>> writeFutures = new ArrayList<>();

            Logger.info(CLASS_NAME, "Calculating sequence statistics and anchor counts (with skip step = " + skip + ")...");
            int effectiveSkip = Math.max(1, skip);

            try (FastaReader reader = new FastaReader(refFastaFile)) {
                for (Fasta fasta : reader) {
                    writeFutures.add(writerExecutor.submit(() -> {
                        KmerIterator it = fasta.kmers(kmerSize, prefixLength, true);
                        long totalKmersSeq = 0;

                        // Use a Set to track distinct unique canonical k-mers in this sequence
                        java.util.Set<Kmer> uniqueKmersSet = new java.util.HashSet<>();

                        long kmerIndexInSeq = 0;
                        while (it.hasNext()) {
                            Kmer kmer = it.next().getKmer();
                            totalKmersSeq++;

                            // Only add to uniqueKmersSet if it matches the skip step (e.g., 0, 5, 10 -> 1st, 6th, 11th)
                            if (kmerIndexInSeq % effectiveSkip == 0) {
                                uniqueKmersSet.add(kmer);
                            }
                            kmerIndexInSeq++;
                        }

                        long uniqueKmersSeq = uniqueKmersSet.size();
                        long uniqueKmersCountEq1Seq = 0;
                        long coveredKmersSeq = 0;
                        long countedReadsSeq = 0;

                        // Evaluate each unique sampled sequence k-mer against the KMC database
                        for (Kmer kmer : uniqueKmersSet) {
                            long kmerIndex = kmc.getIndex(kmer);
                            if (kmerIndex >= 0) {
                                byte[] entry = kmc.getEntry(kmerIndex);
                                int kmerCount = kmc.getCountFromEntry(entry);

                                if (kmerCount == 1) {
                                    uniqueKmersCountEq1Seq++;

                                    int chunkIndex = (int) (kmerIndex / chunkSize);
                                    int indexInChunk = (int) (kmerIndex % chunkSize);
                                    int readCountForKmer = kmerCountChunks[chunkIndex].get(indexInChunk);

                                    countedReadsSeq += readCountForKmer;

                                    if (readCountForKmer > 0) {
                                        coveredKmersSeq++;
                                    }
                                }
                            }
                        }

                        // Format: id, length, totalKmer, uniqueKmer, uniqueKmer with kmcdbcount==1, coveredKmers, countedReads
                        return String.format("%s\t%d\t%d\t%d\t%d\t%d\t%d",
                                fasta.getName(),
                                fasta.getLength(),
                                totalKmersSeq,
                                uniqueKmersSeq,
                                uniqueKmersCountEq1Seq,
                                coveredKmersSeq,
                                countedReadsSeq);
                    }));
                }
            }

            try (BufferedWriter writer = Files.newBufferedWriter(outputFile)) {
                writer.write("HapID\tLength\tTotalKmers\tUniqueKmers\tUniqueKmersSingleCopy\tCoveredKmers\tNumReads\n");

                for (Future<String> future : writeFutures) {
                    try {
                        String resultLine = future.get();
                        writer.write(resultLine + "\n");
                    } catch (Exception e) {
                        Logger.error(CLASS_NAME, "Error retrieving fasta sequence metrics result: " + e.getMessage());
                    }
                }
            }

            writerExecutor.shutdown();
            Logger.info(CLASS_NAME, "Detailed sequence metrics successfully written to: " + outputFile);
        }
    }

    private void processReadBatch(List<Read> batch, KMC kmc, AtomicIntegerArray[] kmerCountChunks,
                                  int chunkSize, AtomicLong totalReadsProcessed, AtomicLong matchedReads,
                                  long reportingInterval, boolean relax) {
        int kmerSize = kmc.getKmerLength();
        int prefixLength = kmc.getPrefixLength();

        long batchCount = batch.size();
        long currentTotal = totalReadsProcessed.addAndGet(batchCount);

        if ((currentTotal / reportingInterval) > ((currentTotal - batchCount) / reportingInterval)) {
            Logger.info(CLASS_NAME, "Processed " + (currentTotal / 1_000_000) + "M reads...");
        }

        long matchedInBatch = 0;

        for (Read read : batch) {
            boolean matchedThisRead = false;
            KmerIterator it = read.kmers(kmerSize, prefixLength, true);

            while (it.hasNext()) {
                Kmer kmer = it.next().getKmer();
                long kmerIndex = kmc.getIndex(kmer);

                if (kmerIndex >= 0) {
                    byte[] entry = kmc.getEntry(kmerIndex);
                    int kmerCount = kmc.getCountFromEntry(entry);

                    if (kmerCount == 1) {
                        int chunkIndex = (int) (kmerIndex / chunkSize);
                        int indexInChunk = (int) (kmerIndex % chunkSize);
                        kmerCountChunks[chunkIndex].incrementAndGet(indexInChunk);

                        matchedThisRead = true;

                        // skip the rest of the read's kmers if relax is true and a matching kmer is found
                        if (relax) {
                            break;
                        }
                    }
                }
            }

            if (matchedThisRead) {
                matchedInBatch++;
            }
        }

        matchedReads.addAndGet(matchedInBatch);
    }

    private Map<Path, String> parseInputFiles(List<Path> files) {
        if (files == null || files.isEmpty()) {
            Logger.error(CLASS_NAME, "Input files list is empty.");
            return null;
        }

        Map<Path, String> resolvedFilesMap = new HashMap<>();

        for (Path file : files) {
            if (!Files.exists(file) || !Files.isRegularFile(file)) {
                Logger.error(CLASS_NAME, "Input file does not exist or is not a regular file: " + file);
                return null;
            }

            String fileName = file.getFileName().toString().toLowerCase();

            // Check if it's a list file
            if (fileName.endsWith(".list")) {
                try {
                    List<String> lines = Files.readAllLines(file);
                    for (String line : lines) {
                        line = line.trim();
                        if (line.isEmpty() || line.startsWith("#")) {
                            continue;
                        }

                        String[] parts = line.split("\\s+");
                        Path readPath = Path.of(parts[0]);
                        if (!readPath.isAbsolute() && file.getParent() != null) {
                            readPath = file.getParent().resolve(readPath);
                        }

                        if (!Files.exists(readPath) || !Files.isRegularFile(readPath)) {
                            Logger.error(CLASS_NAME, "Input read file from list does not exist: " + readPath);
                            return null;
                        }

                        String resolvedFormat = FileUtils.determineInputFormat(readPath.getFileName().toString().toLowerCase());
                        if (resolvedFormat == null) {
                            Logger.error(CLASS_NAME, "Unsupported file format/extension in list for: " + readPath);
                            return null;
                        }
                        resolvedFilesMap.put(readPath, resolvedFormat);
                    }
                } catch (IOException e) {
                    Logger.error(CLASS_NAME, "Error reading file list: " + file + " - " + e.getMessage());
                    return null;
                }
            } else {
                // Direct file check
                String format = FileUtils.determineInputFormat(fileName);
                if (format == null) {
                    Logger.error(CLASS_NAME, "Unsupported file format extension for input file: " + file);
                    return null;
                }
                resolvedFilesMap.put(file, format);
            }
        }

        return resolvedFilesMap;
    }
}