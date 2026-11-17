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

    @Option(names = {"-f", "--format"}, description = "Input reads file format (fq or fa or list)", required = true)
    private String inputFileFormat;

    @Option(names = {"-t", "--threads"}, description = "Number of threads to use", required = false)
    private int nthreads = 1;

    @Option(names = {"-o", "--output"}, description = "Output file path", required = true)
    private String outputFile;

    @Option(names = {"--relax"}, description = "Skip remaining kmers in a read once a matching kmer is found", required = false)
    private boolean relax = false;

    private String CLASS_NAME = ReadCounts.class.getSimpleName();

    public ReadCounts() {}

    public ReadCounts(Path refFastaFile, int kmerSize, List<Path> inputFiles, String inputFileFormat, int nthreads, String outputFile, boolean relax) {
        this.refFastaFile = refFastaFile;
        this.kmerSize = kmerSize;
        this.inputFiles = inputFiles;
        this.inputFileFormat = inputFileFormat;
        this.nthreads = nthreads;
        this.outputFile = outputFile;
        this.relax = relax;
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

        Map<Path, String> validFilesMap = parseInputFiles(inputFiles, inputFileFormat);
        if (validFilesMap == null) {
            Logger.error(CLASS_NAME, "Input files are not valid. Please check the file paths and formats.");
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
                                processReadBatch(batchToProcess, kmc, kmerCountChunks, chunkSize, totalReadsProcessed, reportingInterval, relax);
                                return null;
                            }));
                            batch.clear();
                        }
                    }

                    if (!batch.isEmpty()) {
                        List<Read> batchToProcess = new ArrayList<>(batch);
                        readFutures.add(executor.submit(() -> {
                            processReadBatch(batchToProcess, kmc, kmerCountChunks, chunkSize, totalReadsProcessed, reportingInterval, relax);
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
            Logger.info(CLASS_NAME, "Finished processing all input files sequentially. Total reads processed: " + totalReadsProcessed.get());

            // 1. Set up a thread pool for processing the reference fasta records with correct sequence metrics
            ExecutorService writerExecutor = Executors.newFixedThreadPool(nthreads);
            List<Future<String>> writeFutures = new ArrayList<>();

            Logger.info(CLASS_NAME, "Calculating sequence statistics and anchor counts...");
            try (FastaReader reader = new FastaReader(refFastaFile)) {
                for (Fasta fasta : reader) {
                    writeFutures.add(writerExecutor.submit(() -> {
                        KmerIterator it = fasta.kmers(kmerSize, prefixLength, true);
                        long totalKmersSeq = 0;

                        // Use a Set to track distinct unique canonical k-mers in this sequence
                        java.util.Set<Kmer> uniqueKmersSet = new java.util.HashSet<>();

                        while (it.hasNext()) {
                            Kmer kmer = it.next().getKmer();
                            totalKmersSeq++;
                            uniqueKmersSet.add(kmer); // HashSet deduplicates identical canonical k-mers
                        }

                        long uniqueKmersSeq = uniqueKmersSet.size();
                        long uniqueKmersCountEq1Seq = 0;
                        long countedReadsSeq = 0;

                        // Evaluate each unique sequence k-mer against the KMC database
                        for (Kmer kmer : uniqueKmersSet) {
                            long kmerIndex = kmc.getIndex(kmer);
                            if (kmerIndex >= 0) {
                                byte[] entry = kmc.getEntry(kmerIndex);
                                int kmerCount = kmc.getCountFromEntry(entry);

                                if (kmerCount == 1) {
                                    uniqueKmersCountEq1Seq++;

                                    int chunkIndex = (int) (kmerIndex / chunkSize);
                                    int indexInChunk = (int) (kmerIndex % chunkSize);
                                    countedReadsSeq += kmerCountChunks[chunkIndex].get(indexInChunk);
                                }
                            }
                        }

                        // Format: id, length, totalKmer, uniqueKmer, uniqueKmer with kmcdbcount==1, countedReads
                        return String.format("%s\t%d\t%d\t%d\t%d\t%d",
                                fasta.getName(),
                                fasta.getLength(),
                                totalKmersSeq,
                                uniqueKmersSeq,
                                uniqueKmersCountEq1Seq,
                                countedReadsSeq);
                    }));
                }
            }

            try (BufferedWriter writer = Files.newBufferedWriter(Path.of(outputFile))) {
                // Optional header for clarity
                writer.write("SequenceID\tLength\tTotalKmers\tUniqueKmers\tUniqueKmersSingleCopy\tNumReads\n");

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
                                  int chunkSize, AtomicLong totalReadsProcessed, long reportingInterval, boolean relax) {
        int kmerSize = kmc.getKmerLength();
        int prefixLength = kmc.getPrefixLength();

        long batchCount = batch.size();
        long currentTotal = totalReadsProcessed.addAndGet(batchCount);

        if ((currentTotal / reportingInterval) > ((currentTotal - batchCount) / reportingInterval)) {
            Logger.info(CLASS_NAME, "Processed " + (currentTotal / 1_000_000) + "M reads...");
        }

        for (Read read : batch) {
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

                        // skip the read if relax is true and a matching kmer is found
                        if (relax) {
                            break;
                        }
                    }
                }
            }
        }
    }

    private Map<Path, String> parseInputFiles(List<Path> files, String defaultFormat) {
        if (files == null || files.isEmpty()) {
            Logger.error(CLASS_NAME, "Input files list is empty.");
            return null;
        }

        Map<Path, String> resolvedFilesMap = new HashMap<>();

        if (defaultFormat.equalsIgnoreCase("list")) {
            for (Path listFile : files) {
                if (!Files.exists(listFile)) {
                    Logger.error(CLASS_NAME, "Input list file not found: " + listFile);
                    return null;
                }
                try {
                    List<String> lines = Files.readAllLines(listFile);
                    for (String line : lines) {
                        line = line.trim();
                        if (line.isEmpty() || line.startsWith("#")) {
                            continue;
                        }

                        String[] parts = line.split("\\s+");
                        Path readPath = Path.of(parts[0]);
                        if (!readPath.isAbsolute() && listFile.getParent() != null) {
                            readPath = listFile.getParent().resolve(readPath);
                        }

                        if (FileUtils.isFastaFile(readPath)) {
                            resolvedFilesMap.put(readPath, "fasta");
                        } else if (FileUtils.isFastqFile(readPath)) {
                            resolvedFilesMap.put(readPath, "fastq");
                        } else {
                            Logger.error(CLASS_NAME, "Unsupported file format for input file: " + readPath);
                            return null;
                        }
                    }
                } catch (IOException e) {
                    Logger.error(CLASS_NAME, "Error reading file list: " + listFile + " - " + e.getMessage());
                    return null;
                }
            }
        } else {
            for (Path file : files) {
                if (FileUtils.isFastaFile(file)) {
                    resolvedFilesMap.put(file, "fasta");
                } else if (FileUtils.isFastqFile(file)) {
                    resolvedFilesMap.put(file, "fastq");
                } else {
                    Logger.error(CLASS_NAME, "Unsupported file format for input file: " + file);
                    return null;
                }
            }
        }

        for (Map.Entry<Path, String> entry : resolvedFilesMap.entrySet()) {
            Path file = entry.getKey();
            if (!Files.exists(file) || !Files.isRegularFile(file)) {
                Logger.error(CLASS_NAME, "Input read file does not exist or is invalid: " + file);
                return null;
            }
        }

        return resolvedFilesMap;
    }
}