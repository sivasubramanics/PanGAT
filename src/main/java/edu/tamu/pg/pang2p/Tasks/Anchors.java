package edu.tamu.pg.pang2p.Tasks;

import edu.tamu.pg.pang2p.Data.*;
import edu.tamu.pg.pang2p.Utils.*;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

@Command(
        name = "anchors",
        description = "Extract single copy regions (Anchors) for a single genome with reverse method.",
        mixinStandardHelpOptions = true,
        sortOptions = false
)
public class Anchors implements Callable<Integer>, Runnable {

    private static final String CLASS_NAME = Anchors.class.getSimpleName();

    @Option(names = {"-i", "--input"}, description = "Input FASTA file for the genome", required = true)
    private Path input;

    @Option(names = {"-o", "--output"}, description = "Output FASTA file for extracted single-copy regions", required = true)
    private Path output;

    @Option(names = {"-k", "--kmer-size"}, description = "K-mer size for KMC [default: 31]", defaultValue = "31")
    private int kmerSize;

    @Option(names = {"--min-count"}, description = "Minimum k-mer count cutoff to not consider for anchors construction [default: 2]", defaultValue = "2")
    private int minCount;

    @Option(names = {"--min-length"}, description = "Minimum sequence length for single copy anchors [default: 100]", defaultValue = "100")
    private int minLength;

    @Option(names = {"-t", "--threads"}, description = "Number of threads [default: 2]", defaultValue = "2")
    private int threads;

    @Option(names = {"--max-gap"}, description = "Maximum gap size in bases to merge anchors [default: 50]", defaultValue = "50")
    private int maxGap;
    public Anchors() {}

    public Anchors(Path input, Path output, int kmerSize, int minCount, int minLength, int maxGap, int threads) {
        this.input = input;
        this.output = output;
        this.kmerSize = kmerSize;
        this.minCount = minCount;
        this.minLength = minLength;
        this.maxGap = maxGap;
        this.threads = threads;
    }

    @Override
    public Integer call() throws Exception {
        extractSingleCopyRegions(input, output, kmerSize, minCount, minLength, maxGap, threads);
        return 0;
    }

    @Override
    public void run() {
        try {
            HelperFunctions.printCommandLine(new picocli.CommandLine(this), CLASS_NAME);
            call();
        } catch (Exception e) {
            Logger.error(CLASS_NAME, "Execution failed: " + e.getMessage());
        }
    }

    private void extractSingleCopyRegions(Path input, Path output, int kmerSize, int minCount, int minLength, int maxGap, int threads) throws IOException {
        validateInputs();

        List<Path> tempFastaFiles = new ArrayList<>();
        List<Path> tempBedFiles = new ArrayList<>();
        Path kmcPrefix;
        switch (FileUtils.getExtension(input.toString())) {
            case ".fasta":
                kmcPrefix = Path.of(FileUtils.replaceExtension(input, ".fasta", ""));
                break;
            case ".fa":
                kmcPrefix = Path.of(FileUtils.replaceExtension(input, ".fa", ""));
                break;
            case ".fna":
                kmcPrefix = Path.of(FileUtils.replaceExtension(input, ".fna", ""));
                break;
            default:
                throw new IllegalArgumentException("Output file must have .fasta or .fa or .fna extension.");
        }

        KMCConfig kmcConfig = KMCConfig.builder()
                .inputPath(input)
                .outputPrefix(kmcPrefix)
                .kmerSize(kmerSize)
                .minCount(minCount)
                .threads(threads)
                .maxMemory((int) ((Runtime.getRuntime().maxMemory() / (1024.0 * 1024.0 * 1024.0))+1))
                .build();

        Logger.info(CLASS_NAME, "Running KMC for genome: " + input.getFileName());
        KMCExecutor.run(kmcConfig);

        try (FastaReader fastaReader = new FastaReader(input.toString());
             KMC kmc = new KMC(kmcPrefix, true)) {

            ExecutorService executor = Executors.newFixedThreadPool(threads);
            List<Future<?>> futures = new ArrayList<>();

            for (Fasta fasta : fastaReader) {
                Path tempFastaFile;
                Path tempBedFile;
                switch (FileUtils.getExtension(output.toString())) {
                    case ".fasta":
                        tempFastaFile = Path.of(FileUtils.replaceExtension(output, ".fasta", "." + fasta.getName() + ".tmp"));
                        tempBedFile = Path.of(FileUtils.replaceExtension(output, ".fasta", "." + fasta.getName() + ".tmp.bed"));
                        break;
                    case ".fa":
                        tempFastaFile = Path.of(FileUtils.replaceExtension(output, ".fa", "." + fasta.getName() + ".tmp"));
                        tempBedFile = Path.of(FileUtils.replaceExtension(output, ".fa", "." + fasta.getName() + ".tmp.bed"));
                        break;
                    case ".fna":
                        tempFastaFile = Path.of(FileUtils.replaceExtension(output, ".fna", "." + fasta.getName() + ".tmp"));
                        tempBedFile = Path.of(FileUtils.replaceExtension(output, ".fna", "." + fasta.getName() + ".tmp.bed"));
                        break;
                    default:
                        throw new IllegalArgumentException("Output file must have .fasta or .fa or .fna extension.");
                }
                tempFastaFiles.add(tempFastaFile);
                tempBedFiles.add(tempBedFile);

                futures.add(executor.submit(() -> processChrom(fasta, kmc, tempFastaFile, tempBedFile, minCount, minLength, maxGap)));
            }

            executor.shutdown();
            for (Future<?> future : futures) {
                try {
                    future.get();
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                    Logger.error(CLASS_NAME, "Thread interrupted: " + e.getMessage());
                } catch (ExecutionException e) {
                    Logger.error(CLASS_NAME, "Task execution failed: " + e.getCause().getMessage());
                }
            }
        }

        // Clean up residual progress bar line when all tasks complete
        Logger.clearProgress();

        FileUtils.catFiles(tempFastaFiles, output);
        FileUtils.catFiles(tempBedFiles, Path.of(FileUtils.replaceExtension(output, FileUtils.getExtension(output.toString()), ".bed")));
        FileUtils.rmFiles(tempFastaFiles);
        FileUtils.rmFiles(tempBedFiles);
        FileUtils.rmFile(Path.of(kmcPrefix.toString() + ".kmc_pre"));
        FileUtils.rmFile(Path.of(kmcPrefix.toString() + ".kmc_suf"));
        FileUtils.rmFile(Path.of(kmcPrefix.toString() + ".kmc.log"));
    }

    private void processChrom(Fasta fasta, KMC kmc, Path tempFastaFile, Path tempBedFile, int minCount, int minLength, int maxGap) {
        String chrom = fasta.getName();
        int kmerLen = kmc.getKmerLength();
        KmerIterator it = fasta.kmers(kmerLen, kmc.getPrefixLength(), kmc.isBothStrands());

        long totalKmers = Math.max(1, fasta.getSequence().length() - kmerLen + 1);
        long processedKmers = 0;
        int lastLoggedPercent = -1;

        int start = -1;
        int end = -1;
        int lastKmerPos = -1;

        List<Region> regions = new ArrayList<>();
        String sequence = fasta.getSequence();

        while (it.hasNext()) {
            KmerPosition kp = it.next();
            processedKmers++;

            // Progress logging
            if (processedKmers % 10000 == 0 || processedKmers == totalKmers) {
                int percent = (int) ((processedKmers * 100) / totalKmers);
                if (percent % 5 == 0 && percent != lastLoggedPercent) {
                    lastLoggedPercent = percent;
                    String dots = ".".repeat(percent / 5);
                    Logger.progress(CLASS_NAME, String.format("Processing %s %s%d%%", chrom, dots, percent));
                }
            }

            int count = kmc.getCount(kp.getKmer());
            boolean isSingleCopy = (count < minCount);

            if (isSingleCopy) {
                if (start == -1) {
                    // Start a new region
                    start = kp.getStartPos();
                } else if (kp.getStartPos() > lastKmerPos + 1) {
                    // Invalid k-mer or gap encountered: force-split the single copy region
                    if (end - start >= minLength) {
                        regions.add(new Region(start, end));
                    }
                    start = kp.getStartPos();
                }
                end = kp.getStartPos() + kmerLen;
            } else {
                // High count k-mer: close current region if active
                if (start != -1 && end != -1) {
                    if (end - start >= minLength) {
                        regions.add(new Region(start, end));
                    }
                    start = -1;
                    end = -1;
                }
            }
            lastKmerPos = kp.getStartPos();
        }

        // Handle the last region if it exists
        if (start != -1 && end != -1) {
            if (end - start >= minLength) {
                regions.add(new Region(start, end));
            }
        }

        Logger.info(CLASS_NAME, "Found " + regions.size() + " single copy regions for chromosome " + chrom + " before merging (minLength: " + minLength + ").");

        // Sort regions by start position
        regions.sort(Comparator.comparingInt(r -> r.start));

        // Merge regions separated by less than or equal to maxGap
        List<Region> mergedRegions = mergeRegions(regions, maxGap);

        Logger.info(CLASS_NAME, "Found " + mergedRegions.size() + " single copy regions for chromosome " + chrom + " after merging (maxGap: " + maxGap + ").");

        // Write the final merged regions to the file
        try (FastaWriter writer = new FastaWriter(tempFastaFile.toString(), false);
        BufferedWriter bw = new BufferedWriter(new FileWriter(tempBedFile.toString(), false))) {
            int numRegions = 0;
            for (Region r : mergedRegions) {
                writer.write(new Fasta(numRegions, chrom + "_" + numRegions, sequence.substring(r.start, r.end), chrom + ":" + r.start + "-" + r.end), true);
                bw.write(chrom + "\t" + r.start + "\t" + r.end + "\t" + chrom + "_" + numRegions + "\n");
                numRegions++;
            }
        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Error writing regions for sequence " + fasta.getName() + ": " + e.getMessage());
        }
    }

    private List<Region> mergeRegions(List<Region> regions, int maxGap) {
        List<Region> merged = new ArrayList<>();
        if (regions.isEmpty()) {
            return merged;
        }

        Region current = regions.get(0);
        for (int i = 1; i < regions.size(); i++) {
            Region next = regions.get(i);

            int baseGap = next.start - current.end;

            if (baseGap <= maxGap) {
                current.end = Math.max(current.end, next.end);
            } else {
                merged.add(current);
                current = next;
            }
        }
        merged.add(current);
        return merged;
    }

    // Helper class for region representation
    private static class Region {
        int start;
        int end;

        Region(int start, int end) {
            this.start = start;
            this.end = end;
        }
    }


    private void validateInputs() {
        if (!FileUtils.isFile(input)) {
            Logger.error(CLASS_NAME, "Input file does not exist: " + input);
        }
        if (kmerSize <= 0 || kmerSize % 2 == 0) {
            Logger.error(CLASS_NAME, "K-mer size must be a positive odd integer.");
        }
        if (minCount < 0) {
            Logger.error(CLASS_NAME, "Min count must be a positive integer.");
        }
        if (minLength <= 0) {
            Logger.error(CLASS_NAME, "Min length must be a positive integer.");
        }
        if (threads <= 0) {
            Logger.error(CLASS_NAME, "Number of threads must be a positive integer.");
        }
    }
}