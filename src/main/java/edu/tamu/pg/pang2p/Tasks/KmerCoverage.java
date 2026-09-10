package edu.tamu.pg.pang2p.Tasks;

import edu.tamu.pg.pang2p.Data.*;
import edu.tamu.pg.pang2p.Utils.*;
import picocli.CommandLine;
import picocli.CommandLine.*;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

@Command(name = "kmer-coverage", description = "Calculate k-mer coverage for the given sample", sortOptions = false)
public class KmerCoverage implements Callable<Integer>, Runnable {

    @Option(names = {"-r", "--reference"}, description = "Reference file", required = true)
    private Path referenceFasta;

    @Option(names = {"-s", "--sample"}, description = "Name of the sample", required = true)
    private String sampleName;

    @Option(names = {"-i", "--input"}, description = "Path to the input files (space separated). format is assumed from the file extension.", required = true)
    private Path[] inputFiles;

    @Option(names = {"-o", "--output"}, description = "Path to the output file", required = true)
    private Path outputFile;

    @Option(names = {"-k", "--kmer"}, description = "K-mer size for KMC", required = false, defaultValue = "31")
    private int kmerSize = 31;

    @Option(names = {"-t", "--threads"}, description = "Number of threads to use for processing [2]", required = false)
    private int threads = 2;

    @Option(names = {"-m", "--min_count"}, description = "Minimum count for k-mer coverage calculation [1]", required = false)
    private int minCount = 1;

    @Option(names = {"--memory"}, description = "Should the kmc db loaded into memory", required = false)
    private boolean memory = false;

    @Option(names = {"--keep-kmc"}, description = "Keep the KMC database files after processing", required = false)
    private boolean keepKmc = false;

    private final String CLASS_NAME = this.getClass().getSimpleName();
    private static final int FLUSH_INTERVAL = 5000;

    // Default no-arg constructor required by PicoCLI
    public KmerCoverage() {
    }

    // Parameterized constructor for programmatically calling this task
    public KmerCoverage(Path referenceFasta, String sampleName, Path[] inputFiles, Path outputFile, int kmerSize, int threads, int minCount, boolean memory) {
        this.referenceFasta = referenceFasta;
        this.sampleName = sampleName;
        this.inputFiles = inputFiles;
        this.outputFile = outputFile;
        this.kmerSize = kmerSize;
        this.threads = threads;
        this.minCount = minCount;
        this.memory = memory;
    }

    @Override
    public Integer call() throws Exception {
        HelperFunctions.printCommandLine(new CommandLine(this), CLASS_NAME);
        validateCMD();

        DependencyChecker dependencyChecker = new DependencyChecker();
        dependencyChecker.check();

        calculateCoverage();
        return 0;
    }

    @Override
    public void run() {
        try {
            call();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    /***
     * Main function to execute the KMC database generation and k-mer coverage analysis
     */
    public void calculateCoverage() throws IOException {
        sampleName = cleanSampleName(sampleName);

        Path outputDir = outputFile.toAbsolutePath().getParent();
        if (outputDir != null && !Files.exists(outputDir)) {
            Files.createDirectories(outputDir);
        }

        HelperFunctions.printMemoryUsage("Before KMC run");
        KMCConfig kmcConfig = KMCConfig.builder()
                .inputPaths(inputFiles)
                .outputPrefix(outputDir.resolve(sampleName))
                .kmerSize(kmerSize)
                .threads(threads - 1)
                .minCount(1)
                .build();
        KMCExecutor.run(kmcConfig);
        HelperFunctions.printMemoryUsage("After KMC run");

        FastaIndex fastaIndex = new FastaIndex(referenceFasta);
        int totalSequences = fastaIndex.length();
        if (totalSequences == 0) {
            Logger.info(CLASS_NAME, "No sequences found in reference file: " + referenceFasta);
            fastaIndex.close();
            return;
        }
        HelperFunctions.printMemoryUsage("After FastaIndex");

        Logger.info(CLASS_NAME, "Processing " + totalSequences + " sequences across " + threads + " threads...");

        ExecutorService executor = Executors.newFixedThreadPool(threads);
        ExecutorCompletionService<Void> completionService = new ExecutorCompletionService<>(executor);

        AtomicInteger completedSequences = new AtomicInteger(0);
        AtomicInteger writeCounter = new AtomicInteger(0);

        KMC kmc = new KMC(outputDir.resolve(sampleName).toString(), memory);
        HelperFunctions.printMemoryUsage("After KMC load");

        try (BufferedWriter writer = Files.newBufferedWriter(outputFile);
             FastaReader fastaReader = new FastaReader(referenceFasta)) {

            for (Fasta scr : fastaReader) {
                completionService.submit(() -> {
                    processFastaRecord(scr, kmc, writer, writeCounter);
                    int completed = completedSequences.incrementAndGet();
                    float progress = (float) (completed * 100) / totalSequences;
                    synchronized (System.out) {
                        System.out.printf("\rProgress: %.2f%%", progress);
                    }
                    return null;
                });
            }

            for (int i = 0; i < totalSequences; i++) {
                try {
                    completionService.take().get();
                } catch (Exception e) {
                    Logger.error(CLASS_NAME, "Error during task execution: " + e.getMessage());
                }
            }

            // Clear progress line on console
            System.out.print("\r");
            for (int i = 0; i < 100; i++) {
                System.out.print(" ");
            }
            System.out.print("\r");
            System.out.flush();

            writer.flush();
            Logger.info(CLASS_NAME, "Successfully completed k-mer coverage calculation for " + sampleName);
        } catch (Exception e) {
            Logger.error(CLASS_NAME, "Error during KmerCoverage processing: " + e.getMessage());
        } finally {
            executor.shutdownNow();
        }

        HelperFunctions.printMemoryUsage("After KmerCoverage");

        if (!keepKmc) {
            kmc.rmKMCfiles();
        }
    }

    /**
     * Processes individual FASTA sequence records concurrently.
     */
    private void processFastaRecord(Fasta scr, KMC kmc, BufferedWriter writer, AtomicInteger writeCounter) {
        KmerIterator it = scr.kmers(kmerSize, kmc.getPrefixLength(), kmc.isBothStrands());

        int totalKmers = 0;
        int obsKmers = 0;
        int obsTotalKmers = 0;

        Set<Kmer> walkedKmers = new HashSet<>();

        while (it.hasNext()) {
            KmerPosition kp = it.next();
            Kmer kmer = kp.getKmer();
            totalKmers++;

            if (!walkedKmers.add(kmer)) {
                continue;
            }

            int obsCount = kmc.getCount(kmer);
            if (obsCount >= minCount) {
                obsKmers++;
                obsTotalKmers += obsCount;
            }
        }

        int uniqKmers = walkedKmers.size();
        String resultLine = scr.getName() + "\t" + scr.getLength() + "\t" + uniqKmers + "\t" + totalKmers + "\t" + obsKmers + "\t" + obsTotalKmers + "\n";

        synchronized (writer) {
            try {
                writer.write(resultLine);
                if (writeCounter.incrementAndGet() % FLUSH_INTERVAL == 0) {
                    writer.flush();
                }
            } catch (Exception e) {
                Logger.error(CLASS_NAME, "Failed to write output record for sequence: " + scr.getName());
            }
        }
    }

    /***
     * Validate command line input options
     */
    private void validateCMD() {
        if (threads <= 1) {
            Logger.error(CLASS_NAME, "Number of threads should be greater than 1");
        }
        if (minCount < 1) {
            Logger.error(CLASS_NAME, "Minimum kmer count should be at least 1");
        }
        if (kmerSize <= 0 && kmerSize % 2 == 0) {
            Logger.error(CLASS_NAME, "K-mer size should be greater than 0 and odd");
        }
    }

    private String cleanSampleName(String sampleName) {
        String sanitized = sampleName.replaceAll("[\\\\/:*?\"<>|]", "_");
        if (!sanitized.equals(sampleName)) {
            Logger.warning(CLASS_NAME, "Sample name contains invalid characters, changed to: " + sanitized);
        }
        return sanitized;
    }
}