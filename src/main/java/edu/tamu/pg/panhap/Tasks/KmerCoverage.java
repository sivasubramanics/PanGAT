package edu.tamu.pg.panhap.Tasks;

import edu.tamu.pg.panhap.Data.*;
import edu.tamu.pg.panhap.Utils.*;
import picocli.CommandLine;
import picocli.CommandLine.*;

import java.io.BufferedWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

@Command(name = "KmerCoverage", description = "Calculate k-mer coverage for the given sample", sortOptions = false)
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
    private int kmerSize;

    @Option(names = {"-t", "--threads"}, description = "Number of threads to use for processing", required = false, defaultValue = "2")
    private int threads;

    @Option(names = {"-m", "--min_count"}, description = "Minimum count for k-mer coverage calculation", required = false, defaultValue = "1")
    private int minCount;

    @Option(names = {"--memory"}, description = "Should the kmc db loaded into memory", defaultValue = "false")
    private boolean memory;

    @Option(names = {"--keep-kmc"}, description = "Keep the KMC database files after processing", defaultValue = "false")
    private boolean keepKmc;

    private final String CLASS_NAME = KmerCoverage.class.getSimpleName();
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

        Path outputDir = outputFile.getParent();
        if (outputDir != null && !Files.exists(outputDir)) {
            Files.createDirectories(outputDir);
        }

        KMCConfig kmcConfig = KMCConfig.builder()
                .inputPaths(inputFiles)
                .outputPrefix(outputDir.resolve(sampleName))
                .kmerSize(kmerSize)
                .threads(threads)
                .minCount(1)
                .build();
        KMCExecutor.run(kmcConfig);

        FastaIndex fastaIndex = new FastaIndex(referenceFasta);
        if (fastaIndex.length() == 0) {
            Logger.info(CLASS_NAME, "No sequences found in reference file: " + referenceFasta);
            return 0;
        }

        Logger.info(CLASS_NAME, "Processing " + fastaIndex.length() + " sequences across " + threads + " threads...");

        ExecutorService executor = Executors.newFixedThreadPool(threads);

        // Semaphore to bound memory usage by queuing only up to (2 * threads) tasks at a time
        Semaphore queueBound = new Semaphore(threads * 2);

        AtomicInteger processedSequences = new AtomicInteger(0);
        AtomicInteger writeCounter = new AtomicInteger(0);

        try (BufferedWriter writer = Files.newBufferedWriter(outputFile);
             KMC kmc = new KMC(outputDir.resolve(sampleName).toString(), memory);
             FastaReader fastaReader = new FastaReader(referenceFasta)) {

            for (Fasta scr : fastaReader) {
                queueBound.acquire();

                executor.submit(() -> {
                    try {
                        processFastaRecord(scr, kmc, writer, writeCounter);
                    } finally {
                        queueBound.release();
                        int completed = processedSequences.incrementAndGet();
                        if (completed % 1000 == 0) {
                            double percentage = ((double) completed / fastaIndex.length()) * 100.0;
                            Logger.progress(CLASS_NAME, String.format("Processed %d/%d sequences (%.2f%%)", completed, fastaIndex.length(), percentage));
                        }
                    }
                });
            }

            // Shutdown executor and wait for running tasks to complete
            executor.shutdown();
            if (!executor.awaitTermination(7, TimeUnit.DAYS)) {
                executor.shutdownNow();
            }

            // Final flush to guarantee zero data loss
            writer.flush();
            System.out.print("\r");
            if (!keepKmc) {
                kmc.rmKMCfiles();
            }
            Logger.info(CLASS_NAME, "Successfully completed k-mer coverage calculation for " + sampleName);

        } catch (Exception e) {
            Logger.error(CLASS_NAME, "Error during KmerCoverage processing: " + e.getMessage());
            executor.shutdownNow();
        }

        return 0;
    }

    /**
     * Processes individual FASTA sequence records concurrently.
     */
    private void processFastaRecord(Fasta scr, KMC kmc, BufferedWriter writer, AtomicInteger writeCounter) {
        KmerIterator it = scr.kmers(kmerSize, kmc.getPrefixLength(), kmc.isBothStrands());

        int totalKmers = 0;
        int obsKmers = 0;
        int obsTotalKmers = 0;

        // Sequence-local HashSet (instantiated & garbage collected per task execution)
        Set<Kmer> walkedKmers = new HashSet<>();

        while (it.hasNext()) {
            KmerPosition kp = it.next();
            Kmer kmer = kp.getKmer();
            totalKmers++;

            // Fast presence check & addition in one step
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

        // Thread-safe writing and batch disk flushing
        synchronized (writer) {
            try {
                writer.write(resultLine);
                if (writeCounter.incrementAndGet() % FLUSH_INTERVAL == 0) {
                    writer.flush();
                }
            } catch (Exception e) {
                throw new RuntimeException("Failed to write output record: " + scr.getName(), e);
            }
        }
    }

    @Override
    public void run() {
        try {
            HelperFunctions.printCommandLine(new CommandLine(this), CLASS_NAME);
            call();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }
}