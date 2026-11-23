package edu.tamu.pg.pang2p.Tasks;

import edu.tamu.pg.pang2p.Data.*;
import edu.tamu.pg.pang2p.Utils.*;
import picocli.CommandLine;
import picocli.CommandLine.*;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.atomic.AtomicInteger;

@Command(name = "kmer-coverage", description = "Calculate k-mer coverage for the given sample", sortOptions = false)
public class KmerCoverage implements Callable<Integer>, Runnable {

    @Option(names = {"-r", "--reference"}, description = "Reference file", required = true)
    private Path referenceFasta;

    @Option(names = {"-s", "--sample"}, description = "Name of the sample", required = true)
    private String sampleName;

    @Option(names = {"-i", "--input"}, description = "Input reads file path(s)", required = true, split = ",", arity = "1..*")
    private List<Path> inputFiles;

    @Option(names = {"-o", "--output"}, description = "Path to the output file", required = true)
    private Path outputFile;

    @Option(names = {"-k", "--kmer"}, description = "K-mer size for KMC", required = false, defaultValue = "31")
    private int kmerSize = 31;

    @Option(names = {"-t", "--threads"}, description = "Number of threads to use for processing [2]", required = false)
    private int threads = 2;

    @Option(names = {"-m", "--min_count"}, description = "Minimum count for k-mer coverage calculation [1]", required = false)
    private int minCount = 1;

    @Option(names = {"--skip"}, description = "Only consider every Nth k-mer of each sequence (1 = every k-mer) [1]", required = false)
    private int skip = 1;

    @Option(names = {"--memory"}, description = "Should the kmc db loaded into memory", required = false)
    private boolean memory = false;

    @Option(names = {"--keep-kmc"}, description = "Keep the KMC database files after processing", required = false)
    private boolean keepKmc = false;

    private final String CLASS_NAME = this.getClass().getSimpleName();
    private static final int FLUSH_INTERVAL = 5000;
    private final AtomicInteger lastLoggedPercent = new AtomicInteger(-1);

    // Default no-arg constructor required by PicoCLI
    public KmerCoverage() {
    }

    // Parameterized constructor for programmatically calling this task
    public KmerCoverage(Path referenceFasta, String sampleName, List<Path> inputFiles, Path outputFile, int kmerSize, int threads, int minCount, int skip, boolean memory){
        this(referenceFasta, sampleName, inputFiles, outputFile, kmerSize, threads, minCount, skip, memory, false);
    }

    // Parameterized constructor for programmatically calling this task
    public KmerCoverage(Path referenceFasta, String sampleName, List<Path> inputFiles, Path outputFile, int kmerSize, int threads, int minCount, int skip, boolean memory, boolean keepKmc) {
        this.referenceFasta = referenceFasta;
        this.sampleName = sampleName;
        this.inputFiles = inputFiles;
        this.outputFile = outputFile;
        this.kmerSize = kmerSize;
        this.threads = threads;
        this.minCount = minCount;
        this.memory = memory;
        this.skip = skip;
        this.keepKmc = keepKmc;
    }

    @Override
    public Integer call() throws Exception {
        HelperFunctions.printCommandLine(new CommandLine(this), CLASS_NAME);
        validateCMD();

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

        Logger.info(CLASS_NAME, "Processing " + totalSequences + " sequences across " + threads + " threads...");

        ExecutorService executor = Executors.newFixedThreadPool(threads);

        KMC kmc = new KMC(outputDir.resolve(sampleName).toString(), memory);
        HelperFunctions.printMemoryUsage("After KMC load");

        lastLoggedPercent.set(-1);
        final int maxPending = Math.max(threads * 4, 1);
        Deque<Future<String>> pending = new ArrayDeque<>();
        int written = 0;

        try (BufferedWriter writer = Files.newBufferedWriter(outputFile);
             FastaReader fastaReader = new FastaReader(referenceFasta)) {

            writer.write("#HapID\tLength\tTotalKmers\tUniqueKmers\tObservedKmers\tObservedTotalKmers\n");

            for (Fasta scr : fastaReader) {
                pending.addLast(executor.submit(() -> processFastaRecord(scr, kmc)));

                // Bound memory: once the window is full, write the oldest result (input order)
                while (pending.size() >= maxPending) {
                    writeNext(pending, writer, ++written, totalSequences);
                }
            }

            // Drain whatever is left, still in input order
            while (!pending.isEmpty()) {
                writeNext(pending, writer, ++written, totalSequences);
            }

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
     * Computes the coverage line for one FASTA record. Does no I/O.
     */
    private String processFastaRecord(Fasta scr, KMC kmc) {
        KmerIterator it = scr.kmers(kmerSize, kmc.getPrefixLength(), kmc.isBothStrands());

        int totalKmers = 0;
        int obsKmers = 0;
        int obsTotalKmers = 0;
        long kmerIndex = 0;

        Set<Kmer> walkedKmers = new HashSet<>();
        final int step = Math.max(skip, 1);

        while (it.hasNext()) {
            KmerPosition kp = it.next();

            if (kmerIndex++ % step != 0) {
                continue;
            }
            Kmer kmer = kp.getKmer();
            totalKmers++;

            if (!walkedKmers.add(kmer)) {
                continue;
            }

            int obsCount = kmc.getCount(kmer);
            if (obsCount >= minCount) {
                // if debug is enabled print fasta_name:total_kmers kp.startPos() - kp.startPos() + kmerSize kmer.toString() obsCount
                if (System.getProperty("DEBUG") != null && System.getProperty("DEBUG").equalsIgnoreCase("true")) {
                    Logger.debug(CLASS_NAME, scr.getName() + ":" + totalKmers + " " + kp.getStartPos() + " - " + (kp.getStartPos() + kmerSize) + " " + kmer.toString() + " " + obsCount);
                }
                obsKmers++;
                obsTotalKmers += obsCount;
            }
        }

        int uniqKmers = walkedKmers.size();
        return scr.getName() + "\t" + scr.getLength() + "\t" + uniqKmers + "\t" + totalKmers
                + "\t" + obsKmers + "\t" + obsTotalKmers + "\n";
    }


    /**
     * Waits for the oldest pending task and writes its line, preserving input order.
     * Runs only on the calling thread, so no synchronization is needed.
     */
    private void writeNext(Deque<Future<String>> pending, BufferedWriter writer,
                           int writtenCount, int totalSequences) throws IOException, InterruptedException {
        Future<String> head = pending.pollFirst();
        try {
            writer.write(head.get());
        } catch (ExecutionException e) {
            Logger.error(CLASS_NAME, "Error during task execution: " + e.getCause());
        }

        if (writtenCount % FLUSH_INTERVAL == 0) {
            writer.flush();
        }

        int percent = (int) ((writtenCount * 100L) / totalSequences);
        if (percent % 5 == 0 && percent > lastLoggedPercent.get()) {
            lastLoggedPercent.set(percent);
            Logger.progress(CLASS_NAME, String.format("Processing %s%d%%", ".".repeat(percent / 5), percent));
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
        if (skip < 1) {
            Logger.error(CLASS_NAME, "Skip should be at least 1");
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