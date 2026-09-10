package edu.tamu.pg.pang2p.Tasks;

import edu.tamu.pg.pang2p.Data.*;
import edu.tamu.pg.pang2p.Utils.*;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

@Command(
        name = "anchors",
        description = "Extract single copy regions (Anchors) for a single genome.",
        mixinStandardHelpOptions = true,
        sortOptions = false
)
public class Anchors implements Callable<Integer>, Runnable {

    private static final String CLASS_NAME = Anchors.class.getSimpleName();
    private static final Object PRINT_LOCK = new Object(); // Synchronization lock for terminal output

    @Option(names = {"-i", "--input"}, description = "Input FASTA file for the genome", required = true)
    private Path input;

    @Option(names = {"-o", "--output"}, description = "Output FASTA file for extracted single-copy regions", required = true)
    private Path output;

    @Option(names = {"-k", "--kmer-size"}, description = "K-mer size for KMC [default: 31]", defaultValue = "31")
    private int kmerSize;

    @Option(names = {"--min-count"}, description = "Minimum k-mer count cutoff [default: 1]", defaultValue = "1")
    private int minCount;

    @Option(names = {"--min-length"}, description = "Minimum sequence length for single copy anchors [default: 100]", defaultValue = "100")
    private int minLength;

    @Option(names = {"-t", "--threads"}, description = "Number of threads [default: 2]", defaultValue = "2")
    private int threads;

    public Anchors() {}

    public Anchors(Path input, Path output, int kmerSize, int minCount, int minLength, int threads) {
        this.input = input;
        this.output = output;
        this.kmerSize = kmerSize;
        this.minCount = minCount;
        this.minLength = minLength;
        this.threads = threads;
    }

    @Override
    public Integer call() throws Exception {
        extractSingleCopyRegions(input, output, kmerSize, minCount, minLength, threads);
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

    private void extractSingleCopyRegions(Path input, Path output, int kmerSize, int minCount, int minLength, int threads) throws IOException {
        validateInputs();

        List<Path> tempFiles = new ArrayList<>();
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
                .minCount(0)
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
                Path tempFile;
                switch (FileUtils.getExtension(output.toString())) {
                    case ".fasta":
                        tempFile = Path.of(FileUtils.replaceExtension(output, ".fasta", "." + fasta.getName() + ".tmp"));
                        break;
                    case ".fa":
                        tempFile = Path.of(FileUtils.replaceExtension(output, ".fa", "." + fasta.getName() + ".tmp"));
                        break;
                    case ".fna":
                        tempFile = Path.of(FileUtils.replaceExtension(output, ".fna", "." + fasta.getName() + ".tmp"));
                        break;
                    default:
                        throw new IllegalArgumentException("Output file must have .fasta or .fa extension.");
                }
                tempFiles.add(tempFile);

                futures.add(executor.submit(() -> processChrom(fasta, kmc, tempFile, minCount, minLength)));
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

        FileUtils.catFiles(tempFiles, output);
        FileUtils.rmFiles(tempFiles);
        FileUtils.rmFile(Path.of(kmcPrefix.toString() + ".kmc_pre"));
        FileUtils.rmFile(Path.of(kmcPrefix.toString() + ".kmc_suf"));
        FileUtils.rmFile(Path.of(kmcPrefix.toString() + ".kmc.log"));
    }


    private void processChrom(Fasta fasta, KMC kmc, Path tempFile, int minCount, int minLength) {
        try (FastaWriter writer = new FastaWriter(tempFile.toString())) {
            String chrom = fasta.getName();
            KmerIterator it = fasta.kmers(kmc.getKmerLength(), kmc.getPrefixLength(), kmc.isBothStrands());

            long totalKmers = Math.max(1, fasta.getSequence().length() - kmc.getKmerLength() + 1);
            long processedKmers = 0;
            int lastLoggedPercent = -1;

            int anchorStart = -1;
            int anchorEnd = -1;
            boolean isAnchor = false;
            int numHaptigs = 0;

            while (it.hasNext()) {
                KmerPosition kp = it.next();
                processedKmers++;

                int percent = (int) ((processedKmers * 100) / totalKmers);
                if (percent % 5 == 0 && percent != lastLoggedPercent) {
                    lastLoggedPercent = percent;
                    int dotsCount = percent / 5;
                    String dots = ".".repeat(dotsCount);
                    String progressMsg = String.format("Processing %s %s%d%%", chrom, dots, percent);

                    // Emits progress using Logger
                    Logger.progress(CLASS_NAME, progressMsg);
                }

                int start = kp.getStartPos();
                Kmer kmer = kp.getKmer();
                int count = kmc.getCount(kmer);

                if (count > 0 && count <= minCount) {
                    if (!isAnchor) {
                        anchorStart = start;
                        isAnchor = true;
                    }
                    anchorEnd = start + kmer.getKmerLength();
                } else {
                    if (isAnchor) {
                        int anchorLength = anchorEnd - anchorStart;
                        if (anchorLength >= minLength) {
                            numHaptigs++;
                            String anchorsSeq = fasta.getSubSequence(anchorStart, anchorEnd);
                            writer.write(new Fasta(numHaptigs, chrom + "_" + numHaptigs, anchorsSeq, chrom + ":" + anchorStart + "-" + anchorEnd), true);
                        }
                        isAnchor = false;
                    }
                }
            }

            if (isAnchor) {
                int anchorLength = anchorEnd - anchorStart;
                if (anchorLength >= minLength) {
                    numHaptigs++;
                    String anchorsSeq = fasta.getSubSequence(anchorStart, anchorEnd);
                    writer.write(new Fasta(numHaptigs, chrom + "_" + numHaptigs, anchorsSeq, chrom + ":" + anchorStart + "-" + anchorEnd), true);
                }
            }

            Logger.info(CLASS_NAME, "Completed processing sequence " + fasta.getName() + ". Total anchors found: " + numHaptigs);
        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Error processing sequence " + fasta.getName() + ": " + e.getMessage());
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