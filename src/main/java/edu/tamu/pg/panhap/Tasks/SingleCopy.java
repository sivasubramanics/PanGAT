package edu.tamu.pg.panhap.Tasks;

import edu.tamu.pg.panhap.Data.*;
import edu.tamu.pg.panhap.Utils.Defaults;
import edu.tamu.pg.panhap.Utils.Logger;
import picocli.CommandLine.*;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.atomic.AtomicInteger;

@Command(name = "SingleCopy", description = "Identify single-copy region sequences from multifasta.", sortOptions = false)
public class SingleCopy implements Callable<Integer>, Runnable {

    @Option(names = {"-i", "--input"}, description = "Input multifasta file", required = true)
    private Path input;

    @Option(names = {"-k", "--kmc"}, description = "Path to the KMC database", required = true)
    private String kmcPrefix;

    @Option(names = {"-o", "--output"}, description = "Output fasta file for single-copy regions", required = true)
    private Path output;

    @Option(names = {"--count"}, description = "Max count for kmer for single-copy region extraction", required = false, defaultValue = "1")
    private int count;

    @Option(names = {"--minlength"}, description = "Minimum length for single-copy region extraction", required = false, defaultValue = "100")
    private int minLength;

    @Option(names = {"-t", "--threads"}, description = "Number of threads to use for processing", required = false, defaultValue = "1")
    private int threads;

    private static final String CLASS_NAME = SingleCopy.class.getSimpleName();

    public SingleCopy() {}

    public SingleCopy(Path input, int count, String kmcPrefix, Path output, int threads) {
        this.input = input;
        this.count = count;
        this.kmcPrefix = kmcPrefix;
        this.output = output;
        this.threads = threads;
    }

    @Override
    public Integer call() throws Exception {
        processGenomeParallel(input, kmcPrefix, count, output, threads);
        return 0;
    }

    @Override
    public void run() {
        try {
            call();
        } catch (Exception e) {
            Logger.error(CLASS_NAME, "Execution failed: " + e.getMessage());
        }
    }

    private void processGenomeParallel(Path input, String kmcPrefix, int maxCount, Path output, int numThreads) {
        ForkJoinPool customThreadPool = null;
        try (FastaWriter writer = new FastaWriter(output, false);
             FastaReader genome = new FastaReader(input)) {

            // Read all chromosomes into a list for parallel stream distribution
            List<Fasta> chromosomes = new ArrayList<>();
            for (Fasta f : genome) {
                chromosomes.add(f);
            }

            KMC kmc = new KMC(kmcPrefix, true);

            customThreadPool = new ForkJoinPool(numThreads);

            // Parallel processing across chromosomes
            customThreadPool.submit(() -> chromosomes.parallelStream().forEach(chrom -> {
                // Thread-local KMC reader to avoid locking overhead
                try {
                    processSingleChromosome(chrom, kmc, maxCount, writer);
                } catch (Exception e) {
                    Logger.error(CLASS_NAME, "Error processing chromosome " + chrom.getName() + ": " + e.getMessage());
                }
            })).get();

        } catch (Exception e) {
            Logger.error(CLASS_NAME, "Error processing single-copy regions: " + e.getMessage());
        } finally {
            if (customThreadPool != null) {
                customThreadPool.shutdown();
            }
        }
    }

    private void processSingleChromosome(Fasta chrom, KMC kmc, int maxCount, FastaWriter writer) throws IOException {
        int kSize = kmc.getKmerLength();
        int prefixLength = kmc.getPrefixLength();
        boolean isCanonical = kmc.isBothStrands();

        StringBuilder singleCopySeq = new StringBuilder(1024);

        int regionStart = -1;
        int currentEnd = -1;
        boolean isBeginning = true;
        AtomicInteger numSingleCopy = new AtomicInteger(0);

        KmerIterator it = chrom.kmers(kSize, prefixLength, isCanonical);

        while (it.hasNext()) {
            KmerPosition kp = it.next();
            int start = kp.getStartPos();
            Kmer kmer = kp.getKmer();

            int kmerCount = kmc.getCount(kmer);

            if (kmerCount <= maxCount && kmerCount > 0) {
                if (isBeginning) {
                    regionStart = start;
                    // Extract directly from parent sequence string to avoid array copies
                    singleCopySeq.append(kmer.getOriginalKmer());
                    isBeginning = false;
                } else {
                    // Append only the trailing base directly from sequence
                    singleCopySeq.append(kmer.getOriginalKmer()[kSize - 1]);
                }
                currentEnd = start + kSize;
            } else {
                if (singleCopySeq.length() >= minLength) {
                    writeRegion(writer, chrom.getName(), numSingleCopy.getAndIncrement(),
                            singleCopySeq.toString(), regionStart, currentEnd);
                }
                singleCopySeq.setLength(0);
                isBeginning = true;
                regionStart = -1;
            }
        }

        // Flush remaining region at chromosome end
        if (singleCopySeq.length() >= minLength) {
            writeRegion(writer, chrom.getName(), numSingleCopy.getAndIncrement(),
                    singleCopySeq.toString(), regionStart, currentEnd);
        }

        Logger.info(CLASS_NAME, "Completed " + numSingleCopy.get() + " single copy regions for chromosome " + chrom.getName());
    }

    private synchronized void writeRegion(FastaWriter writer, String chromName, int index,
                                          String seq, int start, int end) throws IOException {
        Fasta scr = new Fasta(index, chromName + "_" + index, seq, start + "-" + end);
        writer.write(scr, true);
    }
}