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

@Command(name = "SingleCopy", description = "Identify single-copy region sequences from multifasta.", sortOptions = false)
public class SingleCopy implements Callable<Integer>, Runnable {
    // input fasta
    @Option(names = {"-i", "--input"}, description = "Input multifasta file", required = true)
    private Path input;

    // kmc db
    @Option(names = {"-k", "--kmc"}, description = "Path to the KMC database", required = true)
    private String kmcPrefix;

    // output fasta
    @Option(names = {"-o", "--output"}, description = "Output fasta file for single-copy regions", required = true)
    private Path output;

    // minimum count
    @Option(names = {"--count"}, description = "Max count for kmer for single-copy region extraction", required = false, defaultValue = "1")
    private int count;

    // min length
    @Option(names = {"--minlength"}, description = "Minimum length for single-copy region extraction", required = false, defaultValue = "100")
    private int minLength;

    // number of threads
    @Option(names = {"-t", "--threads"}, description = "Number of threads to use for processing", required = false, defaultValue = "1")
    private int threads;

    private static final String CLASS_NAME = SingleCopy.class.getSimpleName();
    private static final int dumpSize = Defaults.COUNTER_DUMP_SIZE;

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
        processGenome(input, kmcPrefix, count, output, threads);
        return 0;
    }

    @Override
    public void run() {
        try {
            call();
        }
        catch (Exception e) {
            e.printStackTrace();
        }

    }

    private void processGenome(Path input, String kmcPrefix, int count, Path output, int threads) {
        // Larger buffer: 256KB instead of the default 8KB
        try (FastaWriter writer = new FastaWriter(output, false);
             KMC kmc = new KMC(kmcPrefix, true);
            FastaReader genome = new FastaReader(input)) {

            int kSize = kmc.getKmerLength();
            int prefixLength = kmc.getPrefixLength();
            boolean isCanonical = kmc.isBothStrands();
            for (Fasta chrom:genome) {
                int counter = 0;
                List<Fasta> singleCopyRegions = new ArrayList<>();
                int numSingleCopy = 0;
                StringBuilder singleCopySeq = new StringBuilder();
                boolean isBeginning = true;
                int regionStart = -1;     // track where the region began
                int currentEnd = -1;      // track the latest position seen

                KmerIterator it = chrom.kmers(kSize, prefixLength, isCanonical);

                while (it.hasNext()) {
                    KmerPosition kp = it.next();
                    int start = kp.getStartPos();
                    Kmer kmer = kp.getKmer();

                    int kmerCount = kmc.getCount(kmer);
                    if (kmerCount == 0) {
                        Logger.error(CLASS_NAME, "kmer :" + kmer.toString() + " not found in the kmc database. Make sure the input fasta was used to generate the kmc database.");
                    }
                    if (kmerCount <= count) {
                        if (isBeginning) {
                            regionStart = start;
                            singleCopySeq.append(kmer.getOriginalKmer());
                            isBeginning = false;
                        } else {
                            // Append only the new last char — avoids a full String allocation
                            singleCopySeq.append(kmer.getOriginalKmer()[kSize - 1]);
                        }
                        currentEnd = start + kSize;
                    } else {
                        if (!singleCopySeq.isEmpty()) {
                            if (singleCopySeq.length() >= minLength) {
                                singleCopyRegions.add(new Fasta(numSingleCopy, chrom.getName() + "" + numSingleCopy, singleCopySeq.toString(), regionStart + "-" + currentEnd));
                                if (counter == dumpSize) {
                                    for (Fasta scr : singleCopyRegions) {
                                        writer.write(scr, false);
                                    }
                                    singleCopyRegions.clear();
                                    counter = 0;
                                    Logger.info(CLASS_NAME, "Dumped " + numSingleCopy + " single copy regions for chromsome " + chrom.getName());
                                }
                                counter++;
                                numSingleCopy++;
                            }
                            singleCopySeq.setLength(0);
                        }
                        isBeginning = true;
                        regionStart = -1;
//                        it.skipTo(start + kSize);  // Skip to the next kmer position
                    }
                }

                // ── Flush any region still open at end of chromosome ──
                if (!singleCopySeq.isEmpty()) {
                    if (singleCopySeq.length() >= minLength) {
                        singleCopyRegions.add(new Fasta(numSingleCopy, chrom.getName() + "_" + numSingleCopy, singleCopySeq.toString(), regionStart + "-" + currentEnd));
                        numSingleCopy++;
                    }
                    singleCopySeq.setLength(0);
                }

                if (!singleCopyRegions.isEmpty()) {
                    for (Fasta scr : singleCopyRegions) {
                        writer.write(scr, true);
                    }
                }
                Logger.info(CLASS_NAME, "Dumped " + numSingleCopy + " single copy regions for chromsome " + chrom.getName());
            }
        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Error writing single copy regions: " + e.getMessage());
        }
    }
}
