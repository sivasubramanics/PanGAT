package edu.tamu.pg.panhap.Tasks;

import edu.tamu.pg.panhap.Data.*;
import edu.tamu.pg.panhap.Utils.Defaults;
import edu.tamu.pg.panhap.Utils.FileUtils;
import edu.tamu.pg.panhap.Utils.Logger;
import edu.tamu.pg.panhap.Utils.PipelineSteps;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class ExtractSCR {

    private static final String CLASS_NAME = ExtractSCR.class.getSimpleName();
    private final Path singleCopyDir;
    private final Path kmcDir;
    private final Path genomesDir;
    private static final int dumpSize = Defaults.COUNTER_DUMP_SIZE;
    private final int minCount = Defaults.MIN_COUNT;
    private final int minLength = Defaults.MIN_LENGTH;

    public ExtractSCR(DirectoryStructure dirs) {
        this.kmcDir = dirs.getDirectory(PipelineSteps.KMC);
        this.genomesDir = dirs.getDirectory(PipelineSteps.GENOMES);
        this.singleCopyDir = dirs.getDirectory(PipelineSteps.SINGLE_COPY);
    }

    public void run(GenomeMap genomesMap) throws IOException {
        try {
            run(genomesMap, minCount, minLength, 2);
        } catch (ExecutionException e) {
            throw new RuntimeException(e);
        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        }
    }

    public void run(GenomeMap genomesMap, int minCount, int minLength, int nThreads) throws IOException, ExecutionException, InterruptedException {

        ExecutorService executor = Executors.newFixedThreadPool(nThreads);

        for (String genomeID : genomesMap.getOrgGenomes()) {

            KMC kmc = new KMC(kmcDir.resolve(genomeID).toString(), true);
            FastaReader genome = new FastaReader(genomesDir.resolve(genomeID + ".fasta").toString());

            List<Future<?>> futures = new ArrayList<>();

            for (Fasta chrom : genome) {
                futures.add(executor.submit(() -> {
                    processGenome(chrom, kmc, minCount, minLength);
                    return null;
                }));
            }

            for (Future<?> future : futures) {
                future.get();
            }

//            kmc.dumpKmerTable(kmcDir.resolve(genomeID + "_kmer_table.txt").toString());
            kmc.close();

            try (FastaWriter fastaWriter = new FastaWriter(singleCopyDir.resolve(genomeID + ".fasta").toString())) {
                for (String chromID : genomesMap.getChroms(genomeID)) {
                    Path chromFile = singleCopyDir.resolve(chromID + ".fasta");
                    if (FileUtils.isFile(chromFile)) {
                        try (FastaReader chromReader = new FastaReader(chromFile.toString())) {
                            for (Fasta scr : chromReader) {
                                fastaWriter.write(scr);
                            }
                        } catch (IOException e) {
                            Logger.error(CLASS_NAME, "Error reading single copy regions for chromosome " + chromID + ": " + e.getMessage());
                        }
                        FileUtils.rmFile(chromFile);
                    } else {
                        Logger.info(CLASS_NAME, "No single copy regions found for chromosome " + chromID);
                    }
                }
            } catch (IOException e) {
                Logger.error(CLASS_NAME, "Error writing single copy regions for genome: " + genomeID);
            }
        }
    }

    private void processGenome(Fasta chrom, KMC kmc, int minCount, int minLength) {
        // Larger buffer: 256KB instead of the default 8KB
        try (FastaWriter writer = new FastaWriter(singleCopyDir.resolve(chrom.getName() + ".fasta").toString(), false)) {
            int counter = 0;
            int kSize = kmc.getKmerLength();
            int prefixLength = kmc.getPrefixLength();
            boolean isCanonical = kmc.isBothStrands();
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

                    if (!kmc.isExist(kmer)){
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
                        singleCopyRegions.add(new Fasta(numSingleCopy, chrom.getName() + "" + numSingleCopy, singleCopySeq.toString(), regionStart + "-" + currentEnd));
                        numSingleCopy++;
                    }
                    singleCopySeq.setLength(0);
                }

                if (!singleCopyRegions.isEmpty()) {
                    for (Fasta scr : singleCopyRegions) {
                        writer.write(scr, false);
                    }
                }
                Logger.info(CLASS_NAME, "Dumped " + numSingleCopy + " single copy regions for chromsome " + chrom.getName());
            writer.close();
        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Error writing single copy regions for chromosome " + chrom.getName() + ": " + e.getMessage());
        }
    }
}