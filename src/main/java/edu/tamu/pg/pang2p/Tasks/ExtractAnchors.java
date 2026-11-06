package edu.tamu.pg.pang2p.Tasks;

import edu.tamu.pg.pang2p.Data.*;
import edu.tamu.pg.pang2p.Utils.Defaults;
import edu.tamu.pg.pang2p.Utils.Logger;
import edu.tamu.pg.pang2p.Utils.PipelineSteps;

import java.nio.file.Files;
import java.nio.file.Path;

public class ExtractAnchors {

    private static final String CLASS_NAME = ExtractAnchors.class.getSimpleName();
    private final DirectoryStructure dirs;
    private final Path singleCopyDir;
//    private final Path kmcDir;
    private final Path genomesDir;
    private final boolean isForce;

    private static final int dumpSize = Defaults.COUNTER_DUMP_SIZE;
    private final int minCount = Defaults.MIN_COUNT;
    private final int minLength = Defaults.MIN_LENGTH;

    public ExtractAnchors(DirectoryStructure dirs, boolean isForce) {
        this.dirs = dirs;
//        this.kmcDir = dirs.getDirectory(PipelineSteps.KMC);
        this.genomesDir = dirs.getDirectory(PipelineSteps.GENOMES);
        this.singleCopyDir = dirs.getDirectory(PipelineSteps.ANCHORS);
        this.isForce = isForce;
    }

    public void run(GenomeMap genomesMap, int kmerSize, int minCount, int minLength, int maxGap, int nThreads)
            throws Exception {

        // Check if all output files exist for checkpointing
        if (!isForce && checkOutputsExist(genomesMap)) {
            Logger.info(CLASS_NAME, "Skipping ExtractSCR & KMC step: output files already exist. Use --force to re-run.");
            return;
        }

        // 1. Run KMC
//        RunKMC runKMC = new RunKMC(dirs);
//        runKMC.run(genomesMap, kmerSize, minCount, nThreads);

        // 2. Extract Single Copy Regions
        for (String genomeID : genomesMap.getOrgGenomes()) {
            Path genomeFile = genomesDir.resolve(genomeID + ".fasta");
//            String kmcPrefix = kmcDir.resolve(genomeID).toString();
            Path outputFile = singleCopyDir.resolve(genomeID + ".fasta");

            Logger.info(CLASS_NAME, "Extracting single copy anchors for genome: " + genomeID);

            AnchorsTest anchors = new AnchorsTest(genomeFile, outputFile, kmerSize, minCount, minLength, maxGap, nThreads);
//            Anchors anchors = new Anchors(genomeFile, outputFile, kmerSize, minCount, minLength, maxGap, nThreads);
            int exitCode = anchors.call();

            if (exitCode != 0) {
                Logger.error(CLASS_NAME, "Extraction failed for genome: " + genomeID + " with exit code: " + exitCode);
            }

            // Verify file creation and non-zero size
            if (!Files.exists(outputFile) || Files.size(outputFile) == 0) {
                Logger.error(CLASS_NAME, "Extraction failed or produced empty file for: " + genomeID);
            }
        }

        // 3. Clean up KMC database temporary files safely after ALL extractions finish
        // runKMC.clean();
    }

    private boolean checkOutputsExist(GenomeMap genomesMap) {
        for (String genomeID : genomesMap.getOrgGenomes()) {
            Path outputFile = singleCopyDir.resolve(genomeID + ".fasta");
            if (!Files.exists(outputFile)) {
                return false;
            }
        }
        return true;
    }
}