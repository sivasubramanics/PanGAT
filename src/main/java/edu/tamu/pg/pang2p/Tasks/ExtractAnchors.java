package edu.tamu.pg.pang2p.Tasks;

import edu.tamu.pg.pang2p.Data.*;
import edu.tamu.pg.pang2p.Utils.Logger;
import edu.tamu.pg.pang2p.Utils.PipelineSteps;

import java.nio.file.Files;
import java.nio.file.Path;

public class ExtractAnchors {

    private static final String CLASS_NAME = ExtractAnchors.class.getSimpleName();
    private final DirectoryStructure dirs;
    private final Path singleCopyDir;
    private final Path genomesDir;
    private final boolean isForce;

    public ExtractAnchors(DirectoryStructure dirs, boolean isForce) {
        this.dirs = dirs;
        this.genomesDir = dirs.getDirectory(PipelineSteps.GENOMES);
        this.singleCopyDir = dirs.getDirectory(PipelineSteps.ANCHORS);
        this.isForce = isForce;
    }

    public void run(GenomeMap genomesMap, int kmerSize, int minCount, int minLength, int maxGap, int nThreads, boolean strict)
            throws Exception {

        // Check if all output files exist for checkpointing
        if (!isForce && checkOutputsExist(genomesMap)) {
            Logger.info(CLASS_NAME, "Skipping ExtractSCR & KMC step: output files already exist. Use --force to re-run.");
            return;
        }

        for (String genomeID : genomesMap.getOrgGenomes()) {
            Path genomeFile = genomesDir.resolve(genomeID + ".fasta");
            Path outputFile = singleCopyDir.resolve(genomeID + ".fasta");

            Logger.info(CLASS_NAME, "Extracting single copy anchors for genome: " + genomeID);

            Anchors anchors = new Anchors(genomeFile, outputFile, kmerSize, minCount, minLength, maxGap, nThreads, strict);
            int exitCode = anchors.call();

            if (exitCode != 0) {
                Logger.error(CLASS_NAME, "Extraction failed for genome: " + genomeID + " with exit code: " + exitCode);
            }

            // Verify file creation and non-zero size
            if (!Files.exists(outputFile) || Files.size(outputFile) == 0) {
                Logger.error(CLASS_NAME, "Extraction failed or produced empty file for: " + genomeID);
            }
        }
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