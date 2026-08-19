package edu.tamu.pg.pang2p.Tasks;

import edu.tamu.pg.pang2p.Utils.FileUtils;
import edu.tamu.pg.pang2p.Utils.Logger;
import edu.tamu.pg.pang2p.Utils.PipelineSteps;

import java.nio.file.Files;
import java.nio.file.Path;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class DirectoryStructure {

    private static final String CLASS_NAME = DirectoryStructure.class.getSimpleName();
    private final Path baseDir;

    public DirectoryStructure(Path baseDir) {
        this.baseDir = baseDir;
    }

    /**
     * Creates directories for all pipeline steps.
     *
     * @return List of failed directory paths. If empty, all succeeded.
     */
    public List<Path> createDirectories() {
        List<Path> failedDirs = new ArrayList<>();
        for (PipelineSteps step : PipelineSteps.values()) {
            Path dirPath = baseDir.resolve(step.dirName());
            if (!Files.exists(dirPath)) {
                try {
                    FileUtils.mkdir(dirPath);
                    Logger.info(CLASS_NAME, "Created directory: " + dirPath);
                } catch (IOException e) {
                    Logger.error(CLASS_NAME, "Failed to create directory: " + dirPath + " (" + e.getMessage() + ")");
                    failedDirs.add(dirPath);
                }
            } else {
                Logger.info(CLASS_NAME, "Directory already exists: " + dirPath);
            }
        }
        return failedDirs;
    }

    public Path getBaseDir() {
        return baseDir;
    }

    public Path getDirectory(PipelineSteps step) {
        return baseDir.resolve(step.dirName());
    }

    public Path getSingleCopyReference() {
        return baseDir.resolve(PipelineSteps.SINGLE_COPY.dirName()).resolve("pan_clusters.fasta");
    }

    public Path getSingleCopyRawFile() {
        return baseDir.resolve(PipelineSteps.SINGLE_COPY.dirName()).resolve("all_single_copy_regions.fasta");
    }

    public Path getGenotypeDir() {
        return baseDir.resolve(PipelineSteps.GENOTYPING.dirName());
    }

    public Path getPanClustersRepFile() {
        return baseDir.resolve(PipelineSteps.SINGLE_COPY.dirName()).resolve("pan_clusters_rep.tsv");
    }

    public String getCohortBase() {
        return baseDir.resolve(PipelineSteps.GENOTYPING.dirName()).resolve("cohort").toString();
    }
}