package edu.tamu.pg.pang2p.Tasks;

import edu.tamu.pg.pang2p.Data.GenomeMap;
import edu.tamu.pg.pang2p.Utils.*;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

public class RunKMC {

    private static final String CLASS_NAME = RunKMC.class.getSimpleName();

    private final DirectoryStructure dirs;
    private final int threads;
    private final int kmerSize;
    private final int minCount;

    public RunKMC(DirectoryStructure dirs) {
        this(dirs, 1, Defaults.K_SIZE, Defaults.MIN_COUNT);
    }

    public RunKMC(DirectoryStructure dirs, int threads) {
        this(dirs, threads, Defaults.K_SIZE, Defaults.MIN_COUNT);
    }

    public RunKMC(DirectoryStructure dirs, int threads, int kmerSize, int minCount) {
        this.dirs = dirs;
        this.threads = threads;
        this.kmerSize = kmerSize;
        this.minCount = minCount;
    }

    /**
     * Runs KMC for each genome in the GenomeMap using default/constructor configurations.
     */
    public void run(GenomeMap genomeMap) throws IOException {
        run(genomeMap, this.kmerSize, this.minCount, this.threads);
    }

    /**
     * Core execution method iterating through all genomes in the map.
     * KMCConfig dynamically handles safe memory allocation inside build().
     */
    public void run(GenomeMap genomeMap, int kmerSize, int minCount, int threads) throws IOException {
        if (dirs == null) {
            throw new IllegalStateException("DirectoryStructure must be provided to run KMC with GenomeMap.");
        }

        Path kmcDir = dirs.getDirectory(PipelineSteps.KMC);
        Path genomesDir = dirs.getDirectory(PipelineSteps.GENOMES);

        for (String id : genomeMap.getOrgGenomes()) {
            Path fastaPath = genomesDir.resolve(id + ".fasta");
            Path outputPrefix = kmcDir.resolve(id);

            KMCConfig kmcConfig = KMCConfig.builder()
                    .inputPath(fastaPath)
                    .outputPrefix(outputPrefix)
                    .kmerSize(kmerSize)
                    .minCount(minCount)
                    .threads(threads)
                    .build();

            Logger.info(CLASS_NAME, "Running KMC for genome: " + id);
            KMCExecutor.run(kmcConfig);
        }
    }

    public void clean() {
        Path kmcDir = dirs.getDirectory(PipelineSteps.KMC);
        if (kmcDir != null && Files.exists(kmcDir)) {
            try {
                FileUtils.rmdir(kmcDir);
                Logger.info(CLASS_NAME, "Cleaned up KMC directory: " + kmcDir);
            } catch (IOException e) {
                Logger.error(CLASS_NAME, "Failed to clean KMC directory: " + e.getMessage());
            }
        }
    }
}