package edu.tamu.pg.panhap.Tasks;

import edu.tamu.pg.panhap.Data.GenomeMap;
import edu.tamu.pg.panhap.Utils.*;

import java.io.IOException;
import java.nio.file.Path;

public class RunKMC {

    private static final String CLASS_NAME = RunKMC.class.getSimpleName();

    private final DirectoryStructure dirs;
    private final int threads;
    private final int kmerSize;
    private final int minCount;
    private final double maxMemory;

    public RunKMC(DirectoryStructure dirs) {
        this(dirs, 1, Defaults.K_SIZE, Defaults.MIN_COUNT, Defaults.HEAPSIZE);
    }

    public RunKMC(DirectoryStructure dirs, int threads) {
        this(dirs, threads, Defaults.K_SIZE, Defaults.MIN_COUNT, Defaults.HEAPSIZE);
    }

    public RunKMC(DirectoryStructure dirs, int threads, int kmerSize, int minCount) {
        this(dirs, threads, kmerSize, minCount, Defaults.HEAPSIZE);
    }

    public RunKMC(DirectoryStructure dirs, int threads, int kmerSize, int minCount, double maxMemory) {
        this.dirs = dirs;
        this.threads = threads;
        this.kmerSize = kmerSize;
        this.minCount = minCount;
        this.maxMemory = maxMemory;
    }

    /**
     * Runs KMC for each genome in the GenomeMap using default/constructor configurations.
     */
    public void run(GenomeMap genomeMap) throws IOException {
        run(genomeMap, this.kmerSize, this.minCount, this.threads, this.maxMemory);
    }

    /**
     * Runs KMC for each genome in the GenomeMap with explicit execution parameters.
     */
    public void run(GenomeMap genomeMap, int kmerSize, int minCount, int threads) throws IOException {
        run(genomeMap, kmerSize, minCount, threads, this.maxMemory);
    }

    /**
     * Core execution method iterating through all genomes in the map.
     */
    public void run(GenomeMap genomeMap, int kmerSize, int minCount, int threads, double maxMemory) throws IOException {
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
                    .maxMemory(maxMemory)
                    .build();

            Logger.info(CLASS_NAME, "Running KMC for genome: " + id);
            KMCExecutor.run(kmcConfig);
        }
    }
}