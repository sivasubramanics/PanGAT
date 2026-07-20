package edu.tamu.pg.panhap.Tasks;

import edu.tamu.pg.panhap.Data.GenomeMap;
import edu.tamu.pg.panhap.Utils.*;

import java.io.IOException;
import java.nio.file.Path;

public class RunKMC {

    private static final String CLASS_NAME = RunKMC.class.getSimpleName();

    private final DirectoryStructure dirs;
    private final KMCConfig config;

    public RunKMC(DirectoryStructure dirs, int threads) {
        this(dirs, KMCConfig.builder()
                .threads(threads)
                .build());
    }

    public RunKMC(DirectoryStructure dirs, KMCConfig config) {
        this.dirs = dirs;
        this.config = config;
    }

    public RunKMC(KMCConfig config) {
        this.dirs = null;
        this.config = config;
    }

    public void run(Path fastaPath, Path kmcPrefix) throws IOException {

            String logFile = kmcPrefix + ".kmc.log";

            // kmcPrefix baseDir / tmp
            Path parent = kmcPrefix.getParent();
            Path tmpDir = (parent == null) ? Path.of("tmp") : parent.resolve("tmp");
            Logger.info(CLASS_NAME, "Creating temporary directory for KMC: " + tmpDir);
            FileUtils.mkdir(tmpDir);
        try {
            String command = buildCommand(
                    fastaPath.toString(),
                    kmcPrefix.toString(),
                    tmpDir
            );

            Logger.info(CLASS_NAME, "Running KMC for fasta: " + fastaPath);
            Logger.info(CLASS_NAME, command);

            HelperFunctions.tryExec(command, logFile, logFile);
        }
        finally {
            FileUtils.rmdir(tmpDir);
        }
    }

    public void run(GenomeMap genomeMap) throws IOException {

        // if the maxMemory is less than 2.0GB, return error and exit
        if (config.getMaxMemory() < 2.0) {
            Logger.error(CLASS_NAME, "Max memory must be at least 2.0GB for KMC to run.");
        }

        Path kmcDir = dirs.getDirectory(PipelineSteps.KMC);
        Path genomesDir = dirs.getDirectory(PipelineSteps.GENOMES);

        Path tmpDir = kmcDir.resolve("tmp");

        FileUtils.mkdir(tmpDir);

        try {

            for (String id : genomeMap.getOrgGenomes()) {

                String logFile = kmcDir.resolve(id + ".kmc.log").toString();

                String command = buildCommand(
                        genomesDir.resolve(id + ".fasta").toString(),
                        kmcDir.resolve(id).toString(),
                        tmpDir
                );

                Logger.info(CLASS_NAME, "Running KMC for genome: " + id);
                Logger.info(CLASS_NAME, command);

                HelperFunctions.tryExec(command, logFile, logFile);
            }

        } finally {
            FileUtils.rmdir(tmpDir);
        }
    }

    private String buildCommand(String input,
                                String outputPrefix,
                                Path tmpDir) {

        return String.join(" ",
                "kmc",
                "-k" + config.getKmerSize(),
                "-p" + config.getSignatureLength(),
                "-ci" + config.getMinCount(),
                "-t" + config.getThreads(),
                "-m" + config.getMaxMemory(),
                "-fm",
                input,
                outputPrefix,
                tmpDir.toString()
        );
    }

    public static class KMCConfig {

        private final int kmerSize;
        private final int signatureLength;
        private final int minCount;
        private final double maxMemory;
        private final int threads;

        private KMCConfig(Builder builder) {
            this.kmerSize = builder.kmerSize;
            this.signatureLength = builder.signatureLength;
            this.minCount = builder.minCount;
            this.maxMemory = builder.maxMemory;
            this.threads = builder.threads;
        }

        public int getKmerSize() {
            return kmerSize;
        }

        public int getSignatureLength() {
            return signatureLength;
        }

        public int getMinCount() {
            return minCount;
        }

        public double getMaxMemory() {
            return maxMemory;
        }

        public int getThreads() {
            return threads;
        }

        public static Builder builder() {
            return new Builder();
        }

        public static class Builder {

            private int kmerSize = Defaults.K_SIZE;
            private int signatureLength = Defaults.SIGNATURE_LENGTH;
            private int minCount = Defaults.MIN_COUNT;
            private double maxMemory = Defaults.HEAPSIZE;
            private int threads = 1;

            public Builder kmerSize(int value) {
                this.kmerSize = value;
                return this;
            }

            public Builder signatureLength(int value) {
                this.signatureLength = value;
                return this;
            }

            public Builder minCount(int value) {
                this.minCount = value;
                return this;
            }

            public Builder maxMemory(double value) {
                this.maxMemory = value;
                return this;
            }

            public Builder threads(int value) {
                this.threads = value;
                return this;
            }

            public KMCConfig build() {
                return new KMCConfig(this);
            }
        }
    }
}