package edu.tamu.pg.pang2p.Utils;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class KMCConfig {

    public enum InputFormat {
        FASTA("-fm"),
        FASTQ("-fq"),
        BAM("-fbam"),
        CRAM("-fcram"),
        KMC("KMC");

        private final String flag;

        InputFormat(String flag) {
            this.flag = flag;
        }

        public String getFlag() {
            return flag;
        }

        /**
         * Infers input format based on file extensions (including .gz).
         */
        public static InputFormat fromPath(Path path) {
            String fileName = path.getFileName().toString().toLowerCase();

            // Strip .gz suffix if present to check underlying format
            if (fileName.endsWith(".gz")) {
                fileName = fileName.substring(0, fileName.length() - 3);
            }

            if (fileName.endsWith(".fasta") || fileName.endsWith(".fa") || fileName.endsWith(".fna")) {
                return FASTA;
            } else if (fileName.endsWith(".fastq") || fileName.endsWith(".fq")) {
                return FASTQ;
            } else if (fileName.endsWith(".bam")) {
                return BAM;
            } else if (fileName.endsWith(".cram")) {
                return CRAM;
            } else if (fileName.endsWith(".kmc_suf")) {
                return KMC;
            }

            throw new IllegalArgumentException("Unsupported file extension for path: " + path
                    + ". Supported extensions are: .fasta, .fa, .fna, .fastq, .fq, .bam, .cram, .kmc_suf (and their .gz versions)");
        }
    }

    private final Path[] inputPaths;
    private final Path outputPrefix;
    private final InputFormat inputFormat;
    private final int kmerSize;
    private final int signatureLength;
    private final int minCount;
    private final long maxCount;
    private final double maxMemory;
    private final int threads;

    private KMCConfig(Builder builder) {
        this.inputPaths = builder.inputPaths.toArray(new Path[0]);
        this.outputPrefix = builder.outputPrefix;
        this.inputFormat = builder.inputFormat;
        this.kmerSize = builder.kmerSize;
        this.signatureLength = builder.signatureLength;
        this.minCount = builder.minCount;
        this.maxCount = builder.maxCount;
        this.maxMemory = builder.maxMemory;
        this.threads = builder.threads;
    }

    public Path[] getInputPaths() { return inputPaths.clone(); }
    public Path getOutputPrefix() { return outputPrefix; }
    public InputFormat getInputFormat() { return inputFormat; }
    public int getKmerSize() { return kmerSize; }
    public int getSignatureLength() { return signatureLength; }
    public int getMinCount() { return minCount; }
    public long getMaxCount() { return maxCount; }
    public double getMaxMemory() { return maxMemory; }
    public int getThreads() { return threads; }

    public static Builder builder() {
        return new Builder();
    }

    public static class Builder {
        private static final double AUTO_MEMORY = -1.0;

        private final List<Path> inputPaths = new ArrayList<>();
        private Path outputPrefix;
        private InputFormat inputFormat; // Auto-detected if null
        private int kmerSize = Defaults.K_SIZE;
        private int signatureLength = Defaults.SIGNATURE_LENGTH;
        private int minCount = Defaults.MIN_COUNT;
        private long maxCount = 1_000_000_000L;
        private double maxMemory = AUTO_MEMORY; // Set to auto-calculate by default
        private int threads = 1;

        public Builder inputPath(Path inputPath) {
            if (inputPath != null) {
                this.inputPaths.add(inputPath);
            }
            return this;
        }

        public Builder inputPaths(Path... inputPaths) {
            if (inputPaths != null) {
                this.inputPaths.addAll(Arrays.asList(inputPaths));
            }
            return this;
        }

        public Builder inputPaths(List<Path> inputPaths) {
            if (inputPaths != null) {
                this.inputPaths.addAll(inputPaths);
            }
            return this;
        }

        public Builder outputPrefix(Path outputPrefix) {
            this.outputPrefix = outputPrefix;
            return this;
        }

        public Builder inputFormat(InputFormat inputFormat) {
            this.inputFormat = inputFormat;
            return this;
        }

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

        public Builder maxCount(long value) {
            this.maxCount = value;
            return this;
        }

        /**
         * Sets max memory in GB. Pass a negative value or 0 to let KMCConfig compute a safe dynamic default.
         */
        public Builder maxMemory(double value) {
            this.maxMemory = value;
            return this;
        }

        public Builder threads(int value) {
            this.threads = value;
            return this;
        }

        public KMCConfig build() {
            if (inputPaths.isEmpty()) {
                throw new IllegalArgumentException("At least one input path must be specified.");
            }
            if (outputPrefix == null) {
                throw new IllegalArgumentException("User-defined output prefix must be specified.");
            }

            // Infer input format from the first file if not explicitly set
            if (this.inputFormat == null) {
                this.inputFormat = InputFormat.fromPath(inputPaths.get(0));
            }

            // Ensure all input files match the resolved input format
            for (Path path : inputPaths) {
                InputFormat inferred = InputFormat.fromPath(path);
                if (inferred != this.inputFormat) {
                    throw new IllegalArgumentException("Mixed file formats detected! File '" + path
                            + "' resolves to " + inferred + " but expected " + this.inputFormat);
                }
            }

            // Dynamically calculate memory if default/auto-mode is active
            if (this.maxMemory <= 0) {
                this.maxMemory = calculateSafeMemoryGB();
            }

            return new KMCConfig(this);
        }

        /**
         * Safely estimates memory (in GB) for KMC without over-allocating system/JVM heap resources.
         */
        private double calculateSafeMemoryGB() {
            Runtime runtime = Runtime.getRuntime();

            // Total available headroom before JVM max heap limit
            long maxHeap = runtime.maxMemory();
            long allocatedHeap = runtime.totalMemory();
            long freeHeap = runtime.freeMemory();
            long unallocatedHeap = maxHeap - allocatedHeap;
            long availableHeadroomBytes = freeHeap + unallocatedHeap;

            // Target 30% of remaining available JVM headroom for external KMC execution
            double calculatedGB = (availableHeadroomBytes * 0.30) / (1024.0 * 1024.0 * 1024.0);

            // KMC binary strictly requires a minimum of 2.0 GB (-m2)
            double minMemoryGB = 2.0;
            double maxMemoryGB = 8.0;

            double finalMemory = Math.max(minMemoryGB, Math.min(calculatedGB, maxMemoryGB));

            // Round to 1 decimal place (e.g., 2.5 GB)
            return Math.round(finalMemory * 10.0) / 10.0;
        }
    }
}