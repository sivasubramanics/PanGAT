package edu.tamu.pg.pang2p.Tasks;

import edu.tamu.pg.pang2p.Data.*;
import edu.tamu.pg.pang2p.Utils.*;
import java.io.IOException;
import java.nio.file.DirectoryStream;
import java.nio.file.Path;
import java.nio.file.Files;

public class RunMMseqs {
    private static final String CLASS_NAME = RunMMseqs.class.getSimpleName();
    private final DirectoryStructure dirs;
    private final MMseqsConfig config;

    public RunMMseqs(DirectoryStructure dirs, int threads, boolean isForce) {
        this(dirs, MMseqsConfig.builder().threads(threads).isForce(isForce).build());
    }

    public RunMMseqs(DirectoryStructure dirs, MMseqsConfig config) {
        this.dirs = dirs;
        this.config = config;
    }

    // Example: genomesMap = Map<genomeID, genomeName>
    public void run(GenomeMap genomesMap) throws IOException {
        Path singleCopyDir = dirs.getDirectory(PipelineSteps.SINGLE_COPY);
        Path clusterFile = singleCopyDir.resolve("clusters.tsv");

        // Skip execution if output exists and isForce is false
        if (!config.isForce() && Files.exists(clusterFile)) {
            Logger.info(CLASS_NAME, "Skipping RunMMseqs step: output file already exists ("
                    + clusterFile + "). Use --force to re-run.");
            return;
        }

        // delete all intermediate files from last run.
        cleanFiles(true);

        Logger.info(CLASS_NAME, "Clustering single copy regions using MMseqs2.");
        Path scrFasta = singleCopyDir.resolve("all_single_copy_regions.fasta");

        // Concatenate all single copy regions into one fasta file
        try (FastaWriter mergedFasta = new FastaWriter(scrFasta.toString())) {
            for (String id : genomesMap.getOrgGenomes()) {
                Path scrPath = singleCopyDir.resolve(id + ".fasta");
                if (Files.exists(scrPath)) {
                    try (FastaReader reader = new FastaReader(scrPath.toString())) {
                        for (Fasta fasta : reader) {
                            mergedFasta.write(fasta);
                        }
                    } catch (IOException e) {
                        Logger.error(CLASS_NAME, "Error reading single copy region file: " + e.getMessage());
                    }
                } else {
                    Logger.warning(CLASS_NAME, "Single copy region file not found for genome: " + id);
                }
            }
        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Error concatenating single copy regions: " + e.getMessage());
            return;
        }

        FastaIndex fastaIndex = new FastaIndex(scrFasta.toString());
        fastaIndex.close();

        // CreateDB with mmseqs
        Path dbFile = singleCopyDir.resolve("all_single_copy_regions.db");
        String createDBLog = singleCopyDir.resolve("mmseqs_createDB.log").toString();
        String createDBCmd = buildCreateDBCommand(scrFasta, dbFile);

        Logger.info(CLASS_NAME, "Creating DB for mmseqs clustering. Log file: " + createDBLog);
        Logger.info(CLASS_NAME, "MMseqs2 CMD: " + createDBCmd);
        HelperFunctions.tryExec(createDBCmd, createDBLog, createDBLog);

        // Cluster with mmseqs
        Path tmpDir = singleCopyDir.resolve("tmp");
        if (FileUtils.isDir(tmpDir.toString())) {
            FileUtils.rmdir(tmpDir);
        } else if (FileUtils.isFile(tmpDir.toString())) {
            FileUtils.deleteRecursively(new java.io.File(tmpDir.toString()));
        }
        FileUtils.mkdir(tmpDir);

        Path scrClusters = singleCopyDir.resolve("clusters");
        String mmseqsClusterLog = singleCopyDir.resolve("mmseqs_clustering.log").toString();
        String clusterCmd = buildClusterCommand(dbFile, scrClusters, tmpDir);

        Logger.info(CLASS_NAME, "Clustering single copy regions with MMseqs2. Log file: " + mmseqsClusterLog);
        Logger.info(CLASS_NAME, "MMseqs2 CMD: " + clusterCmd);
        HelperFunctions.tryExec(clusterCmd, mmseqsClusterLog, mmseqsClusterLog);

        // Convert to TSV
        String createTSVLog = singleCopyDir.resolve("mmseqs_clustering.log").toString();
        String createTSVCmd = buildCreateTSVCommand(dbFile, scrClusters, clusterFile);

        Logger.info(CLASS_NAME, "Converting clustering results to TSV. Log file: " + createTSVLog);
        Logger.info(CLASS_NAME, "MMseqs2 CMD: " + createTSVCmd);
        HelperFunctions.tryExec(createTSVCmd, createTSVLog, createTSVLog);

        cleanFiles(false);
        FileUtils.rmdir(tmpDir);
        Logger.info(CLASS_NAME, "Clustering of single copy regions completed.");
    }

    private void cleanFiles(boolean cleanAll) throws IOException {
        Path singleCopyDir = dirs.getDirectory(PipelineSteps.SINGLE_COPY);
        String[] filesToDelete = {
                "clusters.dbtype",
                "clusters.index",
                "all_single_copy_regions.db_h.dbtype",
                "all_single_copy_regions.db.dbtype",
                "all_single_copy_regions.db.lookup",
                "all_single_copy_regions.db_h.index",
                "all_single_copy_regions.db.index",
                "all_single_copy_regions.db",
                "all_single_copy_regions.db.source",
                "all_single_copy_regions.db_h",
        };

        for (String filename : filesToDelete) {
            Path file = singleCopyDir.resolve(filename);
            if (Files.exists(file)) {
                FileUtils.rmFile(file);
            }
        }

        // remove all files with clusters.%d
        try (DirectoryStream<Path> stream = Files.newDirectoryStream(singleCopyDir)) {
            for (Path file : stream) {
                String name = file.getFileName().toString();
                if (name.matches("clusters\\.\\d+")) {
                    FileUtils.rmFile(file);
                }
            }
        }

        Path tmpDir = singleCopyDir.resolve("tmp");
        FileUtils.rmdir(tmpDir);

        if (cleanAll) {
            FileUtils.rmFile(singleCopyDir.resolve("all_single_copy_regions.fasta"));
            FileUtils.rmFile(singleCopyDir.resolve("mmseqs_clustering.log"));
            FileUtils.rmFile(singleCopyDir.resolve("clusters.tsv"));
            FileUtils.rmFile(singleCopyDir.resolve("mmseqs_createDB.log"));
        }
    }

    private String buildCreateDBCommand(Path fasta, Path dbFile) {
        return String.join(" ",
                "mmseqs", "createdb",
                fasta.toString(),
                dbFile.toString(),
                "--threads", String.valueOf(config.getThreads()),
                "--dbtype", String.valueOf(config.getDbType())
        );
    }

    private String buildClusterCommand(Path dbFile, Path clusters, Path tmpDir) {
        return String.join(" ",
                "mmseqs", "cluster",
                dbFile.toString(), clusters.toString(), tmpDir.toString(),
                "--threads", String.valueOf(config.getThreads()),
                "--min-seq-id", String.valueOf(config.getMinSeqIdentity()),
                "-c", String.valueOf(config.getMinCoverage()),
                "--cov-mode", String.valueOf(config.getCovMode())
        );
    }

    private String buildCreateTSVCommand(Path dbFile, Path clusters, Path clusterFile) {
        return String.join(" ",
                "mmseqs", "createtsv",
                dbFile.toString(), dbFile.toString(), clusters.toString(), clusterFile.toString(),
                "--threads", String.valueOf(config.getThreads())
        );
    }

    public static class MMseqsConfig {
        private final double minSeqIdentity;
        private final double minCoverage;
        private final int covMode;
        private final int threads;
        private final int dbType;
        private final boolean isForce;

        private MMseqsConfig(Builder builder) {
            this.minSeqIdentity = builder.minSeqIdentity;
            this.minCoverage = builder.minCoverage;
            this.covMode = builder.covMode;
            this.threads = builder.threads;
            this.dbType = builder.dbType;
            this.isForce = builder.isForce;
        }

        public double getMinSeqIdentity() { return minSeqIdentity; }
        public double getMinCoverage() { return minCoverage; }
        public int getCovMode() { return covMode; }
        public int getThreads() { return threads; }
        public int getDbType() { return dbType; }
        public boolean isForce() { return isForce; }

        public static Builder builder() { return new Builder(); }

        public static class Builder {
            private double minSeqIdentity = Defaults.MIN_SEQ_IDENTITY;
            private double minCoverage = Defaults.MIN_COVERAGE;
            private int covMode = Defaults.COV_MODE;
            private int threads = 1;
            private int dbType = 2;
            private boolean isForce = false;

            public Builder minSeqIdentity(double value) {
                this.minSeqIdentity = value;
                return this;
            }
            public Builder minCoverage(double value) {
                this.minCoverage = value;
                return this;
            }
            public Builder covMode(int value) {
                this.covMode = value;
                return this;
            }
            public Builder threads(int value) {
                this.threads = value;
                return this;
            }
            public Builder dbType(int value) {
                this.dbType = value;
                return this;
            }
            public Builder isForce(boolean value) {
                this.isForce = value;
                return this;
            }
            public MMseqsConfig build() {
                return new MMseqsConfig(this);
            }
        }
    }
}