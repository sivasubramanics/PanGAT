package edu.tamu.pg.panhap.Tasks;

import edu.tamu.pg.panhap.Data.*;
import edu.tamu.pg.panhap.Utils.*;
import java.io.IOException;
import java.nio.file.DirectoryStream;
import java.nio.file.Path;
import java.nio.file.Files;
import java.util.Map;

public class RunMMseqs {
    private static final String CLASS_NAME = RunMMseqs.class.getSimpleName();
    private final DirectoryStructure dirs;
    private final MMseqsConfig config;

    public RunMMseqs(DirectoryStructure dirs, int threads) {
        this(dirs, MMseqsConfig.builder().threads(threads).build());
    }

    public RunMMseqs(DirectoryStructure dirs, MMseqsConfig config) {
        this.dirs = dirs;
        this.config = config;
    }

    // Example: genomesMap = Map<genomeID, genomeName>
    public void run(GenomeMap genomesMap) throws IOException {
        Path singleCopyDir = dirs.getDirectory(PipelineSteps.SINGLE_COPY);
        // delete all intermediate files from last run.
        cleanFiles(true);

        Logger.info(CLASS_NAME, "Clustering single copy regions using MMseqs2.");
        Path scrFasta = singleCopyDir.resolve("all_single_copy_regions.fasta");

        // Concatenate all single copy regions into one fasta file
        try (FastaWriter mergedFasta = new FastaWriter(scrFasta.toString())) {
            for (String id : genomesMap.getGenomes()) {
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
            new FastaIndex(scrFasta.toString());
        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Error concatenating single copy regions: " + e.getMessage());
            return;
        }

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
        Path clusterFile = singleCopyDir.resolve("clusters.tsv");
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
                "clusters.3",
                "clusters.2",
                "clusters.1",
                "clusters.0",
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

    // MMseqsConfig builder is similar to KMCConfig, with appropriate defaults!
    public static class MMseqsConfig {
        private final double minSeqIdentity;
        private final double minCoverage;
        private final int covMode;
        private final int threads;
        private final int dbType;

        private MMseqsConfig(Builder builder) {
            this.minSeqIdentity = builder.minSeqIdentity;
            this.minCoverage = builder.minCoverage;
            this.covMode = builder.covMode;
            this.threads = builder.threads;
            this.dbType = builder.dbType;
        }

        public double getMinSeqIdentity() { return minSeqIdentity; }
        public double getMinCoverage() { return minCoverage; }
        public int getCovMode() { return covMode; }
        public int getThreads() { return threads; }
        public int getDbType() { return dbType; }

        public static Builder builder() { return new Builder(); }

        public static class Builder {
            private double minSeqIdentity = Defaults.MIN_SEQ_IDENTITY;
            private double minCoverage = Defaults.MIN_COVERAGE;
            private int covMode = Defaults.COV_MODE;
            private int threads = 1;
            private int dbType = 2;

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
            public MMseqsConfig build() {
                return new MMseqsConfig(this);
            }
        }
    }
}