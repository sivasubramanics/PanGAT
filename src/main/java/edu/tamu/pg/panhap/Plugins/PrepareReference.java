package edu.tamu.pg.panhap.Plugins;

import edu.tamu.pg.panhap.Data.*;
import edu.tamu.pg.panhap.Tasks.*;
import edu.tamu.pg.panhap.Utils.*;

import picocli.CommandLine;
import picocli.CommandLine.*;
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.*;

@Command(name = "PrepareReferences", description = "Prepare the reference genome for PanHap", sortOptions = false)
public class PrepareReference implements Callable<Integer>, Runnable {

    @Option(names = {"-g", "--genomes"}, description = "Path to the reference genome map file (columns: id, path) in CSV format", required = true)
    private Path genomes;

    @Option(names = {"-b", "--basedir"}, description = "Base directory for the output", required = true)
    private Path baseDir;

    @Option(names = {"-t", "--threads"}, description = "Number of threads to use for processing", required = false, defaultValue = "1")
    private int threads;

    @Option(names = {"-k", "--kmer"}, description = "K-mer size for KMC", required = false, defaultValue = "31")
    private int kmerSize;

    @Option(names = {"-m", "--min_count"}, description = "Minimum count for Single copy region extraction", required = false, defaultValue = "1")
    private int minCount;

    @Option(names = {"--min_seq_identity"}, description = "Minimum sequence identity for MMseqs2 clustering", required = false, defaultValue = "0.9")
    private double minSeqIdentity;

    @Option(names = {"--min_coverage"}, description = "Minimum coverage for MMseqs2 clustering", required = false, defaultValue = "0.8")
    private double minCoverage;

    @Option(names = {"--cov_mode"}, description = "Coverage mode for MMseqs2 clustering", required = false, defaultValue = "0")
    private int covMode;

    @Option(names = {"--min_length"}, description = "Minimum length for single copy regions to be considered", required = false, defaultValue = "100")
    private int minLength;

    private static final String CLASS_NAME = PrepareReference.class.getSimpleName();

    private Path genomesDir;
    private Path kmcDir;
    private Path singleCopyDir;

    double HEAPSIZE = Runtime.getRuntime().maxMemory() / (1024.0 * 1024 * 1024);
    int COUNTER_DUMP_SIZE = 100_000;


    public Integer call() throws Exception {

        // check if necessary tools are accessible
        DependencyChecker dependencyChecker = new DependencyChecker();
        dependencyChecker.check();

        // create directories under baseDir
        DirectoryStructure dirs = new DirectoryStructure(baseDir);
        dirs.createDirectories();

        // prepareReference();
        PrepareData data = new PrepareData(dirs, genomes, threads);
        data.prepare();

        // load genome map from genomesMap.csv
        GenomeMap genomeMap = new GenomeMap(dirs);

        // set commandline options for KMC
        RunKMC.KMCConfig config = RunKMC.KMCConfig.builder()
                .threads(threads)
                .kmerSize(kmerSize)
                .minCount(minCount)
                .build();

        // run KMC for each genome
        RunKMC runKMC = new RunKMC(dirs, config);
        runKMC.run(genomeMap);

        // extract single copy regions for each genome
        ExtractSCR extractSingleCopyRegions = new ExtractSCR(dirs);
        extractSingleCopyRegions.run(genomeMap, minCount, minLength, threads);

        // set commandline options for mmseqs2 clustering
        RunMMseqs.MMseqsConfig mmseqsConfig = RunMMseqs.MMseqsConfig.builder()
                .threads(threads)
                .minSeqIdentity(minSeqIdentity)
                .minCoverage(minCoverage)
                .covMode(covMode)
                .build();
        // run MMseqs2 clustering
        RunMMseqs runMMseqs = new RunMMseqs(dirs, mmseqsConfig);
        runMMseqs.run(genomeMap);

        PanCluster panCluster = new PanCluster(dirs);
        panCluster.run();


        return 0;
    }

    public void run() {
        try {
            HelperFunctions.printCommandLine(new CommandLine(this), CLASS_NAME);
            call();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    /***
     * Prepares the reference genome by
     *      reading the genome map file,
     *      creating the necessary directory structure,
     *      and processing each genome in parallel using the specified number of threads.
     * @throws IOException
     */
//    private void prepareReference() throws IOException {
//
//        // read genome map file
//        Map<String, String> genomesMap = readGenomesMap(genomes);
//
//
//        // fix and copy the genomes to output dir
//        fixGenomes(genomesMap);
//
//        // run KMC for each genome to do kmer counting
//        runKMC(genomesMap);
//
//        // read through the fasta and kmer table and extract single copy regions
//        extractSingleCopyRegions(genomesMap);
//
//        // concatenate and run mmseqs to cluster the single copy regions
//        clusterSCR(genomesMap);
//
//        Logger.info(CLASS_NAME, "Reference genome preparation completed.");
//    }

    private void clusterSCR(Map<String, String> genomesMap) throws IOException {
        Logger.info(CLASS_NAME, "Clustering single copy regions using MMseqs2.");
        String scrFasta = singleCopyDir.resolve("all_single_copy_regions.fasta").toString();

        // concatenate all single copy regions into one fasta file
        try (FastaWriter mergedFasta = new FastaWriter(scrFasta)) {
            for (String id : genomesMap.keySet()) {
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

        // createDB with mmseqs
        String dbFile = singleCopyDir.resolve("all_single_copy_regions.db").toString();
        StringBuilder mmseqsCreateDBCmd = new StringBuilder("mmseqs createdb ");
        mmseqsCreateDBCmd.append(scrFasta).append(" ").append(dbFile).append(" ");
        mmseqsCreateDBCmd.append("--threads").append(" ").append(threads).append(" ");
        mmseqsCreateDBCmd.append("--dbtype").append(" ").append("2");

        String createDBLog = singleCopyDir.resolve("mmseqs_createDB.log").toString();
        Logger.info(CLASS_NAME, "Creating DB for mmseqs clustering. Log file: " + createDBLog);
        Logger.info(CLASS_NAME, "MMseqs2 CMD: " + mmseqsCreateDBCmd.toString());
        HelperFunctions.tryExec(mmseqsCreateDBCmd.toString(), createDBLog, createDBLog);

        // cluster with mmseqs
        String tmpDir = singleCopyDir.resolve("tmp").toString();
        // if tmpDir exists delete
        if (FileUtils.isDir(tmpDir)) {
            FileUtils.rmdir(tmpDir);
        } else if (FileUtils.isFile(tmpDir)) {
            FileUtils.deleteRecursively(new File(tmpDir));
        }

        FileUtils.mkdir(tmpDir);
        String scrClusters = singleCopyDir.resolve("clusters").toString();

        StringBuilder mmseqsClusterCmd = new StringBuilder("mmseqs cluster ");
        mmseqsClusterCmd.append(dbFile).append(" ").append(scrClusters).append(" ").append(tmpDir).append(" ");
        mmseqsClusterCmd.append("--threads").append(" ").append(threads).append(" ");
        mmseqsClusterCmd.append("--min-seq-id").append(" ").append(minSeqIdentity).append(" ");
        mmseqsClusterCmd.append("-c").append(" ").append(minCoverage).append(" ");
        mmseqsClusterCmd.append("--cov-mode").append(" ").append(covMode);

        String mmseqsClusterLog = singleCopyDir.resolve("mmseqs_clustering.log").toString();
        Logger.info(CLASS_NAME, "Clustering single copy regions with MMseqs2. Log file: " + mmseqsClusterLog);
        Logger.info(CLASS_NAME, "MMseqs2 CMD: " + mmseqsClusterCmd.toString());
        HelperFunctions.tryExec(mmseqsClusterCmd.toString(), mmseqsClusterLog, mmseqsClusterLog);

        // convert to tsv
        String clusterFile = singleCopyDir.resolve("clusters.tsv").toString();

        StringBuilder createTSVCmd = new StringBuilder("mmseqs createtsv ");
        createTSVCmd.append(dbFile).append(" ").append(dbFile).append(" ").append(scrClusters).append(" ").append(clusterFile).append(" ");
        createTSVCmd.append("--threads").append(" ").append(threads);

        String createTSVLog = singleCopyDir.resolve("mmseqs_clustering.log").toString();
        Logger.info(CLASS_NAME, "Converting clustering results to TSV. Log file: " + createTSVLog);
        Logger.info(CLASS_NAME, "MMseqs2 CMD: " + createTSVCmd.toString());
        HelperFunctions.tryExec(createTSVCmd.toString(), createTSVLog, createTSVLog);

        FileUtils.rmdir(singleCopyDir.resolve("tmp"));
        Logger.info(CLASS_NAME, "Clustering of single copy regions completed.");
    }

    private void fixGenomes(Map<String, String> genomesMap) throws IOException {
        // create directory structure under baseDir
        // baseDir + "/00.genomes"

    }

    private void extractSingleCopyRegions(Map<String, String> genomesMap) throws IOException {

        FileUtils.mkdir(singleCopyDir);

        for (Map.Entry<String, String> entry : genomesMap.entrySet()) {
            String id = entry.getKey();

            KMC kmc = new KMC(kmcDir.resolve(id).toString(), true);
            FastaReader genome = new FastaReader(genomesDir.resolve(id + ".fasta").toString());
            int k = kmc.getKmerLength();

            // Larger buffer: 256KB instead of the default 8KB
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(singleCopyDir.resolve(id + ".fasta").toString()), 256 * 1024)) {
                int counter = 0;
                List<Fasta> singleCopyRegions = new ArrayList<>();
                int numSingleCopy = 0;
                for (Fasta chrom : genome) {

                    StringBuilder singleCopySeq = new StringBuilder();
                    boolean isBeginning = true;
                    int regionStart = -1;     // track where the region began
                    int currentEnd = -1;      // track the latest position seen

                    KmerIterator it = chrom.kmers(k, kmc.getPrefixLength(), kmc.isBothStrands());

                    while (it.hasNext()) {
                        KmerPosition kp = it.next();
                        int start = kp.getStartPos();
                        Kmer kmer = kp.getKmer();
                        int kmerCount = kmc.getCount(kmer);
//                        int kmerCount = 1;

//                        if (kmerCount == 0){
//                            Logger.error(CLASS_NAME, "Kmer not found in KMC database for genome: " + id + ", chromosome: " + chrom.getName() + ", kmer: " + Arrays.toString(kmer.getOriginalKmer()));
//                        }

                        if (kmerCount > 0 && kmerCount <= minCount) {
                            if (isBeginning) {
                                regionStart = start;
                                singleCopySeq.append(kmer.getOriginalKmer());
                                isBeginning = false;
                            } else {
                                // Append only the new last char — avoids a full String allocation
                                singleCopySeq.append(kmer.getOriginalKmer()[k - 1]);
                            }
                            currentEnd = start + k;

                        } else {
                            if (!singleCopySeq.isEmpty()) {
                                if (singleCopySeq.length() >= minLength) {
                                    singleCopyRegions.add(new Fasta(numSingleCopy, chrom.getName() + ":" + regionStart + "-" + currentEnd, singleCopySeq.toString()));
                                    if (counter == COUNTER_DUMP_SIZE) {
                                        for (Fasta scr : singleCopyRegions) {
                                            writer.write(scr.toFastaString(true));
                                        }
                                        singleCopyRegions.clear();
                                        Logger.info(CLASS_NAME, "Dumped " + numSingleCopy + " single copy regions for genome: " + id);
                                    }
                                    counter++;
                                    numSingleCopy++;
                                }
                                singleCopySeq.setLength(0);
                            }
                            isBeginning = true;
                            regionStart = -1;
                            it.skipTo(start + k);
                        }
                    }

                    // ── Flush any region still open at end of chromosome ──
                    if (!singleCopySeq.isEmpty()) {
                        if (singleCopySeq.length() >= minLength) {
                            singleCopyRegions.add(new Fasta(numSingleCopy, chrom.getName() + ":" + regionStart + "-" + currentEnd, singleCopySeq.toString()));
                            counter++;
                            numSingleCopy++;
                        }
                        singleCopySeq.setLength(0);
                    }

                    if (!singleCopyRegions.isEmpty()) {
                        for (Fasta scr: singleCopyRegions) {
                            writer.write(scr.toFastaString(true));
                        }
                    }
                    Logger.info(CLASS_NAME, "Dumped " + numSingleCopy + " single copy regions for genome: " + id);
                }
                writer.close();
                kmc.close();
                FileUtils.rmdir(kmcDir.resolve(id + ".kmc_pre"));
                FileUtils.rmdir(kmcDir.resolve(id + ".kmc_suf"));
            } catch (IOException e) {
                Logger.error(CLASS_NAME, "Error writing single copy regions for genome: " + id);
                throw e;
            }
        }
    }

    private void runKMC(Map<String, String> genomesMap) throws IOException {

        // create tmp directory for kmc under base_dir / 01.kmc /tmp
        FileUtils.mkdir(kmcDir.resolve("tmp"));


        for  (String id : genomesMap.keySet()) {
            String logKMCpath = kmcDir.resolve(id + ".kmc.log").toString();
            StringBuilder kmcCmd =  new StringBuilder("kmc ");
            kmcCmd.append("-k").append(kmerSize).append(" ");
            kmcCmd.append("-p").append(Defaults.SIGNATURE_LENGTH).append(" ");
            kmcCmd.append("-ci").append("1").append(" ");
            kmcCmd.append("-cx").append(minCount).append(" ");
            kmcCmd.append("-t").append(threads).append(" ");
            kmcCmd.append("-m").append(HEAPSIZE).append(" ");
            kmcCmd.append("-fm").append(" ");
            kmcCmd.append(genomesDir).append(File.separator).append(id).append(".fasta").append(" ");
            kmcCmd.append(kmcDir).append(File.separator).append(id).append(" ");
            kmcCmd.append(kmcDir).append(File.separator).append("tmp");

            Logger.info(CLASS_NAME, "Running KMC for genome: " + id);
            Logger.info(CLASS_NAME, "KMC command: " + kmcCmd.toString());
            HelperFunctions.tryExec(kmcCmd.toString(), logKMCpath, logKMCpath);
        }
        FileUtils.rmdir(kmcDir.resolve("tmp"));
    }





    private void checkDependencies() throws IOException {
        List<String> requiredTools = List.of(
                "kmc",
                "mmseqs"
        );

        List<String> missingTools = new ArrayList<>();

        for (String tool : requiredTools) {
            if (!HelperFunctions.isInstalled(tool)) {
                missingTools.add(tool);
            }
        }

        if (!missingTools.isEmpty()) {
            String message = "Missing required tools: " + String.join(", ", missingTools)
                    + ". Please install them and ensure they are available in your system PATH.";

            Logger.error(CLASS_NAME, message);
            throw new IOException(message);
        }
    }
}


