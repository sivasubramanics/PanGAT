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

    @Option(names = {"--cov_mode"}, description = "Coverage mode for MMseqs2 clustering", required = false, defaultValue = "1")
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
//        PrepareData data = new PrepareData(dirs, genomes, threads);
//        data.prepare();

        // load genome map from genomesMap.csv
        GenomeMap genomeMap = new GenomeMap(dirs);

        // set commandline options for KMC
//        RunKMC.KMCConfig config = RunKMC.KMCConfig.builder()
//                .threads(threads)
//                .kmerSize(kmerSize)
//                .minCount(minCount)
//                .build();

        // run KMC for each genome
//        RunKMC runKMC = new RunKMC(dirs, config);
//        runKMC.run(genomeMap);

        // extract single copy regions for each genome
//        ExtractSCR extractSingleCopyRegions = new ExtractSCR(dirs);
//        extractSingleCopyRegions.run(genomeMap, minCount, minLength, threads);

        // set commandline options for mmseqs2 clustering
//        RunMMseqs.MMseqsConfig mmseqsConfig = RunMMseqs.MMseqsConfig.builder()
//                .threads(threads)
//                .minSeqIdentity(minSeqIdentity)
//                .minCoverage(minCoverage)
//                .covMode(covMode)
//                .build();
//        // run MMseqs2 clustering
//        RunMMseqs runMMseqs = new RunMMseqs(dirs, mmseqsConfig);
//        runMMseqs.run(genomeMap);

        PanCluster panCluster = new PanCluster(dirs);
        panCluster.run(genomeMap);


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
}


