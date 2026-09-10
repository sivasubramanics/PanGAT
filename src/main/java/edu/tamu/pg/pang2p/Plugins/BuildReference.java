package edu.tamu.pg.pang2p.Plugins;

import edu.tamu.pg.pang2p.Data.*;
import edu.tamu.pg.pang2p.Tasks.*;
import edu.tamu.pg.pang2p.Utils.*;

import picocli.CommandLine;
import picocli.CommandLine.*;

import java.nio.file.Path;
import java.util.concurrent.*;

@Command(name = "build-ref", description = "Prepare and Build pan anchor reference for PanG2P analysis", sortOptions = false)
public class BuildReference implements Callable<Integer>, Runnable {

    @Option(names = {"-g", "--genomes"}, description = "Path to the reference genome map file (columns: id, path) in CSV format", required = true)
    private Path genomes;

    @Option(names = {"-b", "--basedir"}, description = "Base directory for the output", required = true)
    private Path baseDir;

    @Option(names = {"-t", "--threads"}, description = "Number of threads to use for processing", required = false, defaultValue = "2")
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

    @Option(names = {"--min_length"}, description = "Minimum length for single copy regions to be considered", required = false, defaultValue = "200")
    private int minLength;

    @Option(names = {"--force"}, description = "Force overwrite of existing files", required = false)
    private boolean isForce;

    private static final String CLASS_NAME = BuildReference.class.getSimpleName();

    public Integer call() throws Exception {

        // Check if necessary tools are accessible
        DependencyChecker dependencyChecker = new DependencyChecker();
        dependencyChecker.check();

        // Create directories under baseDir
        DirectoryStructure dirs = new DirectoryStructure(baseDir);
        dirs.createDirectories();

        // Prepare raw reference genome data
        PrepareData data = new PrepareData(dirs, genomes, threads, isForce);
        data.prepare();

        // Load genome map from genomesMap.csv
        GenomeMap genomeMap = new GenomeMap(dirs);

        // Extract single copy regions (Includes KMC step & automatic cleanup)
        ExtractAnchors extractAnchors = new ExtractAnchors(dirs, isForce);
        extractAnchors.run(genomeMap, kmerSize, minCount, minLength, threads);

        // Set commandline options for mmseqs2 clustering
        RunMMseqs.MMseqsConfig mmseqsConfig = RunMMseqs.MMseqsConfig.builder()
                .threads(threads)
                .minSeqIdentity(minSeqIdentity)
                .minCoverage(minCoverage)
                .covMode(covMode)
                .isForce(isForce)
                .build();

        RunMMseqs runMMseqs = new RunMMseqs(dirs, mmseqsConfig);
        runMMseqs.run(genomeMap);

        PanCluster panCluster = new PanCluster(dirs, threads, isForce);
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