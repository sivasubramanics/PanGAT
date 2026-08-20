package edu.tamu.pg.pang2p.Plugins;

import edu.tamu.pg.pang2p.Tasks.Cohort;
import edu.tamu.pg.pang2p.Tasks.DependencyChecker;
import edu.tamu.pg.pang2p.Tasks.DirectoryStructure;
import edu.tamu.pg.pang2p.Tasks.KmerCoverage;
import edu.tamu.pg.pang2p.Utils.HelperFunctions;
import edu.tamu.pg.pang2p.Utils.Logger;
import picocli.CommandLine;
import picocli.CommandLine.*;

import java.io.BufferedReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;

@Command(name = "Genotype", description = "Genotype the samples based on the reference pan haplotypes", sortOptions = false)
public class Genotype implements Callable<Integer> , Runnable {
    @Option(names = {"-b", "--basedir"}, description = "Base directory for the output", required = true)
    private Path baseDir;

    @Option(names = {"-s", "--samples"}, description = "Path to the samples file", required = true)
    private Path samples;

    @Option(names = {"-k", "--kmer"}, description = "K-mer size for KMC", required = false, defaultValue = "31")
    private int kmerSize;

    @Option(names = {"-t", "--threads"}, description = "Number of threads to use for processing", required = false, defaultValue = "2")
    private int threads;

    @Option(names = {"-m", "--min_count"}, description = "Minimum count for k-mer coverage calculation", required = false, defaultValue = "1")
    private int minCount;

    @Option(names = {"--memory"}, description = "Should the kmc db loaded into memory", defaultValue = "false")
    private boolean memory;

    @Option(names = {"--keep-kmc"}, description = "Keep the KMC database files after processing", defaultValue = "false")
    private boolean keepKmc;

    private final String CLASS_NAME = Genotype.class.getSimpleName();

    public Genotype() {}

    public Genotype(Path baseDir, Path samples, int threads) {
        this.baseDir = baseDir;
        this.samples = samples;
        this.threads = threads;
    }

    @Override
    public Integer call() throws Exception {

        // create directories under baseDir
        DirectoryStructure dirs = new DirectoryStructure(baseDir);
        dirs.createDirectories();

        HashMap<String, Set<Path>> samplePaths = new HashMap<>();
        try (BufferedReader reader = Files.newBufferedReader(samples)) {
            String line;
            int lineNumber = 0;
            while ((line = reader.readLine()) != null) {
                String[] parts = line.split(",");
                if (lineNumber == 0) {
                    if (parts[0].equals("id")) {
                        lineNumber++;
                        continue;
                    } else {
                        Logger.error(CLASS_NAME, "Invalid samples file format. Header should be 'id,path'");
                    }
                }
                if (parts.length == 2){
                    String sampleName = parts[0];
                    String genotypeFilePath = parts[1];
                    samplePaths.computeIfAbsent(sampleName, k -> new HashSet<>()).add(Path.of(genotypeFilePath));
                }
                else {
                    Logger.error(CLASS_NAME, "Invalid samples file format. File should have only two columns: sample_name and genotype_file_path, separated by a comma. header should be id,path");
                }
            }
        } catch (Exception e) {
            Logger.error(CLASS_NAME, "Error reading samples file: " + e.getMessage());
            return 1;
        }

        Path referenceFile = dirs.getSingleCopyReference();

        for (Map.Entry<String, Set<Path>> entry : samplePaths.entrySet()) {
            String sampleName = entry.getKey();
            Set<Path> genotypeFiles = entry.getValue();
            Path outCovFile = dirs.getGenotypeDir().resolve(sampleName + ".cov.tsv");
            KmerCoverage kmerCoverage = new KmerCoverage(referenceFile, sampleName, genotypeFiles.toArray(new Path[0]), outCovFile, kmerSize, threads, minCount, memory);
            kmerCoverage.call();
        }

        String[] samples = samplePaths.keySet().toArray(new String[0]);
        Path[] covFiles = new Path[samples.length];
        for (int i = 0; i < samples.length; i++) {
            covFiles[i] = dirs.getGenotypeDir().resolve(samples[i] + ".cov.tsv");
        }

        Cohort cohort = new Cohort(dirs.getPanClustersRepFile(), covFiles, samples, dirs.getCohortBase(), true);
        cohort.run();

        return 0;
    }

    @Override
    public void run() {
        try {
            HelperFunctions.printCommandLine(new CommandLine(this), CLASS_NAME);
            call();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }


}
