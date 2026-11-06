package edu.tamu.pg.pang2p.Plugins;

import edu.tamu.pg.pang2p.Tasks.Cohort;
import edu.tamu.pg.pang2p.Tasks.DirectoryStructure;
import edu.tamu.pg.pang2p.Tasks.KmerCoverage;
import edu.tamu.pg.pang2p.Utils.HelperFunctions;
import edu.tamu.pg.pang2p.Utils.Logger;
import picocli.CommandLine;
import picocli.CommandLine.*;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.Callable;

@Command(name = "genotype", description = "Genotype the samples based on the reference pan anchors", sortOptions = false)
public class Genotype implements Callable<Integer>, Runnable {

    @Option(names = {"-b", "--basedir"}, description = "Base directory for the output", required = true)
    private Path baseDir;

    @Option(names = {"-s", "--samples"}, description = "Path to the samples CSV file (format: id,path)", required = true)
    private Path samplesFile;

    @Option(names = {"-k", "--kmer"}, description = "K-mer size for KMC [31]", required = false, defaultValue = "31")
    private int kmerSize = 31;

    @Option(names = {"-t", "--threads"}, description = "Number of threads to use for processing [2]", required = false, defaultValue = "2")
    private int threads = 2;

    @Option(names = {"-m", "--min_count"}, description = "Minimum count for k-mer coverage calculation [1]", required = false, defaultValue = "1")
    private int minCount = 1;

    @Option(names = {"--memory"}, description = "Load the KMC database into memory [false]", required = false, defaultValue = "false")
    private boolean memory = false;

    @Option(names = {"--keep-kmc"}, description = "Keep the KMC database files after processing [false]", required = false, defaultValue = "false")
    private boolean keepKmc = false;

    private final String CLASS_NAME = this.getClass().getSimpleName();

    public Genotype() {
    }

    public Genotype(Path baseDir, Path samplesFile, int threads) {
        this.baseDir = baseDir;
        this.samplesFile = samplesFile;
        this.threads = threads;
    }

    @Override
    public Integer call() throws Exception {
        HelperFunctions.printCommandLine(new CommandLine(this), CLASS_NAME);
        validateInputs();

        executeGenotyping();
        return 0;
    }

    @Override
    public void run() {
        try {
            call();
        } catch (Exception e) {
            Logger.error(CLASS_NAME, "Execution failed: " + e.getMessage());
            throw new RuntimeException(e);
        }
    }

    /**
     * Executes the main workflow for sample k-mer coverage calculation and cohort aggregation.
     */
    private void executeGenotyping() throws IOException {
        DirectoryStructure dirs = new DirectoryStructure(baseDir);
        dirs.createDirectories();

        Map<String, Set<Path>> sampleMap = parseSamplesFile(samplesFile);
        if (sampleMap.isEmpty()) {
            Logger.error(CLASS_NAME, "No valid sample entries found in: " + samplesFile);
            return;
        }

        Path referenceFile = dirs.getSingleCopyReference();
        List<String> sampleNames = new ArrayList<>(sampleMap.keySet());
        List<Path> coverageFiles = new ArrayList<>();

        for (String sampleName : sampleNames) {
            Set<Path> inputPaths = sampleMap.get(sampleName);
            Path outputCovFile = dirs.getGenotypeDir().resolve(sampleName + ".cov.tsv");
            coverageFiles.add(outputCovFile);

            Logger.info(CLASS_NAME, "Calculating k-mer coverage for sample: " + sampleName);

            KmerCoverage kmerCoverage = new KmerCoverage(
                    referenceFile,
                    sampleName,
                    inputPaths.toArray(new Path[0]),
                    outputCovFile,
                    kmerSize,
                    threads,
                    minCount,
                    memory
            );

            try {
                kmerCoverage.call();
            } catch (Exception e) {
                Logger.error(CLASS_NAME, "Failed to compute coverage for sample: " + sampleName + " - " + e.getMessage());
            }
        }

        Logger.info(CLASS_NAME, "Aggregating coverage across cohort samples...");
        Cohort cohort = new Cohort(
                dirs.getPanClustersRepFile(),
                coverageFiles.toArray(new Path[0]),
                sampleNames.toArray(new String[0]),
                dirs.getCohortBase(),
                true
        );
        cohort.run();

        Logger.info(CLASS_NAME, "Genotyping completed successfully.");
    }

    /**
     * Parses the CSV file containing sample mapping (id,path).
     */
    private Map<String, Set<Path>> parseSamplesFile(Path filePath) throws IOException {
        Map<String, Set<Path>> sampleMap = new LinkedHashMap<>();

        try (BufferedReader reader = Files.newBufferedReader(filePath)) {
            String line;
            int lineNumber = 0;

            while ((line = reader.readLine()) != null) {
                lineNumber++;
                line = line.trim();
                if (line.isEmpty()) {
                    continue;
                }

                String[] parts = line.split(",");
                if (lineNumber == 1 && parts[0].equalsIgnoreCase("id")) {
                    continue; // Skip header line
                }

                if (parts.length == 2) {
                    String sampleId = parts[0].trim();
                    Path inputPath = Path.of(parts[1].trim());
                    sampleMap.computeIfAbsent(sampleId, k -> new HashSet<>()).add(inputPath);
                } else {
                    Logger.warning(CLASS_NAME, "Skipping malformed line " + lineNumber + " in " + filePath + ": " + line);
                }
            }
        }

        return sampleMap;
    }

    /**
     * Validates command-line parameters and file paths.
     */
    private void validateInputs() {
        if (!Files.exists(samplesFile)) {
            Logger.error(CLASS_NAME, "No samples file found: " + samplesFile);
        }
        if (threads < 1) {
            Logger.error(CLASS_NAME, "Number of threads must be at least 1.");
        }
        if (kmerSize <= 0 || kmerSize % 2 == 0) {
            Logger.error(CLASS_NAME, "K-mer size must be a positive odd integer.");
        }
        if (minCount < 1) {
            Logger.error(CLASS_NAME, "Minimum k-mer count must be at least 1.");
        }
    }
}