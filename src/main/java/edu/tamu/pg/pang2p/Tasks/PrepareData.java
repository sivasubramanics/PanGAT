package edu.tamu.pg.pang2p.Tasks;

import edu.tamu.pg.pang2p.Data.*;
import edu.tamu.pg.pang2p.Utils.*;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.*;

public class PrepareData {
    private static final String CLASS_NAME = PrepareData.class.getSimpleName();
    private final Path genomesFile;
    private final int nThreads;
    private final DirectoryStructure dirs;

    public PrepareData(DirectoryStructure dirs, Path genomesFile, int nThreads) {
        this.dirs = Objects.requireNonNull(dirs, "DirectoryStructure cannot be null");
        this.genomesFile = Objects.requireNonNull(genomesFile, "Genomes file cannot be null");
        if (!FileUtils.isFile(genomesFile)) {
            throw new IllegalArgumentException("Genomes file does not exist: " + genomesFile);
        }
        this.nThreads = nThreads > 0 ? nThreads : 1;
    }

    public void prepare() throws IOException {
        Path genomesDir = dirs.getDirectory(PipelineSteps.GENOMES);
        if (!FileUtils.isDir(genomesDir)) {
            Logger.error(CLASS_NAME, "Genome directory does not exist: " + genomesDir);
        }

        Map<String, Path> genomes = readGenomesFile(genomesFile);

        Logger.info(CLASS_NAME, "Preparing reference genomes from map file: " + genomesFile);

        ExecutorService executor = Executors.newFixedThreadPool(nThreads);
        List<Future<List<String[]>>> futures = new ArrayList<>();

        int index = 0;
        for (Map.Entry<String, Path> entry : genomes.entrySet()) {
            String genomeId = entry.getKey();
            Path genomePath = entry.getValue();
            Path destPath = genomesDir.resolve(genomeId + ".fasta");

            int finalIndex = index;
            futures.add(executor.submit(() -> processGenome(destPath, genomePath, finalIndex, genomeId)));
            index++;
        }

        executor.shutdown();

        int success = 0, failed = 0;
        List<String[]> allProcessedData = new ArrayList<>();

        for (Future<List<String[]>> future : futures) {
            try {
                allProcessedData.addAll(future.get());
                success++;
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                throw new IOException("Genome processing interrupted.", e);
            } catch (ExecutionException e) {
                Logger.error(CLASS_NAME, "Unexpected error while processing genomes: " + e.getCause());
                failed++;
            }
        }

        Logger.info(CLASS_NAME, String.format(
                "Read %d entries from genomes map file. Success: %d, Failed: %d.",
                genomes.size(), success, failed
        ));

        if (failed > 0) {
            Logger.error(CLASS_NAME, "Some genomes failed to process. Check logs.");
        }

        Path genomesMapFile = dirs.getBaseDir().resolve("genomesMap.csv");
        try (BufferedWriter writer = Files.newBufferedWriter(genomesMapFile)) {
            for (String[] data : allProcessedData) {
                writer.write(String.join(",", data));
                writer.newLine();
            }
        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Error writing genomes map file: " + e.getMessage());
        }
    }

    private Map<String, Path> readGenomesFile(Path genomesFile) throws IOException {
        Map<String, Path> genomesMap = new LinkedHashMap<>();
        Logger.info(CLASS_NAME, "Reading genomes map from: " + genomesFile);

        try (BufferedReader reader = new BufferedReader(new FileReader(genomesFile.toString()))) {
            String header = reader.readLine();
            if (header == null) {
                throw new IOException("Genomes map file is empty: " + genomesFile);
            }

            String[] headerParts = header.split(",");
            if (headerParts.length != 2
                    || !headerParts[0].trim().equalsIgnoreCase("id")
                    || !headerParts[1].trim().equalsIgnoreCase("path")) {
                throw new IOException("Expected header: id,path but got: " + header);
            }

            String line;
            int lineNumber = 2;
            while ((line = reader.readLine()) != null) {
                String[] parts = line.split(",", 2);
                if (parts.length != 2) {
                    Logger.error(CLASS_NAME, String.format(
                            "Invalid line %d in genomes map file: '%s'. Skipping.", lineNumber, line));
                    lineNumber++;
                    continue;
                }

                String genomeId = parts[0].trim();
                String genomePathStr = parts[1].trim();

                if (genomeId.isEmpty() || genomePathStr.isEmpty()) {
                    Logger.warning(CLASS_NAME,
                            String.format("Empty ID or path on line %d. Skipping entry.", lineNumber));
                    lineNumber++;
                    continue;
                }

                Path genomePath = Path.of(genomePathStr);
                if (!FileUtils.isFile(genomePath)) {
                    Logger.error(CLASS_NAME, String.format(
                            "Invalid genome path on line %d: %s. File does not exist. Skipping.", lineNumber, genomePath));
                    lineNumber++;
                    continue;
                }

                if (genomesMap.containsKey(genomeId)) {
                    Logger.warning(CLASS_NAME,
                            String.format("Duplicate entry for ID '%s' on line %d. Overwriting previous value.", genomeId, lineNumber));
                }
                genomesMap.put(genomeId, genomePath);
                lineNumber++;
            }
        }
        return genomesMap;
    }

    /**
     * Reads, rewrites and indexes a genome FASTA. Throws IOException on error.
     */
    private List<String[]> processGenome(Path destPath, Path genomePath, int index, String genomeId) throws IOException {
        Logger.info(CLASS_NAME, String.format(
                "Processing genome [%s]: %s => %s", genomeId, genomePath, destPath));
        String genomeCode = HelperFunctions.encode(index, Defaults.CHROM_CODE_LENGTH);
        int chromNum = 0;
        List<String[]> processedData = new ArrayList<>();
        try (
                FastaWriter fw = new FastaWriter(destPath);
                FastaReader freader = new FastaReader(genomePath)
        ) {
            for (Fasta fa : freader) {
                String chromCode = HelperFunctions.encode(chromNum, Defaults.CHROM_CODE_LENGTH);
                String newID = genomeCode + chromCode;
                String newDescription = String.format("genome=%s,chrom=%s", genomeId, fa.getName());
                processedData.add(new String[]{newID, genomeId, fa.getName()});
                fa.setName(newID);
                fa.setDescription(newDescription);
                fw.write(fa, false);
                chromNum++;
            }
        } catch (IOException e) {
            Logger.error(CLASS_NAME, String.format("Error copying genome [%s]: %s to %s: %s", genomeId, genomePath, destPath, e.getMessage()));
            throw e;
        }
        // create fasta index for each fasta file
        try {
            new FastaIndex(destPath);
        } catch (Exception e) {
            Logger.error(CLASS_NAME, String.format("Error creating FASTA index for genome [%s]: %s: %s", genomeId, destPath, e.getMessage()));
            throw new IOException("Error indexing FASTA for genome [" + genomeId + "]", e);
        }
        return processedData;
    }
}