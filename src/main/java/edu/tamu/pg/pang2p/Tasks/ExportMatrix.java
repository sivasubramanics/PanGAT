package edu.tamu.pg.pang2p.Tasks;

import edu.tamu.pg.pang2p.Utils.Defaults;
import edu.tamu.pg.pang2p.Utils.FileUtils;
import edu.tamu.pg.pang2p.Utils.HelperFunctions;
import edu.tamu.pg.pang2p.Utils.Logger;
import picocli.CommandLine.*;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.Callable;

@Command(name = "ExportMatrix", description = "Export normalized coverage matrix and discrete genotype matrix from mapped_kmers matrix", sortOptions = false)
public class ExportMatrix implements Runnable, Callable<Integer> {

    @Option(names = {"-p", "--panmatrix"}, description = "Path to input pan-genome matrix file.", required = true)
    private Path panMatrixFile;

    @Option(names = {"-o", "--outprefix"}, description = "Prefix for output files.", required = true)
    private String outPrefix;

    @Option(names = {"-g", "--genotypes"}, description = "File with list of genotype IDs to include in the output.", required = false)
    private Path genotypeIdsFile;

    private final String CLASS_NAME = ExportMatrix.class.getSimpleName();

    public ExportMatrix() {
    }

    public ExportMatrix(Path panMatrixFile, String outPrefix, Path genotypeIdsFile) {
        this.panMatrixFile = panMatrixFile;
        this.outPrefix = outPrefix;
        this.genotypeIdsFile = genotypeIdsFile;
    }

    @Override
    public Integer call() throws Exception {
        Set<String> targetGenotypeIDs = new HashSet<>();
        if (genotypeIdsFile != null) {
            targetGenotypeIDs = getGenotypeIDs(genotypeIdsFile);
        }

        Path tempDir = Files.createTempDirectory("pang2p_transpose_");
        List<BufferedWriter> tempWriters = new ArrayList<>();
        List<Path> tempFiles = new ArrayList<>();
        List<String> selectedSampleNames = new ArrayList<>();
        List<String> markerIDs = new ArrayList<>();

        String genoFile = outPrefix + "_genotypes.tsv";
        String mapFile = outPrefix + "_map.tsv";

        try (BufferedReader reader = FileUtils.reader(panMatrixFile);
             BufferedWriter mapWriter = FileUtils.writer(mapFile)) {

            // Write header for the map file
            mapWriter.write("Name\tChromosome\tPosition");
            mapWriter.newLine();

            String line;

            while ((line = reader.readLine()) != null) {
                if (line.startsWith("##")) continue;

                List<String> fields = splitByTab(line);

                if (line.startsWith("#")) {
                    if (fields.size() < 14) {
                        Logger.error(CLASS_NAME, "Invalid header format in input file.");
                        return 1;
                    }

                    // Create disk buffers for each genotype/sample
                    for (int i = 13; i < fields.size(); i++) {
                        String sample = fields.get(i);
                        if (targetGenotypeIDs.isEmpty() || targetGenotypeIDs.contains(sample)) {
                            selectedSampleNames.add(sample);
                            Path tempFile = tempDir.resolve("sample_" + (i - 13) + ".tmp");
                            tempFiles.add(tempFile);
                            BufferedWriter bw = Files.newBufferedWriter(tempFile);
                            bw.write(sample); // Write Sample ID as row head
                            tempWriters.add(bw);
                        } else {
                            tempWriters.add(null); // Ignored genotype marker slot
                        }
                    }
                    continue;
                }

                if (tempWriters.isEmpty()) {
                    Logger.error(CLASS_NAME, "Data row encountered before header row (#).");
                    return 1;
                }

                String markerID = fields.get(0);
                int chrom = HelperFunctions.alphaToInt(fields.get(3), Defaults.CHROM_CODE_LENGTH) + 1;
                String pos = fields.get(4);

                markerIDs.add(markerID);

                // Export marker details to the map file stream
                mapWriter.write(markerID + "\t" + chrom + "\t" + pos);
                mapWriter.newLine();

                // Stream genotype calls directly to individual temporary files
                for (int i = 13; i < fields.size(); i++) {
                    BufferedWriter bw = tempWriters.get(i - 13);
                    if (bw != null) {
                        String fieldVal = fields.get(i);
                        int colonPos = fieldVal.indexOf(':');
                        String genotypeVal = (colonPos != -1) ? fieldVal.substring(0, colonPos) : fieldVal;

                        bw.write('\t');
                        bw.write(genotypeVal);
                    }
                }
            }

            // Close all sample buffer streams
            for (BufferedWriter bw : tempWriters) {
                if (bw != null) bw.close();
            }

            Logger.info(CLASS_NAME, "Writing final genotype matrix file...");

            try (BufferedWriter finalWriter = FileUtils.writer(genoFile)) {
                // Header row: Genotype followed by Marker IDs
                finalWriter.write("taxa\t" + String.join("\t", markerIDs));
                finalWriter.newLine();

                // Append each sample's transposed row to final output file
                for (Path tempFile : tempFiles) {
                    try (BufferedReader tempReader = Files.newBufferedReader(tempFile)) {
                        String tempLine;
                        while ((tempLine = tempReader.readLine()) != null) {
                            finalWriter.write(tempLine);
                            finalWriter.newLine();
                        }
                    }
                }
            }

            Logger.info(CLASS_NAME, "Export completed successfully. Map and genotype matrices created.");

        } catch (Exception e) {
            Logger.error(CLASS_NAME, "Error processing matrix: " + e.getMessage());
            return 1;
        } finally {
            // Clean up temporary files
            for (Path tempFile : tempFiles) {
                Files.deleteIfExists(tempFile);
            }
            Files.deleteIfExists(tempDir);
        }

        return 0;
    }

    private static List<String> splitByTab(String line) {
        List<String> list = new ArrayList<>();
        int start = 0;
        int end;
        while ((end = line.indexOf('\t', start)) != -1) {
            list.add(line.substring(start, end));
            start = end + 1;
        }
        list.add(line.substring(start));
        return list;
    }

    @Override
    public void run() {
        try {
            int exitCode = call();
            if (exitCode != 0) System.exit(exitCode);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    private Set<String> getGenotypeIDs(Path genotypeIdsFile) throws IOException {
        Set<String> genotypeIDs = new HashSet<>();
        try (BufferedReader reader = FileUtils.reader(genotypeIdsFile)) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (!line.trim().isEmpty()) {
                    genotypeIDs.add(line.trim());
                }
            }
        } catch (Exception e) {
            Logger.error(CLASS_NAME, "Error reading genotype IDs file: " + e.getMessage());
            throw e;
        }
        return genotypeIDs;
    }
}