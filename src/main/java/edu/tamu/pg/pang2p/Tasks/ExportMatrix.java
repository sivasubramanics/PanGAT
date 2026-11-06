package edu.tamu.pg.pang2p.Tasks;

import edu.tamu.pg.pang2p.Utils.Defaults;
import edu.tamu.pg.pang2p.Utils.FileUtils;
import edu.tamu.pg.pang2p.Utils.HelperFunctions;
import edu.tamu.pg.pang2p.Utils.Logger;
import picocli.CommandLine;
import picocli.CommandLine.*;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.Callable;

@Command(name = "export-matrix", description = "Export normalized coverage matrix and discrete genotype matrix from mapped_kmers matrix", sortOptions = false)
public class ExportMatrix implements Runnable, Callable<Integer> {

    @Option(names = {"-p", "--panmatrix"}, description = "Path to input pan-genome matrix file.", required = true)
    private Path panMatrixFile;

    @Option(names = {"-o", "--outprefix"}, description = "Prefix for output files.", required = true)
    private String outPrefix;

    @Option(names = {"-i", "--include"}, description = "File with list of genotype IDs to include in the output.", required = false)
    private Path includeGenotypeIdsFile;

    @Option(names = {"-e", "--exclude"}, description = "File with list of genotype IDs to exclude from the output.", required = false)
    private Path excludeGenotypeIdsFile;

    @Option(names = {"--min-maf"}, description = "Minimum minor allele frequency for filtering markers [0.05]", required = false, defaultValue = "0.05")
    private double minMAF;

    @Option(names = {"--max-missing-marker"}, description = "Maximum missing data proportion for filtering markers [0.1]", required = false, defaultValue = "0.1")
    private double maxMissingMarker;

    @Option(names = {"--max-heterozygosity"}, description = "Maximum heterozygosity for filtering markers [0.5]", required = false, defaultValue = "0.5")
    private double maxHeterozygosity;

    @Option(names = {"--cv"}, description = "CV cutoff for filtering markers [0.3]", required = false, defaultValue = "0.3")
    private double cvCutoff;

    @Option(names = {"--gapit"}, description = "Convert output to GAPIT format (requires R and GAPIT installed).", required = false, defaultValue = "false")
    private boolean gapit = false;

    private final String CLASS_NAME = ExportMatrix.class.getSimpleName();

    public ExportMatrix() {
    }

    public ExportMatrix(Path panMatrixFile, String outPrefix) {
        this.panMatrixFile = panMatrixFile;
        this.outPrefix = outPrefix;
    }

    @Override
    public Integer call() throws Exception {
        HelperFunctions.printCommandLine(new CommandLine(this), CLASS_NAME);
        validateInputs();

        try (BufferedReader reader = FileUtils.reader(panMatrixFile);
             BufferedWriter writer = FileUtils.writer(outPrefix + ".genotype.tsv")) {

            String line;
            Set<String> includeGenotypes = loadGenotypeSet(includeGenotypeIdsFile);
            Set<String> excludeGenotypes = loadGenotypeSet(excludeGenotypeIdsFile);

            if (includeGenotypes != null && excludeGenotypes != null) {
                Set<String> overlap = new HashSet<>(excludeGenotypes);
                overlap.retainAll(includeGenotypes);
                for (String genotypeId : overlap) {
                    Logger.warning(CLASS_NAME, "Genotype ID " + genotypeId + " is present in both include and exclude lists. It will be included.");
                }
                excludeGenotypes.removeAll(overlap);
            }

            int[] targetIndices = null;
            int totalSamplesInHeader = 0;
            int totalMarkers = 0;
            int failedMAF = 0;
            int failedMissingMarker = 0;
            int failedHeterozygosity = 0;
            int failedCV = 0;
            int filteredMarkers = 0;

            long globalMissing = 0, globalHomRef = 0, globalHet = 0, globalHomAlt = 0;

            // Track seen positions to handle duplicates across genomes/anchors
            Set<String> seenCoordinates = new HashSet<>();

            StringBuilder outLine = new StringBuilder(1024);

            Logger.info(CLASS_NAME, "Reading input file " + panMatrixFile + "...");
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("##")) continue;

                if (line.startsWith("#")) {
                    String[] header = line.split("\t");
                    totalSamplesInHeader = header.length - 13;

                    List<Integer> indexList = new ArrayList<>();
                    for (int i = 13; i < header.length; i++) {
                        String genotypeId = header[i];
                        if (includeGenotypes != null && !includeGenotypes.contains(genotypeId)) continue;
                        if (excludeGenotypes != null && excludeGenotypes.contains(genotypeId)) continue;
                        indexList.add(i);
                    }

                    targetIndices = indexList.stream().mapToInt(Integer::intValue).toArray();

                    // Write header
                    outLine.append("ID\tCHR\tSTART\tEND");
                    for (int idx : targetIndices) {
                        outLine.append("\t").append(header[idx]);
                    }
                    outLine.append("\n");
                    writer.write(outLine.toString());
                    outLine.setLength(0);
                    continue;
                }

                String[] tokens = line.split("\t");
                if (tokens.length < 2 || targetIndices == null) continue;

                totalMarkers++;

                boolean cvFailed = Double.parseDouble(tokens[12]) < cvCutoff;
                if (cvFailed) {
                    failedCV++;
                }

                MarkerStats stats = new MarkerStats(tokens, targetIndices);

                if (stats.missingProportion() > maxMissingMarker) {
                    failedMissingMarker++;
                    stats.setGood(false);
                }
                if (stats.maf() < minMAF) {
                    failedMAF++;
                    stats.setGood(false);
                }
                if (stats.heterozygosity() > maxHeterozygosity) {
                    failedHeterozygosity++;
                    stats.setGood(false);
                }

                if (cvFailed || !stats.isGood()) {
                    continue;
                }

                filteredMarkers++;
                globalMissing += stats.alleleCounts[0];
                globalHomRef += stats.alleleCounts[1];
                globalHet += stats.alleleCounts[2];
                globalHomAlt += stats.alleleCounts[3];

                // Resolve duplicate CHR + START positions
                String chr = String.valueOf(HelperFunctions.alphaToInt(tokens[3], Defaults.CHROM_CODE_LENGTH) + 1);
                long startPos = Long.parseLong(tokens[4]);
                long endPos = Long.parseLong(tokens[5]);
                long lengthDiff = endPos - startPos; // preserve original marker length

                String coordKey = chr + ":" + startPos;
                while (seenCoordinates.contains(coordKey)) {
                    startPos++;
                    coordKey = chr + ":" + startPos;
                }
                seenCoordinates.add(coordKey);

                // Recalculate end to match the shifted start if you want to keep width constant
                endPos = startPos + lengthDiff;

                outLine.append(tokens[0]).append("\t")  // ID
                        .append(chr).append("\t")       // CHR
                        .append(startPos).append("\t")  // Start (adjusted if duplicate)
                        .append(endPos);                // End (adjusted if duplicate)
                for (int gt : stats.gtData) {
                    outLine.append("\t").append(gt);
                }
                outLine.append("\n");
                writer.write(outLine.toString());
                outLine.setLength(0);
            }

            // Summary Logs
            int failedTotal = totalMarkers - filteredMarkers;
            double retainPct = totalMarkers > 0 ? (filteredMarkers * 100.0 / totalMarkers) : 0.0;

            Logger.info(CLASS_NAME, "-------------------- Sample Summary --------------------");
            Logger.info(CLASS_NAME, String.format("%-35s: %10d", "Total samples", totalSamplesInHeader));
            Logger.info(CLASS_NAME, String.format("%-35s: %10d", "Retained samples", targetIndices != null ? targetIndices.length : 0));

            Logger.info(CLASS_NAME, "-------------------- Marker Filters --------------------");
            Logger.info(CLASS_NAME, String.format("%-35s: %10d", "Total markers input", totalMarkers));
            Logger.info(CLASS_NAME, String.format("%-35s: %10d", String.format("Failed CV (< %.2f)", cvCutoff), failedCV));
            Logger.info(CLASS_NAME, String.format("%-35s: %10d", String.format("Failed missing (> %.2f)", maxMissingMarker), failedMissingMarker));
            Logger.info(CLASS_NAME, String.format("%-35s: %10d", String.format("Failed MAF (< %.2f)", minMAF), failedMAF));
            Logger.info(CLASS_NAME, String.format("%-35s: %10d", String.format("Failed het (> %.2f)", maxHeterozygosity), failedHeterozygosity));
            Logger.info(CLASS_NAME, String.format("%-35s: %10d", "Total markers removed", failedTotal));
            Logger.info(CLASS_NAME, String.format("%-35s: %10d (%.2f%%)", "Passed markers", filteredMarkers, retainPct));

            Logger.info(CLASS_NAME, "------------------ Genotype Distribution ---------------");
            long totalCalls = globalHomRef + globalHet + globalHomAlt + globalMissing;
            Logger.info(CLASS_NAME, String.format("%-35s: %10d (%.2f%%)", "Hom Ref (0)", globalHomRef, totalCalls > 0 ? (globalHomRef * 100.0 / totalCalls) : 0));
            Logger.info(CLASS_NAME, String.format("%-35s: %10d (%.2f%%)", "Het (1)", globalHet, totalCalls > 0 ? (globalHet * 100.0 / totalCalls) : 0));
            Logger.info(CLASS_NAME, String.format("%-35s: %10d (%.2f%%)", "Hom Alt (2)", globalHomAlt, totalCalls > 0 ? (globalHomAlt * 100.0 / totalCalls) : 0));
            Logger.info(CLASS_NAME, String.format("%-35s: %10d (%.2f%%)", "Missing (-1)", globalMissing, totalCalls > 0 ? (globalMissing * 100.0 / totalCalls) : 0));
            Logger.info(CLASS_NAME, "--------------------------------------------------------");

            writer.flush();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }

        // Run GAPIT conversion outside try-with-resources to guarantee file handles are closed
        if (gapit) {
            convertGTtoGAPIT();
        }
        return 0;
    }

    @Override
    public void run() {
        try {
            call();
        } catch (Exception e) {
            Logger.error(CLASS_NAME, "Error during task execution: " + e.getMessage());
            e.printStackTrace();
        }
    }

    private Set<String> loadGenotypeSet(Path path) throws IOException {
        if (path != null && Files.exists(path)) {
            return new HashSet<>(Arrays.asList(FileUtils.fileToArray(path)));
        }
        return null;
    }

    private void validateInputs() {
        if (!Files.exists(panMatrixFile)) {
            Logger.error(CLASS_NAME, "Pan-genome matrix file does not exist: " + panMatrixFile);
        }
        if (outPrefix == null || outPrefix.isEmpty()) {
            Logger.error(CLASS_NAME, "Output prefix is null or empty: " + outPrefix);
        }
    }

    private static class MarkerStats {
        final int[] gtData;
        final int[] alleleCounts = new int[4]; // [missing, hom_ref, het, hom_alt]
        final double maf;
        final double missingProportion;
        final double heterozygosity;
        boolean isGood = true;

        MarkerStats(String[] tokens, int[] targetIndices) {
            this.gtData = new int[targetIndices.length];

            for (int i = 0; i < targetIndices.length; i++) {
                String val = tokens[targetIndices[i]];
                int colonIdx = val.indexOf(':');
                String gtStr = (colonIdx != -1) ? val.substring(0, colonIdx) : val;
                int gt = Integer.parseInt(gtStr);
                this.gtData[i] = gt;

                switch (gt) {
                    case -1 -> alleleCounts[0]++;
                    case 0 -> alleleCounts[1]++;
                    case 1 -> alleleCounts[2]++;
                    case 2 -> alleleCounts[3]++;
                    default -> throw new IllegalArgumentException("Unexpected genotype value: " + gt);
                }
            }

            int total = gtData.length;
            int nonMissingCount = total - alleleCounts[0];

            this.missingProportion = (double) alleleCounts[0] / total;
            this.maf = nonMissingCount == 0 ? 0.0 : Math.min((double) alleleCounts[1] / nonMissingCount, (double) alleleCounts[3] / nonMissingCount);
            this.heterozygosity = nonMissingCount == 0 ? 0.0 : (double) alleleCounts[2] / nonMissingCount;
        }

        void setGood(boolean good) {
            this.isGood = good;
        }

        boolean isGood() {
            return isGood;
        }

        double maf() {
            return maf;
        }

        double missingProportion() {
            return missingProportion;
        }

        double heterozygosity() {
            return heterozygosity;
        }
    }

    private void convertGTtoGAPIT() {
        Logger.info(CLASS_NAME, "Converting genotype matrix to GAPIT format via temporary RScript...");

        File tempScript = null;
        try {
            tempScript = File.createTempFile("convertGTmatrixToRdata_", ".R");

            try (BufferedWriter writer = new BufferedWriter(new FileWriter(tempScript))) {
                writer.write(String.format(
                        "gt <- read.table('%s.genotype.tsv', header = TRUE, check.names = FALSE)\n" +
                                "\n" +
                                "# 1. Extract Genotype Map (GM)\n" +
                                "GM <- data.frame(\n" +
                                "  Name = gt$ID,\n" +
                                "  Chromosome = gt$CHR,\n" +
                                "  Position = gt$START,\n" +
                                "  stringsAsFactors = FALSE\n" +
                                ")\n" +
                                "\n" +
                                "# 2. Extract Genotype Data (GD)\n" +
                                "gt_matrix <- as.matrix(gt[, 5:ncol(gt)])\n" +
                                "rownames(gt_matrix) <- gt$ID\n" +
                                "GD_transposed <- t(gt_matrix)\n" +
                                "\n" +
                                "GD <- data.frame(\n" +
                                "  Taxa = rownames(GD_transposed),\n" +
                                "  GD_transposed,\n" +
                                "  row.names = NULL,\n" +
                                "  check.names = FALSE,\n" +
                                "  stringsAsFactors = FALSE\n" +
                                ")\n" +
                                "\n" +
                                "# 3. Save output files\n" +
                                "save(GD, GM, file = '%s.GAPIT.RData')\n" +
                                "# write.table(GM, '%s.map.tsv', sep = '\\t', quote = FALSE, row.names = FALSE)\n" +
                                "# write.table(GD, '%s.gt.tsv', sep = '\\t', quote = FALSE, row.names = FALSE)\n",
                        outPrefix, outPrefix, outPrefix, outPrefix
                ));
            }

            ProcessBuilder pb = new ProcessBuilder("Rscript", tempScript.getAbsolutePath());
            pb.redirectErrorStream(true);
            Process proc = pb.start();

            try (BufferedReader rReader = new BufferedReader(new InputStreamReader(proc.getInputStream()))) {
                String rLine;
                while ((rLine = rReader.readLine()) != null) {
                    Logger.info(CLASS_NAME, "[Rscript] " + rLine);
                }
            }

            int exitCode = proc.waitFor();
            if (exitCode == 0) {
                Logger.info(CLASS_NAME, "Successfully generated GAPIT file with GD and GM objects (.GAPIT.RData)");
            } else {
                Logger.error(CLASS_NAME, "Rscript exited with non-zero code: " + exitCode);
            }

        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Failed to create or write temporary R script: " + e.getMessage());
        } catch (InterruptedException e) {
            Logger.error(CLASS_NAME, "Rscript execution was interrupted: " + e.getMessage());
            Thread.currentThread().interrupt();
        } finally {
            if (tempScript != null && tempScript.exists()) {
                tempScript.delete();
            }
        }
    }
}