package edu.tamu.pg.pang2p.Data;

import edu.tamu.pg.pang2p.Utils.FileUtils;
import edu.tamu.pg.pang2p.Utils.KMCConfig;
import edu.tamu.pg.pang2p.Utils.KMCExecutor;
import edu.tamu.pg.pang2p.Utils.Logger;

import java.io.*;
import java.nio.file.Path;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

public class PanHap {
    private final String CLASS_NAME = PanHap.class.getSimpleName();

    ConcurrentHashMap<String, HapData> hapData = new ConcurrentHashMap<>();
    String[] samples = new String[0];
    double[] sampleSums = new double[0];

    private String version = "0.0.1";
    private final List<String> commandLines = new ArrayList<>();
    private final Map<String, String> parameters = new LinkedHashMap<>();

    public PanHap(Path panClustersRepFile) {
        this.samples = new String[0];
        try (BufferedReader reader = FileUtils.reader(panClustersRepFile)) {
            String line;
            int lineNumber = 0;
            while ((line = reader.readLine()) != null) {
                if (line.trim().isEmpty()) continue;
                if (lineNumber == 0) {
                    lineNumber++;
                    continue;
                }
                String[] parts = line.split("\t");
                String clusterId = parts[0];
                String hapID = parts[1];
                String newGenomeID = parts[2];
                String newChromID = parts[3];
                String orgGenomeID = parts[4];
                String orgChromID = parts[5];
                int length = Integer.parseInt(parts[6]);
                int start = Integer.parseInt(parts[7]);
                int end = Integer.parseInt(parts[8]);
                int totalKmers = Integer.parseInt(parts[9]);
                int uniqueKmers = Integer.parseInt(parts[10]);
                int globalUniqueKmers = Integer.parseInt(parts[11]);
                hapData.computeIfAbsent(hapID, k -> new HapData(clusterId, hapID, newGenomeID, newChromID, orgGenomeID, orgChromID, start, end, length, totalKmers, uniqueKmers, globalUniqueKmers));
            }
            Logger.info(CLASS_NAME, "Loaded " + hapData.size() + " hap data");
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }


    public PanHap(Path panClustersRepFile, Path mapKmersFile) {
        this(panClustersRepFile);

        try (BufferedReader reader = FileUtils.reader(mapKmersFile)) {
            Logger.info(CLASS_NAME, "Loading pan-matrix from: " + mapKmersFile);
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.trim().isEmpty()) continue;

                if (line.startsWith("##")) {
                    parseMetadataHeader(line);
                    continue;
                }

                if (line.startsWith("hapID") || line.startsWith("# hapID")) {
                    String[] parts = line.replace("# ", "").split("\t");
                    if (parts.length > 15) {
                        samples = Arrays.copyOfRange(parts, 15, parts.length);
                        sampleSums = new double[samples.length];
                    }
                    continue;
                }

                String[] parts = line.split("\t");
                String hapID = parts[0];
                int totalKmers = Integer.parseInt(parts[7]);
                int uniqueKmers = Integer.parseInt(parts[8]);
                int globalUniqueKmers = Integer.parseInt(parts[9]);

                int sampleCount = parts.length - 15;
                int[] genotype = new int[sampleCount];
                double[] gtScaled = new double[sampleCount];
                double[] normalizedCounts = new double[sampleCount];
                int[] obsCounts = new int[sampleCount];
                int[] obsTotalCounts = new int[sampleCount];


                for (int i = 0; i < sampleCount; i++) {
                    String[] attrs = parts[i + 15].split(":");

                    genotype[i] = Integer.parseInt(attrs[0]);
                    gtScaled[i] = Double.parseDouble(attrs[1]);
                    normalizedCounts[i] = Double.parseDouble(attrs[2]);
                    obsCounts[i] = Integer.parseInt(attrs[3]);
                    obsTotalCounts[i] = Integer.parseInt(attrs[4]);

                    sampleSums[i] += obsCounts[i];
                }

                HapData hap = hapData.get(hapID);
                if (hap != null) {
                    hap.setTotKmers(totalKmers);
                    hap.setUniqueKmers(uniqueKmers);
                    hap.setGlobalUniqueKmers(globalUniqueKmers);
                    hap.setGenotype(genotype);
                    hap.setGtScaled(gtScaled);
                    hap.setNormalizedCounts(normalizedCounts);
                    hap.setObsCounts(obsCounts);
                    hap.setObsTotalCounts(obsTotalCounts);
                }
            }
            Logger.info(CLASS_NAME, "Loaded " + hapData.size() + " hap data");
            Logger.info(CLASS_NAME, "Loaded " + samples.length + " samples");
        } catch (IOException e) {
            throw new RuntimeException("Error reading matrix file: " + mapKmersFile, e);
        }
    }

    // --- Metadata Accessors & Mutators ---

    public String getVersion() {
        return version;
    }

    public void setVersion(String version) {
        this.version = version;
    }

    public List<String> getCommandLines() {
        return commandLines;
    }

    public void addCommandLine(String commandLine) {
        if (commandLine != null && !commandLine.trim().isEmpty()) {
            this.commandLines.add(commandLine);
        }
    }

    public void addParameter(String key, String value) {
        this.parameters.put(key, value);
    }

    public Map<String, String> getParameters() {
        return parameters;
    }

    private void parseMetadataHeader(String line) {
        String content = line.substring(2);
        if (content.startsWith("fileformat=")) {
            // Optional handling for format checks
        } else if (content.startsWith("version=")) {
            this.version = content.split("=", 2)[1];
        } else if (content.startsWith("commandline=")) {
            // Append header commandlines as they are read from existing files
            addCommandLine(content.split("=", 2)[1]);
        } else if (content.startsWith("parameter=")) {
            String[] paramParts = content.substring("parameter=<".length(), content.length() - 1).split(",");
            String key = "", value = "";
            for (String p : paramParts) {
                String[] kv = p.split("=");
                if (kv[0].equals("ID")) key = kv[1];
                if (kv[0].equals("Value")) value = kv[1];
            }
            if (!key.isEmpty()) this.parameters.put(key, value);
        }
    }

    public void addCountFile(String sampleID, Path countFile) {
        Logger.info(CLASS_NAME, "Adding sample: " + sampleID + " from file: " + countFile);

        List<String> sampleNames = (samples != null) ? new ArrayList<>(Arrays.asList(samples)) : new ArrayList<>();

        if (sampleNames.contains(sampleID)) {
            Logger.info(CLASS_NAME, "Warning: Sample " + sampleID + " already exists. Skipping.");
            return;
        }
        sampleNames.add(sampleID);
        samples = sampleNames.toArray(new String[0]);
        sampleSums = Arrays.copyOf(sampleSums, samples.length);

        int sampleIndex = samples.length - 1;

        try (BufferedReader reader = FileUtils.reader(countFile)) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.trim().isEmpty()) continue;
                // skip first line if it starts with HapID or # HapID
                if (line.startsWith("HapID") || line.startsWith("# HapID")) {
                    continue;
                }
                String[] parts = line.split("\t");
                String hapID = parts[0];
                int length = Integer.parseInt(parts[1]);
                int totalKmers = Integer.parseInt(parts[2]);
                int uniqueKmers = Integer.parseInt(parts[3]);
                int globalUniqueKmers = Integer.parseInt(parts[4]);
                int countedKmers = Integer.parseInt(parts[5]);
                int numReads = Integer.parseInt(parts[6]);

                sampleSums[sampleIndex] += numReads;
                HapData hap = hapData.get(hapID);
                if (hap != null) {
                    if (hap.getTotKmers() == -1) {
                        hap.setTotKmers(totalKmers);
                        hap.setUniqueKmers(uniqueKmers);
                        hap.setGlobalUniqueKmers(globalUniqueKmers);
                    } else if (hap.getTotKmers() != totalKmers) {
                        Logger.error(CLASS_NAME, "HapID " + hapID + " total kmers mismatch. Expected: " + hap.getTotKmers() + ", Found: " + totalKmers);
                    }

                    if (hap.getLength() != length) {
                        Logger.error(CLASS_NAME, "HapID " + hapID + " length mismatch. Expected: " + hap.getLength() + ", Found: " + length);
                    }
                    hap.addSample(countedKmers, numReads);
                } else {
                    Logger.error(CLASS_NAME, "HapID " + hapID + " not found in panClustersRepFile.");
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public void normalizePanHap() {
        Logger.info(CLASS_NAME, "Normalizing pan-hap data for " + samples.length + " samples.");
        for (HapData hap : hapData.values()) {
            int[] obsCounts = hap.getObsCounts();
            double[] normalizedCounts = new double[samples.length];
            if (obsCounts != null) {
                for (int i = 0; i < obsCounts.length; i++) {
                    double normalizedCount = (sampleSums[i] > 0) ? (obsCounts[i] / sampleSums[i]) * 1_000_000_000.0 : 0.0;
                    normalizedCounts[i] = normalizedCount;
                }
            }
            hap.setNormalizedCounts(normalizedCounts);
            hap.calcMetrics();
        }
    }

    /**
     * Returns hapData values sorted by newChromID and start.
     */
    public List<HapData> getSortedHapData() {
        this.normalizePanHap();
        List<HapData> sortedList = new ArrayList<>(hapData.values());
        Collections.sort(sortedList);
        return sortedList;
    }

    /**
     * Writes the complete Pan-Genome Genotype and Coverage Matrix with VCF-style metadata headers.
     * Output cell format: genotype:cv:normalizedCount:obsCount:obsTotalCount
     *
     * Recommended output extension: .panmat or .tsv
     */
    public void writePanMatrix(Path panMatrixFile, double cutoff) {
        Logger.info(CLASS_NAME, "Writing pan-genome matrix to: " + panMatrixFile);
        try (BufferedWriter writer = FileUtils.writer(panMatrixFile)) {
            // Write VCF-style metadata headers
            writer.write("##fileformat=\"PanG2P v0.1.0\"");
            writer.newLine();
            writer.write("##fileDate=" + LocalDateTime.now().format(DateTimeFormatter.ofPattern("yyyyMMdd_HHmmss")));
            writer.newLine();
            writer.write("##version=" + version);
            writer.newLine();
            for (String cmd : commandLines) {
                writer.write("##commandline=" + cmd);
                writer.newLine();
            }
            for (Map.Entry<String, String> entry : parameters.entrySet()) {
                writer.write("##parameter=<ID=" + entry.getKey() + ",Value=" + entry.getValue() + ">");
                writer.newLine();
            }
            writer.write("##FORMAT=<ID=GT:GTScaled:NORM:OBS:TOTAL,Number=5,Description=\"Genotype:GenotypeScaled:NormalizedCount:MappedKmers:MappedTotalKmers\">");
            writer.newLine();

            // Column Headers
            writer.write("# hapID\tGenomeID\tChromID\tSimChromID\tStart\tEnd\tLength\tTotKmers\tUniqueKmers\tGlobalUniqueKmers\tMinRPB\tMaxRPB\tMeanRPB\tSdRPB\tCV");
            for (String sample : samples) {
                writer.write("\t" + sample);
            }
            writer.newLine();

            // Data Rows
            for (HapData hap : getSortedHapData()) {
                writer.write(hap.toString());
                writer.newLine();
            }
            Logger.info(CLASS_NAME, "Writing pan-genome matrix. Hap data written: " + hapData.size() + ", Samples: " + samples.length);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public void writePanMatrix(Path panMatrixFile) {
        writePanMatrix(panMatrixFile, 0.0);
    }

    public int getHapNum() {
        return hapData.size();
    }
}