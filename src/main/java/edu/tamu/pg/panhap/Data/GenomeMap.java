package edu.tamu.pg.panhap.Data;

import edu.tamu.pg.panhap.Tasks.DirectoryStructure;
import edu.tamu.pg.panhap.Utils.Defaults;
import edu.tamu.pg.panhap.Utils.Logger;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

/**
 * The GenomeMap class encapsulates loading genome and chromosome mappings from a file,
 * providing easy access by new id, and keeping the GenomeInfo data structure as a static inner class.
 */
public class GenomeMap {
    private final DirectoryStructure dirs;
    private final Map<String, GenomeInfo> genomeInfoMap = new HashMap<>();
    private static final String CLASS_NAME = GenomeMap.class.getSimpleName();
    private final Path genomeFile;
    private final Set<String> genomes = new LinkedHashSet<>();
    private final Map<String, List<String>> genomeToChromosomes = new LinkedHashMap<>();
    // Add other maps if you want lookup by oldGenomeID, etc.

    public GenomeMap(DirectoryStructure dirs) throws IOException {
        this.dirs = dirs;
        this.genomeFile = dirs.getBaseDir().resolve("genomesMap.csv");
        loadFromFile(genomeFile);
        Logger.info(CLASS_NAME, "Loaded genome map from file: " + genomeFile);
        Logger.info(CLASS_NAME, "Total entries loaded: " + genomeInfoMap.size());
    }

    private void loadFromFile(Path genomeFile) throws IOException {
        try (BufferedReader reader = Files.newBufferedReader(genomeFile)) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.trim().isEmpty() || line.startsWith("#")) continue;

                String[] tokens = line.split(",");
                if (tokens.length < 3) continue;

                String newID = tokens[0].trim();
                String genomeID = tokens[1].trim();
                String oldID = tokens[2].trim();

                GenomeInfo info = new GenomeInfo(newID, genomeID, oldID);
                genomeInfoMap.put(newID, info);

                if (!genomes.contains(genomeID)) {
                    genomes.add(genomeID);
                }

                genomeToChromosomes
                        .computeIfAbsent(genomeID, k -> new ArrayList<>())
                        .add(newID);
            }
        }
    }

    public List<String> getChromosomes(String genomeID) {
        return genomeToChromosomes.getOrDefault(genomeID, Collections.emptyList());
    }

    public Map<String, List<String>> getGenomeToChromosomes() {
        return Collections.unmodifiableMap(genomeToChromosomes);
    }

    public GenomeInfo getInfo(String newID) {
        return genomeInfoMap.get(newID);
    }

    public List<String> getGenomes() {
        return new ArrayList<>(genomes);
    }

    // --- Static inner class holding the info ---
    public static class GenomeInfo {
        private final String newID;         // e.g. aaabbb
        private final String GenomeID; // e.g. Lakhan
        private final String oldID;  // e.g. Chr01

        public GenomeInfo(String newID, String genomeID, String oldID) {
            this.newID = newID;
            this.GenomeID = genomeID;
            this.oldID = oldID;
        }

        public String getNewID() { return newID; }
        public String getGenomeCode() { return newID.substring(0, Defaults.CHROM_CODE_LENGTH); }
        public String getChromCode() { return newID.substring(Defaults.CHROM_CODE_LENGTH, Defaults.CHROM_CODE_LENGTH + Defaults.CHROM_CODE_LENGTH); }
        public String getGenomeID() { return GenomeID; }
        public String getOldID() { return oldID; }
    }
}   