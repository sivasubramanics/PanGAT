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
    private final Map<String, GenomeInfo> newIdMap = new HashMap<>();
    private static final String CLASS_NAME = GenomeMap.class.getSimpleName();
    // Add other maps if you want lookup by oldGenomeName, etc.

    public GenomeMap(DirectoryStructure dirs, Path genomeFile) throws IOException {
        this.dirs = dirs;
        loadFromFile(genomeFile);
        Logger.info(CLASS_NAME, "Loaded genome map from file: " + genomeFile);
        Logger.info(CLASS_NAME, "Total entries loaded: " + newIdMap.size());
    }

    private void loadFromFile(Path genomeFile) throws IOException {
        try (BufferedReader reader = Files.newBufferedReader(genomeFile)) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.trim().isEmpty() || line.startsWith("#")) continue;
                String[] tokens = line.split(",");
                if (tokens.length < 3) continue;
                String newId = tokens[0].trim();
                String oldGenomeName = tokens[1].trim();
                String oldChromName = tokens[2].trim();
                GenomeInfo info = new GenomeInfo(newId, oldGenomeName, oldChromName);
                newIdMap.put(newId, info);
            }
        }
    }

    public GenomeInfo getByNewId(String newId) {
        return newIdMap.get(newId);
    }

    // --- Static inner class holding the info ---
    public static class GenomeInfo {
        private final String newId;         // e.g. aaabbb
        private final String oldGenomeName; // e.g. Lakhan
        private final String oldChromName;  // e.g. Chr01

        public GenomeInfo(String newId, String oldGenomeName, String oldChromName) {
            this.newId = newId;
            this.oldGenomeName = oldGenomeName;
            this.oldChromName = oldChromName;
        }

        public String getNewId() { return newId; }
        public String getGenomeCode() { return newId.substring(0, Defaults.CHROM_CODE_LENGTH); }
        public String getChromCode() { return newId.substring(Defaults.CHROM_CODE_LENGTH, Defaults.CHROM_CODE_LENGTH + Defaults.CHROM_CODE_LENGTH); }
        public String getOldGenomeName() { return oldGenomeName; }
        public String getOldChromName() { return oldChromName; }
    }
}