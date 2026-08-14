package edu.tamu.pg.pang2p.Data;

import edu.tamu.pg.pang2p.Tasks.DirectoryStructure;
import edu.tamu.pg.pang2p.Utils.Defaults;
import edu.tamu.pg.pang2p.Utils.Logger;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

/**
 * GenomeMap provides lookup for genome and chromosome codes and their original names.
 * Supports lookup by chrom, genome, orgChrom, orgGenome, etc.
 */
public class GenomeMap {
    private final Path genomeFile;

    // Lookup maps
    private final Map<String, GenomeInfo> chromToInfo = new HashMap<>(); // chrom => info
    private final Map<String, List<GenomeInfo>> orgGenomeToInfos = new LinkedHashMap<>(); // orgGenome => info
    private final Map<String, GenomeInfo> orgChromToInfo = new HashMap<>(); // orgChrom => info
    private final Set<String> orgGenomes = new LinkedHashSet<>(); // keeps insertion order

    public GenomeMap(DirectoryStructure dirs) throws IOException {
        this.genomeFile = dirs.getBaseDir().resolve("genomesMap.csv");
        loadFromFile();
        Logger.info(getClass().getSimpleName(), "Loaded genome map from: " + genomeFile);
        Logger.info(getClass().getSimpleName(), "Total genome entries: " + chromToInfo.size());
    }

    private void loadFromFile() throws IOException {
        try (BufferedReader reader = Files.newBufferedReader(genomeFile)) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) continue;
                String[] tokens = line.split(",");
                if (tokens.length < 3) continue;
                String chrom = tokens[0].trim();    // aabb
                String orgGenome = tokens[1].trim(); // Lakhan
                String orgChrom = tokens[2].trim();   // Chr01
                String genome = chrom.substring(0, Defaults.CHROM_CODE_LENGTH);

                GenomeInfo info = new GenomeInfo(chrom, orgChrom, genome, orgGenome);

                chromToInfo.put(chrom, info);
                orgGenomes.add(orgGenome);

                orgGenomeToInfos.computeIfAbsent(orgGenome, k -> new ArrayList<>()).add(info);

                orgChromToInfo.put(orgChrom, info); // Lookup by orgChrom name
            }
        }
    }

    /** Lookup by new chrom code (e.g. aabb) */
    public GenomeInfo getByChrom(String chrom) {
        return chromToInfo.get(chrom);
    }

    /** Lookup all GenomeInfo entries for an orgGenome name (e.g. Lakhan) */
    public List<GenomeInfo> getByOrgGenome(String orgGenome) {
        return orgGenomeToInfos.getOrDefault(orgGenome, Collections.emptyList());
    }

    /** Lookup GenomeInfo by original chromosome name */
    public GenomeInfo getByOrgChrom(String orgChrom) {
        return orgChromToInfo.get(orgChrom);
    }

    /** Get all orgGenome names (ordered) */
    public List<String> getOrgGenomes() {
        return new ArrayList<>(orgGenomes);
    }

    /** Get all GenomeInfo objects */
    public List<GenomeInfo> getAllEntries() {
        return new ArrayList<>(chromToInfo.values());
    }

    /**
     * Returns the list of new chrom codes (e.g. aabb, aabc) for a genome code (e.g. aa) or orgGenome (e.g. Lakhan)
     */
    public List<String> getChroms(String genomeIdentifier) {
        // 1. Lookup by orgGenome
        List<GenomeInfo> infos = orgGenomeToInfos.getOrDefault(genomeIdentifier, Collections.emptyList());
        if (!infos.isEmpty()) {
            List<String> chroms = new ArrayList<>(infos.size());
            for (GenomeInfo info : infos) chroms.add(info.getChrom());
            return chroms;
        }
        // 2. Lookup by new genome code (e.g. aa)
        List<String> chroms = new ArrayList<>();
        for (GenomeInfo info : chromToInfo.values()) {
            if (info.getGenome().equals(genomeIdentifier)) {
                chroms.add(info.getChrom());
            }
        }
        return chroms;
    }

    /**
     * Returns the list of orgChroms (e.g. Chr01, Chr02, ...) for a genome code (e.g. aa) or orgGenome (e.g. Lakhan)
     */
    public List<String> getOrgChroms(String genomeIdentifier) {
        // 1. Lookup by orgGenome
        List<GenomeInfo> infos = orgGenomeToInfos.getOrDefault(genomeIdentifier, Collections.emptyList());
        if (!infos.isEmpty()) {
            List<String> orgChroms = new ArrayList<>(infos.size());
            for (GenomeInfo info : infos) orgChroms.add(info.getOrgChrom());
            return orgChroms;
        }
        // 2. Lookup by new genome code (e.g. aa)
        List<String> orgChroms = new ArrayList<>();
        for (GenomeInfo info : chromToInfo.values()) {
            if (info.getGenome().equals(genomeIdentifier)) {
                orgChroms.add(info.getOrgChrom());
            }
        }
        return orgChroms;
    }

    // --- Static inner data class ---
    public static class GenomeInfo {
        private final String chrom;       // new chrom code, e.g. aabb
        private final String orgChrom;    // original chrom name, e.g. Chr01
        private final String genome;      // new genome code, e.g. aa
        private final String orgGenome;   // original genome name, e.g. Lakhan

        public GenomeInfo(String chrom, String orgChrom, String genome, String orgGenome) {
            this.chrom = chrom;
            this.orgChrom = orgChrom;
            this.genome = genome;
            this.orgGenome = orgGenome;
        }

        public String getChrom()      { return chrom; }
        public String getOrgChrom()   { return orgChrom; }
        public String getGenome()     { return genome; }
        public String getOrgGenome()  { return orgGenome; }

        @Override
        public String toString() {
            return String.format("GenomeInfo{chrom='%s', orgChrom='%s', genome='%s', orgGenome='%s'}",
                    chrom, orgChrom, genome, orgGenome);
        }
    }
}