package edu.tamu.pg.panhap.Tasks;

import edu.tamu.pg.panhap.Data.*;
import edu.tamu.pg.panhap.Utils.*;
import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;
import org.jgrapht.alg.connectivity.ConnectivityInspector;

import java.nio.file.*;
import java.io.*;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

public class PanCluster {
    private final Path clustersFile;
    private final Path singleCopyDir;
    private final Path scrFile;
    private int nThreads;
    private final String CLASS_NAME = PanCluster.class.getSimpleName();

    // Accept DirectoryStructure and retrieve file paths internally
    public PanCluster(DirectoryStructure dirs) {
        this.singleCopyDir = dirs.getDirectory(PipelineSteps.SINGLE_COPY);
        this.clustersFile = singleCopyDir.resolve("clusters.tsv");
        this.scrFile = dirs.getSingleCopyReference();
        this.nThreads = 2;
    }

    public PanCluster(DirectoryStructure dirs, int nThreads) {
        this(dirs);
        this.nThreads = nThreads;
    }

    public void run(GenomeMap genomeMap) throws IOException {

        Graph<String, DefaultEdge> graph = new SimpleGraph<>(DefaultEdge.class);
        try (BufferedReader br = new BufferedReader(new FileReader(clustersFile.toFile()))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] fields = line.split("\t");
                graph.addVertex(fields[0]);
                graph.addVertex(fields[1]);
                if (!fields[0].equals(fields[1])) {
                    graph.addEdge(fields[0], fields[1]);
                }
            }
        }

        Map<String, SeqMeta> seqMetaMap = new HashMap<>();
        FastaIndex fastaIndex = new FastaIndex(scrFile);

        for (String seqName : fastaIndex.getSequenceNames()) {
            if (!graph.containsVertex(seqName)) {
                graph.addVertex(seqName);
            }
            FastaIndexEntry entry = fastaIndex.getEntry(seqName);
            SeqMeta seqMeta = new SeqMeta(seqName, genomeMap);
            seqMeta.setLength(entry.getLength());

            String desc = entry.getDescription();
            int start = -1, end = -1;
            if (desc != null && desc.contains("-")) {
                String[] parts = desc.trim().split("-");
                try {
                    start = Integer.parseInt(parts[0]);
                    end = Integer.parseInt(parts[1]);
                } catch (Exception e) {
                    start = end = -1;
                }
            }
            seqMeta.setStart(start);
            seqMeta.setEnd(end);

            seqMetaMap.put(seqName, seqMeta);
        }

        ConnectivityInspector<String, DefaultEdge> inspector = new ConnectivityInspector<>(graph);
        List<Set<String>> graphClusters = inspector.connectedSets();

        List<Cluster> clusters = new ArrayList<>();
        int idx = 0;
        for (Set<String> clusterSeqIds : graphClusters) {
            List<SeqMeta> memberSeqs = clusterSeqIds.stream()
                    .map(seqMetaMap::get)
                    .collect(Collectors.toList());
            clusters.add(new Cluster("clstr" + idx, memberSeqs));
            idx++;
        }

        Path repOut = singleCopyDir.resolve("pan_clusters_rep.tsv");
        Path membOut = singleCopyDir.resolve("pan_clusters_members.tsv");
        Path bedOut = singleCopyDir.resolve("pan_clusters.bed");

        try (
                BufferedWriter repBw = Files.newBufferedWriter(repOut);
                BufferedWriter membBw = Files.newBufferedWriter(membOut);
                BufferedWriter bedBw = Files.newBufferedWriter(bedOut)
        ) {
            // 1. Rep header
            repBw.write(String.join(Defaults.FIELD_SEPERATOR,
                    "clusterId",
                    "repId",
                    "genome",
                    "chrom",
                    "orgGenome",
                    "orgChrom",
                    "length",
                    "start",
                    "end",
                    "size",
                    "ngenotypes",
                    "minLength",
                    "maxLength",
                    "meanLength",
                    "sdLength",
                    "cvLength",
                    "mainChrom",
                    "mainChromFrac",
                    "nChroms"));
            repBw.newLine();

            // 2. Members header
            membBw.write(String.join(Defaults.FIELD_SEPERATOR,
                    "clusterId",
                    "seqId",
                    "genome",
                    "chrom",
                    "orgGenome",
                    "orgChrom",
                    "length",
                    "start",
                    "end"));
            membBw.newLine();

            // 3. Loop clusters
            for (Cluster c : clusters) {
                // Write rep
                repBw.write(c.repSummaryTSV());
                repBw.newLine();

                // Write members
                for (String line : c.membersSummaryTSV()) {
                    membBw.write(line);
                    membBw.newLine();
                }

                // Write BED for rep
                bedBw.write(c.repBEDLine());
                bedBw.newLine();
            }
        }

        SortBED sorter = new SortBED(bedOut);
        sorter.sort(nThreads);

        // Read BED sequence IDs (preserve order, and set for quick lookup)
        List<String> panScrIDs = new ArrayList<>();
        Set<String> panScrSet = new HashSet<>();
        try (BufferedReader br = new BufferedReader(new FileReader(bedOut.toFile()))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] fields = line.split("\t");
                if (fields.length < 4) continue;
                String seqId = fields[3];
                panScrIDs.add(seqId);
                panScrSet.add(seqId);
            }
        }


        Map<String, Fasta> fastaMap = new HashMap<>();
        try (FastaReader scrReader = new FastaReader(scrFile)) {
            for (Fasta scr : scrReader) {
                if (panScrSet.contains(scr.getName())) {
                    fastaMap.put(scr.getName(), scr);
                }
            }
        }

        Path panClustersFasta = singleCopyDir.resolve("pan_clusters.fasta");
        try (FastaWriter fastaWriter = new FastaWriter(panClustersFasta.toString())) {
            for (String seqId : panScrIDs) {
                Fasta scr = fastaMap.get(seqId);
                if (scr != null) {
                    fastaWriter.write(scr);
                } else {
                    Logger.error(CLASS_NAME, "Sequence not found in fasta file: " + seqId);
                }
            }
        }

        Path kmcPrefix = singleCopyDir.resolve("pan_clusters");
        KMCConfig kmcConfig = KMCConfig.builder()
                .inputPath(panClustersFasta)
                .outputPrefix(kmcPrefix)
                .threads(nThreads)
                .kmerSize(51)
                .minCount(1)
                .build();

        Logger.info(CLASS_NAME, "Running KMC on pan_clusters fasta...");
        KMCExecutor.run(kmcConfig);
    }

    /**
     * Contains sequence metadata, parsed from sequence ID and a GenomeMap.
     */
    public static class SeqMeta implements Comparable<SeqMeta> {
        private final String id;
        private final String chrom;
        private final String orgChrom;
        private final String genome;
        private final String orgGenome;
        private int length = -1;
        private int start = -1;
        private int end = -1;
        private boolean isForward = true;

        /**
         * Constructs SeqMeta from sequence ID and genomeMap.
         * ID should have at least 2*CHROM_CODE_LENGTH characters.
         */
        public SeqMeta(String id, GenomeMap genomeMap) {
            if (id == null) throw new IllegalArgumentException("ID cannot be null");
            int chromLen = Defaults.CHROM_CODE_LENGTH;
            if (id.length() < 2 * chromLen) {
                throw new IllegalArgumentException("ID length is too short for expected genome/chrom extraction");
            }
            this.id = id;
            this.genome = id.substring(0, chromLen);
            this.chrom = id.substring(0, 2 * chromLen);

            GenomeMap.GenomeInfo genomeInfo = genomeMap.getByChrom(this.chrom);
            if (genomeInfo == null) {
                this.orgGenome = this.genome;
                this.orgChrom = this.chrom;
            } else {
                this.orgChrom = genomeInfo.getOrgChrom();
                this.orgGenome = genomeInfo.getOrgGenome();
            }
        }

        // --- Getters ---
        public String getId()        { return id; }
        public String getChrom()     { return chrom; }
        public String getChromCode() { return chrom.substring(Defaults.CHROM_CODE_LENGTH); }
        public String getOrgChrom()  { return orgChrom; }
        public String getGenome()    { return genome; }
        public String getOrgGenome() { return orgGenome; }
        public int getLength()       { return length; }
        public int getStart()        { return start; }
        public int getEnd()          { return end; }
        public boolean isForward()   { return isForward; }

        // --- Setters ---
        public void setLength(int length)   { this.length = length; }
        public void setStart(int start)     { this.start = start; }
        public void setEnd(int end)         { this.end = end; }
        public void setDirection(boolean isForward) { this.isForward = isForward; }

        public String toTSV(boolean isStrand) {
            if (isStrand) {
                return String.format("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s",
                        id, genome, chrom, orgGenome, orgChrom, length, start, end, isForward
                );
            } else {
                return String.format("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d",
                        id, genome, chrom, orgGenome, orgChrom, length, start, end
                );
            }
        }

        @Override
        public String toString() {
            return String.format(
                    "SeqMeta{id='%s', genome='%s', chrom='%s', orgGenome='%s', orgChrom='%s', length=%d, start=%d, end=%d, forward=%s}",
                    id, genome, chrom, orgGenome, orgChrom, length, start, end, isForward
            );
        }

        @Override
        public int compareTo(SeqMeta o) {
            // sort based on length
            return Integer.compare(this.length, o.length);
        }
    }

    public class Cluster {
        // Fields
        private final String clusterId;
        private final List<SeqMeta> members;  // All member sequences (sorted by length descending)
        private final SeqMeta rep;            // Longest sequence
        private final int size;
        private final int ngenotypes;
        private final int minLength;
        private final int maxLength;
        private final double meanLength;
        private final double sdLength;
        private final double cvLength;
        private final String mainChrom;       // Most common chromosome (excluding "chrUn")
        private final double mainChromFraction;
        private final int nChroms;            // Number of unique chromosomes (excluding "chrUn")
        private final Map<String, List<SeqMeta>> membersByGenotype;  // genotype -> list of SeqMeta

        // Constructor
        public Cluster(String clusterId, List<SeqMeta> members) {
            this.clusterId = clusterId;
            if (members == null || members.isEmpty())
                throw new IllegalArgumentException("Cluster must have at least one member");
            // Defensive copy and sort by descending length
            this.members = new ArrayList<>(members);
            this.members.sort(Comparator.comparingInt(SeqMeta::getLength).reversed());
            this.rep = this.members.get(0);
            this.size = this.members.size();

            // Unique genotypes
            Set<String> genotypes = this.members.stream()
                    .map(SeqMeta::getGenome)
                    .collect(Collectors.toSet());
            this.ngenotypes = genotypes.size();

            // Length statistics
            List<Integer> lengths = this.members.stream()
                    .map(SeqMeta::getLength)
                    .collect(Collectors.toList());
            this.minLength = Collections.min(lengths);
            this.maxLength = Collections.max(lengths);
            this.meanLength = lengths.stream().mapToInt(i -> i).average().orElse(0);
            // Calculate standard deviation
            double mean = this.meanLength;
            this.sdLength = Math.sqrt(
                    lengths.stream()
                            .mapToDouble(l -> Math.pow(l - mean, 2))
                            .average().orElse(0)
            );
            this.cvLength = mean == 0 ? 0 : this.sdLength / mean;

            // Chromosome stats (ignore "chrUn")
            Map<String, Integer> chromCounts = new HashMap<>();
            int chromValidTotal = 0;
            for (SeqMeta seq : this.members) {
                String chrom = seq.getChrom();
                if (chrom != null && !chrom.equals("chrUn")) {
                    chromCounts.put(chrom, chromCounts.getOrDefault(chrom, 0) + 1);
                    chromValidTotal++;
                }
            }
            String majorChrom = "NA";
            int maxCount = 0;
            for (Map.Entry<String, Integer> entry : chromCounts.entrySet()) {
                if (entry.getValue() > maxCount) {
                    maxCount = entry.getValue();
                    majorChrom = entry.getKey();
                }
            }
            this.mainChrom = majorChrom;
            this.nChroms = chromCounts.size();
            this.mainChromFraction = chromValidTotal == 0 ? 0 : ((double) maxCount) / chromValidTotal;

            // Per-genotype lists
            Map<String, List<SeqMeta>> byGeno = new HashMap<>();
            for (SeqMeta seq : this.members) {
                String geno = seq.getGenome();
                byGeno.computeIfAbsent(geno, k -> new ArrayList<>()).add(seq);
            }
            this.membersByGenotype = Collections.unmodifiableMap(byGeno);
        }

        // --- Getters ---
        public String getClusterId()               { return clusterId; }
        public List<SeqMeta> getMembers()          { return Collections.unmodifiableList(members); }
        public SeqMeta getRep()                    { return rep; }
        public int getSize()                       { return size; }
        public int getNgenotypes()                 { return ngenotypes; }
        public int getMinLength()                  { return minLength; }
        public int getMaxLength()                  { return maxLength; }
        public double getMeanLength()              { return meanLength; }
        public double getSdLength()                { return sdLength; }
        public double getCvLength()                { return cvLength; }
        public String getMainChrom()               { return mainChrom; }
        public double getMainChromFraction()       { return mainChromFraction; }
        public int getNChroms()                    { return nChroms; }
        public Map<String, List<SeqMeta>> getMembersByGenotype() { return membersByGenotype; }

        public String repSummaryTSV() {
            SeqMeta r = getRep();
            return String.join("\t",
                    clusterId,
                    r.getId(),
                    r.getGenome(),
                    r.getChrom(),
                    r.getOrgGenome(),
                    r.getOrgChrom(),
                    Integer.toString(r.getLength()),
                    Integer.toString(r.getStart()),
                    Integer.toString(r.getEnd()),
                    Integer.toString(getSize()),
                    Integer.toString(getNgenotypes()),
                    Integer.toString(getMinLength()),
                    Integer.toString(getMaxLength()),
                    HelperFunctions.floatFormat(getMeanLength()),
                    HelperFunctions.floatFormat(getSdLength()),
                    HelperFunctions.floatFormat(getCvLength()),
                    getMainChrom(),
                    HelperFunctions.floatFormat(getMainChromFraction()),
                    Integer.toString(getNChroms())
            );
        }

        public List<String> membersSummaryTSV() {
            List<String> lines = new ArrayList<>();
            for (SeqMeta m : getMembers()) {
                lines.add(String.join("\t",
                        clusterId,
                        m.getId(),
                        m.getGenome(),
                        m.getChrom(),
                        m.getOrgGenome(),
                        m.getOrgChrom(),
                        Integer.toString(m.getLength()),
                        Integer.toString(m.getStart()),
                        Integer.toString(m.getEnd())
                ));
            }
            return lines;
        }

        public String repBEDLine() {
            SeqMeta r = getRep();
            // BED: chrom, start-1 (0-based BED), end, name
            return String.join("\t",
                    r.getChromCode(),
                    Integer.toString(Math.max(0, r.getStart() - 1)),   // start in BED is 0-based
                    Integer.toString(r.getEnd()),
                    r.getId()
            );
        }
    }
}