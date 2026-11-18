package edu.tamu.pg.pang2p.Tasks;

import edu.tamu.pg.pang2p.Data.*;
import edu.tamu.pg.pang2p.Utils.*;
import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;
import org.jgrapht.alg.connectivity.ConnectivityInspector;
import picocli.CommandLine.*;

import java.nio.file.*;
import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;


@Command(name = "pan-clusters", description = "Parse mmseqs2 clustering output to generate members and representative files.")
public class PanCluster implements Runnable, Callable<Integer> {
    @Option(names = {"-s", "--single-copy-fasta"}, description = "Path to the all_single_copy_regions.fasta file", required = true)
    private Path scrFile;

    @Option(names = {"-c", "--clusters"}, description = "Path to the mmseqs2 clusters.tsv file", required = true)
    private Path clustersFile;

    @Option(names = {"-g", "--genome-map"}, description = "Path to the genome map file", required = true)
    private Path genomeMapFile;

    @Option(names = {"-o", "--outprefix"}, description = "Output prefix for the generated files", required = true)
    private Path outPrefix;

    @Option(names = {"-k", "--kmer"}, description = "K-mer size for Pan Cluster KMC", required = false, defaultValue = "51")
    private int kmerSize;

    @Option(names = {"-t", "--threads"}, description = "Number of threads to use for processing", required = false, defaultValue = "2")
    private int nThreads;

    @Option(names = {"--memory"}, description = "Use memory mode for processing", required = false, defaultValue = "false")
    private boolean memory = false;

    @Option(names = {"--force"}, description = "Force re-execution of the step", required = false, defaultValue = "false")
    private boolean isForce = false;

    private final String CLASS_NAME = PanCluster.class.getSimpleName();
    private final AtomicInteger lastLoggedPercent = new AtomicInteger(-1);

    public PanCluster() {}

    public PanCluster(Path scrFile, Path clustersFile, Path genomeMapFile, Path outPrefix, int kmerSize, int nThreads, boolean isForce) {
        this.scrFile = scrFile;
        this.clustersFile = clustersFile;
        this.genomeMapFile = genomeMapFile;
        this.outPrefix = outPrefix;
        this.kmerSize = kmerSize;
        this.nThreads = nThreads;
        this.isForce = isForce;
    }

    public void run() {
        try {
            call();
        } catch (Exception e) {
            Logger.error(CLASS_NAME, "Error during PanCluster execution: " + e.getMessage());
            e.printStackTrace();
        }
    }

    @Override
    public Integer call() throws Exception {
        GenomeMap genomeMap = new GenomeMap(genomeMapFile);
        if (FileUtils.isDir(outPrefix)) {
            outPrefix = outPrefix.resolve("pan_clusters");
        }
        Path panClustersFasta = Path.of(outPrefix + ".fasta");
        Path repOut = Path.of(outPrefix + "_rep.tsv");
        Path membOut = Path.of(outPrefix + "_members.tsv");
        Path bedOut = Path.of(outPrefix + ".bed");
        Path kmcDB = Path.of(outPrefix + ".k" + kmerSize);

        if (!isForce && Files.exists(panClustersFasta) && Files.exists(repOut)
                && Files.exists(membOut) && Files.exists(bedOut)) {
            Logger.info(CLASS_NAME, "Skipping PanCluster step: output files already exist. Use --force to re-run.");
            return 0;
        }

        Logger.info(CLASS_NAME, "Running PanCluster step.");

        // 1. Build Graph from clustering output
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

        // 2. Parse sequence metadata
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
            if (desc != null && desc.contains(":") && desc.contains("-")) {
                try {
                    String[] chromAndRange = desc.trim().split(":");
                    if (chromAndRange.length == 2) {
                        String[] coordinates = chromAndRange[1].split("-");
                        if (coordinates.length == 2) {
                            start = Integer.parseInt(coordinates[0]);
                            end = Integer.parseInt(coordinates[1]);
                        }
                    }
                } catch (Exception e) {
                    start = end = -1;
                }
            }
            seqMeta.setStart(start);
            seqMeta.setEnd(end);
            seqMetaMap.put(seqName, seqMeta);
        }

        // 3. Identify connected components and determine representative IDs
        ConnectivityInspector<String, DefaultEdge> inspector = new ConnectivityInspector<>(graph);
        List<Set<String>> graphClusters = inspector.connectedSets();

        List<String> repIds = new ArrayList<>();
        List<List<SeqMeta>> validClusterMembersList = new ArrayList<>();

        for (Set<String> clusterSeqIds : graphClusters) {
            List<SeqMeta> memberSeqs = clusterSeqIds.stream()
                    .map(seqMetaMap::get)
                    .filter(Objects::nonNull) // Prevent NullPointerException if ID is missing from fasta
                    .collect(Collectors.toList());

            if (memberSeqs.isEmpty()) {
                continue;
            }

            memberSeqs.sort(Comparator.comparingInt(SeqMeta::getLength).reversed());
            repIds.add(memberSeqs.get(0).getId());
            validClusterMembersList.add(memberSeqs);
        }

        // 4. Write representative sequences to FASTA
        Set<String> repIdSet = new HashSet<>(repIds);

        try (FastaWriter fastaWriter = new FastaWriter(panClustersFasta);
             FastaReader fastaReader = new FastaReader(scrFile)) {
            for (Fasta in : fastaReader) {
                if (repIdSet.contains(in.getName())) {
                    fastaWriter.write(in);
                }
            }
        }

        // 5. Run KMC on representative FASTA
        KMCConfig kmcConfig = KMCConfig.builder()
                .inputPaths(panClustersFasta)
                .outputPrefix(kmcDB)
                .kmerSize(kmerSize)
                .threads(nThreads)
                .minCount(1)
                .build();
        KMCExecutor.run(kmcConfig);

        long totalKmersSum = 0;
        Map<String, Long> seqTotalKmersMap = new ConcurrentHashMap<>();
        Map<String, Long> seqUniqueKmersMap = new ConcurrentHashMap<>();
        Map<String, Long> seqGlobalUniqueKmersMap = new ConcurrentHashMap<>();

        // 6. Multithreaded processing with progress tracking
        ExecutorService executor = Executors.newFixedThreadPool(nThreads);
        ExecutorCompletionService<Void> completionService = new ExecutorCompletionService<>(executor);

        List<Fasta> fastasToProcess = new ArrayList<>();
        try (FastaReader reader = new FastaReader(panClustersFasta)) {
            for (Fasta fasta : reader) {
                fastasToProcess.add(fasta);
            }
        }

        int totalSequences = fastasToProcess.size();
        AtomicInteger completedSequences = new AtomicInteger(0);
        lastLoggedPercent.set(-1);

        try (KMC kmc = new KMC(kmcDB, memory)) {
            int prefixLength = kmc.getPrefixLength();

            for (Fasta fasta : fastasToProcess) {
                completionService.submit(() -> {
                    KmerIterator it = fasta.kmers(kmerSize, prefixLength, true);

                    Set<Kmer> sequenceUniqueKmers = new HashSet<>();
                    long sequenceTotalKmers = 0;
                    long sequenceGlobalUniqueKmers = 0;

                    while (it.hasNext()) {
                        Kmer kmer = it.next().getKmer();
                        sequenceTotalKmers++;
                        sequenceUniqueKmers.add(kmer);

                        if (kmc.getCount(kmer) == 1) {
                            sequenceGlobalUniqueKmers++;
                        }
                    }

                    seqTotalKmersMap.put(fasta.getName(), sequenceTotalKmers);
                    seqUniqueKmersMap.put(fasta.getName(), (long) sequenceUniqueKmers.size());
                    seqGlobalUniqueKmersMap.put(fasta.getName(), sequenceGlobalUniqueKmers);

                    // Track progress updates
                    int completed = completedSequences.incrementAndGet();
                    int percent = (int) ((completed * 100L) / totalSequences);

                    if (percent % 5 == 0) {
                        synchronized (lastLoggedPercent) {
                            if (percent > lastLoggedPercent.get()) {
                                lastLoggedPercent.set(percent);
                                int dotsCount = percent / 5;
                                String dots = ".".repeat(dotsCount);
                                String progressMsg = String.format("Processing K-mers %s%d%%", dots, percent);
                                Logger.progress(CLASS_NAME, progressMsg);
                            }
                        }
                    }
                    return null;
                });
            }


            for (int i = 0; i < totalSequences; i++) {
                try {
                    completionService.take().get();
                } catch (Exception e) {
                    Logger.error(CLASS_NAME, "Error during parallel K-mer processing task: " + e.getMessage());
                }
            }
        } catch (Exception e) {
            throw new RuntimeException("Error during parallel K-mer processing", e);
        } finally {
            executor.shutdownNow();
        }

        for (long count : seqTotalKmersMap.values()) {
            totalKmersSum += count;
        }

        Logger.info(CLASS_NAME, "Summary -> Total K-mers processed: " + totalKmersSum);

        // 7. Instantiate final Cluster objects using collected k-mer maps
        List<Cluster> clusters = new ArrayList<>();
        int idx = 0;
        for (List<SeqMeta> memberSeqs : validClusterMembersList) {
            String repId = memberSeqs.get(0).getId();

            long repTotalKmers = seqTotalKmersMap.getOrDefault(repId, 0L);
            long repUniqueKmers = seqUniqueKmersMap.getOrDefault(repId, 0L);
            long repGlobalUniqueKmers = seqGlobalUniqueKmersMap.getOrDefault(repId, 0L);

            clusters.add(new Cluster("clstr" + idx, memberSeqs, repTotalKmers, repUniqueKmers, repGlobalUniqueKmers));
            idx++;
        }

        // 8. Write outputs (rep.tsv, members.tsv, bed)
        try (
                BufferedWriter repBw = Files.newBufferedWriter(repOut);
                BufferedWriter membBw = Files.newBufferedWriter(membOut);
                BufferedWriter bedBw = Files.newBufferedWriter(bedOut)
        ) {
            repBw.write(String.join(Defaults.FIELD_SEPERATOR,
                    "clusterId", "repId", "genome", "chrom", "orgGenome", "orgChrom",
                    "length", "start", "end", "totalKmers", "uniqueKmers", "globalUniqueKmers", "size", "ngenotypes", "minLength",
                    "maxLength", "meanLength", "sdLength", "cvLength", "mainChrom",
                    "mainChromFraction", "nChroms"));
            repBw.newLine();

            membBw.write(String.join(Defaults.FIELD_SEPERATOR,
                    "clusterId", "seqId", "genome", "chrom", "orgGenome", "orgChrom",
                    "length", "start", "end"));
            membBw.newLine();

            for (Cluster c : clusters) {
                repBw.write(c.repSummaryTSV());
                repBw.newLine();

                for (String line : c.membersSummaryTSV()) {
                    membBw.write(line);
                    membBw.newLine();
                }

                bedBw.write(c.repBEDLine());
                bedBw.newLine();
            }
        }

        SortBED sorter = new SortBED(bedOut);
        sorter.sort(nThreads);

        Logger.info(CLASS_NAME, "PanCluster step completed successfully.");
        return 0;
    }

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

        public void setLength(int length)   { this.length = length; }
        public void setStart(int start)     { this.start = start; }
        public void setEnd(int end)         { this.end = end; }
        public void setDirection(boolean isForward) { this.isForward = isForward; }

        @Override
        public int compareTo(SeqMeta o) {
            return Integer.compare(this.length, o.length);
        }
    }

    public class Cluster {
        private final String clusterId;
        private final List<SeqMeta> members;
        private final SeqMeta rep;
        private final int size;
        private final int ngenotypes;
        private final int minLength;
        private final int maxLength;
        private final double meanLength;
        private final double sdLength;
        private final double cvLength;
        private final String mainChrom;
        private final double mainChromFraction;
        private final int nChroms;
        private final Map<String, List<SeqMeta>> membersByGenotype;
        private final long totalKmers;
        private final long uniqueKmers;
        private final long globalUniqueKmers;

        public Cluster(String clusterId, List<SeqMeta> members, long totalKmers, long uniqueKmers, long globalUniqueKmers) {
            this.clusterId = clusterId;
            if (members == null || members.isEmpty())
                throw new IllegalArgumentException("Cluster must have at least one member");
            this.members = new ArrayList<>(members);
            this.members.sort(Comparator.comparingInt(SeqMeta::getLength).reversed());
            this.rep = this.members.get(0);
            this.size = this.members.size();
            this.totalKmers = totalKmers;
            this.uniqueKmers = uniqueKmers;
            this.globalUniqueKmers = globalUniqueKmers;

            Set<String> genotypes = this.members.stream()
                    .map(SeqMeta::getGenome)
                    .collect(Collectors.toSet());
            this.ngenotypes = genotypes.size();

            List<Integer> lengths = this.members.stream()
                    .map(SeqMeta::getLength)
                    .collect(Collectors.toList());
            this.minLength = Collections.min(lengths);
            this.maxLength = Collections.max(lengths);
            this.meanLength = lengths.stream().mapToInt(i -> i).average().orElse(0);
            double mean = this.meanLength;
            this.sdLength = Math.sqrt(
                    lengths.stream()
                            .mapToDouble(l -> Math.pow(l - mean, 2))
                            .average().orElse(0)
            );
            this.cvLength = mean == 0 ? 0 : this.sdLength / mean;

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

            Map<String, List<SeqMeta>> byGeno = new HashMap<>();
            for (SeqMeta seq : this.members) {
                String geno = seq.getGenome();
                byGeno.computeIfAbsent(geno, k -> new ArrayList<>()).add(seq);
            }
            this.membersByGenotype = Collections.unmodifiableMap(byGeno);
        }

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
                    Long.toString(totalKmers),
                    Long.toString(uniqueKmers),
                    Long.toString(globalUniqueKmers),
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
            return String.join("\t",
                    r.getChromCode(),
                    Integer.toString(Math.max(0, r.getStart() - 1)),
                    Integer.toString(r.getEnd()),
                    r.getId()
            );
        }

        public String getClusterId() { return clusterId; }
        public List<SeqMeta> getMembers() { return Collections.unmodifiableList(members); }
        public SeqMeta getRep() { return rep; }
        public int getSize() { return size; }
        public int getNgenotypes() { return ngenotypes; }
        public int getMinLength() { return minLength; }
        public int getMaxLength() { return maxLength; }
        public double getMeanLength() { return meanLength; }
        public double getSdLength() { return sdLength; }
        public double getCvLength() { return cvLength; }
        public String getMainChrom() { return mainChrom; }
        public double getMainChromFraction() { return mainChromFraction; }
        public int getNChroms() { return nChroms; }
        public long getTotalKmers() { return totalKmers; }
        public long getUniqueKmers() { return uniqueKmers; }
        public long getGlobalUniqueKmers() { return globalUniqueKmers; }
    }
}