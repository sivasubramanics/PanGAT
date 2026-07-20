package edu.tamu.pg.panhap.Tasks;

import edu.tamu.pg.panhap.Data.Fasta;
import edu.tamu.pg.panhap.Data.FastaIndex;
import edu.tamu.pg.panhap.Data.FastaReader;
import edu.tamu.pg.panhap.Data.GenomeMap;
import edu.tamu.pg.panhap.Utils.Defaults;
import edu.tamu.pg.panhap.Utils.PipelineSteps;
import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;
import org.jgrapht.alg.connectivity.ConnectivityInspector;

import java.nio.file.*;
import java.io.*;
import java.util.*;

public class PanCluster {
    private final Path clustersFile;
    private final Path singleCopyDir;
    private final Path scrFile;

    // Accept DirectoryStructure and retrieve file paths internally
    public PanCluster(DirectoryStructure dirs) {
        this.singleCopyDir = dirs.getDirectory(PipelineSteps.SINGLE_COPY);
        this.clustersFile = singleCopyDir.resolve("clusters.tsv");
        this.scrFile = singleCopyDir.resolve("all_single_copy_regions.fasta");
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
        FastaReader fastaReader = new FastaReader(scrFile);
        for (Fasta scr : fastaReader) {
            if (!graph.containsVertex(scr.getName())) {
                graph.addVertex(scr.getName());
            }
            SeqMeta seqMeta = new SeqMeta(scr.getName(), genomeMap);
            seqMeta.setLength(scr.getLength());
            // description is start-end
            int start = scr.getDescription().indexOf('-') > 0 ? Integer.parseInt(scr.getDescription().split("-")[0]) : -1;
            int end = scr.getDescription().indexOf('-') > 0 ? Integer.parseInt(scr.getDescription().split("-")[1]) : -1;
            seqMeta.setStart(start);
            seqMeta.setEnd(end);
            seqMetaMap.put(scr.getName(), seqMeta);
        }
        fastaReader.close();

        ConnectivityInspector<String, DefaultEdge> inspector = new ConnectivityInspector<>(graph);
        List<Set<String>> clusters = inspector.connectedSets();

        int clusterID = 0;
        Map<String, List<SeqMeta>> clusterMap = new LinkedHashMap<>();
        for (Set<String> cluster: clusters) {
            Set<SeqMeta> seqMetas = new HashSet<>();
            for (String seq: cluster) {
                seqMetas.add(seqMetaMap.get(seq));
            }
            // sort seqMetas
            List<SeqMeta> sortedSeqMetas = new ArrayList<>(seqMetas);

            // descending by length
            sortedSeqMetas.sort(Comparator.comparingInt(SeqMeta::getLength).reversed());

            clusterMap.put("cluster_" + clusterID, sortedSeqMetas);
            clusterID++;
        }

        for (String clusterName: clusterMap.keySet()) {
            for (SeqMeta seqMeta: clusterMap.get(clusterName)) {
                System.out.println(clusterName + " : " + seqMeta.toString());
            }
        }

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
}