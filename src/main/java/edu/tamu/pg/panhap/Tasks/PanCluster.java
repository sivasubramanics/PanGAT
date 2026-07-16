package edu.tamu.pg.panhap.Tasks;

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
    private final Path faiFile;
    private final Path outDir;

    // Accept DirectoryStructure and retrieve file paths internally
    public PanCluster(DirectoryStructure dirs) {
        Path singleCopyDir = dirs.getDirectory(PipelineSteps.SINGLE_COPY);
        this.clustersFile = singleCopyDir.resolve("clusters.tsv");
        this.faiFile = singleCopyDir.resolve("all_single_copy_regions.faidx");
        this.outDir = singleCopyDir; // outputs go to same directory, change if needed
    }

    // Simple run() method (no arguments)
    public void run() throws IOException {
        Files.createDirectories(outDir);

        List<Edge> edges = readEdges(clustersFile);
        Graph<String, DefaultEdge> graph = new SimpleGraph<>(DefaultEdge.class);
        for (Edge e : edges) {
            graph.addVertex(e.seq1);
            graph.addVertex(e.seq2);
            graph.addEdge(e.seq1, e.seq2);
        }

        Map<String, SeqMeta> seqMeta = readFai(faiFile);

        // Add singletons
        Set<String> allSeq = new HashSet<>(seqMeta.keySet());
        Set<String> graphSeq = graph.vertexSet();
        for (String seq : allSeq) {
            if (!graphSeq.contains(seq))
                graph.addVertex(seq);
        }

        ConnectivityInspector<String, DefaultEdge> inspector = new ConnectivityInspector<>(graph);
        List<Set<String>> clusters = inspector.connectedSets();

        Map<String, String> seqToCluster = new HashMap<>();
        int clusterIdx = 1;
        for (Set<String> cluster : clusters) {
            String clusterId = "cluster_" + clusterIdx++;
            for (String seq : cluster)
                seqToCluster.put(seq, clusterId);
        }

        try (BufferedWriter writer = Files.newBufferedWriter(outDir.resolve("pan_cluster_all_member.txt"))) {
            writer.write("seq\tseq.length\tgenotype\tchr\tstart\tend\tcomponents\n");
            for (String seq : seqMeta.keySet()) {
                SeqMeta meta = seqMeta.get(seq);
                writer.write(String.join("\t",
                        seq, String.valueOf(meta.length),
                        meta.genotype, meta.chr,
                        String.valueOf(meta.start), String.valueOf(meta.end),
                        seqToCluster.get(seq)) + "\n");
            }
        }

        Map<String, String> repSeq = new HashMap<>();
        Map<String, Integer> repLen = new HashMap<>();
        for (Set<String> cluster : clusters) {
            String clusterId = seqToCluster.get(cluster.iterator().next());
            String rep = null;
            int maxLen = -1;
            for (String seq : cluster) {
                int len = seqMeta.get(seq).length;
                if (len > maxLen) {
                    rep = seq;
                    maxLen = len;
                }
            }
            repSeq.put(clusterId, rep);
            repLen.put(clusterId, maxLen);
        }

        try (BufferedWriter writer = Files.newBufferedWriter(outDir.resolve("pan_cluster_mmseq_rep.txt"))) {
            writer.write("components\trep\trep.length\n");
            for (String clusterId : repSeq.keySet()) {
                writer.write(clusterId + "\t" + repSeq.get(clusterId) + "\t" + repLen.get(clusterId) + "\n");
            }
        }

        try (BufferedWriter writer = Files.newBufferedWriter(outDir.resolve("pan_cluster_bed_rep.txt"))) {
            writer.write("chr\tstart\tend\trep\trep.genotype\n");
            for (String clusterId : repSeq.keySet()) {
                String rep = repSeq.get(clusterId);
                SeqMeta meta = seqMeta.get(rep);
                writer.write(String.join("\t",
                        meta.chr, String.valueOf(meta.start), String.valueOf(meta.end), rep, meta.genotype) + "\n");
            }
        }
    }

    private static class Edge {
        String seq1, seq2;
        Edge(String seq1, String seq2) { this.seq1 = seq1; this.seq2 = seq2; }
    }

    private static class SeqMeta {
        String seq;
        int length, start, end;
        String genotype, chr;
        SeqMeta(String seq, int length) {
            this.seq = seq;
            this.length = length;
            genotype = seq.contains("|") ? seq.split("\\|")[0] : "";
            String rest = seq.contains("|") ? seq.split("\\|")[1] : seq;
            chr = rest.contains(":") ? rest.split(":")[0] : "";
            String[] coords = rest.split(":");
            if (coords.length > 1) {
                String[] range = coords[1].split("-");
                if (range.length == 2) {
                    start = Integer.parseInt(range[0]);
                    end = Integer.parseInt(range[1]);
                }
            }
        }
    }

    private List<Edge> readEdges(Path clustersFile) throws IOException {
        List<Edge> edges = new ArrayList<>();
        try (BufferedReader reader = Files.newBufferedReader(clustersFile)) {
            String line;
            while ((line = reader.readLine()) != null) {
                String[] parts = line.trim().split("\t");
                if (parts.length >= 2)
                    edges.add(new Edge(parts[0], parts[1]));
            }
        }
        return edges;
    }

    private Map<String, SeqMeta> readFai(Path faiFile) throws IOException {
        Map<String, SeqMeta> seqMeta = new HashMap<>();
        try (BufferedReader reader = Files.newBufferedReader(faiFile)) {
            String line;
            while ((line = reader.readLine()) != null) {
                String[] parts = line.trim().split("\t");
                if (parts.length >= 2) {
                    String seq = parts[0];
                    int len = Integer.parseInt(parts[1]);
                    seqMeta.put(seq, new SeqMeta(seq, len));
                }
            }
        }
        return seqMeta;
    }
}