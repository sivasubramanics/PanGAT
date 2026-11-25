package edu.tamu.pg.pang2p.Data;

import edu.tamu.pg.pang2p.Utils.*;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.*;
import java.util.concurrent.*;

public class PanHap {
    private static final String CLASS_NAME = PanHap.class.getSimpleName();
    private static final int FIXED_COLS = 15;          // columns before the sample columns
    private static final double SCALE = 1_000_000_000.0;

    private final Map<String, HapData> hapData = new ConcurrentHashMap<>();
    private final List<String> samples = new ArrayList<>();
    private final Map<String, Integer> sampleIndex = new HashMap<>();
    private double[] sampleSums = new double[0];

    private String version = "0.0.1";
    private final List<String> commandLines = new ArrayList<>();
    private final Params parameters = new Params();
    private boolean minCountSet = false;

    private PanHap() {}

    // ------------------------------------------------------------------
    // Factories
    // ------------------------------------------------------------------

    /** Start a new PanHap from a pan-cluster representative file and k-mer counts. */
    public PanHap(Path repFile, Path panClusterFastaFile, int kmerSize, int skipSize, int threads) {
        this.readRepFile(repFile);
        try {
            this.addOrUpdateKmerCounts(panClusterFastaFile, kmerSize, skipSize, threads);
        } catch (IOException e) {
            throw new RuntimeException("Error counting k-mers for pan-cluster fasta: " + panClusterFastaFile, e);
        }
        parameters.set("kmerSize", kmerSize);
        parameters.set("skipSize", skipSize);
    }

    /** Resume from a previously written matrix. */
    public PanHap(Path panHapFile) {
        this.readMatrix(panHapFile);
    }

    // ------------------------------------------------------------------
    // Growing
    // ------------------------------------------------------------------

    /** Add a sample from a per-sample count file. Returns this for chaining. */
    public PanHap add(Path covFile) {
        String localSampleName = null;
        int localKmerSize = 0;
        int localMinCount = 0;
        int localSkipSize = 0;
        int localNumHaps = 0;
        int sampleIdx = -1;

        try (BufferedReader reader = FileUtils.reader(covFile)) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.trim().isEmpty()) continue;

                if (line.startsWith("##parameter=<") && line.endsWith(">")) {
                    String body = line.substring("##parameter=<".length(), line.length() - 1);
                    String key = "", value = "";
                    for (String kv : body.split(",")) {
                        String[] t = kv.split("=", 2);
                        if (t.length < 2) continue;
                        if (t[0].equals("ID")) key = t[1];
                        else if (t[0].equals("Value")) value = t[1];
                    }
                    switch (key) {
                        case "sampleName" -> localSampleName = value;
                        case "kmerSize" -> localKmerSize = Integer.parseInt(value);
                        case "minCount" -> localMinCount = Integer.parseInt(value);
                        case "skipSize" -> localSkipSize = Integer.parseInt(value);
                        case "numHaps" -> localNumHaps = Integer.parseInt(value);
                        default -> Logger.error(CLASS_NAME, "Unknown parameter in count file: " + key);
                    }
                    continue;
                }

                if (line.startsWith("HapID") || line.startsWith("# HapID")) {
                    sampleIdx = applyAndValidateCovParams(localSampleName, localKmerSize, localMinCount, localSkipSize, localNumHaps);
                    continue;
                }

                if (line.startsWith("#")) continue;
                if (sampleIdx < 0) {
                    sampleIdx = applyAndValidateCovParams(localSampleName, localKmerSize, localMinCount, localSkipSize, localNumHaps);
                }

                String[] p = line.split("\t");
                String hapID = p[0];
                int obsKmers = Integer.parseInt(p[4]);
                int obsTotalKmers = Integer.parseInt(p[5]);

                HapData hap = hapData.get(hapID);
                if (hap == null) {
                    Logger.error(CLASS_NAME, "HapID " + hapID + " not found in rep file.");
                    continue;
                }
                hap.addSample(obsKmers, obsTotalKmers);
                sampleSums[sampleIdx] += obsKmers;
            }
        } catch (IOException e) {
            throw new RuntimeException("Error reading coverage file: " + covFile, e);
        }
        return this;
    }

    private int applyAndValidateCovParams(
            String localSampleName,
            int localKmerSize,
            int localMinCount,
            int localSkipSize,
            int localNumHaps
    ) {
        if (parameters.get("kmerSize") != 0 && parameters.get("kmerSize") != localKmerSize) {
            Logger.error(CLASS_NAME, "K-mer size mismatch: expected " + parameters.get("kmerSize") + ", found " + localKmerSize);
        }
        if (parameters.get("skipSize") != 0 && parameters.get("skipSize") != localSkipSize) {
            Logger.error(CLASS_NAME, "Skip size mismatch: expected " + parameters.get("skipSize") + ", found " + localSkipSize);
        }
        if (parameters.get("numHaps") != 0 && parameters.get("numHaps") != localNumHaps) {
            Logger.error(CLASS_NAME, "Number of haplotypes mismatch: expected " + parameters.get("numHaps") + ", found " + localNumHaps);
        }

        if (samples.isEmpty() && !minCountSet) {
            parameters.set("minCount", localMinCount);
            minCountSet = true;
        } else if (parameters.get("minCount") != localMinCount) {
            Logger.error(CLASS_NAME, "Min count mismatch: expected " + parameters.get("minCount") + ", found " + localMinCount);
        }

        return registerSample(localSampleName);
    }

    // ------------------------------------------------------------------
    // Reading internals
    // ------------------------------------------------------------------

    private void readRepFile(Path repFile) {
        Logger.info(CLASS_NAME, "Loading pan-cluster representative file: " + repFile);
        try (BufferedReader reader = FileUtils.reader(repFile)) {
            String line;
            boolean header = true;
            while ((line = reader.readLine()) != null) {
                if (line.trim().isEmpty()) continue;
                if (header) { header = false; continue; }
                String[] p = line.split("\t");
                String simChromID = getSimChromID(p[3]);
                hapData.computeIfAbsent(p[1], k -> new HapData(
                        p[1], p[4], p[5], simChromID,
                        Integer.parseInt(p[7]),    // start
                        Integer.parseInt(p[8]),    // end
                        Integer.parseInt(p[6])     // length
                ));
            }
            parameters.set("numHaps", hapData.size());
            Logger.info(CLASS_NAME, "Loaded " + parameters.get("numHaps") + " hap data");
        } catch (IOException e) {
            throw new RuntimeException("Error reading rep file: " + repFile, e);
        }
    }

    public PanHap addOrUpdateKmerCounts(Path panClusterFastaFile, int kmerSize, int skipSize, int nThreads) throws IOException {
        Logger.info(CLASS_NAME, "Counting k-mers for pan-cluster fasta: " + panClusterFastaFile);
        Path kmcDB = Path.of(panClusterFastaFile + ".k" + kmerSize);
        KMCConfig kmcConfig = KMCConfig.builder()
                .inputPaths(panClusterFastaFile)
                .outputPrefix(kmcDB)
                .kmerSize(kmerSize)
                .threads(nThreads)
                .minCount(1)
                .build();
        KMCExecutor.run(kmcConfig);

        ExecutorService executor = Executors.newFixedThreadPool(nThreads);
        Deque<Future<?>> pending = new ArrayDeque<>();
        int lastLoggedPercent = 0;

        try (FastaReader reader = new FastaReader(panClusterFastaFile);
             KMC kmc = new KMC(kmcDB, true)) {

            int prefixLength = kmc.getPrefixLength();
            for (Fasta fasta : reader) {
                pending.addLast(executor.submit(() -> countKmers(fasta, kmc, kmerSize, skipSize, prefixLength)));
                while (!pending.isEmpty() && pending.peekFirst().isDone()) {
                    pending.pollFirst().get();
                }
            }

            int totalTasks = pending.size();
            int completed = 0;
            while (!pending.isEmpty()) {
                pending.pollFirst().get();
                completed++;
                int percent = (int) ((completed * 100L) / Math.max(1, totalTasks));
                if (percent % 5 == 0 && percent > lastLoggedPercent) {
                    lastLoggedPercent = percent;
                    Logger.progress(
                            CLASS_NAME,
                            String.format("Counting k-mers %s%d%%", ".".repeat(percent / 5), percent)
                    );
                }
            }

        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Interrupted while counting k-mers", e);
        } catch (ExecutionException e) {
            throw new RuntimeException("K-mer counting failed", e.getCause());
        } finally {
            executor.shutdownNow();
        }

        parameters.set("kmerSize", kmerSize);
        parameters.set("skipSize", skipSize);
        return this;
    }

    private void readMatrix(Path matrixFile) {
        Logger.info(CLASS_NAME, "Loading pan-matrix from: " + matrixFile);
        try (BufferedReader reader = FileUtils.reader(matrixFile)) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.trim().isEmpty()) continue;
                if (line.startsWith("##")) { parseMetadata(line.substring(2)); continue; }
                if (line.startsWith("hapID") || line.startsWith("# hapID")) {
                    String[] p = line.replace("# ", "").split("\t");
                    for (int i = FIXED_COLS; i < p.length; i++) registerSample(p[i]);
                    continue;
                }
                parseMatrixRow(line.split("\t"));
            }
            parameters.set("numHaps", hapData.size());
            parameters.set("numSamples", samples.size());
            Logger.info(CLASS_NAME, "Loaded " + parameters.get("numHaps") + " hap data for " + parameters.get("numSamples") + " samples");
        } catch (IOException e) {
            throw new RuntimeException("Error reading matrix file: " + matrixFile, e);
        }
    }

    private void parseMatrixRow(String[] p) {
        int n = p.length - FIXED_COLS;
        int[] gt = new int[n], obs = new int[n], obsTotal = new int[n];
        double[] gtScaled = new double[n], norm = new double[n];

        for (int i = 0; i < n; i++) {
            String[] a = p[i + FIXED_COLS].split(":");
            gt[i] = Integer.parseInt(a[0]);
            gtScaled[i] = Double.parseDouble(a[1]);
            norm[i] = Double.parseDouble(a[2]);
            obs[i] = Integer.parseInt(a[3]);
            obsTotal[i] = Integer.parseInt(a[4]);
            sampleSums[i] += obs[i];
        }

        hapData.computeIfAbsent(p[0], k -> new HapData(
                p[0], p[1], p[2], p[3],
                Integer.parseInt(p[4]),    // start
                Integer.parseInt(p[5]),    // end
                Integer.parseInt(p[6])     // length
        ));
        HapData hap = hapData.get(p[0]);
        hap.setTotKmers(Integer.parseInt(p[7]));
        hap.setUniqueKmers(Integer.parseInt(p[8]));
        hap.setGlobalUniqueKmers(Integer.parseInt(p[9]));
        hap.setGenotype(gt);
        hap.setGtScaled(gtScaled);
        hap.setNormalizedCounts(norm);
        hap.setObsCounts(obs);
        hap.setObsTotalCounts(obsTotal);
    }

    private int registerSample(String id) {
        if (id == null || id.trim().isEmpty()) {
            throw new IllegalArgumentException("Sample name is null or empty. Sample names must be non-empty strings.");
        }
        if (sampleIndex.containsKey(id)) {
            throw new IllegalArgumentException("Sample name " + id + " already exists. Duplicate sample names are not allowed.");
        }
        int idx = samples.size();
        samples.add(id);
        sampleIndex.put(id, idx);
        sampleSums = Arrays.copyOf(sampleSums, samples.size());
        parameters.set("numSamples", samples.size());
        return idx;
    }

    private void parseMetadata(String content) {
        if (content.startsWith("version=")) {
            version = content.split("=", 2)[1];
        } else if (content.startsWith("commandline=")) {
            addCommandLine(content.split("=", 2)[1]);
        } else if (content.startsWith("parameter=<") && content.endsWith(">")) {
            String body = content.substring("parameter=<".length(), content.length() - 1);
            String key = "", value = "";
            for (String kv : body.split(",")) {
                String[] t = kv.split("=", 2);
                if (t.length < 2) continue;
                if (t[0].equals("ID")) key = t[1];
                else if (t[0].equals("Value")) value = t[1];
            }
            if (!key.isEmpty()) {
                parameters.set(key, Integer.parseInt(value));
                if ("minCount".equals(key)) {
                    minCountSet = true;
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // Normalization / output
    // ------------------------------------------------------------------

    public PanHap normalize() {
        Logger.info(CLASS_NAME, "Normalizing pan-hap data for " + samples.size() + " samples.");
        for (HapData hap : hapData.values()) {
            int[] obs = hap.getObsCounts();
            double[] norm = new double[samples.size()];
            if (obs != null) {
                for (int i = 0; i < obs.length; i++) {
                    norm[i] = sampleSums[i] > 0 ? (obs[i] / sampleSums[i]) * SCALE : 0.0;
                }
            }
            hap.setNormalizedCounts(norm);
            hap.calcMetrics();
        }
        return this;
    }

    public List<HapData> getSortedHapData() {
        List<HapData> list = new ArrayList<>(hapData.values());
        Collections.sort(list);
        return list;
    }

    private void countKmers(Fasta fasta, KMC kmc, int kmerSize, int skipSize, int prefixLength) {
        HapData hap = hapData.get(fasta.getName());
        if (hap == null) {
            Logger.error(CLASS_NAME, "Sequence " + fasta.getName() + " not found in rep file.");
            return;
        }

        KmerIterator it = fasta.kmers(kmerSize, prefixLength, true);
        Set<Kmer> unique = new HashSet<>();
        int total = 0, globalUnique = 0;
        int kmerIndex = 0;

        while (it.hasNext()) {
            Kmer kmer = it.next().getKmer();
            if (skipSize > 1 && kmerIndex % skipSize != 0) {
                kmerIndex++;
                continue;
            }
            kmerIndex++;
            total++;
            unique.add(kmer);
            if (kmc.getCount(kmer) == 1) globalUnique++;
        }

        hap.setTotKmers(total);
        hap.setUniqueKmers(unique.size());
        hap.setGlobalUniqueKmers(globalUnique);
    }

    public PanHap write(Path out) {
        normalize();
        Logger.info(CLASS_NAME, "Writing pan-genome matrix to: " + out);
        try (BufferedWriter w = FileUtils.writer(out)) {
            writeHeader(w);
            for (HapData hap : getSortedHapData()) {
                w.write(hap.toString());
                w.newLine();
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return this;
    }

    private void writeHeader(BufferedWriter w) throws IOException {
        parameters.set("numHaps", hapData.size());
        parameters.set("numSamples", samples.size());

        w.write("##fileformat=\"PanG2P v0.1.0\"\n");
        w.write("##fileDate=" + LocalDateTime.now().format(DateTimeFormatter.ofPattern("yyyyMMdd_HHmmss")) + "\n");
        w.write("##version=" + version + "\n");
        for (String cmd : commandLines) w.write("##commandline=" + cmd + "\n");
        for (Map.Entry<String, Integer> e : parameters.entries()) {
            w.write("##parameter=<ID=" + e.getKey() + ",Value=" + e.getValue() + ">\n");
        }
        w.write("##FORMAT=<ID=GT:GTScaled:NORM:OBS:TOTAL,Number=5,Description=\"Genotype:GenotypeScaled:NormalizedCount:MappedKmers:MappedTotalKmers\">\n");
        w.write("# hapID\tGenomeID\tChromID\tSimChromID\tStart\tEnd\tLength\tTotKmers\tUniqueKmers\tGlobalUniqueKmers\tMinRPB\tMaxRPB\tMeanRPB\tSdRPB\tCV");
        for (String s : samples) w.write("\t" + s);
        w.newLine();
    }

    // ------------------------------------------------------------------
    // Accessors
    // ------------------------------------------------------------------

    public int getHapNum() { return hapData.size(); }
    public int getSampleNum() { return samples.size(); }
    public int getKmerSize() { return parameters.get("kmerSize"); }
    public int getSkipSize() { return parameters.get("skipSize"); }
    public int getMinCount() { return parameters.get("minCount"); }
    public int getNumHaps() { return parameters.get("numHaps"); }
    public int getNumSamples() { return parameters.get("numSamples"); }
    public List<String> getSamples() { return Collections.unmodifiableList(samples); }
    public String getVersion() { return version; }
    public PanHap setVersion(String v) { this.version = v; return this; }

    public PanHap addCommandLine(String cmd) {
        if (cmd != null && !cmd.trim().isEmpty()) commandLines.add(cmd);
        return this;
    }

    public String getSimChromID(String newChromID) {
        if (newChromID != null && newChromID.length() >= Defaults.CHROM_CODE_LENGTH) {
            return newChromID.substring(Defaults.CHROM_CODE_LENGTH);
        }
        return newChromID;
    }

    private static final class Params {
        int kmerSize;
        int minCount;
        int skipSize;
        int numHaps;
        int numSamples;

        void apply(String key, String value) {
            switch (key) {
                case "kmerSize" -> this.kmerSize = Integer.parseInt(value);
                case "minCount" -> this.minCount = Integer.parseInt(value);
                case "skipSize" -> this.skipSize = Integer.parseInt(value);
                case "numHaps" -> this.numHaps = Integer.parseInt(value);
                case "numSamples" -> this.numSamples = Integer.parseInt(value);
                default -> Logger.error(CLASS_NAME, "Unknown parameter: " + key);
            }
        }

        void apply(String key, int value) {
            switch (key) {
                case "kmerSize" -> this.kmerSize = value;
                case "minCount" -> this.minCount = value;
                case "skipSize" -> this.skipSize = value;
                case "numHaps" -> this.numHaps = value;
                case "numSamples" -> this.numSamples = value;
                default -> Logger.error(CLASS_NAME, "Unknown parameter: " + key);
            }
        }

        public void set(String key, int value) {
            apply(key, value);
        }

        public int get(String key) {
            return switch (key) {
                case "kmerSize" -> kmerSize;
                case "minCount" -> minCount;
                case "skipSize" -> skipSize;
                case "numHaps" -> numHaps;
                case "numSamples" -> numSamples;
                default -> 0;
            };
        }

        public Set<Map.Entry<String, Integer>> entries() {
            Map<String, Integer> map = new LinkedHashMap<>();
            map.put("kmerSize", kmerSize);
            map.put("minCount", minCount);
            map.put("skipSize", skipSize);
            map.put("numHaps", numHaps);
            map.put("numSamples", numSamples);
            return map.entrySet();
        }
    }
}