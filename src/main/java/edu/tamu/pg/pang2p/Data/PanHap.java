package edu.tamu.pg.pang2p.Data;

import edu.tamu.pg.pang2p.Utils.FileUtils;
import edu.tamu.pg.pang2p.Utils.KMCConfig;
import edu.tamu.pg.pang2p.Utils.KMCExecutor;
import edu.tamu.pg.pang2p.Utils.Logger;

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
    private int kmerSize = 0;
    private int skipSize = 0;
    private int minCount = 0;
    private int numHaps = 0;
    private int numSamples = 0;

    private final Map<String, HapData> hapData = new ConcurrentHashMap<>();
    private final List<String> samples = new ArrayList<>();
    private final Map<String, Integer> sampleIndex = new HashMap<>();
    private double[] sampleSums = new double[0];

    private String version = "0.0.1";
    private final List<String> commandLines = new ArrayList<>();
    private final Map<String, String> parameters = new LinkedHashMap<>();

    private PanHap() {}

    // ------------------------------------------------------------------
    // Factories
    // ------------------------------------------------------------------

    /** Start a new PanHap from a pan-cluster representative file (no samples yet). */
    public static PanHap init(Path repFile) {
        PanHap ph = new PanHap();
        ph.readRepFile(repFile);
        return ph;
    }

    /** Resume from a previously written matrix. The rep file supplies hap definitions. */
    public static PanHap init(Path panHapFile, Path repFile) {
        PanHap ph = init(repFile);
        ph.readMatrix(panHapFile);
        return ph;
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
        try (BufferedReader reader = FileUtils.reader(covFile)) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.trim().isEmpty()) continue;
                if (line.startsWith("##")) {
                    if (line.startsWith("parameter=<") && line.endsWith(">")) {
                        String body = line.substring("parameter=<".length(), line.length() - 1);
                        String key = "", value = "";
                        for (String kv : body.split(",")) {
                            String[] t = kv.split("=", 2);
                            if (t.length < 2) continue;
                            if (t[0].equals("ID")) key = t[1];
                            else if (t[0].equals("Value")) value = t[1];
                        }
                        if (key.equals("sampleName")) {
                            localSampleName = value;
                        } else if (key.equals("kmerSize")) {
                            localKmerSize = Integer.parseInt(value);
                        } else if (key.equals("minCount")) {
                            localMinCount = Integer.parseInt(value);
                        } else if (key.equals("skipSize")) {
                            localSkipSize = Integer.parseInt(value);
                        } else if (key.equals("numHaps")) {
                            // ignore these, they are not needed for adding a sample
                            localNumHaps = Integer.parseInt(value);
                        } else {
                            Logger.error(CLASS_NAME, "Unknown parameter in count file: " + key);
                        }
                    }
                }
                if (line.startsWith("HapID") || line.startsWith("# HapID")) {
                    if (kmerSize != localKmerSize) {
                        Logger.error(CLASS_NAME, "K-mer size mismatch: expected " + kmerSize + ", found " + localKmerSize);
                    }
                    if (minCount != localMinCount) {
                        Logger.error(CLASS_NAME, "Min count mismatch: expected " + minCount + ", found " + localMinCount);
                    }
                    if (skipSize != localSkipSize) {
                        Logger.error(CLASS_NAME, "Skip size mismatch: expected " + skipSize + ", found " + localSkipSize);
                    }
                    if (numHaps != localNumHaps) {
                        Logger.error(CLASS_NAME, "Number of haplotypes mismatch: expected " + numHaps + ", found " + localNumHaps);
                    }
                    registerSample(localSampleName);
                    continue;
                }

                String[] p = line.split("\t");
                String hapID = p[0];
                int length = Integer.parseInt(p[1]);
                int totalKmers = Integer.parseInt(p[2]);
                int uniqueKmers = Integer.parseInt(p[3]);
                int obsKmers = Integer.parseInt(p[4]);
                int obsTotalKmers = Integer.parseInt(p[5]);

                HapData hap = hapData.get(hapID);
                if (hap == null) {
                    Logger.error(CLASS_NAME, "HapID " + hapID + " not found in rep file.");
                    continue;
                }
                hap.addSample(obsKmers, obsTotalKmers);
            }
        } catch (IOException e) {
            throw new RuntimeException("Error reading coverage file: " + covFile, e);
        }
        return this;
    }

        // ------------------------------------------------------------------
    // Reading internals
    // ------------------------------------------------------------------

    private void readRepFile(Path repFile) {
        try (BufferedReader reader = FileUtils.reader(repFile)) {
            String line;
            boolean header = true;
            while ((line = reader.readLine()) != null) {
                if (line.trim().isEmpty()) continue;
                if (header) { header = false; continue; }
                String[] p = line.split("\t");
                hapData.computeIfAbsent(p[1], k -> new HapData(
                        p[0], p[1], p[2], p[3], p[4], p[5],
                        Integer.parseInt(p[7]),    // start
                        Integer.parseInt(p[8]),    // end
                        Integer.parseInt(p[6]),    // length
                        Integer.parseInt(p[9]),    // totalKmers
                        Integer.parseInt(p[10]),   // uniqueKmers
                        Integer.parseInt(p[11]))); // globalUniqueKmers
            }
            numHaps = hapData.size();
            Logger.info(CLASS_NAME, "Loaded " + numHaps + " hap data");
        } catch (IOException e) {
            throw new RuntimeException("Error reading rep file: " + repFile, e);
        }
    }

    public PanHap addOrUpdateKmerCounts(Path panClusterFastaFile, int kmerSize, int stepSize, int nThreads) throws IOException {
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
        try (FastaReader reader = new FastaReader(panClusterFastaFile);
             KMC kmc = new KMC(kmcDB, true)) {

            int prefixLength = kmc.getPrefixLength();
            List<Future<?>> futures = new ArrayList<>();
            for (Fasta fasta : reader) {
                futures.add(executor.submit(() -> countKmers(fasta, kmc, kmerSize, prefixLength)));
            }
            // Wait for every task BEFORE the try block closes kmc/reader,
            // and surface any exception thrown inside a worker.
            for (Future<?> f : futures) f.get();

        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Interrupted while counting k-mers", e);
        } catch (ExecutionException e) {
            throw new RuntimeException("K-mer counting failed", e.getCause());
        } finally {
            executor.shutdownNow();
        }

        // Stored as parameters so they are written to the matrix header
        // and restored when the matrix is reloaded.
        parameters.put("kmerSize", String.valueOf(kmerSize));
        parameters.put("minCount", String.valueOf(minCount));
        parameters.put("skipSize", String.valueOf(skipSize));
        parameters.put("numHaps", String.valueOf(numHaps));
        parameters.put("numSamples", String.valueOf(numSamples));
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
            numHaps = hapData.size();
            numSamples = samples.size();
            Logger.info(CLASS_NAME, "Loaded " + numHaps + " haps, " + numSamples + " samples");
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

        HapData hap = hapData.get(p[0]);
        if (hap == null) {
            Logger.error(CLASS_NAME, "HapID " + p[0] + " in matrix not found in rep file.");
            return;
        }
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
        int idx = samples.size();
        if (sampleIndex.containsKey(id)) {
            Logger.error(CLASS_NAME, "Sample " + id + " already exists in the pan-matrix. Sample names must be unique.");
        }
        samples.add(id);
        sampleIndex.put(id, idx);
        sampleSums = Arrays.copyOf(sampleSums, samples.size());
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
            if (!key.isEmpty()) parameters.put(key, value);
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

    private void countKmers(Fasta fasta, KMC kmc, int kmerSize, int prefixLength) {
        HapData hap = hapData.get(fasta.getName());
        if (hap == null) {
            Logger.error(CLASS_NAME, "Sequence " + fasta.getName() + " not found in rep file.");
            return;
        }

        KmerIterator it = fasta.kmers(kmerSize, prefixLength, true);
        Set<Kmer> unique = new HashSet<>();
        int total = 0, globalUnique = 0;

        while (it.hasNext()) {
            Kmer kmer = it.next().getKmer();
            total++;
            unique.add(kmer);
            if (kmc.getCount(kmer) == 1) globalUnique++;
        }

        hap.setTotKmers(total);
        hap.setUniqueKmers(unique.size());
        hap.setGlobalUniqueKmers(globalUnique);
    }

    public PanHap write(Path out) {
        normalize(); // explicit, instead of hiding inside getSortedHapData()
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
        w.write("##fileformat=\"PanG2P v0.1.0\"\n");
        w.write("##fileDate=" + LocalDateTime.now().format(DateTimeFormatter.ofPattern("yyyyMMdd_HHmmss")) + "\n");
        w.write("##version=" + version + "\n");
        for (String cmd : commandLines) w.write("##commandline=" + cmd + "\n");
        for (Map.Entry<String, String> e : parameters.entrySet())
            w.write("##parameter=<ID=" + e.getKey() + ",Value=" + e.getValue() + ">\n");
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
    public int getKmerSize() { return kmerSize; }
    public int getSkipSize() { return skipSize; }
    public int getMinCount() { return minCount; }
    public int getNumHaps() { return numHaps; }
    public int getNumSamples() { return numSamples; }
    public List<String> getSamples() { return Collections.unmodifiableList(samples); }
    public String getVersion() { return version; }
    public PanHap setVersion(String v) { this.version = v; return this; }
    public List<String> getCommandLines() { return commandLines; }
    public Map<String, String> getParameters() { return parameters; }

    public PanHap addCommandLine(String cmd) {
        if (cmd != null && !cmd.trim().isEmpty()) commandLines.add(cmd);
        return this;
    }

    public PanHap addParameter(String key, String value) {
        parameters.put(key, value);
        return this;
    }
}