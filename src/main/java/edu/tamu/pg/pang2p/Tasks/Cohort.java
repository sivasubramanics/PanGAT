package edu.tamu.pg.pang2p.Tasks;

import edu.tamu.pg.pang2p.Utils.Logger;
import picocli.CommandLine.*;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

@Command(name = "Cohort", description = "Compile raw mapped k-mer counts across sample coverage files into mapped_kmers.tsv.gz", sortOptions = false)
public class Cohort implements Runnable, Callable<Integer> {

    @Option(names = {"-i", "--input"}, description = "Path to input sample coverage files.", required = true, arity = "1..*", split = "\\s+")
    private Path[] input;

    @Option(names = {"-o", "--outprefix"}, description = "Prefix for output files.", required = true)
    private String outPrefix;

    @Option(names = {"-f", "--force"}, description = "Force overwrite of existing mapped_kmers file instead of updating.")
    private boolean force;

    private final String CLASS_NAME = Cohort.class.getSimpleName();

    public Cohort() {
    }

    public Cohort(Path[] input, String outPrefix, boolean force) {
        this.input = input;
        this.outPrefix = outPrefix;
        this.force = force;
    }

    private static class KmerData {
        String accId;
        String len;
        String totKmers;
        List<Double> rawCounts = new ArrayList<>();

        KmerData(String accId, String len, String totKmers) {
            this.accId = accId;
            this.len = len;
            this.totKmers = totKmers;
        }

        void ensureCapacity(int numSamples) {
            while (rawCounts.size() < numSamples) {
                rawCounts.add(0.0);
            }
        }
    }

    @Override
    public Integer call() throws Exception {
        if (input == null || input.length == 0) {
            Logger.error(CLASS_NAME, "Error: No input sample files provided.");
            return 1;
        }

        List<String> sampleNames = new ArrayList<>();
        Map<String, KmerData> kmerMap = new LinkedHashMap<>();
        Path mappedKmersPath = Path.of(outPrefix + ".mapped_kmers.tsv.gz");

        // 1. Check if existing mapped_kmers matrix exists
        if (Files.exists(mappedKmersPath)) {
            if (force) {
                Logger.info(CLASS_NAME, "Force flag (-f) set. Overwriting existing matrix: " + mappedKmersPath);
            } else {
                Logger.info(CLASS_NAME, "Existing mapped kmers matrix found. Loading for update: " + mappedKmersPath);
                try (BufferedReader reader = getBufferedReader(mappedKmersPath)) {
                    String header = reader.readLine();
                    if (header == null) {
                        throw new IllegalArgumentException("Existing mapped_kmers file is empty.");
                    }

                    String[] cols = header.split("\t");
                    for (int i = 3; i < cols.length; i++) {
                        sampleNames.add(cols[i]);
                    }

                    String line;
                    while ((line = reader.readLine()) != null) {
                        if (line.trim().isEmpty()) continue;
                        String[] parts = line.split("\t");
                        String accId = parts[0];
                        String len = parts[1];
                        String totKmers = parts[2];

                        KmerData var = new KmerData(accId, len, totKmers);
                        for (int i = 3; i < parts.length; i++) {
                            var.rawCounts.add(Double.parseDouble(parts[i]));
                        }
                        kmerMap.put(accId, var);
                    }
                }
                Logger.info(CLASS_NAME, "Loaded " + kmerMap.size() + " existing variants across " + sampleNames.size() + " samples.");
            }
        }

        // 2. Read new input coverage files
        Logger.info(CLASS_NAME, "Loading input coverage files...");

        for (Path file : input) {
            String sampleName = extractSampleName(file);

            if (sampleNames.contains(sampleName)) {
                Logger.info(CLASS_NAME, "Warning: Duplicate sample name '" + sampleName + "'. Skipping " + file);
                continue;
            }

            sampleNames.add(sampleName);
            int currentSampleIdx = sampleNames.size() - 1;
            Logger.info(CLASS_NAME, "Processing sample (" + sampleNames.size() + "): " + sampleName + " from " + file.getFileName());

            try (BufferedReader reader = getBufferedReader(file)) {
                String line;
                while ((line = reader.readLine()) != null) {
                    if (line.trim().isEmpty() || line.startsWith("#")) {
                        continue;
                    }

                    String[] parts = line.split("\t");
                    if (parts.length < 5) continue;

                    String accId = parts[0];
                    String len = parts[1];
                    String totKmers = parts[2];
                    double mappedKmers = Double.parseDouble(parts[4]);

                    KmerData var = kmerMap.computeIfAbsent(accId, k -> new KmerData(accId, len, totKmers));
                    var.ensureCapacity(currentSampleIdx);
                    var.rawCounts.add(mappedKmers);
                }
            }

            for (KmerData var : kmerMap.values()) {
                var.ensureCapacity(sampleNames.size());
            }
        }

        int numSamples = sampleNames.size();
        List<KmerData> masterList = new ArrayList<>(kmerMap.values());

        // 3. Write output mapped_kmers.tsv.gz
        Logger.info(CLASS_NAME, "Writing mapped_kmers matrix -> " + mappedKmersPath);

        try (BufferedWriter writer = getBufferedWriter(mappedKmersPath)) {
            writer.write("acc.id\tlen\ttot_kmers");
            for (String sample : sampleNames) {
                writer.write("\t" + sample);
            }
            writer.write("\n");

            for (KmerData var : masterList) {
                writer.write(var.accId + "\t" + var.len + "\t" + var.totKmers);
                for (int s = 0; s < numSamples; s++) {
                    writer.write(String.format(Locale.US, "\t%.0f", var.rawCounts.get(s)));
                }
                writer.write("\n");
            }
        }

        Logger.info(CLASS_NAME, "Aggregated " + masterList.size() + " variants across " + numSamples + " samples.");
        return 0;
    }

    @Override
    public void run() {
        try {
           call();
        } catch (Exception e) {
            Logger.error(CLASS_NAME, e.getMessage());
        }
    }

    private BufferedReader getBufferedReader(Path path) throws Exception {
        InputStream is = Files.newInputStream(path);
        if (path.getFileName().toString().endsWith(".gz")) is = new GZIPInputStream(is);
        return new BufferedReader(new InputStreamReader(is));
    }

    private BufferedWriter getBufferedWriter(Path path) throws Exception {
        OutputStream os = Files.newOutputStream(path);
        if (path.getFileName().toString().endsWith(".gz")) os = new GZIPOutputStream(os);
        return new BufferedWriter(new OutputStreamWriter(os));
    }

    private String extractSampleName(Path file) {
        String name = file.getFileName().toString();
        if (name.endsWith(".gz")) name = name.substring(0, name.length() - 3);
        int dotIndex = name.indexOf('.');
        return (dotIndex > 0) ? name.substring(0, dotIndex) : name;
    }
}