package edu.tamu.pg.pang2p.Tasks;

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
import edu.tamu.pg.pang2p.Utils.Logger;

@Command(name = "ExportMatrix", description = "Export normalized coverage matrix and discrete genotype matrix from mapped_kmers matrix", sortOptions = false)
public class ExportMatrix implements Runnable, Callable<Integer> {

    @Option(names = {"-k", "--kmers"}, description = "Path to input mapped_kmers.tsv.gz file.", required = true)
    private Path kmersFile;

    @Option(names = {"-o", "--outprefix"}, description = "Prefix for output files.", required = true)
    private String outPrefix;

    @Option(names = {"-v", "--cv-cutoff"}, description = "CV cutoff filter threshold (default: 0.3)", defaultValue = "0.3")
    private double cvCutoff;

    private final String CLASS_NAME = ExportMatrix.class.getSimpleName();

    public ExportMatrix() {
    }

    public ExportMatrix(Path kmersFile, String outPrefix, double cvCutoff) {
        this.kmersFile = kmersFile;
        this.outPrefix = outPrefix;
        this.cvCutoff = cvCutoff;
    }

    private static class VariantRecord {
        String accId;
        String len;
        String totKmers;
        double[] rawCounts;
        double[] rpbValues;
        double meanRpb;
        double sdRpb;
        double cv;
        int[] dosages;

        VariantRecord(String accId, String len, String totKmers, int numSamples) {
            this.accId = accId;
            this.len = len;
            this.totKmers = totKmers;
            this.rawCounts = new double[numSamples];
            this.rpbValues = new double[numSamples];
            this.dosages = new int[numSamples];
        }
    }

    @Override
    public Integer call() throws Exception {
        if (!Files.exists(kmersFile)) {
            Logger.error(CLASS_NAME, "Input file does not exist.");
            return 1;
        }

        List<String> sampleNames = new ArrayList<>();
        List<VariantRecord> variants = new ArrayList<>();

        // 1. Load mapped_kmers matrix
        Logger.info(CLASS_NAME, "Loading mapped kmers matrix.");
        try (BufferedReader reader = getBufferedReader(kmersFile)) {
            String header = reader.readLine();
            if (header == null) {
                throw new IllegalArgumentException("Input kmers file is empty.");
            }

            String[] cols = header.split("\t");
            for (int i = 3; i < cols.length; i++) {
                sampleNames.add(cols[i]);
            }

            int numSamples = sampleNames.size();
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.trim().isEmpty()) continue;
                String[] parts = line.split("\t");

                VariantRecord var = new VariantRecord(parts[0], parts[1], parts[2], numSamples);
                for (int s = 0; s < numSamples; s++) {
                    var.rawCounts[s] = Double.parseDouble(parts[3 + s]);
                }
                variants.add(var);
            }
        }

        int numSamples = sampleNames.size();
        Logger.info(CLASS_NAME, "Loaded " + variants.size() + " variants across " + numSamples + " samples.");

        // 2. Compute column sums (total mapped kmers per sample)
        double[] colSums = new double[numSamples];
        for (VariantRecord var : variants) {
            for (int s = 0; s < numSamples; s++) {
                colSums[s] += var.rawCounts[s];
            }
        }

        // 3. Normalize RPB, compute summary stats, and calculate dosages
        Logger.info(CLASS_NAME, "Computing RPB normalization and dosage genotypes...");
        List<VariantRecord> filteredVariants = new ArrayList<>();

        for (VariantRecord var : variants) {
            double sumRpb = 0.0;
            for (int s = 0; s < numSamples; s++) {
                var.rpbValues[s] = (colSums[s] > 0) ? (var.rawCounts[s] / colSums[s]) * 1e9 : 0.0;
                sumRpb += var.rpbValues[s];
            }

            var.meanRpb = sumRpb / numSamples;
            double varianceSum = 0.0;
            for (int s = 0; s < numSamples; s++) {
                varianceSum += Math.pow(var.rpbValues[s] - var.meanRpb, 2);
            }
            var.sdRpb = Math.sqrt(varianceSum / numSamples);
            var.cv = (var.meanRpb > 0) ? (var.sdRpb / var.meanRpb) : 0.0;

            // Discrete dosage rounding (0, 1, 2)
            if (var.cv >= cvCutoff) {
                double minRpb = Double.MAX_VALUE;
                double maxRpb = -Double.MAX_VALUE;

                for (double r : var.rpbValues) {
                    if (r < minRpb) minRpb = r;
                    if (r > maxRpb) maxRpb = r;
                }

                for (int s = 0; s < numSamples; s++) {
                    if (maxRpb > minRpb) {
                        double scaled = 2.0 * (var.rpbValues[s] - minRpb) / (maxRpb - minRpb);
                        int dosage = (int) Math.round(scaled);
                        if (dosage < 0) dosage = 0;
                        if (dosage > 2) dosage = 2;
                        var.dosages[s] = dosage;
                    } else {
                        var.dosages[s] = 0;
                    }
                }
                filteredVariants.add(var);
            }
        }

        // 4. Write output files
        Path covMatrixPath = Path.of(outPrefix + ".cov_matrix.tsv.gz");
        Path genotypePath  = Path.of(outPrefix + ".genotype_matrix.tsv.gz");

        // Output cov_matrix.tsv.gz
        Logger.info(CLASS_NAME, "Writing cov_matrix -> " + covMatrixPath);
        try (BufferedWriter writer = getBufferedWriter(covMatrixPath)) {
            writer.write("acc.id\tlen\ttot_kmers");
            for (String sample : sampleNames) {
                writer.write("\t" + sample);
            }
            writer.write("\tmean_rpb\tsd_rpb\tcv\n");

            for (VariantRecord var : variants) {
                writer.write(var.accId + "\t" + var.len + "\t" + var.totKmers);
                for (int s = 0; s < numSamples; s++) {
                    writer.write(String.format(Locale.US, "\t%.4f", var.rpbValues[s]));
                }
                writer.write(String.format(Locale.US, "\t%.4f\t%.4f\t%.4f\n", var.meanRpb, var.sdRpb, var.cv));
            }
        }

        // Output genotype_matrix.tsv.gz
        Logger.info(CLASS_NAME, "Writing genotype_matrix (0:2 rounded) -> " + genotypePath);
        try (BufferedWriter writer = getBufferedWriter(genotypePath)) {
            writer.write("acc.id");
            for (String sample : sampleNames) {
                writer.write("\t" + sample);
            }
            writer.write("\n");

            for (VariantRecord var : filteredVariants) {
                writer.write(var.accId);
                for (int s = 0; s < numSamples; s++) {
                    writer.write("\t" + var.dosages[s]);
                }
                writer.write("\n");
            }
        }

        Logger.info(CLASS_NAME, "Exported cov_matrix (" + variants.size() +
                " variants) and genotype_matrix (" + filteredVariants.size() + " variants with CV >= " + cvCutoff + ").");
        return 0;
    }

    @Override
    public void run() {
        try {
            int exitCode = call();
            if (exitCode != 0) System.exit(exitCode);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
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
}