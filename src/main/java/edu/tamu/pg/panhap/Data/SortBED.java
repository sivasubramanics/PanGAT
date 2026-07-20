package edu.tamu.pg.panhap.Data;

import edu.tamu.pg.panhap.Utils.Defaults;
import edu.tamu.pg.panhap.Utils.Logger;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.ForkJoinPool;

public class SortBED {
    private final Path bedFile;
    private final String CLASS_NAME = SortBED.class.getSimpleName();

    public SortBED(Path bedFile) {
        this.bedFile = bedFile;
    }

    static class BEDEntry {
        String line;
        String chr;
        int start;
        int end;

        BEDEntry(String line) {
            this.line = line;
            String[] f = line.split("\t");
            if (f.length < 3)
                throw new IllegalArgumentException("BED line must have at least 3 columns: " + line);
            this.chr = f[0];
            this.start = Integer.parseInt(f[1]);
            this.end = Integer.parseInt(f[2]);
        }
    }

    public void sort() throws IOException {
        sort(Defaults.THREADS);
    }

    public void sort(int nThreads) throws IOException {

        Logger.info(CLASS_NAME, "Sorting BED file: " + bedFile + " using " + nThreads + " threads.");

        List<String> header = new ArrayList<>();
        List<BEDEntry> entries = new ArrayList<>();
        // Read the file: preserve header lines (starting with #), parse BED lines
        try (BufferedReader br = Files.newBufferedReader(bedFile)) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.trim().isEmpty()) continue;
                if (line.startsWith("#")) {
                    header.add(line);
                } else {
                    entries.add(new BEDEntry(line));
                }
            }
        }

        // Sorting
        Comparator<BEDEntry> bedComp = Comparator
                .comparing((BEDEntry e) -> e.chr, Comparator.naturalOrder())
                .thenComparingInt(e -> e.start)
                .thenComparingInt(e -> e.end);

        BEDEntry[] arr = entries.toArray(new BEDEntry[0]);
        ForkJoinPool customThreadPool = new ForkJoinPool(nThreads);
        try {
            customThreadPool.submit(() -> Arrays.parallelSort(arr, bedComp)).get();
        } catch (Exception e) {
            throw new RuntimeException(e);
        } finally {
            customThreadPool.shutdown();
        }

        // Write results: temp file and replace original atomically
        Path tmp = Files.createTempFile(bedFile.getParent(), "sortbed_", ".tmp");
        try (BufferedWriter bw = Files.newBufferedWriter(tmp)) {
            for (String h : header) {
                bw.write(h); bw.newLine();
            }
            for (BEDEntry be : arr) {
                bw.write(be.line);
                bw.newLine();
            }
        }
        Files.move(tmp, bedFile, StandardCopyOption.REPLACE_EXISTING);
    }
}