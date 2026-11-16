package edu.tamu.pg.pang2p.Data;

import edu.tamu.pg.pang2p.Utils.Logger;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.zip.GZIPInputStream;

public class FastqReader implements Iterable<Fastq>, AutoCloseable {
    private final Path fastqPath;
    private FastqIterator iterator;

    private final String CLASS_NAME = this.getClass().getSimpleName();
    private static final int BUFFER_SIZE = 65536; // 64 KB buffer for efficient HPC/GPFS reads

    public FastqReader(Path fastqPath) {
        this.fastqPath = fastqPath;
        Logger.info(CLASS_NAME, "Reading fastq file: " + fastqPath);
    }

    public FastqReader(String fastqFile) {
        this(Path.of(fastqFile));
    }

    @Override
    public Iterator<Fastq> iterator() {
        if (iterator == null) {
            iterator = new FastqIterator();
        }
        return iterator;
    }

    @Override
    public void close() throws IOException {
        if (iterator != null) {
            iterator.close();
        }
    }

    private BufferedReader openReader(Path path) throws IOException {
        InputStream fileStream = new FileInputStream(path.toFile());
        InputStream bufferedStream = new BufferedInputStream(fileStream, BUFFER_SIZE);

        // Detect gzip by extension OR magic bytes (0x1f 0x8b)
        boolean isGzip = path.toString().toLowerCase().endsWith(".gz");
        if (!isGzip) {
            bufferedStream.mark(2);
            int b1 = bufferedStream.read();
            int b2 = bufferedStream.read();
            bufferedStream.reset();
            isGzip = (b1 == 0x1f && b2 == 0x8b);
        }

        InputStream finalStream = isGzip
                ? new GZIPInputStream(bufferedStream, BUFFER_SIZE)
                : bufferedStream;

        return new BufferedReader(new InputStreamReader(finalStream, StandardCharsets.UTF_8), BUFFER_SIZE);
    }

    private class FastqIterator implements Iterator<Fastq>, AutoCloseable {
        private String currentHeader;
        private int sequenceId;
        private BufferedReader reader;
        private final StringBuilder seqBuilder = new StringBuilder();
        private final StringBuilder qualBuilder = new StringBuilder();

        public FastqIterator() {
            try {
                reader = openReader(fastqPath);
                currentHeader = nextNonEmptyLine();
                sequenceId = 0;
            } catch (IOException e) {
                Logger.error(CLASS_NAME, "Error opening file: " + fastqPath + " (" + e.getMessage() + ")");
            }
        }

        private String nextNonEmptyLine() throws IOException {
            String line;
            while ((line = reader.readLine()) != null) {
                if (!line.isEmpty()) {
                    return line;
                }
            }
            return null;
        }

        @Override
        public boolean hasNext() {
            return currentHeader != null;
        }

        @Override
        public Fastq next() {
            if (!hasNext()) {
                throw new NoSuchElementException();
            }

            if (!currentHeader.startsWith("@")) {
                Logger.error(CLASS_NAME, "Invalid FASTQ header (expected '@'): " + currentHeader);
                throw new IllegalStateException("Invalid FASTQ header: " + currentHeader);
            }

            // Parse @name description (matching FastaReader header splitting)
            String rawHeader = currentHeader.substring(1);
            int spaceIdx = -1;
            for (int i = 0; i < rawHeader.length(); i++) {
                if (Character.isWhitespace(rawHeader.charAt(i))) {
                    spaceIdx = i;
                    break;
                }
            }
            String name = (spaceIdx == -1) ? rawHeader : rawHeader.substring(0, spaceIdx);
            String description = (spaceIdx == -1) ? "" : rawHeader.substring(spaceIdx + 1).trim();

            seqBuilder.setLength(0);
            qualBuilder.setLength(0);

            try {
                // 1. Read sequence line(s) until '+' separator (supports both standard 4-line and wrapped FASTQ)
                String line;
                while ((line = reader.readLine()) != null && !line.startsWith("+")) {
                    seqBuilder.append(line.trim());
                }

                if (line == null) {
                    Logger.error(CLASS_NAME, "Unexpected EOF before '+' separator for record: " + name);
                    throw new IllegalStateException("Truncated FASTQ record: " + name);
                }

                // 2. Read quality line(s) until quality length matches sequence length
                //    (Critical: quality strings can start with '@', so we match length rather than checking startsWith("@"))
                int targetLen = seqBuilder.length();
                while (qualBuilder.length() < targetLen && (line = reader.readLine()) != null) {
                    qualBuilder.append(line.trim());
                }

                if (qualBuilder.length() != targetLen) {
                    Logger.error(CLASS_NAME, "Sequence/quality length mismatch for record: " + name
                            + " (seq=" + targetLen + ", qual=" + qualBuilder.length() + ")");
                    throw new IllegalStateException("Corrupt FASTQ record (length mismatch): " + name);
                }

                // 3. Advance to next record header
                currentHeader = nextNonEmptyLine();

            } catch (IOException e) {
                Logger.error(CLASS_NAME, "Error reading FASTQ data: " + e.getMessage());
                currentHeader = null;
            }

            return new Fastq(sequenceId++, name, seqBuilder.toString(), qualBuilder.toString(), description);
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("Remove operation is not supported.");
        }

        @Override
        public void close() throws IOException {
            if (reader != null) {
                reader.close();
            }
        }
    }
}