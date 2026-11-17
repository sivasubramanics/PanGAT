package edu.tamu.pg.pang2p.Data;

import edu.tamu.pg.pang2p.Utils.Logger;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.zip.GZIPInputStream;

public class FastqReader implements Iterable<Fastq>, AutoCloseable {
    private final Path fastqPath;
    private FastqIterator iterator;

    private final String CLASS_NAME = this.getClass().getSimpleName();
    private static final int BUFFER_SIZE = 128 * 1024; // 128 KB optimized for HPC parallel filesystems (Lustre/GPFS)

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

    private InputStream openStream(Path path) throws IOException {
        InputStream fileStream = new FileInputStream(path.toFile());
        BufferedInputStream bufferedStream = new BufferedInputStream(fileStream, BUFFER_SIZE);

        boolean isGzip = path.toString().toLowerCase().endsWith(".gz");
        if (!isGzip) {
            bufferedStream.mark(2);
            int b1 = bufferedStream.read();
            int b2 = bufferedStream.read();
            bufferedStream.reset();
            isGzip = (b1 == 0x1f && b2 == 0x8b);
        }

        return isGzip
                ? new GZIPInputStream(bufferedStream, BUFFER_SIZE)
                : bufferedStream;
    }

    private class FastqIterator implements Iterator<Fastq>, AutoCloseable {
        private final InputStream inputStream;
        private final byte[] readBuffer = new byte[BUFFER_SIZE];
        private int bufferLength = 0;
        private int bufferOffset = 0;
        private boolean EOF = false;
        private int sequenceId = 0;

        private String nextHeader = null;
        private final ByteArrayOutputStream lineCache = new ByteArrayOutputStream(512);

        public FastqIterator() {
            InputStream tempStream = null;
            try {
                tempStream = openStream(fastqPath);
            } catch (IOException e) {
                Logger.error(CLASS_NAME, "Error opening file: " + fastqPath + " (" + e.getMessage() + ")");
            }
            this.inputStream = tempStream;
            try {
                advanceToNextHeader();
            } catch (IOException e) {
                Logger.error(CLASS_NAME, "Error initializing FASTQ iterator: " + e.getMessage());
                nextHeader = null;
            }
        }

        private String readLineString() throws IOException {
            lineCache.reset();
            if (EOF) return null;

            while (true) {
                if (bufferOffset >= bufferLength) {
                    bufferLength = inputStream.read(readBuffer);
                    if (bufferLength == -1) {
                        EOF = true;
                        break;
                    }
                    bufferOffset = 0;
                }

                int start = bufferOffset;
                while (bufferOffset < bufferLength) {
                    byte b = readBuffer[bufferOffset++];
                    if (b == '\n') {
                        int length = (bufferOffset - 1) - start;
                        if (length > 0 && readBuffer[start + length - 1] == '\r') {
                            length--;
                        }
                        lineCache.write(readBuffer, start, length);
                        return lineCache.toString(StandardCharsets.US_ASCII);
                    }
                }
                lineCache.write(readBuffer, start, bufferLength - start);
            }

            return lineCache.size() > 0 ? lineCache.toString(StandardCharsets.US_ASCII) : null;
        }

        private void advanceToNextHeader() throws IOException {
            String line;
            while ((line = readLineString()) != null) {
                line = line.trim();
                if (!line.isEmpty()) {
                    if (line.startsWith("@")) {
                        nextHeader = line;
                        return;
                    }
                }
            }
            nextHeader = null;
        }

        @Override
        public boolean hasNext() {
            return nextHeader != null;
        }

        @Override
        public Fastq next() {
            if (!hasNext()) {
                throw new NoSuchElementException();
            }

            String currentHeader = nextHeader;
            String rawHeader = currentHeader.substring(1);
            int spaceIdx = rawHeader.indexOf(' ');
            String name = (spaceIdx == -1) ? rawHeader : rawHeader.substring(0, spaceIdx);
            String description = (spaceIdx == -1) ? "" : rawHeader.substring(spaceIdx + 1).trim();

            StringBuilder seqBuilder = new StringBuilder(150);
            StringBuilder qualBuilder = new StringBuilder(150);

            try {
                // 1. Read sequence line(s) until '+' separator
                String line;
                while ((line = readLineString()) != null) {
                    if (line.startsWith("+")) {
                        break;
                    }
                    seqBuilder.append(line.trim());
                }

                if (line == null) {
                    throw new IllegalStateException("Unexpected EOF before '+' separator for record: " + name);
                }

                // 2. Read quality line(s) until length matches sequence length
                int targetLen = seqBuilder.length();
                while (qualBuilder.length() < targetLen && (line = readLineString()) != null) {
                    qualBuilder.append(line.trim());
                }

                if (qualBuilder.length() != targetLen) {
                    throw new IllegalStateException("Corrupt FASTQ record length mismatch: " + name);
                }

                // 3. Preload the next header for the subsequent hasNext() check
                advanceToNextHeader();

            } catch (IOException e) {
                Logger.error(CLASS_NAME, "Error parsing FASTQ record: " + e.getMessage());
                nextHeader = null;
                throw new RuntimeException(e);
            }

            return new Fastq(sequenceId++, name, seqBuilder.toString(), qualBuilder.toString(), description);
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("Remove operation is not supported.");
        }

        @Override
        public void close() throws IOException {
            if (inputStream != null) {
                inputStream.close();
            }
        }
    }
}