package edu.tamu.pg.pang2p.Data;

import edu.tamu.pg.pang2p.Utils.Logger;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.zip.GZIPInputStream;

public class FastaReader implements Iterable<Fasta>, AutoCloseable {
    private final Path fastaPath;
    private FastaIterator iterator;

    private final String CLASS_NAME = this.getClass().getSimpleName();
    private static final int BUFFER_SIZE = 128 * 1024; // 128 KB optimized for HPC parallel filesystems (Lustre/GPFS)

    public FastaReader(Path fastaPath) {
        this.fastaPath = fastaPath;
        Logger.info(CLASS_NAME, "Reading fasta file: " + fastaPath);
    }

    public FastaReader(String fastaFile) {
        this(Path.of(fastaFile));
    }

    @Override
    public Iterator<Fasta> iterator() {
        if (iterator == null) {
            iterator = new FastaIterator();
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

    private class FastaIterator implements Iterator<Fasta>, AutoCloseable {
        private final InputStream inputStream;
        private final byte[] readBuffer = new byte[BUFFER_SIZE];
        private int bufferLength = 0;
        private int bufferOffset = 0;
        private boolean EOF = false;
        private int sequenceId = 0;

        private String nextHeader = null;
        private final ByteArrayOutputStream lineCache = new ByteArrayOutputStream(512);

        public FastaIterator() {
            InputStream tempStream = null;
            try {
                tempStream = openStream(fastaPath);
            } catch (IOException e) {
                Logger.error(CLASS_NAME, "Error opening file: " + fastaPath + " (" + e.getMessage() + ")");
            }
            this.inputStream = tempStream;
            try {
                advanceToNextHeader();
            } catch (IOException e) {
                Logger.error(CLASS_NAME, "Error initializing FASTA iterator: " + e.getMessage());
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
                    if (line.startsWith(">")) {
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
        public Fasta next() {
            if (!hasNext()) {
                throw new NoSuchElementException();
            }

            String currentHeader = nextHeader;
            String headerContent = currentHeader.startsWith(">") ? currentHeader.substring(1) : currentHeader;

            // Fast whitespace split without regex overhead
            int spaceIdx = -1;
            for (int i = 0; i < headerContent.length(); i++) {
                if (Character.isWhitespace(headerContent.charAt(i))) {
                    spaceIdx = i;
                    break;
                }
            }
            String name = (spaceIdx == -1) ? headerContent : headerContent.substring(0, spaceIdx);
            String description = (spaceIdx == -1) ? "" : headerContent.substring(spaceIdx + 1).trim();

            StringBuilder sequenceBuilder = new StringBuilder(1024);

            try {
                String line;
                while ((line = readLineString()) != null) {
                    line = line.trim();
                    if (line.startsWith(">")) {
                        nextHeader = line; // Cache next header for future hasNext()
                        return new Fasta(sequenceId++, name, sequenceBuilder.toString(), description);
                    }
                    if (!line.isEmpty()) {
                        sequenceBuilder.append(line);
                    }
                }

                // If we reach here, we hit EOF
                nextHeader = null;

            } catch (IOException e) {
                Logger.error(CLASS_NAME, "Error parsing FASTA record: " + e.getMessage());
                nextHeader = null;
                throw new RuntimeException(e);
            }

            return new Fasta(sequenceId++, name, sequenceBuilder.toString(), description);
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