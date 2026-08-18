package edu.tamu.pg.pang2p.Data;

import edu.tamu.pg.pang2p.Utils.FileUtils;
import edu.tamu.pg.pang2p.Utils.Logger;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Iterator;
import java.util.NoSuchElementException;

public class FastaReader implements Iterable<Fasta>, AutoCloseable {
    private final Path fastaPath;
    private FastaIterator iterator;

    private final String CLASS_NAME = this.getClass().getSimpleName();

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

    private class FastaIterator implements Iterator<Fasta>, AutoCloseable {
        private String currentLine;
        private int sequence_id;
        private BufferedReader reader;
        private StringBuilder sequenceBuilder;

        public FastaIterator() {
            try {
                reader = FileUtils.reader(fastaPath);
                currentLine = reader.readLine();
                sequence_id = 0;
                sequenceBuilder = new StringBuilder();
            } catch (Exception e) {
                Logger.error(CLASS_NAME, "Error opening file: " + fastaPath);
            }
        }

        @Override
        public boolean hasNext() {
            return currentLine != null;
        }

        @Override
        public Fasta next() {
            if (!hasNext()) {
                throw new NoSuchElementException();
            }

            // Parse the header line (assume starts with ">")
            String header = currentLine.startsWith(">") ? currentLine.substring(1) : currentLine;
            String name = header.split("\\s+")[0];
            String description = header.length() > name.length() ? header.substring(name.length()).trim() : "";
            sequenceBuilder.setLength(0);

            try {
                // Read sequence lines until the next header or end of file
                while ((currentLine = reader.readLine()) != null && !currentLine.startsWith(">")) {
                    sequenceBuilder.append(currentLine.trim());
                }
            } catch (IOException e) {
                Logger.error(CLASS_NAME, "Error reading sequence data: " + e.getMessage());
            }

            return new Fasta(sequence_id++, name, sequenceBuilder.toString(), description);
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