package edu.tamu.pg.pang2p.Data;

import edu.tamu.pg.pang2p.Utils.Logger;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Iterator;
import java.util.NoSuchElementException;

public class FastaReader implements Iterable<Fasta>, AutoCloseable {
    private final String fastaFile;
    private FastaIterator iterator;

    private final String CLASS_NAME = this.getClass().getSimpleName();

    public FastaReader(Path fastaPath) {
        this(fastaPath.toString());
    }

    public FastaReader(String fastaFile) {
        this.fastaFile = fastaFile;
        Logger.info(CLASS_NAME, "Reading fasta file: " + fastaFile);
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
        private BufferedReader reader;
        private String currentLine;
        private StringBuilder sequenceBuilder;
        private int sequence_id;

        public FastaIterator() {
            try {
                reader = new BufferedReader(new FileReader(fastaFile));
                currentLine = reader.readLine();
                sequence_id = 0;
                sequenceBuilder = new StringBuilder();
            } catch (IOException e) {
                Logger.error(CLASS_NAME, "Error opening file: " + fastaFile);
                throw new RuntimeException("Error opening file: " + fastaFile, e);
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
            String header = currentLine.substring(1);
            String name = header.split(" ")[0];
            String description = header.substring(name.length()).trim();
            sequenceBuilder.setLength(0);

            try {
                // Read sequence lines until the next header or end of file
                while ((currentLine = reader.readLine()) != null && !currentLine.startsWith(">")) {
                    sequenceBuilder.append(currentLine.trim());
                }
            } catch (IOException e) {
                throw new RuntimeException("Error reading file: " + fastaFile, e);
            }

            // Return a new Fasta object (assuming Fasta class implements Sequence)
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
