package edu.tamu.pg.pang2p.Data;

import edu.tamu.pg.pang2p.Utils.Logger;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Iterator;
import java.util.NoSuchElementException;

public class SeqReader implements Iterable<Read>, AutoCloseable {
    private final Path filePath;
    private final String format; // "fasta" or "fastq"
    private SeqIterator iterator;

    private final String CLASS_NAME = this.getClass().getSimpleName();

    public SeqReader(Path filePath, String format) {
        this.filePath = filePath;
        this.format = format.toLowerCase();
        Logger.info(CLASS_NAME, "Reading " + this.format + " file: " + filePath);
    }

    public SeqReader(String filePath, String format) {
        this(Path.of(filePath), format);
    }

    @Override
    public Iterator<Read> iterator() {
        if (iterator == null) {
            iterator = new SeqIterator();
        }
        return iterator;
    }

    @Override
    public void close() throws IOException {
        if (iterator != null) {
            iterator.close();
        }
    }

    private class SeqIterator implements Iterator<Read>, AutoCloseable {
        private FastaReader fastaReader;
        private FastqReader fastqReader;
        private Iterator<Fasta> fastaIterator;
        private Iterator<Fastq> fastqIterator;
        private boolean isFastqFormat;

        public SeqIterator() {
            try {
                if (format.contains("fastq") || format.equals("fq")) {
                    isFastqFormat = true;
                    fastqReader = new FastqReader(filePath);
                    fastqIterator = fastqReader.iterator();
                } else if (format.contains("fasta") || format.equals("fa") || format.equals("fna")) {
                    isFastqFormat = false;
                    fastaReader = new FastaReader(filePath);
                    fastaIterator = fastaReader.iterator();
                } else {
                    Logger.error(CLASS_NAME, "Unknown sequence format specified: " + format);
                    isFastqFormat = false;
                }
            } catch (Exception e) {
                Logger.error(CLASS_NAME, "Error initializing reader for: " + filePath + " - " + e.getMessage());
            }
        }

        @Override
        public boolean hasNext() {
            if (isFastqFormat) {
                return fastqIterator != null && fastqIterator.hasNext();
            } else {
                return fastaIterator != null && fastaIterator.hasNext();
            }
        }

        @Override
        public Read next() {
            if (!hasNext()) {
                throw new NoSuchElementException();
            }

            if (isFastqFormat) {
                Fastq fq = fastqIterator.next();
                return new Read(fq);
            } else {
                Fasta fa = fastaIterator.next();
                return new Read(fa);
            }
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("Remove operation is not supported.");
        }

        @Override
        public void close() throws IOException {
            // Close the underlying resource-owning readers safely
            if (fastqReader != null) {
                fastqReader.close();
            }
            if (fastaReader != null) {
                fastaReader.close();
            }
        }
    }
}