package edu.tamu.pg.pang2p.Data;

import edu.tamu.pg.pang2p.Utils.Logger;
import edu.tamu.pg.pang2p.Utils.Defaults;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;

public class FastaWriter implements AutoCloseable {
    private final BufferedWriter writer;

    private final String CLASS_NAME = this.getClass().getSimpleName();

    /**
     * Creates a FASTA writer.
     *
     * @param fastaFile Output FASTA file.
     * @param append    Whether to append to an existing file.
     */
    public FastaWriter(String fastaFile, boolean append) {
        try {
            Logger.debug(CLASS_NAME, "Writing fasta file: " + fastaFile);
            writer = new BufferedWriter(new FileWriter(fastaFile, append));
        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Error opening file: " + fastaFile);
            throw new RuntimeException("Error opening file: " + fastaFile, e);
        }
    }

    public FastaWriter(Path fastaPath) {
        this(fastaPath.toString(), false);
    }

    public FastaWriter(Path fastaPath, boolean append) {
        this(fastaPath.toString(), append);
    }

    /**
     * Creates a FASTA writer that overwrites the file.
     */
    public FastaWriter(String fastaFile) {
        this(fastaFile, false);
    }

    public void write(Fasta fasta, boolean isFold){
        if (isFold) {
            this.write(fasta, Defaults.FASTA_FOLD_LENGTH);
        }
        else {
            this.write(fasta, fasta.getSequence().length());
        }
    }


    public void write(Fasta fasta) {
        write(fasta, Defaults.FASTA_FOLD_LENGTH);
    }

    /**
     * Writes a FASTA record.
     */
    public void write(Fasta fasta, int foldLenth){
        try {
            writer.write(">");
            writer.write(fasta.getName());

            if (fasta.getDescription() != null && !fasta.getDescription().isEmpty()) {
                writer.write(" ");
                writer.write(fasta.getDescription());
            }

            writer.newLine();

            // Write sequence (wrapped to 80 characters)
            String sequence = fasta.getSequence();

            for (int i = 0; i < sequence.length(); i += foldLenth) {
                int end = Math.min(i + foldLenth, sequence.length());
                writer.write(sequence, i, end - i);
                writer.newLine();
            }
        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Error writing FASTA record.");
            throw new RuntimeException("Error writing FASTA record.", e);
        }
    }

    /**
     * Flushes buffered data.
     */
    public void flush() {
        try {
            writer.flush();
        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Error flushing FASTA writer.");
            throw new RuntimeException("Error flushing FASTA writer.", e);
        }
    }

    @Override
    public void close() throws IOException {
        writer.close();
    }
}