package edu.tamu.pg.pang2p.Data;

public class Read {
    private final Fasta fasta; // Underlying wrapper (Fastq can be represented or converted)
    private final String quality; // Null if it's a FASTA read
    private final boolean isFastq;

    // Constructor for FASTA
    public Read(Fasta fasta) {
        this.fasta = fasta;
        this.quality = null;
        this.isFastq = false;
    }

    // Constructor for FASTQ
    public Read(Fastq fastq) {
        this.fasta = fastq.toFasta();
        this.quality = fastq.getQuality();
        this.isFastq = true;
    }

    public int getId() {
        return fasta.getId();
    }

    public String getName() {
        return fasta.getName();
    }

    public String getDescription() {
        return fasta.getDescription();
    }

    public String getSequence() {
        return fasta.getSequence();
    }

    public boolean hasQuality() {
        return isFastq;
    }

    public String getQuality() {
        return quality;
    }

    public int length() {
        return fasta.getLength();
    }

    public Fasta toFasta() {
        return fasta;
    }

    /**
     * Streamlines kmer generation directly from the unified Read object
     */
    public KmerIterator kmers(int kmerLength, int prefixLength, boolean useCanonical) {
        return fasta.kmers(kmerLength, prefixLength, useCanonical);
    }

    @Override
    public String toString() {
        return fasta.toString();
    }
}