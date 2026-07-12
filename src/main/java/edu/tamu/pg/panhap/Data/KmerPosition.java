package edu.tamu.pg.panhap.Data;

public class KmerPosition {

    private final int startPos;
    private final Kmer kmer;

    public KmerPosition(int startPos, Kmer kmer) {
        this.startPos = startPos;
        this.kmer = kmer;
    }

    public int getStartPos() {
        return startPos;
    }

    public Kmer getKmer() {
        return kmer;
    }
}