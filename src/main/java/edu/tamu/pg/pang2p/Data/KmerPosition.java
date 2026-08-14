package edu.tamu.pg.pang2p.Data;

import java.io.Serializable;

/**
 * Immutable pair of a k-mer and its start position within a sequence.
 */
public final class KmerPosition implements Serializable {

    private static final long serialVersionUID = 1L;

    private final int startPos;
    private final Kmer kmer;

    public KmerPosition(int startPos, Kmer kmer) {
        this.startPos = startPos;
        this.kmer = kmer;
    }

    /**
     * Returns the 0-based start position of this k-mer.
     */
    public int getStartPos() {
        return startPos;
    }

    /**
     * Returns the k-mer.
     */
    public Kmer getKmer() {
        return kmer;
    }

    @Override
    public String toString() {
        return "KmerPosition{" +
                "startPos=" + startPos +
                ", kmer=" + kmer +
                '}';
    }

    @Override
    public int hashCode() {
        int result = Integer.hashCode(startPos);
        result = 31 * result + (kmer != null ? kmer.hashCode() : 0);
        return result;
    }

    @Override
    public boolean equals(Object obj) {

        if (this == obj)
            return true;

        if (!(obj instanceof KmerPosition))
            return false;

        KmerPosition other = (KmerPosition) obj;

        if (startPos != other.startPos)
            return false;

        if (kmer == null)
            return other.kmer == null;

        return kmer.equals(other.kmer);
    }
}