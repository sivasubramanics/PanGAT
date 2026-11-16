package edu.tamu.pg.pang2p.Data;

import java.util.Iterator;

public class Fastq {
    private final int id;
    private final String name;
    private final String description;
    private final String sequence;
    private final String quality;

    public Fastq(int id, String name, String sequence, String quality, String description) {
        this.id = id;
        this.name = name;
        this.sequence = sequence;
        this.quality = quality;
        this.description = description;
    }

    public int getId() {
        return id;
    }

    public String getName() {
        return name;
    }

    public String getDescription() {
        return description;
    }

    public String getSequence() {
        return sequence;
    }

    public String getQuality() {
        return quality;
    }

    public int length() {
        return sequence.length();
    }

    @Override
    public String toString() {
        String header = description.isEmpty() ? "@" + name : "@" + name + " " + description;
        return header + "\n" + sequence + "\n+\n" + quality;
    }

    public Fasta toFasta() {
        return new Fasta(id, name, sequence, description);
    }

    public KmerIterator kmers(int k, int prefixLength, boolean useCanonical) {
        return this.toFasta().kmers(k, prefixLength, useCanonical);
    }
}