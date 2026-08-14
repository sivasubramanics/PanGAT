package edu.tamu.pg.pang2p.Data;

import org.jetbrains.annotations.NotNull;

/***
 * Class representing an entry in a FASTA index file.
 * Contains information about the sequence name, offset, length, line bases, line width, and sequence ID.
 */
public class FastaIndexEntry implements Comparable<FastaIndexEntry> {
    private final String sequenceName;
    private final long offset;
    private final int length;
    private final int lineBases;
    private final int lineWidth;
    private final int seqId;
    private String description;

    public FastaIndexEntry(int seqId, String sequenceName, int length, long offset, int lineBases, int lineWidth, String description) {
        this.seqId = seqId;
        this.sequenceName = sequenceName;
        this.offset = offset;
        this.length = length;
        this.lineBases = lineBases;
        this.lineWidth = lineWidth;
        this.description = description;
    }

    public FastaIndexEntry(int id, String name, int length, long offset, int lineBases, int lineWidth) {
        this(id, name, length, offset, lineBases, lineWidth, null); // Default: description = null
    }

    public String getSequenceName() {
        return sequenceName;
    }

    public long getOffset() {
        return offset;
    }

    public int getLength() {
        return length;
    }

    public int getLineBases() {
        return lineBases;
    }

    public int getLineWidth() {
        return lineWidth;
    }

    public void setDescription(String desc) { this.description = desc; }
    public String getDescription() { return description; }

    @Override
    public String toString() {
        return sequenceName + "\t" + length + "\t" + offset + "\t" + lineBases + "\t" + lineWidth;
    }

    @Override
    public int compareTo(@NotNull FastaIndexEntry o) {
        return Integer.compare(this.seqId, o.seqId);
    }

    public int getSeqId() {
        return seqId;
    }
}
//EOF