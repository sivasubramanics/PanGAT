package edu.tamu.pg.pang2p.Data;

import java.util.Iterator;
import java.util.NoSuchElementException;

public class KmerIterator implements Iterator<KmerPosition> {

    private final String sequence;
    private final int sequenceLength;
    private final int kmerLength;
    private final int prefixLength;
    private final boolean useCanonical;

    private final char[] kmerChars;

    private int index;
    private int validStart;

    private int validStretchLength;
    private KmerPosition next;

    public KmerIterator(String sequence, int kmerLength, int prefixLength, boolean useCanonical) {
        this.sequence = sequence;
        this.sequenceLength = sequence.length();
        this.kmerLength = kmerLength;
        this.prefixLength = prefixLength;
        this.useCanonical = useCanonical;
        this.kmerChars = new char[kmerLength];

        rewind();
    }

    @Override
    public boolean hasNext() {
        if (next == null) {
            next = computeNext();
        }
        return next != null;
    }

    @Override
    public KmerPosition next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        KmerPosition result = next;
        next = null;
        return result;
    }

    public void skipTo(int startPos) {
        if (startPos < 0) startPos = 0;
        if (startPos > sequenceLength) startPos = sequenceLength;

        this.index = startPos;
        this.validStart = startPos;
        this.validStretchLength = 0;
        this.next = null;
    }

    public void rewind() {
        skipTo(0);
    }

    public int getCurrentPosition() {
        return index;
    }

    public String getSequence() {
        return sequence;
    }

//    private KmerPosition computeNext() {
//        final int kLength = this.kmerLength;
//        final int seqLen = this.sequenceLength;
//        final char[] chars = this.kmerChars;
//
//        while (index < seqLen) {
//            char base = (char) (sequence.charAt(index) & 0xDF);
//
//            if (base != 'A' && base != 'C' && base != 'G' && base != 'T') {
//                validStart = index + 1;
//                validStretchLength = 0;
//                index++;
//                continue;
//            }
//
//            if (validStretchLength < kLength) {
//                chars[validStretchLength] = base;
//                validStretchLength++;
//
//                if (validStretchLength == kLength) {
//                    Kmer kmer = new Kmer(chars, prefixLength, useCanonical);
//                    KmerPosition result = new KmerPosition(validStart, kmer);
//                    index++;
//                    return result;
//                }
//            } else {
//                int writeIndex = (index - validStart) % kLength;
//                chars[writeIndex] = base;
//
//                char[] linearKmer = new char[kLength];
//                int part2 = kLength - 1 - writeIndex;
//
//                System.arraycopy(chars, writeIndex + 1, linearKmer, 0, part2);
//                System.arraycopy(chars, 0, linearKmer, part2, writeIndex + 1);
//
//                Kmer kmer = new Kmer(linearKmer, prefixLength, useCanonical);
//                KmerPosition result = new KmerPosition(index - kLength + 1, kmer);
//                index++;
//                return result;
//            }
//            index++;
//        }
//        return null;
//    }
    private KmerPosition computeNext() {
        final int kLength = this.kmerLength;
        final int seqLen = this.sequenceLength;
        final char[] chars = this.kmerChars;

        while (index < seqLen) {
            char base = (char) (sequence.charAt(index) & 0xDF);

            if (base != 'A' && base != 'C' && base != 'G' && base != 'T') {
                index++;
                validStart = index;
                validStretchLength = 0;
                continue;
            }

            if (validStretchLength < kLength) {
                chars[validStretchLength] = base;
                validStretchLength++;

                if (validStretchLength == kLength) {
                    Kmer kmer = new Kmer(chars, prefixLength, useCanonical);
                    KmerPosition result = new KmerPosition(validStart, kmer);
                    index++;
                    return result;
                }
            } else {
                int writeIndex = (index - validStart) % kLength;
                chars[writeIndex] = base;

                char[] linearKmer = new char[kLength];
                int part2 = kLength - 1 - writeIndex;

                System.arraycopy(chars, writeIndex + 1, linearKmer, 0, part2);
                System.arraycopy(chars, 0, linearKmer, part2, writeIndex + 1);

                Kmer kmer = new Kmer(linearKmer, prefixLength, useCanonical);
                KmerPosition result = new KmerPosition(index - kLength + 1, kmer);
                index++;
                return result;
            }
            index++;
        }
        return null;
    }
}