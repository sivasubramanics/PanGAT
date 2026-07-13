package edu.tamu.pg.panhap.Data;

import java.util.Arrays;
import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Iterator over all valid kmers in a DNA sequence.
 *
 * Invalid bases (anything other than A/C/G/T) break a kmer stretch.
 * The iterator lazily generates kmers as they are requested.
 *
 * This iterator also supports skipping directly to an arbitrary sequence
 * position without iterating through all intermediate kmers.
 */
public class KmerIterator implements Iterator<KmerPosition> {

    private final String sequence;
    private final int sequenceLength;
    private final int kmerLength;
    private final int prefixLength;
    private final boolean useCanonical;

    private final char[] kmerChars;

    /**
     * Current sequence index being processed.
     */
    private int index;

    /**
     * Beginning of the current valid ATGC stretch.
     */
    private int validStart;

    /**
     * Cached next value for hasNext().
     */
    private KmerPosition next;

    public KmerIterator(String sequence,
                        int kmerLength,
                        int prefixLength,
                        boolean useCanonical) {

        this.sequence = sequence;
        this.sequenceLength = sequence.length();
        this.kmerLength = kmerLength;
        this.prefixLength = prefixLength;
        this.useCanonical = useCanonical;

        this.kmerChars = new char[kmerLength];

        this.index = 0;
        this.validStart = 0;
        this.next = null;
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

    /**
     * Skip directly to the requested sequence position.
     *
     * The next call to next() will produce the first valid kmer beginning
     * at or after this position.
     */
    public void skipTo(int startPos) {

        if (startPos < 0) {
            startPos = 0;
        }

        if (startPos > sequenceLength) {
            startPos = sequenceLength;
        }

        index = startPos;
        validStart = startPos;

        Arrays.fill(kmerChars, '\0');
        next = null;
    }

    /**
     * Restart iteration from the beginning.
     */
    public void rewind() {
        skipTo(0);
    }

    /**
     * Returns the current sequence index.
     */
    public int getCurrentPosition() {
        return index;
    }

    /**
     * Returns the underlying sequence.
     */
    public String getSequence() {
        return sequence;
    }

    private KmerPosition computeNext() {

        while (index < sequenceLength) {

            char base = Character.toUpperCase(sequence.charAt(index));

            if (!isValidBase(base)) {

                validStart = index + 1;

                Arrays.fill(kmerChars, '\0');

                index++;

                continue;
            }

            int offset = index - validStart;

            if (offset < kmerLength) {

                kmerChars[offset] = base;

                if (offset == kmerLength - 1) {

                    Kmer kmer =
                            new Kmer(kmerChars, prefixLength, useCanonical);

                    KmerPosition result =
                            new KmerPosition(validStart, kmer);

                    index++;

                    return result;
                }

            } else {

                System.arraycopy(
                        kmerChars,
                        1,
                        kmerChars,
                        0,
                        kmerLength - 1);

                kmerChars[kmerLength - 1] = base;

                Kmer kmer =
                        new Kmer(kmerChars, prefixLength, useCanonical);

                KmerPosition result =
                        new KmerPosition(index - kmerLength + 1, kmer);

                index++;

                return result;
            }

            index++;
        }

        return null;
    }

    private boolean isValidBase(char base) {

        return base == 'A'
                || base == 'C'
                || base == 'G'
                || base == 'T';
    }
}