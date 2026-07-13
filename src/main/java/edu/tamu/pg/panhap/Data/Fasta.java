package edu.tamu.pg.panhap.Data;

import java.io.BufferedWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;
import java.util.Iterator;
import java.util.NoSuchElementException;

import edu.tamu.pg.panhap.Utils.HelperFunctions;
import edu.tamu.pg.panhap.Utils.Defaults;

/***
 * This class represents a Fasta object with thread-safe operations.
 */
public class Fasta implements Serializable {
    // sync lock for thread safety
    private final Object lock = new Object();

    private int id;
    private String name;
    private String sequence;
    private String description;

    // no-arg constructor
    public Fasta() {
    }

    public Fasta(int id, String name, String sequence) {
        this.id = id;
        this.name = name;
        this.sequence = sequence;
        this.description = "";
    }

    public Fasta(int id, String name, String sequence, String description) {
        this.id = id;
        this.name = name;
        this.sequence = sequence;
        this.description = description;
    }

    /***
     * This method will write the fasta file to the given writer
     */
    public void writeFasta(BufferedWriter writer) {
        synchronized (lock) {
            try {
                writer.write(">" + name + " " + description + "\n");
                writer.write(HelperFunctions.fold_seq(sequence, Defaults.FASTA_FOLD_LENGTH));
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }



    /***
     * Get fasta formatted string for the sequence
     */
    public String toFastaString(boolean fold) {
        synchronized (lock) {
            if (fold) {
                return ">" + name + " " + description + "\n" + HelperFunctions.fold_seq(sequence, Defaults.FASTA_FOLD_LENGTH);
            } else {
                return ">" + name + " " + description + "\n" + sequence + "\n";
            }
        }
    }

    /***
     * Get the length of the sequence
     */
    public int getLength() {
        synchronized (lock) {
            return sequence.length();
        }
    }
    /***
     * Convert the sequence to uppercase
     */
    public String getUppercaseSequence() {
        synchronized (lock) {
            return sequence.toUpperCase();
        }
    }

    /**
     * Returns a lazy iterator over all valid kmers in this FASTA sequence.
     *
     * Invalid bases (anything other than A/C/G/T) break kmer generation.
     * The returned iterator works on a snapshot of the sequence, so iteration
     * is safe even if the Fasta object is modified later.
     *
     * @param kmerLength     length of each kmer
     * @param prefixLength   prefix length used by Kmer encoding
     * @param useCanonical   whether to use canonical strand representation
     * @return KmerIterator
     */
    public KmerIterator kmers(int kmerLength,
                              int prefixLength,
                              boolean useCanonical) {

        final String sequenceSnapshot;

        synchronized (lock) {
            sequenceSnapshot = this.sequence;
        }

        return new KmerIterator(
                sequenceSnapshot,
                kmerLength,
                prefixLength,
                useCanonical);
    }

    /***
     * This method will return the kmers of the sequence as a list (NOT UNIQUE KMERS)
     */
    public List<Kmer> getKmersList(int kmerLength, int prefixLength, boolean useCanonical) {
        synchronized (lock) {
            List<Kmer> mers = new ArrayList<>();
            char[] kmerChars = new char[kmerLength];
            Kmer kmer = null;

            int validStart = 0;
            for (int i = 0; i < sequence.length(); i++) {
                char base = Character.toUpperCase(sequence.charAt(i));
                if (!isValidBase(base)) {
                    validStart = i + 1;
                    kmer = null;
                    kmerChars = new char[kmerLength];
                    continue;
                }

                int offset = i - validStart;
                if (offset < kmerLength) {
                    kmerChars[offset] = base;
                    if (offset == kmerLength - 1) {
                        kmer = new Kmer(kmerChars, prefixLength, useCanonical);
                    }
                } else if (kmer != null) {
                    // TODO:
                    // below line is to be implemented in a better way to improve the performance
//                    kmer = kmer.insertBase(base);
                    System.arraycopy(kmerChars, 1, kmerChars, 0, kmerLength - 1);
                    kmerChars[kmerLength - 1] = base;
                    kmer = new Kmer(kmerChars, prefixLength, useCanonical);
                }

                if (kmer != null) {
                    mers.add(kmer);
                }
            }
            return mers;
        }
    }

    /***
     * This method will return true if the base is valid
     */
    private boolean isValidBase(char base) {
        return base == 'A' || base == 'C' || base == 'G' || base == 'T';
    }

    /***
     * Get number of effective ATGC bases in the sequence that can be used to extract kmers
     * if < kmerLength ATGC is flagged with N, then it will skip
     */
    public int getEffectiveATGCCount(int kmerLength) {
        synchronized (lock) {
            int count = 0;
            int stretchLength = 0; // Track the length of consecutive ATGC stretch

            for (int i = 0; i < sequence.length(); i++) {
                char base = Character.toUpperCase(sequence.charAt(i));

                // check if the base is one of A, T, G, C
                if (base == 'A' || base == 'T' || base == 'G' || base == 'C') {
                    stretchLength++;
                } else {
                    // end of a stretch: check if it's long enough
                    if (stretchLength >= kmerLength) {
                        count += stretchLength;
                    }
                    stretchLength = 0;
                }
            }

            // handle the last stretch, if the sequence ends with ATGC
            if (stretchLength >= kmerLength) {
                count += stretchLength;
            }

            return count;
        }
    }

    @Override
    public String toString() {
        synchronized (lock) {
            return ">" + name + " " + description + "\n" + HelperFunctions.fold_seq(sequence, Defaults.FASTA_FOLD_LENGTH);
        }
    }

    public int getId() {
        synchronized (lock) {
            return id;
        }
    }

    public String getName() {
        synchronized (lock) {
            return name;
        }
    }

    public String getSequence() {
        synchronized (lock) {
            return sequence;
        }
    }

    public String getDescription() {
        synchronized (lock) {
            return description;
        }
    }

    public void setName(String name) {
        synchronized (lock) {
            this.name = name;
        }
    }

    public void setSequence(String sequence) {
        synchronized (lock) {
            this.sequence = sequence;
        }
    }

    public void setId(int id) {
        synchronized (lock) {
            this.id = id;
        }
    }

    public void setDescription(String description) {
        synchronized (lock) {
            this.description = description;
        }
    }
}
//EOF