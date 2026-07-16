package edu.tamu.pg.panhap.Utils;

import java.util.List;

public final class Defaults {
    public static final int FASTA_FOLD_LENGTH = 80;
    public static final int SIGNATURE_LENGTH = 9;
    public static final List<String> TOOLS = List.of(
            "kmc",
            "mmseqs"
        );
    public static final int CHROM_CODE_LENGTH = 2;
    public static final int K_SIZE = 31;
    public static final int MIN_COUNT = 1;
    public static final double HEAPSIZE = Runtime.getRuntime().maxMemory() / (1024.0 * 1024.0 * 1024.0); // in GB
    public static final int COUNTER_DUMP_SIZE = 10_000;
    public static final int MIN_LENGTH = 100;
    public static final double MIN_SEQ_IDENTITY = 0.9;
    public static final double MIN_COVERAGE = 0.8;
    public static final int COV_MODE = 1;


    private Defaults() {}
}
