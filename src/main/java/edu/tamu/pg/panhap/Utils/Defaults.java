package edu.tamu.pg.panhap.Utils;

import java.util.List;

public final class Defaults {
    public static final int FASTA_FOLD_LENGTH = 60;
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

    private Defaults() {}
}
