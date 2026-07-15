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

    private Defaults() {}
}
