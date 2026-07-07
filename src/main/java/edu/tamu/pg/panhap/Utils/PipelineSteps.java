package edu.tamu.pg.panhap.Utils;

public enum PipelineSteps {
    PREP_GENOME("00.genomes"),
    SINGLE_COPY("01.single_copy_regions");

    private final String dirName;

    PipelineSteps(String dirName) {
        this.dirName = dirName;
    }

    public String dirName(){
        return dirName;
    }
}
