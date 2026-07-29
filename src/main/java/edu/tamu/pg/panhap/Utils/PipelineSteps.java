package edu.tamu.pg.panhap.Utils;

public enum PipelineSteps {
    GENOMES("00.genomes"),
    KMC("01.kmc"),
    SINGLE_COPY("02.single_copy"),
    GENOTYPING("03.genotyping");

    private final String dirName;

    PipelineSteps(String dirName) {
        this.dirName = dirName;
    }

    public String dirName(){
        return dirName;
    }
}
