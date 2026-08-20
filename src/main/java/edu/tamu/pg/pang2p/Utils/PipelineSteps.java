package edu.tamu.pg.pang2p.Utils;

public enum PipelineSteps {
    GENOMES("PrepReference/genomes"),
    KMC("PrepReference/kmc"),
    SINGLE_COPY("PrepReference/HapCatalog"),
    GENOTYPING("Genotyping");

    private final String dirName;

    PipelineSteps(String dirName) {
        this.dirName = dirName;
    }

    public String dirName(){
        return dirName;
    }
}
