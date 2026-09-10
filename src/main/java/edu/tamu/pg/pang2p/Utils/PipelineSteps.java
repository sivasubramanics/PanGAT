package edu.tamu.pg.pang2p.Utils;

public enum PipelineSteps {
    GENOMES("00_prep_ref/00_genomes"),
    ANCHORS("00_prep_ref/01_pan_anchors"),
    GENOTYPING("01_genotyping");

    private final String dirName;

    PipelineSteps(String dirName) {
        this.dirName = dirName;
    }

    public String dirName(){
        return dirName;
    }
}
