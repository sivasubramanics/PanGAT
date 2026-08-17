package edu.tamu.pg.pang2p.Data;

import java.util.Arrays;

public class HapData {

    private String clusterID;
    private String hapID;
    private String newGenomeID;
    private String newChromID;
    private String orgGenomeID;
    private String orgChromID;
    private int start;
    private int end;
    private int length;
    private int totKmers;

    // Primitive arrays indexed by sample position
    private int[] obsCounts;
    private int[] obsTotalCounts;

    public HapData(String clusterID, String hapID, String newGenomeID, String newChromID,
                   String orgGenomeID, String orgChromID, int start, int end, int length) {
        this.clusterID = clusterID;
        this.hapID = hapID;
        this.newGenomeID = newGenomeID;
        this.newChromID = newChromID;
        this.orgGenomeID = orgGenomeID;
        this.orgChromID = orgChromID;
        this.start = start;
        this.end = end;
        this.length = length;
        this.totKmers = 0;
    }

    public void initCounts(int nSamples) {
        this.obsCounts = new int[nSamples];
        this.obsTotalCounts = new int[nSamples];
    }

    public void addSample(int obsCount, int obsTotalCount) {
        int[] tmpObs = new int[this.obsCounts.length + 1];
        int[] tmpObsTotals = new int[this.obsTotalCounts.length + 1];
        for (int i = 0; i < this.obsCounts.length; i++) {
            tmpObs[i] = this.obsCounts[i];
            tmpObsTotals[i] = this.obsTotalCounts[i];
        }
        this.obsCounts = tmpObs;
        this.obsTotalCounts = tmpObsTotals;
    }

    // Getters and Setters
    public String getClusterID() { return clusterID; }
    public String getHapID() { return hapID; }
    public String getNewGenomeID() { return newGenomeID; }
    public String getNewChromID() { return newChromID; }
    public String getOrgGenomeID() { return orgGenomeID; }
    public String getOrgChromID() { return orgChromID; }
    public int getStart() { return start; }
    public int getEnd() { return end; }
    public int getLength() { return length; }
    public int getTotKmers() { return totKmers; }

    public void setClusterID(String clusterID) { this.clusterID = clusterID; }
    public void setHapID(String hapID) { this.hapID = hapID; }
    public void setNewGenomeID(String newGenomeID) { this.newGenomeID = newGenomeID; }
    public void setNewChromID(String newChromID) { this.newChromID = newChromID; }
    public void setOrgGenomeID(String orgGenomeID) { this.orgGenomeID = orgGenomeID; }
    public void setOrgChromID(String orgChromID) { this.orgChromID = orgChromID; }
    public void setStart(int start) { this.start = start; }
    public void setEnd(int end) { this.end = end; }
    public void setLength(int length) { this.length = length; }
    public void setTotKmers(int totKmers) { this.totKmers = totKmers; }
}