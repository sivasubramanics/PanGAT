package edu.tamu.pg.pang2p.Data;

import edu.tamu.pg.pang2p.Utils.Defaults;

import java.util.Arrays;

public class HapData implements Comparable<HapData> {

    private String clusterID;
    private String hapID;
    private String newGenomeID;
    private String newChromID;
    private String orgGenomeID;
    private String orgChromID;
    private String simChromID;
    private int start;
    private int end;
    private int length;
    private int totKmers;
    private int uniqueKmers;
    private int globalUniqueKmers;

    private double minRpb;
    private double maxRpb;
    private double meanRPB;
    private double sdRPB;
    private double cv;

    private int[] genotype = new int[0];
    private double[] gtScaled = new double[0];
    private double[] normalizedCounts = new double[0];
    private int[] obsCounts = new int[0];
    private int[] obsTotalCounts = new int[0];

    public HapData(String clusterID, String hapID, String newGenomeID, String newChromID,
                   String orgGenomeID, String orgChromID, int start, int end, int length, int totKmers, int uniqueKmers, int globalUniqueKmers) {
        this.clusterID = clusterID;
        this.hapID = hapID;
        this.newGenomeID = newGenomeID;
        setNewChromID(newChromID);
        this.orgGenomeID = orgGenomeID;
        this.orgChromID = orgChromID;
        this.start = start;
        this.end = end;
        this.length = length;
        this.totKmers = totKmers;
        this.uniqueKmers = uniqueKmers;
        this.globalUniqueKmers = globalUniqueKmers;
    }

    public void initCounts(int nSamples) {
        this.genotype = new int[nSamples];
        Arrays.fill(this.genotype, -1);
        this.gtScaled = new double[nSamples];
        this.normalizedCounts = new double[nSamples];
        this.obsCounts = new int[nSamples];
        this.obsTotalCounts = new int[nSamples];
    }

    public synchronized void addSample(int obsCount, int obsTotalCount) {
        int currentLen = (this.obsCounts != null) ? this.obsCounts.length : 0;
        int newLen = currentLen + 1;

        // Resize all sample-indexed arrays
        this.obsCounts = (this.obsCounts != null) ? Arrays.copyOf(this.obsCounts, newLen) : new int[1];
        this.obsTotalCounts = (this.obsTotalCounts != null) ? Arrays.copyOf(this.obsTotalCounts, newLen) : new int[1];
        this.genotype = (this.genotype != null) ? Arrays.copyOf(this.genotype, newLen) : new int[1];
        this.gtScaled = (this.gtScaled != null) ? Arrays.copyOf(this.gtScaled, newLen) : new double[1];
        this.normalizedCounts = (this.normalizedCounts != null) ? Arrays.copyOf(this.normalizedCounts, newLen) : new double[1];
        this.obsCounts = (this.obsCounts != null) ? Arrays.copyOf(this.obsCounts, newLen) : new int[1];
        this.obsTotalCounts = (this.obsTotalCounts != null) ? Arrays.copyOf(this.obsTotalCounts, newLen) : new int[1];

        // Assign observed counts and initialize genotype to -1
        this.genotype[currentLen] = -1;
        this.gtScaled[currentLen] = 0.0;
        this.normalizedCounts[currentLen] = 0.0;
        this.obsCounts[currentLen] = obsCount;
        this.obsTotalCounts[currentLen] = obsTotalCount;
    }

    @Override
    public int compareTo(HapData other) {
        if (other == null) {
            return 1;
        }

        if (this.simChromID == null && other.simChromID != null) return -1;
        if (this.simChromID != null && other.simChromID == null) return 1;

        int chromComparison = 0;
        if (this.simChromID != null && other.simChromID != null) {
            chromComparison = this.simChromID.compareTo(other.simChromID);
        }

        if (chromComparison != 0) {
            return chromComparison;
        }

        return Integer.compare(this.start, other.start);
    }

    public void calcMetrics() {
        if (this.normalizedCounts == null || this.normalizedCounts.length == 0) {
            this.minRpb = 0.0;
            this.maxRpb = 0.0;
            this.meanRPB = 0.0;
            this.sdRPB = 0.0;
            this.cv = 0.0;
            this.gtScaled = new double[0];
            return;
        }

        int len = this.normalizedCounts.length;
        if (this.gtScaled == null || this.gtScaled.length != len) {
            this.gtScaled = new double[len];
        }

        double cumNormalizedCounts = 0.0;
        this.minRpb = this.normalizedCounts[0];
        this.maxRpb = this.normalizedCounts[0];

        for (int i = 0; i < len; i++) {
            double val = this.normalizedCounts[i];
            cumNormalizedCounts += val;
            if (val < minRpb) minRpb = val;
            if (val > maxRpb) maxRpb = val;
        }

        meanRPB = cumNormalizedCounts / len;

        double sumSquaredDiffs = 0.0;
        for (double val : this.normalizedCounts) {
            double diff = val - meanRPB;
            sumSquaredDiffs += diff * diff;
        }
        sdRPB = Math.sqrt(sumSquaredDiffs / len);

        cv = (meanRPB != 0.0) ? (sdRPB / meanRPB) : 0.0;

        double range = maxRpb - minRpb;
//        for (int i = 0; i < len; i++) {
//            gtScaled[i] = (range > 0.0) ? (2.0 * (this.normalizedCounts[i] - minRpb) / range) : 0.0;
//        }
        for (int i = 0; i < len; i++) {
            double scaled = (range > 0.0)
                    ? 2.0 * (normalizedCounts[i] - minRpb) / range
                    : 0.0;

            gtScaled[i] = scaled;
            genotype[i] = (int) Math.round(gtScaled[i]);
        }
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
    public int getUniqueKmers() { return uniqueKmers; }
    public int getGlobalUniqueKmers() { return globalUniqueKmers; }
    public int[] getObsCounts() { return obsCounts; }
    public int[] getObsTotalCounts() { return obsTotalCounts; }
    public int[] getGenotype() { return genotype; }
    public double[] getGtScaled() { return gtScaled; }
    public double[] getNormalizedCounts() { return normalizedCounts; }

    public void setClusterID(String clusterID) { this.clusterID = clusterID; }
    public void setHapID(String hapID) { this.hapID = hapID; }
    public void setNewGenomeID(String newGenomeID) { this.newGenomeID = newGenomeID; }

    public void setNewChromID(String newChromID) {
        this.newChromID = newChromID;
        if (newChromID != null && newChromID.length() >= Defaults.CHROM_CODE_LENGTH) {
            this.simChromID = newChromID.substring(Defaults.CHROM_CODE_LENGTH);
        } else {
            this.simChromID = newChromID;
        }
    }

    public void setOrgGenomeID(String orgGenomeID) { this.orgGenomeID = orgGenomeID; }
    public void setOrgChromID(String orgChromID) { this.orgChromID = orgChromID; }
    public void setStart(int start) { this.start = start; }
    public void setEnd(int end) { this.end = end; }
    public void setLength(int length) { this.length = length; }
    public void setTotKmers(int totKmers) { this.totKmers = totKmers; }
    public void setUniqueKmers(int uniqueKmers) { this.uniqueKmers = uniqueKmers; }
    public void setGlobalUniqueKmers(int globalUniqueKmers) { this.globalUniqueKmers = globalUniqueKmers; }
    public void setObsCounts(int[] obsCounts) { this.obsCounts = (obsCounts != null) ? obsCounts : new int[0]; }
    public void setObsTotalCounts(int[] obsTotalCounts) { this.obsTotalCounts = (obsTotalCounts != null) ? obsTotalCounts : new int[0]; }
    public void setGenotype(int[] genotype) { this.genotype = (genotype != null) ? genotype : new int[0]; }
    public void setGtScaled(double[] gtScaled) { this.gtScaled = (gtScaled != null) ? gtScaled : new double[0]; }
    public void setNormalizedCounts(double[] normalizedCounts) { this.normalizedCounts = (normalizedCounts != null) ? normalizedCounts : new double[0]; }

    private double roundTo2Decimals(double val) {
        return Math.round(val * 100.0) / 100.0;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append(hapID).append("\t");
        sb.append(orgGenomeID).append("\t");
        sb.append(orgChromID).append("\t");
        sb.append(simChromID).append("\t");
        sb.append(start).append("\t");
        sb.append(end).append("\t");
        sb.append(length).append("\t");
        sb.append(totKmers).append("\t");
        sb.append(uniqueKmers).append("\t");
        sb.append(globalUniqueKmers).append("\t");
        sb.append(roundTo2Decimals(minRpb)).append("\t");
        sb.append(roundTo2Decimals(maxRpb)).append("\t");
        sb.append(roundTo2Decimals(meanRPB)).append("\t");
        sb.append(roundTo2Decimals(sdRPB)).append("\t");
        sb.append(roundTo2Decimals(cv)).append("\t");

        int len = (this.normalizedCounts != null) ? this.normalizedCounts.length : 0;
        for (int i = 0; i < len; i++) {
            sb.append((i < genotype.length) ? this.genotype[i] : -1).append(":");

            double gtVal = (i < gtScaled.length) ? this.gtScaled[i] : 0.0;
            sb.append(roundTo2Decimals(gtVal)).append(":");
            sb.append(roundTo2Decimals(this.normalizedCounts[i])).append(":");
            sb.append((i < obsCounts.length) ? this.obsCounts[i] : 0).append(":");
            sb.append((i < obsTotalCounts.length) ? this.obsTotalCounts[i] : 0);

            if (i < len - 1) {
                sb.append("\t");
            }
        }
        return sb.toString();
    }
}