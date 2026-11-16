package edu.tamu.pg.pang2p.Tasks;

import picocli.CommandLine.*;

import java.util.concurrent.Callable;

@Command(name = "counts-reads", description = "Count reads to each anchors based on the kmers", mixinStandardHelpOptions = true)
public class ReadCounts implements Callable<Integer>, Runnable {
    @Option(names = {"-r", "--reference"}, description = "Input reference fasta file path", required = true)
    private String refFastaFile;

    @Option(names = {"-k", "--kmer"}, description = "Kmer size", required = true)
    private int kmerSize;

    @Option(names = {"-i", "--input"}, description = "Input reads file path", required = true)
    private String inputReadsFile;

    @Option(names = {"-f", "--format"}, description = "Input reads file format (fastq or fasta or list)", required = true)
    private String inputReadsFormat;

    @Option(names = {"-t", "--threads"}, description = "Number of threads to use", required = false)
    private int nthreads = 1;

    @Option(names = {"-o", "--output"}, description = "Output file path", required = true)
    private String outputFile;

    @Override
    public void run() {
        try {
            call();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Override
    public Integer call() throws Exception {
        return 0;
    }
}
