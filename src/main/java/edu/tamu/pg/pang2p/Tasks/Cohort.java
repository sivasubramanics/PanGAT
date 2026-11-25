package edu.tamu.pg.pang2p.Tasks;

import edu.tamu.pg.pang2p.Data.PanHap;
import edu.tamu.pg.pang2p.Utils.FileUtils;
import edu.tamu.pg.pang2p.Utils.HelperFunctions;
import edu.tamu.pg.pang2p.Utils.Logger;
import picocli.CommandLine;
import picocli.CommandLine.*;

import java.nio.file.Path;
import java.util.concurrent.Callable;

@Command(name = "cohort", description = "Compile raw mapped k-mer counts across sample coverage files into cohort.pan_matrix.tsv.gz", sortOptions = false)
public class Cohort implements Runnable, Callable<Integer> {

    @Option(names = {"-r", "--rep"}, description = "Path to pan_clusters_rep.tsv file.", required = true)
    private Path panclustersRep;

    @Option(names = {"-f", "--fasta"}, description = "Path to pan_clusters_rep.fasta file.", required = true)
    private Path panclustersFasta;

    @Option(names = {"-i", "--input"}, description = "Path to input sample coverage files.", required = true, arity = "1..*", split = "\\s+")
    private Path[] input;

    @Option(names = {"-o", "--outprefix"}, description = "Prefix for output files.", required = true)
    private String outPrefix;

    @Option(names = {"-k", "--kmer"}, description = "K-mer size for counting (useful only for --fasta).", defaultValue = "51")
    private int kmerSize;

    @Option(names = {"--skip"}, description = "Skip size for k-mer counting (useful only for --fasta).", defaultValue = "1")
    private int skipSize;

    @Option(names = {"--threads"}, description = "Number of threads to use for k-mer counting (useful only for --fasta).", defaultValue = "1")
    private int threads;

    @Option(names = {"--force"}, description = "Force overwrite of existing mapped_kmers file instead of updating.", defaultValue = "false")
    private boolean force;

    private final String CLASS_NAME = Cohort.class.getSimpleName();
    private PanHap panHap;

    public Cohort() {
    }

    public Cohort(Path panclustersRep, Path panclustersFasta, int kmerSize, int skipSize, int threads, Path[] input, String[] samples, String outPrefix, boolean force) {
        this.panclustersRep = panclustersRep;
        this.panclustersFasta = panclustersFasta;
        this.kmerSize = kmerSize;
        this.skipSize = skipSize;
        this.threads = threads;
        this.input = input;
        this.outPrefix = outPrefix;
        this.force = force;
    }

    @Override
    public Integer call() throws Exception {

        if (!FileUtils.isFile(panclustersRep)) {
            Logger.error(CLASS_NAME, "Pan clusters rep file does not exist.");
            return 1;
        }

        Path panMatrixFile = Path.of(outPrefix + ".pan_matrix.tsv.gz");

        if (!FileUtils.isFile(panMatrixFile) || force) {
            Logger.info(CLASS_NAME, "Creating new pan matrix file.");
            panHap = new PanHap(panclustersRep, panclustersFasta, kmerSize, skipSize, threads);
        } else {
            Logger.info(CLASS_NAME, "Pan matrix file already exists. Samples will be added to the existing file. Use --force to overwrite.");
            panHap = new PanHap(panMatrixFile);
        }

        for (int i = 0; i < input.length; i++) {
            panHap.add(input[i]);
        }

        // Record metadata after samples are added so the counts are accurate.
        panHap.addCommandLine(HelperFunctions.getCommandLine());

        panHap.write(panMatrixFile);

        return 0;
    }

    @Override
    public void run() {
        try {
            HelperFunctions.printCommandLine(new CommandLine(this), CLASS_NAME);
            call();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }
}