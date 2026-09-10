package edu.tamu.pg.pang2p.Tasks;

import edu.tamu.pg.pang2p.Data.PanHap;
import edu.tamu.pg.pang2p.Utils.HelperFunctions;
import edu.tamu.pg.pang2p.Utils.Logger;
import picocli.CommandLine;
import picocli.CommandLine.*;
import picocli.CommandLine.Model.CommandSpec; // Import CommandSpec

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.Callable;

@Command(name = "cohort", description = "Compile raw mapped k-mer counts across sample coverage files into cohort.pan_matrix.tsv.gz", sortOptions = false)
public class Cohort implements Runnable, Callable<Integer> {

    @Spec
    private CommandSpec spec; // 1. Inject Picocli CommandSpec

    @Option(names = {"-p", "--panclusters"}, description = "Path to pan_clusters_rep.tsv file.", required = true)
    private Path panclustersFile;

    @Option(names = {"-i", "--input"}, description = "Path to input sample coverage files.", required = true, arity = "1..*", split = "\\s+")
    private Path[] input;

    @Option(names = {"-s", "--samples"}, description = "List of sample names.", required = true, arity = "1..*", split = "\\s+")
    private String[] samples;

    @Option(names = {"-o", "--outprefix"}, description = "Prefix for output files.", required = true)
    private String outPrefix;

    @Option(names = {"-f", "--force"}, description = "Force overwrite of existing mapped_kmers file instead of updating.", defaultValue = "false")
    private boolean force;

    @Option(names = {"--cv-cutoff"}, description = "Coefficient of variation cutoff.", defaultValue = "0.0")
    private double cvCutoff;

    private final String CLASS_NAME = Cohort.class.getSimpleName();
    private PanHap panHap;

    public Cohort() {
    }

    public Cohort(Path panclustersFile, Path[] input, String[] samples, String outPrefix, boolean force) {
        this.panclustersFile = panclustersFile;
        this.input = input;
        this.samples = samples;
        this.outPrefix = outPrefix;
        this.force = force;
    }

    @Override
    public Integer call() throws Exception {
        if (!Files.exists(panclustersFile)) {
            Logger.error(CLASS_NAME, "Pan clusters file does not exist.");
            return 1;
        }

        if (input.length != samples.length) {
            Logger.error(CLASS_NAME, "Number of input files must match number of sample names.");
            return 1;
        }

        Path panMatrixFile = Path.of(outPrefix + ".pan_matrix.tsv.gz");
        if (Files.exists(panMatrixFile) && !force) {
            Logger.info(CLASS_NAME, "Pan matrix file already exists. Samples will be added to the existing file. Use --force to overwrite.");
            panHap = new PanHap(panclustersFile, panMatrixFile);
        } else {
            Logger.info(CLASS_NAME, "Creating new pan matrix file.");
            panHap = new PanHap(panclustersFile);
        }

        // 2. Set the command line string on PanHap
        panHap.addCommandLine(HelperFunctions.getCommandLine());

        // Optional: Add parameters to metadata
        panHap.addParameter("cvCutoff", String.valueOf(cvCutoff));

        panHap.addParameter("numHaps", String.valueOf(panHap.getHapNum()));

        for (int i = 0; i < input.length; i++) {
            panHap.addSample(samples[i], input[i]);
        }

        panHap.writePanMatrix(panMatrixFile, cvCutoff);

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