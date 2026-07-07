package edu.tamu.pg.panhap.Plugins;


import edu.tamu.pg.panhap.Utils.FileUtils;
import edu.tamu.pg.panhap.Utils.HelperFunctions;
import edu.tamu.pg.panhap.Utils.Logger;
import picocli.CommandLine;
import picocli.CommandLine.*;

import java.util.concurrent.Callable;

@Command(name = "PrepareReferences", description = "Prepare the reference genome for PanHap", sortOptions = false)
public class PrepareReference implements Callable<Integer>, Runnable {

    @Option(names = {"-g", "--genomes"}, description = "Path to the reference genome map file (columns: id, path) in CSV format", required = true)
    private String genomes;

    @Option(names = {"-b", "--basedir"}, description = "Base directory for the output", required = true)
    private String basedir;

    private static final String CLASS_NAME = PrepareReference.class.getSimpleName();
    public Integer call() throws Exception {
        // Implement the logic for preparing the reference here
        prepareReference();
        return 0;
    }

    public void run() {
        try {
            HelperFunctions.printCommandLine(new CommandLine(this), CLASS_NAME);
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    private void prepareReference() {

        Logger.info(CLASS_NAME, "Preparing reference genome from map file: " + genomes);
        Logger.info(CLASS_NAME, "Base directory for output: " + basedir);

        // create directory structure under baseDir
        // baseDir + "/00.genomes"
        FileUtils.mkdir(basedir + "/00.genomes");

        // Add the logic to read the genomes map file and prepare the reference genome
        // This may include reading the CSV file, validating paths, and performing necessary operations

        Logger.info(CLASS_NAME, "Reference genome preparation completed.");

    }
}


