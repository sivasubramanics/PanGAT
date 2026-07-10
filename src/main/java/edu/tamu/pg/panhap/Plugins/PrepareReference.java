package edu.tamu.pg.panhap.Plugins;

import edu.tamu.pg.panhap.Data.*;
import edu.tamu.pg.panhap.Utils.*;
import picocli.CommandLine;
import picocli.CommandLine.*;

import java.io.*;
import java.util.HashMap;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

@Command(name = "PrepareReferences", description = "Prepare the reference genome for PanHap", sortOptions = false)
public class PrepareReference implements Callable<Integer>, Runnable {

    @Option(names = {"-g", "--genomes"}, description = "Path to the reference genome map file (columns: id, path) in CSV format", required = true)
    private String genomes;

    @Option(names = {"-b", "--basedir"}, description = "Base directory for the output", required = true)
    private String basedir;

    @Option(names = {"-t", "--threads"}, description = "Number of threads to use for processing", required = false, defaultValue = "1")
    private int threads;

    private static final String CLASS_NAME = PrepareReference.class.getSimpleName();
    public Integer call() throws Exception {
        prepareReference();
        return 0;
    }

    public void run() {
        try {
            HelperFunctions.printCommandLine(new CommandLine(this), CLASS_NAME);
            call();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    /***
     * Prepares the reference genome by
     *      reading the genome map file,
     *      creating the necessary directory structure,
     *      and processing each genome in parallel using the specified number of threads.
     * @throws IOException
     */
    private void prepareReference() throws IOException {

        // read genome map file
        HashMap<String, String> genomesMap = readGenomesMap(genomes);
        Logger.info(CLASS_NAME, "Read " + genomesMap.size() + " entries from genomes map file.");
        // create directory structure under baseDir
        // baseDir + "/00.genomes"

        String GENOMES_DIR = basedir + File.separator + PipelineSteps.PREP_GENOME.dirName();
        Logger.info(CLASS_NAME, "Preparing reference genome from map file: " + genomes);
        FileUtils.mkdir(GENOMES_DIR);
        // read the genome fasta and rewrite into 00.genomes directory. This is mainly to fix the genome fasta in the format we can work with.
        // sometimes the genomes fasta header might contain tabs instead of spaces, also it can be single [long] line fasta instead of multi line.
        ExecutorService executor = Executors.newFixedThreadPool(threads);
        for (String id : genomesMap.keySet()) {
            String genomePath = genomesMap.get(id);
            String destPath = GENOMES_DIR + File.separator + id + ".fasta";
            executor.submit(() -> {
                try {
                    processGenome(destPath, genomePath);
                } catch (IOException e) {
                    Logger.error(CLASS_NAME, "Error processing genome: " + genomePath);
                    throw new RuntimeException(e);
                }
            });
        }
        executor.shutdown();
        while (!executor.isTerminated()) {
            // Wait for all tasks to complete
        }
        Logger.info(CLASS_NAME, "Reference genome preparation completed.");
    }

    /***
     * Processes a genome by reading it from the given path and writing it to the destination path in a standardized format.
     * It also creates a FASTA index for the processed genome.
     * @param destPath
     * @param genomePath
     * @throws IOException
     */
    private static void processGenome(String destPath, String genomePath) throws IOException {
        Logger.info(CLASS_NAME, "Processing genome: " + genomePath);
        try (FastaWriter fw = new FastaWriter(destPath); FastaReader freader = new FastaReader(genomePath)) {
            for (Fasta fa: freader) {
                fw.write(fa);
            }
        }
        catch (IOException e) {
            Logger.error(CLASS_NAME, "Error copying genome: " + genomePath + " to " + destPath);
            throw e;
        }
        // create fasta index for each fasta file
        new FastaIndex(destPath);
    }

    /***
     * Reads the genomes map file and returns a HashMap of genome IDs to their corresponding paths.
     * @param genomesMapFile
     * @return
     * @throws IOException
     */
    private HashMap<String, String> readGenomesMap(String genomesMapFile) throws IOException {
        HashMap<String, String> genomesMap = new HashMap<>();
        Logger.info(CLASS_NAME, "Reading genomes map from: " + genomesMapFile);
        try (BufferedReader reader = new BufferedReader(new FileReader(genomesMapFile))) {

            // Validate header
            String header = reader.readLine();
            if (header == null) {
                throw new IOException("Genomes map file is empty.");
            }

            String[] headerParts = header.split(",");
            if (headerParts.length != 2
                    || !headerParts[0].trim().equalsIgnoreCase("id")
                    || !headerParts[1].trim().equalsIgnoreCase("path")) {

                Logger.error(CLASS_NAME, "Invalid header in genomes map file: " + header);
                throw new IOException("Expected header: id,path");
            }

            String line;
            int lineNumber = 2; // First data line

            while ((line = reader.readLine()) != null) {

                String[] parts = line.split(",", 2);

                if (parts.length != 2) {
                    Logger.error(CLASS_NAME,
                            "Invalid line " + lineNumber + " in genomes map file: " + line);
                    lineNumber++;
                    continue;
                }

                String id = parts[0].trim();
                String path = parts[1].trim();

                if (genomesMap.containsKey(id)) {
                    Logger.warning(CLASS_NAME,
                            "Duplicate entry for ID '" + id + "' on line " + lineNumber + ". Overwriting previous value.");
                }
                genomesMap.put(id, path);
                lineNumber++;
            }
        }
        return genomesMap;
    }
}


