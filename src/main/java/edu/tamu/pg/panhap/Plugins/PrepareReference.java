package edu.tamu.pg.panhap.Plugins;

import edu.tamu.pg.panhap.Data.*;
import edu.tamu.pg.panhap.Utils.*;

import picocli.CommandLine;
import picocli.CommandLine.*;
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.ExecutionException;

@Command(name = "PrepareReferences", description = "Prepare the reference genome for PanHap", sortOptions = false)
public class PrepareReference implements Callable<Integer>, Runnable {

    @Option(names = {"-g", "--genomes"}, description = "Path to the reference genome map file (columns: id, path) in CSV format", required = true)
    private String genomes;

    @Option(names = {"-b", "--basedir"}, description = "Base directory for the output", required = true)
    private String baseDir;

    @Option(names = {"-t", "--threads"}, description = "Number of threads to use for processing", required = false, defaultValue = "1")
    private int threads;

    @Option(names = {"-k", "--kmer"}, description = "K-mer size for KMC", required = false, defaultValue = "31")
    private int kmerSize;

    @Option(names = {"-m", "--min_count"}, description = "Minimum count for Single copy region extraction", required = false, defaultValue = "1")
    private int minCount;

    private static final String CLASS_NAME = PrepareReference.class.getSimpleName();

    private Path genomesDir;
    private Path kmcDir;
    private Path singleCopyDir;

    double HEAPSIZE = Runtime.getRuntime().maxMemory() / (1024.0 * 1024 * 1024);
    int COUNTER_DUMP_SIZE = 100_000;


    public Integer call() throws Exception {
        // if heap size is less than 2GB, KMC won't run
        if (HEAPSIZE < 2.2) {
            Logger.warning(CLASS_NAME, "Max heap size is less than 2.2GB. KMC may fail due to insufficient memory.");
        }

        // check if necessary tools are accessible
        checkDependencies();

        // remove trailing "/" from the base_dir
        while (baseDir.endsWith(File.separator)) {
            baseDir = baseDir.substring(0, baseDir.length() - File.separator.length());
        }
        genomesDir = Paths.get(baseDir, PipelineSteps.GENOMES.dirName());
        kmcDir = Paths.get(baseDir, PipelineSteps.KMC.dirName());
        singleCopyDir = Paths.get(baseDir, PipelineSteps.SINGLE_COPY.dirName());
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
        Map<String, String> genomesMap = readGenomesMap(genomes);
        Logger.info(CLASS_NAME, "Read " + genomesMap.size() + " entries from genomes map file.");

        // fix and copy the genomes to output dir
        fixGenomes(genomesMap);

        // run KMC for each genome to do kmer counting
        runKMC(genomesMap);

        extractSingleCopyRegions(genomesMap);

        Logger.info(CLASS_NAME, "Reference genome preparation completed.");
    }

    private void fixGenomes(Map<String, String> genomesMap) throws IOException {
        // create directory structure under baseDir
        // baseDir + "/00.genomes"
        Logger.info(CLASS_NAME, "Preparing reference genome from map file: " + genomes);
        FileUtils.mkdir(genomesDir);
        // read the genome fasta and rewrite into 00.genomes directory. This is mainly to fix the genome fasta in the format we can work with.
        // sometimes the genomes fasta header might contain tabs instead of spaces, also it can be single [long] line fasta instead of multi line.
        ExecutorService executor = Executors.newFixedThreadPool(threads);
        List<Future<?>> futures = new ArrayList<>();

        try {
            for (Map.Entry<String, String> entry : genomesMap.entrySet()) {
                String id = entry.getKey();
                String genomePath = entry.getValue();
                Path destPath = genomesDir.resolve(id + ".fasta");

                futures.add(executor.submit(() -> {
                    processGenome(destPath.toString(), genomePath);
                    return null;
                }));
            }

            // Wait for all tasks and propagate failures.
            for (Future<?> future : futures) {
                try {
                    future.get();
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                    throw new IOException("Genome processing interrupted.", e);
                } catch (ExecutionException e) {
                    Throwable cause = e.getCause();

                    if (cause instanceof IOException io) {
                        throw io;
                    }

                    throw new IOException("Unexpected error while processing genomes.", cause);
                }
            }
        } finally {
            executor.shutdownNow();
        }
    }

    private void extractSingleCopyRegions(Map<String, String> genomesMap) throws IOException {

        FileUtils.mkdir(singleCopyDir);

        for (Map.Entry<String, String> entry : genomesMap.entrySet()) {
            String id = entry.getKey();

            KMC kmc = new KMC(kmcDir.resolve(id).toString());
            FastaReader genome = new FastaReader(genomesDir.resolve(id + ".fasta").toString());
            int k = kmc.getKmerLength();

            // Larger buffer: 256KB instead of the default 8KB
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(singleCopyDir.resolve(id + ".fasta").toString()), 256 * 1024)) {
                int counter = 0;
                List<Fasta> singleCopyRegions = new ArrayList<>();
                int numSingleCopy = 0;
                for (Fasta chrom : genome) {

                    StringBuilder singleCopySeq = new StringBuilder();
                    boolean isBeginning = true;
                    int regionStart = -1;     // track where the region began
                    int currentEnd = -1;      // track the latest position seen

                    KmerIterator it = chrom.kmers(k, kmc.getPrefixLength(), kmc.isBothStrands());

                    while (it.hasNext()) {
                        KmerPosition kp = it.next();
                        int start = kp.getStartPos();
                        Kmer kmer = kp.getKmer();
                        int kmerCount = kmc.getCount(kmer);

//                        if (kmerCount == 0){
//                            Logger.error(CLASS_NAME, "Kmer not found in KMC database for genome: " + id + ", chromosome: " + chrom.getName() + ", kmer: " + Arrays.toString(kmer.getOriginalKmer()));
//                        }

                        if (kmerCount > 0 && kmerCount <= minCount) {
                            if (isBeginning) {
                                regionStart = start;
                                singleCopySeq.append(kmer.getOriginalKmer());
                                isBeginning = false;
                            } else {
                                // Append only the new last char — avoids a full String allocation
                                singleCopySeq.append(kmer.getOriginalKmer()[k - 1]);
                            }
                            currentEnd = start + k;

                        } else {
                            if (!singleCopySeq.isEmpty()) {
                                singleCopyRegions.add(new Fasta(numSingleCopy, chrom.getName() + ":" + regionStart + "-" + currentEnd, singleCopySeq.toString()));
                                if (counter == COUNTER_DUMP_SIZE) {
                                    for (Fasta scr: singleCopyRegions) {
                                        writer.write(scr.toFastaString(true));
                                    }
                                    singleCopyRegions.clear();
                                    Logger.info(CLASS_NAME, "Dumped " + numSingleCopy + " single copy regions for genome: " + id);
                                    counter = 0;
                                }
                                counter++;
                                numSingleCopy++;
                                singleCopySeq.setLength(0);
                            }
                            isBeginning = true;
                            regionStart = -1;
                            it.skipTo(start + k);
                        }
                    }

                    // ── Flush any region still open at end of chromosome ──
                    if (!singleCopySeq.isEmpty()) {
                        singleCopyRegions.add(new Fasta(numSingleCopy, chrom.getName() + ":" + regionStart + "-" + currentEnd, singleCopySeq.toString()));
                        for (Fasta scr: singleCopyRegions) {
                            writer.write(scr.toFastaString(true));
                        }
                        counter++;
                        numSingleCopy++;
                        singleCopySeq.setLength(0);
                    }
                    Logger.info(CLASS_NAME, "Dumped " + numSingleCopy + " single copy regions for genome: " + id);
                }

                kmc.close();    // once per genome, still inside try so writer flushes first

            } catch (IOException e) {
                Logger.error(CLASS_NAME, "Error writing single copy regions for genome: " + id);
                throw e;
            }
        }
    }

    private void runKMC(Map<String, String> genomesMap) throws IOException {

        // create tmp directory for kmc under base_dir / 01.kmc /tmp
        FileUtils.mkdir(kmcDir.resolve("tmp"));


        for  (String id : genomesMap.keySet()) {
            String logKMCpath = kmcDir.resolve(id + ".kmc.log").toString();
            StringBuilder kmcCmd =  new StringBuilder("kmc ");
            kmcCmd.append("-k").append(kmerSize).append(" ");
            kmcCmd.append("-p").append(Defaults.SIGNATURE_LENGTH).append(" ");
            kmcCmd.append("-ci").append("1").append(" ");
            kmcCmd.append("-cx").append(minCount).append(" ");
            kmcCmd.append("-t").append(threads).append(" ");
            kmcCmd.append("-m").append(HEAPSIZE).append(" ");
            kmcCmd.append("-fm").append(" ");
            kmcCmd.append(genomesDir).append(File.separator).append(id).append(".fasta").append(" ");
            kmcCmd.append(kmcDir).append(File.separator).append(id).append(" ");
            kmcCmd.append(kmcDir).append(File.separator).append("tmp");

            Logger.info(CLASS_NAME, "Running KMC for genome: " + id);
            Logger.info(CLASS_NAME, "KMC command: " + kmcCmd.toString());
            HelperFunctions.tryExec(kmcCmd.toString(), logKMCpath, logKMCpath);
        }
        FileUtils.rmdir(kmcDir.resolve("tmp"));
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
    private Map<String, String> readGenomesMap(String genomesMapFile) throws IOException {
        Map<String, String> genomesMap = new HashMap<>();
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

                // if the number of columns are not 2, throw an error
                if (parts.length != 2) {
                    Logger.error(CLASS_NAME,
                            "Invalid line " + lineNumber + " in genomes map file: " + line);
                    lineNumber++;
                    continue;
                }

                String id = parts[0].trim();
                String path = parts[1].trim();

                // if id or path is empty skipt the entry
                if (id.isEmpty() || path.isEmpty()) {
                    Logger.warning(CLASS_NAME,
                            "Empty ID or path on line " + lineNumber + ". Skipping entry.");
                    lineNumber++;
                    continue;
                }

                // if the ids are repeated, skip
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

    private void checkDependencies() throws IOException {
        List<String> requiredTools = List.of(
                "kmc",
                "mmseqs"
        );

        List<String> missingTools = new ArrayList<>();

        for (String tool : requiredTools) {
            if (!HelperFunctions.isInstalled(tool)) {
                missingTools.add(tool);
            }
        }

        if (!missingTools.isEmpty()) {
            String message = "Missing required tools: " + String.join(", ", missingTools)
                    + ". Please install them and ensure they are available in your system PATH.";

            Logger.error(CLASS_NAME, message);
            throw new IOException(message);
        }
    }
}


