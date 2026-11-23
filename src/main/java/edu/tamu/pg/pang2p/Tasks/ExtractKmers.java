package edu.tamu.pg.pang2p.Tasks;

import edu.tamu.pg.pang2p.Data.*;
import picocli.CommandLine.*;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

@Command(name = "extract-kmers", description = "Extract kmers from a multi-fasta file", sortOptions = false)
public class ExtractKmers implements Runnable, Callable<Integer> {
    @Option(names = {"-i", "--input"}, description = "Path to the multi-fasta file", required = true)
    private Path fastaFile;

    @Option(names = {"-k", "--kmer"}, description = "K-mer size for extraction [31]", required = false, defaultValue = "31")
    private int kmerSize = 31;

    @Option(names = {"-s", "--skip"}, description = "Process every `skip`-th k-mer for extraction [1]", required = false, defaultValue = "1")
    private int skip = 1;

    @Option(names = {"-o", "--output"}, description = "Output fasta file path for extracted kmers", required = true)
    private Path outputFile;

    @Option(names = {"-t", "--threads"}, description = "Number of threads to use for processing [2]", required = false, defaultValue = "2")
    private int threads = 2;

    @Option(names = {"-c", "--canonical"}, description = "Use canonical kmers (default: false)", required = false, defaultValue = "false")
    private boolean canonical = false;

    public ExtractKmers() {
    }

    public ExtractKmers(Path fastaFile, int kmerSize, Path outputFile) {
        this.fastaFile = fastaFile;
        this.kmerSize = kmerSize;
        this.outputFile = outputFile;
    }

    @Override
    public void run() {
        try{
            call();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    @Override
    public Integer call() throws Exception {
        try {
            extractKmers();
            return 0;
        } catch (Exception e) {
            throw new Exception("Error occurred while extracting kmers", e);
        }
    }

    private void extractKmers() throws Exception {
        int effectiveSkip = Math.max(skip, 1);
        int numThreads = Math.max(1, Runtime.getRuntime().availableProcessors());
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);
        List<Future<List<Fasta>>> futures = new ArrayList<>();

        try (FastaReader reader = new FastaReader(fastaFile);
             FastaWriter writer = new FastaWriter(outputFile)) {

            for (Fasta fasta : reader) {
                futures.add(executor.submit(() -> extractForRecord(fasta, effectiveSkip)));
            }

            // Iterating in submission order guarantees output matches input FASTA order
            for (Future<List<Fasta>> future : futures) {
                List<Fasta> records = future.get();
                for (Fasta record : records) {
                    writer.write(record);
                }
            }
        } catch (ExecutionException e) {
            Throwable cause = e.getCause() != null ? e.getCause() : e;
            throw new Exception("Error occurred while extracting kmers from fasta file", cause);
        } catch (Exception e) {
            throw new Exception("Error occurred while extracting kmers from fasta file", e);
        } finally {
            executor.shutdownNow();
        }
    }

    private List<Fasta> extractForRecord(Fasta fasta, int effectiveSkip) {
        KmerIterator it = fasta.kmers(kmerSize, 5, canonical);
        String baseName = fasta.getName();

        List<Fasta> results = new ArrayList<>();
        int kmerIndex = 0;
        int kmerNum = 0;

        while (it.hasNext()) {
            KmerPosition kp = it.next();
            int start = kp.getStartPos();

            if (kmerIndex % effectiveSkip != 0) {
                kmerIndex++;
                continue;
            }
            kmerIndex++;
            kmerNum++;
            results.add(new Fasta(1, baseName + ":" + kmerNum, kp.getKmer().toString(), baseName + ":" + start + "-" + (start + kmerSize)));
        }
        return results;
    }
}
