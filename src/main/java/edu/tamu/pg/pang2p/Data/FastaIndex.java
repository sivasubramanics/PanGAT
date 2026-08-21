package edu.tamu.pg.pang2p.Data;

import edu.tamu.pg.pang2p.Utils.HelperFunctions;
import edu.tamu.pg.pang2p.Utils.Logger;

import java.io.*;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

/**
 * FASTA index for random-access (samtools faidx style).
 * Memory maps the whole file once. Scales to millions of sequences.
 */
public class FastaIndex implements AutoCloseable {
    private final Map<String, FastaIndexEntry> index;
    private final List<String> sequenceNames;
    private MappedByteBuffer fastaBuffer;
    private FileChannel fastaFileChannel;
    private File fastaFile;
    private final Object bufferLock = new Object();
    private final String CLASS_NAME = this.getClass().getSimpleName();
    private static final Set<Character> VALID_BASES = Set.of(
            'A','C','G','T','Y','R','W','S','M','K','H','B','V','D','N',
            'a','c','g','t','y','r','w','s','m','k','h','b','v','d','n'
    );

    public FastaIndex(Path fastaFilePath) throws IOException {
        this(fastaFilePath.toAbsolutePath().toString());
    }
    public FastaIndex(File fastaFile) throws IOException {
        this(fastaFile.getAbsolutePath());
    }

    public FastaIndex(String fastaFilePath) throws IOException {
        this.fastaFile = new File(fastaFilePath);
        String faiFilePath = fastaFilePath + ".faidx";
        File indexFile = new File(faiFilePath);

        if (!indexFile.exists() || HelperFunctions.isOlder(indexFile, fastaFile)) {
            Logger.info(CLASS_NAME, "Generating/Updating index file: " + faiFilePath);
            generateIndexFile(fastaFile, indexFile);
        } else {
            Logger.info(CLASS_NAME, "Using existing index file: " + faiFilePath);
        }

        // Load sorted index
        this.index = sortByValue(Collections.unmodifiableMap(loadIndex(indexFile)));
        this.sequenceNames = Collections.unmodifiableList(new ArrayList<>(index.keySet()));
    }

    /**
     * Explicitly loads or memory-maps the FASTA file.
     * Can be called manually, or will be called automatically when sequence extraction is requested.
     */
    public void loadToMemory() throws IOException {
        synchronized (bufferLock) {
            if (fastaFileChannel != null) {
                return; // Already loaded/mapped
            }

            try {
                RandomAccessFile raf = new RandomAccessFile(this.fastaFile, "r");
                this.fastaFileChannel = raf.getChannel();
                long fileSize = fastaFileChannel.size();

                if (fileSize > Integer.MAX_VALUE) {
                    Logger.info(CLASS_NAME, "Large FASTA file detected (>2GB), using RandomAccessFile (not memory-mapped).");
                    this.fastaBuffer = null;
                } else {
                    Logger.info(CLASS_NAME, "Memory mapping entire fasta file (size: " + fileSize + ")");
                    this.fastaBuffer = fastaFileChannel.map(FileChannel.MapMode.READ_ONLY, 0, fileSize);
                }
            } catch (IOException e) {
                Logger.error(CLASS_NAME, "Error memory-mapping fasta file: " + e.getMessage());
                throw e;
            }
        }
    }

    private Map<String, FastaIndexEntry> loadIndex(File indexFile) throws IOException {
        Map<String, FastaIndexEntry> tempIndex = new ConcurrentHashMap<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(indexFile))) {
            String line;
            int seqId = 0;
            while ((line = reader.readLine()) != null) {
                String[] fields = line.split("\t", 6);
                String desc = (fields.length >= 6 && !fields[5].isEmpty()) ? fields[5] : null;
                FastaIndexEntry entry = new FastaIndexEntry(seqId,
                        fields[0],
                        Integer.parseInt(fields[1]),
                        Long.parseLong(fields[2]),
                        Integer.parseInt(fields[3]),
                        Integer.parseInt(fields[4]),
                        desc
                );
                seqId++;
                if (tempIndex.putIfAbsent(entry.getSequenceName(), entry) != null) {
                    Logger.error(CLASS_NAME, "Duplicate sequence name in index: " + entry.getSequenceName());
                }
            }
        }
        return tempIndex;
    }

    public FastaIndexEntry getEntry(String sequenceName) {
        return index.get(sequenceName);
    }

    public String getSequence(String sequenceName) {
        FastaIndexEntry entry = getEntry(sequenceName);
        if (entry == null) {
            Logger.error(CLASS_NAME, "Sequence not found in index: " + sequenceName);
            return null;
        }
        return getSequence(sequenceName, 0, entry.getLength());
    }

    /**
     * Efficient random-access for subsequence.
     * Auto-calls loadToMemory() if the sequence buffer is not open yet.
     */
    public String getSequence(String name, int start, int length) {
        FastaIndexEntry entry = getEntry(name);
        if (entry == null) {
            Logger.error(CLASS_NAME, "Sequence not found in index: " + name);
            return null;
        }
        int end = start + length;
        if (start < 0 || end > entry.getLength() || start >= end) {
            Logger.error(CLASS_NAME, "Invalid range: " + start + "-" + end + " for sequence: " + name);
            return null;
        }

        StringBuilder sequence = new StringBuilder(end - start);

        synchronized (bufferLock) {
            try {
                // Ensure file channel / buffer is ready
                if (fastaFileChannel == null) {
                    loadToMemory();
                }

                // get index info
                long seqOffset = entry.getOffset();
                int lineBases = entry.getLineBases();
                int lineWidth = entry.getLineWidth();

                // calculate where in the file this region starts
                int startLine = start / lineBases;
                int startLineBaseIndex = start % lineBases;
                long startOffsetInFile = seqOffset + ((long) startLine * lineWidth) + startLineBaseIndex;

                int basesToExtract = end - start;
                long curPosInFile = startOffsetInFile;
                int curBaseInLine = startLineBaseIndex;

                while (basesToExtract > 0) {
                    int remainingBasesInLine = lineBases - curBaseInLine;
                    int basesToRead = Math.min(basesToExtract, remainingBasesInLine);

                    if (fastaBuffer != null) {
                        fastaBuffer.position((int) curPosInFile);
                        for (int i = 0; i < basesToRead; i++) {
                            sequence.append((char) fastaBuffer.get());
                        }
                    } else { // fallback for >2GB files
                        byte[] chunk = new byte[basesToRead];
                        fastaFileChannel.position(curPosInFile);
                        int read = fastaFileChannel.read(java.nio.ByteBuffer.wrap(chunk));
                        for (int i = 0; i < read; i++) {
                            sequence.append((char) chunk[i]);
                        }
                    }
                    curPosInFile += basesToRead + (lineWidth - lineBases);
                    basesToExtract -= basesToRead;
                    curBaseInLine = 0;
                }
            } catch (Exception e) {
                Logger.error(CLASS_NAME, "Error reading sequence: " + e.getMessage());
                return null;
            }
        }

        return sequence.toString();
    }

    public boolean containsSequence(String sequenceName) {
        return index.containsKey(sequenceName);
    }

    @Override
    public void close() throws IOException {
        synchronized (bufferLock) {
            // 1. Explicitly unmap native memory buffer
            if (fastaBuffer != null) {
                unmap(fastaBuffer);
                fastaBuffer = null;
            }

            // 2. Close the underlying file channel
            if (fastaFileChannel != null && fastaFileChannel.isOpen()) {
                fastaFileChannel.close();
                fastaFileChannel = null;
            }

            // 3. Nullify large internal index data structures
            if (index instanceof Map) {
                try {
                    // Handle immutable wrapping cleanly if needed
                    index.clear();
                } catch (UnsupportedOperationException ignored) {}
            }
        }
    }

    /**
     * Helper method to immediately release MappedByteBuffer off-heap memory.
     */
    private static void unmap(MappedByteBuffer buffer) {
        if (buffer == null) return;
        try {
            // Direct byte buffer cleaner invocation via reflection
            java.lang.reflect.Method cleanerMethod = buffer.getClass().getMethod("cleaner");
            cleanerMethod.setAccessible(true);
            Object cleaner = cleanerMethod.invoke(buffer);
            if (cleaner != null) {
                java.lang.reflect.Method cleanMethod = cleaner.getClass().getMethod("clean");
                cleanMethod.setAccessible(true);
                cleanMethod.invoke(cleaner);
            }
        } catch (Exception e) {
            // Direct buffer cleanup fallback for modern Java runtimes (JDK 9+)
            try {
                Class<?> unsafeClass = Class.forName("sun.misc.Unsafe");
                java.lang.reflect.Field unsafeField = unsafeClass.getDeclaredField("theUnsafe");
                unsafeField.setAccessible(true);
                Object unsafe = unsafeField.get(null);
                java.lang.reflect.Method invokeCleaner = unsafeClass.getMethod("invokeCleaner", java.nio.ByteBuffer.class);
                invokeCleaner.invoke(unsafe, buffer);
            } catch (Exception ignored) {
                // If reflection is restricted, normal GC will clean up the buffer later
                Logger.error("FastaIndex", "Unable to unmap MappedByteBuffer. It will be cleaned up by GC later.");
            }
        }
    }

    public Fasta getFasta(String tr) {
        synchronized (bufferLock) {
            int seqId = sequenceNames.indexOf(tr);
            FastaIndexEntry entry = getEntry(tr);
            if (entry == null) return null;
            String desc = entry.getDescription();
            return new Fasta(seqId, tr, getSequence(tr), desc);
        }
    }

    public int length() {
        return sequenceNames.size();
    }

    public List<String> getSequenceNames() {
        return sequenceNames;
    }

    public int getSequenceLength(String name) {
        FastaIndexEntry entry = getEntry(name);
        if (entry == null) {
            Logger.error(CLASS_NAME, "Sequence not found in index: " + name);
            return -1;
        }
        return entry.getLength();
    }

    /**
     * Index generation -- unchanged (does not load entire FASTA into memory).
     */
    private void generateIndexFile(File fastaFile, File indexFile) throws IOException {
        if (HelperFunctions.isCompressed(fastaFile)) {
            Logger.error(CLASS_NAME, "Fasta file is compressed. Please decompress before indexing: " + fastaFile);
            throw new IllegalArgumentException("Fasta file is compressed. Please decompress before indexing: " + fastaFile);
        }
        Set<String> sequenceNames = new HashSet<>();
        try (BufferedReader reader = new BufferedReader(
                new InputStreamReader(new FileInputStream(fastaFile), java.nio.charset.StandardCharsets.ISO_8859_1));
             BufferedWriter writer = new BufferedWriter(new FileWriter(indexFile))) {
            String line;
            int lineNumber = 0;
            long offset = 0;
            String currentName = null;
            String currentDesc = null;
            int lineBases = 0;
            int lineWidth = 0;
            int seqLength = 0;
            long sequenceStartOffset = 0;
            boolean inSequence = false;
            while ((line = reader.readLine()) != null) {
                lineNumber++;
                int originalLength = line.length();
                String trimmedLine = line.trim();
                if (offset == 0 && !trimmedLine.startsWith(">")) {
                    Logger.error(CLASS_NAME, "Invalid fasta file: " + fastaFile);
                }
                if (trimmedLine.startsWith(">")) {
                    // WRITE previous record if any
                    if (inSequence) {
                        if (currentDesc != null && !currentDesc.isEmpty()) {
                            writer.write(currentName + "\t" + seqLength + "\t" + sequenceStartOffset
                                    + "\t" + lineBases + "\t" + lineWidth + "\t" + currentDesc + "\n");
                        } else {
                            writer.write(currentName + "\t" + seqLength + "\t" + sequenceStartOffset
                                    + "\t" + lineBases + "\t" + lineWidth + "\n");
                        }
                    }
                    // Extract name & desc robustly
                    String[] tokens = trimmedLine.substring(1).split("\\s+", 2);
                    currentName = tokens[0];
                    currentDesc = (tokens.length > 1) ? tokens[1] : null;
                    offset += originalLength + System.lineSeparator().length();
                    sequenceStartOffset = offset;
                    inSequence = true;
                    seqLength = 0;
                    if (sequenceNames.contains(currentName)) {
                        Logger.error(CLASS_NAME, "Duplicate sequence name in fasta file: " + currentName + " at line " + lineNumber);
                    } else {
                        sequenceNames.add(currentName);
                    }
                } else {
                    if (seqLength == 0 && !trimmedLine.isEmpty()) {
                        lineBases = originalLength;
                        lineWidth = originalLength + System.lineSeparator().length();
                    }
                    seqLength += originalLength;
                    offset += originalLength + System.lineSeparator().length();
                }
            }
            // Final record
            if (inSequence && currentName != null) {
                if (currentDesc != null && !currentDesc.isEmpty()) {
                    writer.write(currentName + "\t" + seqLength + "\t" + sequenceStartOffset
                            + "\t" + lineBases + "\t" + lineWidth + "\t" + currentDesc + "\n");
                } else {
                    writer.write(currentName + "\t" + seqLength + "\t" + sequenceStartOffset
                            + "\t" + lineBases + "\t" + lineWidth + "\n");
                }
            }
        }
    }


    /**
     * Sort index
     */
    public static Map<String, FastaIndexEntry> sortByValue(Map<String, FastaIndexEntry> map) {
        List<Map.Entry<String, FastaIndexEntry>> entryList = new ArrayList<>(map.entrySet());
        entryList.sort(Map.Entry.comparingByValue());
        Map<String, FastaIndexEntry> sortedMap = new LinkedHashMap<>();
        for (Map.Entry<String, FastaIndexEntry> entry : entryList) {
            sortedMap.put(entry.getKey(), entry.getValue());
        }
        return sortedMap;
    }
}
//EOF