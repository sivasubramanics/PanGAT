package edu.tamu.pg.pang2p.Utils;

import java.io.*;
import java.nio.channels.FileChannel;
import java.nio.file.*;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class FileUtils {

    private static final String CLASS_NAME = FileUtils.class.getSimpleName();
    private static final int BUFFER_SIZE = 64 * 1024;

    /***
     * Create directory on the file system recursively (something like mkdir -p dirName)
     *
     * @param dir
     * @return
     * @throws IOException
     */
    public static boolean mkdir(Path dir) throws IOException {
        Files.createDirectories(dir);
        return true;
    }

    public static boolean mkdir(String dirName) throws IOException {
        return mkdir(Path.of(dirName));
    }

    public static boolean mkdir(File dir) throws IOException {
        return mkdir(dir.toPath());
    }

    public static void rmdir(Path dir) throws IOException {
        deleteRecursively(dir.toFile());
    }

    public static void rmdir(String dirName) throws IOException {
        rmdir(Path.of(dirName));
    }

    public static void rmdir(File dir) throws IOException {
        rmdir(dir.toPath());
    }


    public static boolean cleanTmpDir(File tmpDir) {
        if (tmpDir == null || !tmpDir.exists() || !tmpDir.isDirectory()) {
            return false;
        }

        try {
            Files.walkFileTree(tmpDir.toPath(), new SimpleFileVisitor<>() {
                @Override
                public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
                    Files.delete(file);
                    return FileVisitResult.CONTINUE;
                }

                @Override
                public FileVisitResult postVisitDirectory(Path dir, IOException exc) throws IOException {
                    if (!dir.equals(tmpDir.toPath())) {
                        Files.delete(dir);
                    }
                    return FileVisitResult.CONTINUE;
                }
            });
            return true;
        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Failed to clean tmp directory: " + e.getMessage());
            return false;
        }
    }

    public static String getParentDirName(File file) {
        File parent = file.getParentFile();
        return parent != null ? parent.getName() : null;
    }

    public static String getParentDirFullPath(File file) {
        File parent = file.getParentFile();
        return parent != null ? parent.getAbsolutePath() : null;
    }

    public static List<File> getChildrenWithPrefix(File dir, String prefix) {
        List<File> result = new ArrayList<>();
        if (dir != null && dir.isDirectory()) {
            File[] children = dir.listFiles();
            if (children != null) {
                for (File child : children) {
                    if (child.getName().startsWith(prefix)) {
                        result.add(child);
                    }
                }
            }
        }
        return result;
    }

    public static boolean deleteRecursively(Path path) throws IOException {
        return deleteRecursively(path.toFile());
    }

    public static boolean deleteRecursively(File file) {
        try {
            // if it doesn't exist skip
            if (!file.isDirectory() || !file.exists() || !file.canRead()) {
                return true;
            }
            Files.walkFileTree(file.toPath(), new SimpleFileVisitor<>() {
                @Override
                public FileVisitResult visitFile(Path path, BasicFileAttributes attrs) throws IOException {
                    Files.delete(path);
                    return FileVisitResult.CONTINUE;
                }

                @Override
                public FileVisitResult postVisitDirectory(Path dir, IOException exc) throws IOException {
                    Files.delete(dir);
                    return FileVisitResult.CONTINUE;
                }
            });
            return true;
        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Failed to delete recursively: " + e.getMessage());
            return false;
        }
    }

    public static String replaceExtension(Path file, String oldExtension, String newExtension) {
        String path = file.toString();
        int lastIndex = path.lastIndexOf(oldExtension);
        if (lastIndex != -1) {
            return path.substring(0, lastIndex) + newExtension + path.substring(lastIndex + oldExtension.length());
        } else {
            return path; // return original if extension not found
        }
    }

    public static File replaceExtension(File file, String oldExtension, String newExtension) {
        String path = file.getAbsolutePath();
        int lastIndex = path.lastIndexOf(oldExtension);
        if (lastIndex != -1) {
            String newPath = path.substring(0, lastIndex) + newExtension + path.substring(lastIndex + oldExtension.length());
            return new File(newPath);
        } else {
            return file; // return original if extension not found
        }
    }

    public static void copyFile(String sourcePath, String targetPath) {
        File sourceFile = new File(sourcePath);
        File targetFile = new File(targetPath);
        copyFile(sourceFile, targetFile);
    }

    public static void copyFile(File clusterFile, File targetFile) {
        try {
            if (!targetFile.exists()) {
                Files.copy(clusterFile.toPath(), targetFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
            } else {
                Logger.info(CLASS_NAME, "Target file already exists: " + targetFile.getAbsolutePath() + ", replacing it.");
                Files.deleteIfExists(targetFile.toPath());
                Files.copy(clusterFile.toPath(), targetFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
            }
        } catch (IOException e){
            Logger.error(CLASS_NAME, "Failed to copy file: " + e.getMessage());
        }
    }

    public static boolean isDir(String tmpDir) {
        return tmpDir != null && new File(tmpDir).isDirectory();
    }

    public static boolean isDir(Path path) {
        return isDir(path.toString());
    }

    public static boolean isFile(String tmpDir) {
        return tmpDir != null && new File(tmpDir).isFile();
    }

    public static boolean isFile(Path path) {
        return isFile(path.toString());
    }

    public static void rmFile(Path file) {
        try {
            Files.deleteIfExists(file);
        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Failed to delete file: " + e.getMessage());
        }
    }

    public static BufferedReader reader(Path path) throws IOException {
        InputStream is = Files.newInputStream(path);
        if (path.getFileName().toString().endsWith(".gz")) {
            is = new GZIPInputStream(is);
        }
        return new BufferedReader(new InputStreamReader(is));
    }

    public static BufferedReader reader(File file) throws IOException {
        return reader(file.toPath());
    }

    public static BufferedReader reader(String filePath) throws IOException {
        return reader(Path.of(filePath));
    }

    /**
     * Creates a BufferedWriter for a file, automatically handling gzipped files (.gz).
     */
    public static BufferedWriter writer(Path path) throws IOException {
        OutputStream os = Files.newOutputStream(path);
        if (path.getFileName().toString().endsWith(".gz")) {
            os = new GZIPOutputStream(os);
        }
        return new BufferedWriter(new OutputStreamWriter(os));
    }

    public static BufferedWriter writer(File file) throws IOException {
        return writer(file.toPath());
    }

    public static BufferedWriter writer(String filePath) throws IOException {
        return writer(Path.of(filePath));
    }

    public static String[] fileToArray(Path filePath) throws IOException {
        List<String> lines = Files.readAllLines(filePath);
        return lines.stream().map(String::trim).filter(line -> !line.isEmpty()).toArray(String[]::new);
    }

    public static void catFiles(String[] inFiles, String outFile) throws IOException {
        if (inFiles == null || inFiles.length == 0) {
            return;
        }

        Path outputPath = Path.of(outFile);
        boolean isGz = outFile.endsWith(".gz");

        if (!isGz) {
            // Zero-Copy OS transfer mode for uncompressed files
            try (FileChannel outChannel = FileChannel.open(outputPath,
                    StandardOpenOption.CREATE,
                    StandardOpenOption.WRITE,
                    StandardOpenOption.TRUNCATE_EXISTING)) {

                for (String inFile : inFiles) {
                    Path inPath = Path.of(inFile);
                    if (!Files.exists(inPath)) continue;

                    try (FileChannel inChannel = FileChannel.open(inPath, StandardOpenOption.READ)) {
                        long transferred = 0;
                        long size = inChannel.size();
                        while (transferred < size) {
                            transferred += inChannel.transferTo(transferred, size - transferred, outChannel);
                        }
                    }
                }
            }
        } else {
            // Optimized stream mode for GZIP streams
            try (OutputStream os = new BufferedOutputStream(new GZIPOutputStream(Files.newOutputStream(outputPath)), BUFFER_SIZE)) {
                byte[] buffer = new byte[BUFFER_SIZE];
                for (String inFile : inFiles) {
                    Path inPath = Path.of(inFile);
                    if (!Files.exists(inPath)) continue;

                    try (InputStream is = new BufferedInputStream(Files.newInputStream(inPath), BUFFER_SIZE)) {
                        int bytesRead;
                        while ((bytesRead = is.read(buffer)) != -1) {
                            os.write(buffer, 0, bytesRead);
                        }
                    }
                }
            }
        }
    }

    public static void catFiles(List<Path> inFiles, Path outFile) throws IOException {
        String[] inPaths = inFiles.stream().map(Path::toString).toArray(String[]::new);
        catFiles(inPaths, outFile.toString());
    }

    public static void rmFiles(List<Path> tempFiles) {
        for (Path tempFile : tempFiles) {
            try {
                Files.deleteIfExists(tempFile);
            } catch (IOException e) {
                // Ignore
            }
        }
    }

    public static String getExtension(String string) {
        int lastDot = string.lastIndexOf('.');
        return lastDot == -1 ? "" : string.substring(lastDot);
    }

    public static boolean isFastaFile(String refFastaFile) {
        String lower = refFastaFile.toLowerCase();
        return lower.endsWith(".fasta") || lower.endsWith(".fa") || lower.endsWith(".fna") ||
                lower.endsWith(".fasta.gz") || lower.endsWith(".fa.gz") || lower.endsWith(".fna.gz");
    }

    public static boolean isFastaFile(Path refFastaFile) {
        return isFastaFile(refFastaFile.toString());
    }

    public static Path replaceExtension(String inFile, String s) {
        String newFileName = inFile.replaceAll("\\.[^.]+$", s);
        return Path.of(newFileName);
    }

    public static Path replaceExtension(Path inFile, String s) {
        String pathStr = inFile.toString();
        String lower = pathStr.toLowerCase();

        // Check for multi-part or standard fasta extensions and replace them cleanly
        if (lower.endsWith(".fasta.gz")) {
            pathStr = pathStr.substring(0, pathStr.length() - 9);
        } else if (lower.endsWith(".fasta")) {
            pathStr = pathStr.substring(0, pathStr.length() - 6);
        } else if (lower.endsWith(".fna")) {
            pathStr = pathStr.substring(0, pathStr.length() - 4);
        } else if (lower.endsWith(".fa")) {
            pathStr = pathStr.substring(0, pathStr.length() - 3);
        } else {
            // Fallback to standard last-extension replacement
            pathStr = pathStr.replaceAll("\\.[^.]+$", "");
        }

        return Path.of(pathStr + s);
    }

    public static boolean isFastqFile(String file) {
        String lower = file.toLowerCase();
        return lower.endsWith(".fastq") || lower.endsWith(".fq") || lower.endsWith(".fastq.gz") || lower.endsWith(".fq.gz");
    }

    public static boolean isFastqFile(Path file) {
        return isFastqFile(file.toString());
    }

    /**
     * Infers the input file format ("fasta", "fastq", or "list") based on file extensions.
     * Returns null if the format/extension is unsupported.
     */
    public static String determineInputFormat(String fileName) {
        if (fileName == null) return null;
        String lower = fileName.toLowerCase();

        if (lower.endsWith(".fa") || lower.endsWith(".fasta") || lower.endsWith(".fna") ||
                lower.endsWith(".fa.gz") || lower.endsWith(".fasta.gz") || lower.endsWith(".fna.gz")) {
            return "fasta";
        } else if (lower.endsWith(".fq") || lower.endsWith(".fastq") ||
                lower.endsWith(".fq.gz") || lower.endsWith(".fastq.gz")) {
            return "fastq";
        } else if (lower.endsWith(".list")) {
            return "list";
        }
        return null;
    }

    public static String determineInputFormat(Path filePath) {
        if (filePath == null || filePath.getFileName() == null) return null;
        return determineInputFormat(filePath.getFileName().toString());
    }
}