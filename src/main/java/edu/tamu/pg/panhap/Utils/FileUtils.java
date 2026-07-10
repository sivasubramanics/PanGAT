package edu.tamu.pg.panhap.Utils;

import java.io.File;
import java.io.IOException;
import java.nio.file.*;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.ArrayList;
import java.util.List;

public class FileUtils {

    private static final String CLASS_NAME = FileUtils.class.getSimpleName();

    /***
     * Create directory on the file system recursively (something like mkdir -p dirName)
     * @param dirName
     * @return
     * @throws IOException
     */
    public static boolean mkdir(String dirName) throws IOException {
        File dir = new File(dirName);
        if (!dir.exists()) {
            return dir.mkdirs();
        }
        return true;
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

    public static boolean deleteRecursively(File file) {
        try {
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
}