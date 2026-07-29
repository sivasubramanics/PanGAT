package edu.tamu.pg.panhap.Utils;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;

public class KMCExecutor {

    private static final String CLASS_NAME = KMCExecutor.class.getSimpleName();

    public static void run(KMCConfig config) throws IOException {
        if (config.getMaxMemory() < 2.0) {
            String errorMsg = "Max memory must be at least 2.0GB for KMC to run. Current: " + config.getMaxMemory();
            Logger.error(CLASS_NAME, errorMsg);
            throw new IllegalArgumentException(errorMsg);
        }

        Path parent = config.getOutputPrefix().getParent();
        Path tmpDir = (parent == null) ? Path.of("tmp") : parent.resolve("tmp");
        FileUtils.mkdir(tmpDir);

        String logFile = config.getOutputPrefix() + ".kmc.log";
        Path listFile = null;

        try {
            String inputArgument;
            Path[] inputs = config.getInputPaths();

            if (inputs.length == 1) {
                inputArgument = inputs[0].toString();
            } else {
                // Multi-input handling: create intermediate input list
                listFile = tmpDir.resolve("input_files.list");
                List<String> lines = Arrays.stream(inputs).map(Path::toString).toList();
                Files.write(listFile, lines);
                inputArgument = "@" + listFile;
            }

            String command = buildCommand(config, inputArgument, tmpDir);

            Logger.info(CLASS_NAME, "Executing KMC for " + inputs.length + " input file(s)...");
            Logger.info(CLASS_NAME, command);

            HelperFunctions.tryExec(command, logFile, logFile);

        } finally {
            if (listFile != null && Files.exists(listFile)) {
                Files.deleteIfExists(listFile);
            }
            FileUtils.rmdir(tmpDir);
        }
    }

    private static String buildCommand(KMCConfig config, String inputArg, Path tmpDir) {
        return String.join(" ",
                "kmc",
                "-k" + config.getKmerSize(),
                "-p" + config.getSignatureLength(),
                "-ci" + config.getMinCount(),
                "-cx" + config.getMaxCount(),
                "-t" + config.getThreads(),
                "-m" + (int) config.getMaxMemory(),
                config.getInputFormat().getFlag(),
                inputArg,
                config.getOutputPrefix().toString(),
                tmpDir.toString()
        );
    }
}