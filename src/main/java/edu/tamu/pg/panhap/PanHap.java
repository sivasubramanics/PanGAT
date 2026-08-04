package edu.tamu.pg.panhap;

import edu.tamu.pg.panhap.Tasks.KmerCoverage;
import edu.tamu.pg.panhap.Tasks.SingleCopy;
import picocli.CommandLine;
import picocli.CommandLine.*;
import edu.tamu.pg.panhap.Plugins.*;
import edu.tamu.pg.panhap.Utils.*;


import java.util.Locale;


@Command(name = "panhap", description = "PANgenome HAPlotype catalog and genotyping toolkit", mixinStandardHelpOptions = true, versionProvider = Configs.VersionProvider.class,
        subcommands = {
            PrepareReference.class,
                Genotype.class,
                KmerCoverage.class,
                SingleCopy.class
        })
public class PanHap {
    public static void main(String[] args) {
        Locale.setDefault(Locale.US);
        final PanHap panHap = new PanHap();
        final CommandLine cmd = new CommandLine(panHap);

        boolean isHelpOrVersion = false;

        for (String arg : args) {
            if (arg.equals("--help") || arg.equals("-h") || arg.equals("--version") || arg.equals("-V")) {
                isHelpOrVersion = true;
                break;
            }
        }

        long startTime = isHelpOrVersion ? 0 : System.currentTimeMillis();
        int exitCode = cmd.execute(args);
        long endTime = System.currentTimeMillis();

        if (!isHelpOrVersion) {
            String executionTime = String.format("%02d:%02d:%02d",
                    (endTime - startTime) / 3600000,
                    ((endTime - startTime) % 3600000) / 60000,
                    ((endTime - startTime) % 60000) / 1000);

            switch (exitCode) {
                case 0:
                    Logger.info(PanHap.class.getSimpleName(), "Total execution time: " + executionTime);
                    break;
                case 2:
                    break;
                default:
                    Logger.error(PanHap.class.getSimpleName(), "Execution failed");
                    break;
            }
        }

        System.exit(exitCode);
        System.out.println("Hello World!");
    }
}