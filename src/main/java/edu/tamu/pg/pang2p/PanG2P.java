package edu.tamu.pg.pang2p;

import edu.tamu.pg.pang2p.Tasks.Cohort;
import edu.tamu.pg.pang2p.Tasks.ExportMatrix;
import edu.tamu.pg.pang2p.Tasks.KmerCoverage;
import edu.tamu.pg.pang2p.Tasks.Anchors;
import picocli.CommandLine;
import picocli.CommandLine.*;
import edu.tamu.pg.pang2p.Plugins.*;
import edu.tamu.pg.pang2p.Utils.*;

import java.util.Locale;

@Command(name = "pang2p", description = "pangenome-informed genome-to-phenome framework for trait mapping and prediction", mixinStandardHelpOptions = true, versionProvider = Configs.VersionProvider.class,
        subcommands = {
                // Main Pipeline
                BuildReference.class,
                Genotype.class,
                ExportMatrix.class,

                // Utility Commands
                Anchors.class,
                KmerCoverage.class,
                Cohort.class,
        })
public class PanG2P {
    public static void main(String[] args) {
        Locale.setDefault(Locale.US);
        final PanG2P panG2P = new PanG2P();
        final CommandLine cmd = new CommandLine(panG2P);

        // Custom section renderer for grouped subcommands
        cmd.getHelpSectionMap().put(
                CommandLine.Model.UsageMessageSpec.SECTION_KEY_COMMAND_LIST,
                help ->
                        renderCommands(
                                cmd,
                                "Pipeline Commands:",
                                "build-ref",
                                "genotype",
                                "export-matrix"
                        ) + renderCommands(
                                        cmd,
                                        "Utility Commands:",
                                        "anchors",
                                        "kmer-coverage",
                                        "cohort"
                                )
        );

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
                    Logger.info(PanG2P.class.getSimpleName(), "Total execution time: " + executionTime);
                    break;
                case 2:
                    break;
                default:
                    Logger.error(PanG2P.class.getSimpleName(), "Execution failed");
                    break;
            }
        }

        System.exit(exitCode);
    }

    private static String renderCommands(
            CommandLine cmd,
            String heading,
            String... commandNames) {

        StringBuilder sb = new StringBuilder();

        sb.append(System.lineSeparator())
                .append(heading)
                .append(System.lineSeparator());

        for (String name : commandNames) {
            CommandLine command = cmd.getSubcommands().get(name);

            if (command != null) {
                String description =
                        command.getCommandSpec()
                                .usageMessage()
                                .description()[0];

                sb.append(String.format(
                        "  %-15s %s%n",
                        name,
                        description
                ));
            }
        }

        return sb.toString();
    }
}