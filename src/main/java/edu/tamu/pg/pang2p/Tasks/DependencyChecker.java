package edu.tamu.pg.pang2p.Tasks;

import edu.tamu.pg.pang2p.Utils.Defaults;
import edu.tamu.pg.pang2p.Utils.HelperFunctions;
import edu.tamu.pg.pang2p.Utils.Logger;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

public class DependencyChecker {

    private static final String CLASS_NAME = DependencyChecker.class.getSimpleName();
    private final List<String> requiredTools;

    public DependencyChecker() {
        this(Defaults.TOOLS);
    }

    public DependencyChecker(List<String> requiredTools) {
        Objects.requireNonNull(requiredTools, "requiredTools must not be null");
        this.requiredTools = List.copyOf(requiredTools);
    }

    /**
     * Checks all configured required tools.
     */
    public void check() throws IOException {
        checkTools(this.requiredTools);
    }

    /**
     * Checks a single tool.
     */
    public void check(String tool) throws IOException {
        Objects.requireNonNull(tool, "tool must not be null");
        checkTools(List.of(tool));
    }

    /**
     * Core validation logic for checking a collection of tools.
     */
    public void checkTools(List<String> tools) throws IOException {
        if (tools == null || tools.isEmpty()) {
            return;
        }

        List<String> missingTools = new ArrayList<>();

        for (String tool : tools) {
            if (HelperFunctions.isInstalled(tool)) {
                Logger.info(CLASS_NAME, "Tool \"" + tool + "\" is installed.");
            } else {
                Logger.info(CLASS_NAME, "Tool \"" + tool + "\" is NOT installed.");
                missingTools.add(tool);
            }
        }

        if (!missingTools.isEmpty()) {
            boolean single = missingTools.size() == 1;
            String message = String.format(
                    "Missing required tool%s: %s. Please install %s and ensure %s available in your system PATH.",
                    single ? "" : "s",
                    String.join(", ", missingTools),
                    single ? "it" : "them",
                    single ? "it is" : "they are"
            );

            Logger.error(CLASS_NAME, message);
            throw new IOException(message);
        }

        Logger.info(CLASS_NAME, (tools.size() == 1 ? "Tool is" : "All required tools are") + " installed.");
    }
}