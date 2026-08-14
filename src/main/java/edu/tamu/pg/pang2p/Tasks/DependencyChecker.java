package edu.tamu.pg.pang2p.Tasks;

import edu.tamu.pg.pang2p.Utils.Defaults;
import edu.tamu.pg.pang2p.Utils.HelperFunctions;
import edu.tamu.pg.pang2p.Utils.Logger;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class DependencyChecker {

    private final List<String> requiredTools;

    private static final String  CLASS_NAME = DependencyChecker.class.getSimpleName();

    public DependencyChecker() {
        this.requiredTools = Defaults.TOOLS;
    }

    public DependencyChecker(List<String> requiredTools) {
        this.requiredTools = new ArrayList<>(requiredTools);
    }

    public void check() throws IOException {
        List<String> missingTools = new ArrayList<>();

        for (String tool : requiredTools) {
            if (!HelperFunctions.isInstalled(tool)) {
                missingTools.add(tool);
                Logger.info(CLASS_NAME, "Tool \"" + tool + "\" is NOT installed.");
            } else {
                Logger.info(CLASS_NAME, "Tool \"" + tool + "\" is installed.");
            }
        }

        if (!missingTools.isEmpty()) {
            String message = "Missing required tools: " + String.join(", ", missingTools)
                    + ". Please install them and ensure they are available in your system PATH.";
            Logger.error(CLASS_NAME, message);
            throw new IOException(message);
        } else {
            Logger.info(CLASS_NAME, "All required tools are installed.");
        }
    }
}