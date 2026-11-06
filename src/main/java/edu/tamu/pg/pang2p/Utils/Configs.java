package edu.tamu.pg.pang2p.Utils;

import picocli.CommandLine;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;
/***
 * Enum to store the configuration values for the KCF file format
 */
public enum Configs {
    VERSION(getVersion()),
    SOURCE("PanG2P");

    private final String value;

    Configs(String value) {
        this.value = value;
    }

    public String getValue() {
        return value;
    }

    public static String getVersion() {
        try (InputStream in = Configs.class.getResourceAsStream("/version.properties")) {
            Properties props = new Properties();
            props.load(in);
            return props.getProperty("version", "UNKNOWN");
        } catch (IOException e) {
            return "UNKNOWN";
        }
    }

    public static class VersionProvider implements CommandLine.IVersionProvider {
        @Override
        public String[] getVersion() {
            return new String[]{ Configs.getVersion() };
        }
    }

}
//EOF