package edu.tamu.pg.pang2p.Utils;

import java.text.SimpleDateFormat;
import java.util.Date;

/***
 * Logger class to log messages to the console with different log levels
 */
public class Logger {
    public enum LogLevel {
        INFO, WARNING, ERROR, DEBUG, UNKNOWN, PROGRESS
    }

    // Synchronized to prevent multi-threaded log interleaving
    public static synchronized void log(LogLevel logLevel, String className, String message) {
        if (className.length() < 20) {
            className = String.format("%-20s", className);
        }

        String timeString = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss:SSS").format(new Date());
        String logPrefix = timeString + " - " + logLevelToString(logLevel) + " - " + className + " - ";

        if (logLevel == LogLevel.DEBUG && System.getProperty("DEBUG") == null) {
            return;
        }

        if (logLevel == LogLevel.PROGRESS) {
            // \r returns to start of line, \033[K clears everything to the right before writing
            System.out.print("\r\033[K" + logPrefix + message);
            System.out.flush();
        } else {
            // Clear inline progress line first, then print standard newline log
            System.out.println("\r\033[K" + logPrefix + message);
        }

        if (logLevel == LogLevel.ERROR) {
            System.exit(1);
        }
    }

    /**
     * Clear the current inline terminal line manually if needed
     */
    public static synchronized void clearProgress() {
        System.out.print("\r\033[K");
        System.out.flush();
    }

    private static String logLevelToString(LogLevel logLevel) {
        return switch (logLevel) {
            case INFO -> "INFO    ";
            case WARNING -> "WARNING ";
            case ERROR -> "ERROR   ";
            case DEBUG -> "DEBUG   ";
            case PROGRESS -> "PROGRESS";
            default -> "UNKNOWN ";
        };
    }

    public static void info(String className, String message) {
        log(LogLevel.INFO, className, message);
    }

    public static void warning(String className, String message) {
        log(LogLevel.WARNING, className, message);
    }

    public static void error(String className, String message) {
        log(LogLevel.ERROR, className, message);
    }

    public static void progress(String className, String message) {
        log(LogLevel.PROGRESS, className, message);
    }

    public static void debug(String className, String message) {
        log(LogLevel.DEBUG, className, message);
    }
}