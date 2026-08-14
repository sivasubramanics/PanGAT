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

    public static void log(LogLevel logLevel, String className, String message) {
        // assuming max className is 20 characters long
        if (className.length() < 20) {
            className = String.format("%-20s", className);
        }

        String timeString = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss:SSS").format(new Date());
        String logPrefix = timeString + " - " + logLevelToString(logLevel) + " - " + className + " - ";

        if (logLevel == LogLevel.DEBUG && System.getProperty("DEBUG") == null) {
            return;
        }

        if (logLevel == LogLevel.PROGRESS) {
            System.out.printf("\r%s", logPrefix + message);
        } else {
            System.out.println(logPrefix + message);
        }

        if (logLevel == LogLevel.ERROR) {
            System.exit(1);
        }
    }

    /**
     * Converts the log level to a string representation
     */
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

    /***
     * Log a message with the INFO log level
     */
    public static void info(String className, String message) {
        log(LogLevel.INFO, className, message);
    }

    /***
     * Log a message with the WARNING log level
     */
    public static void warning(String className, String message) {
        log(LogLevel.WARNING, className, message);
    }

    /***
     * Log a message with the ERROR log level
     */
    public static void error(String className, String message) {
        log(LogLevel.ERROR, className, message);
    }

    /***
     * print progress with inplacement texts
     */
    public static void progress(String className, String message) {
        log(LogLevel.PROGRESS, className, message);
    }

    /***
     * Log a message with the DEBUG log level
     */
    public static void debug(String className, String message) {
        log(LogLevel.DEBUG, className, message);
    }
}
//EOF