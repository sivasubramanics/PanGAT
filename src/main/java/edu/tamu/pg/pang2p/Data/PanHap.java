package edu.tamu.pg.pang2p.Data;

import java.io.*;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;

public class PanHap {

    ConcurrentHashMap<String, HapData> hapData = new ConcurrentHashMap<>();
    String[] samples;

    public PanHap(Path panClustersRepFile){
        List<HapData> hapDataList = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(panClustersRepFile.toFile()))){
            String line;
            int lineNumber = 0;
            while ((line = reader.readLine()) != null){
                if (line.trim().isEmpty()) continue;
                if (lineNumber == 0){
                    lineNumber++;
                    continue;
                }
                String[] parts = line.split("\t");
                String clusterId = parts[0];
                String hapID = parts[1];
                String newGenomeID = parts[2];
                String newChromID = parts[3];
                String orgGenomeID = parts[4];
                String orgChromID = parts[5];
                int length = Integer.parseInt(parts[6]);
                int start = Integer.parseInt(parts[7]);
                int end = Integer.parseInt(parts[8]);
                hapData.computeIfAbsent(hapID, new HapData((clusterId, hapID, newGenomeID, newChromID, orgGenomeID, orgChromID, start, end, length));
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}