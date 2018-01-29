package proteomics.TheoSeq;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.*;
import java.util.regex.*;

public class DbTool {

    private static final Logger logger = LoggerFactory.getLogger(DbTool.class);

    private Map<String, String> pro_seq_map = new HashMap<>();
    private Map<String, String> pro_annotate_map = new HashMap<>();

    public DbTool(String db_name, String databaseType) {
        String id = "";
        String annotate;
        StringBuilder seq = new StringBuilder(99999);

        boolean new_pro = true;

        Pattern header_pattern;
        if (databaseType.contentEquals("TAIR")) {
            header_pattern = Pattern.compile("^>([^\\s]+)[\\s|]+(.+)$");
        } else if (databaseType.contentEquals("UniProt") || databaseType.contentEquals("SwissProt")) {
            header_pattern = Pattern.compile("^>[^|]+\\|(.+)\\|(.+)$");
        } else if (databaseType.contentEquals("Others")) {
            header_pattern = Pattern.compile("^>(.+)$");
        }
        else {
            header_pattern = null;
            logger.error("Incorrect database type ({}) in the parameter file.", databaseType);
            System.exit(1);
        }

        try (BufferedReader db_reader = new BufferedReader(new FileReader(db_name))) {
            String line;
            while ((line = db_reader.readLine()) != null) {
                line = line.trim();
                Matcher head_matcher = header_pattern.matcher(line);
                if (head_matcher.matches()) {
                    // This line is a header
                    if (!new_pro) {
                        // This isn't the first protein
                        pro_seq_map.put(id, seq.toString());
                    }
                    id = head_matcher.group(1).trim();
                    if (databaseType.contentEquals("Others")) {
                        annotate = id;
                    } else {
                        annotate = head_matcher.group(2).trim();
                    }
                    pro_annotate_map.put(id, annotate);
                    new_pro = true;
                } else if (!line.isEmpty()) {
                    // This line is a body
                    if (new_pro) {
                        seq = new StringBuilder(99999);
                        seq.append(line);
                        new_pro = false;
                    } else {
                        seq.append(line);
                    }
                }
            }
            // Last protein
            pro_seq_map.put(id, seq.toString());
        } catch (IOException | PatternSyntaxException ex) {
            logger.error(ex.toString());
            ex.printStackTrace();
            System.exit(1);
        }
    }

    public Map<String, String> getProSeqMap() {
        return pro_seq_map;
    }

    public Map<String, String> getProAnnotateMap() {
        return pro_annotate_map;
    }

    public Set<Integer> findPeptideLocation(String proteinId, String peptide) throws NullPointerException {
        peptide = peptide.trim().replaceAll("[^A-Z]+", "");
        Set<Integer> output = new HashSet<>();
        int idx = pro_seq_map.get(proteinId).indexOf(peptide);
        while (idx >= 0) {
            output.add(idx);
            idx = pro_seq_map.get(proteinId).indexOf(peptide, idx + 1);
        }
        if (!output.isEmpty()) {
            return output;
        } else {
            throw new NullPointerException(String.format(Locale.US, "Cannot find the peptide %s from the protein %s.", peptide, proteinId));
        }
    }

    public static Set<String> reduceProteinIdSet(Set<String> input) {
        if (input.size() == 1) {
            return input;
        } else {
            Map<String, Integer> tempMap = new HashMap<>();
            for (String s : input) {
                String[] tempArray = s.split("\\.");
                if (tempMap.containsKey(tempArray[0])) {
                    if (tempMap.get(tempArray[0]) > Integer.valueOf(tempArray[1])) {
                        tempMap.put(tempArray[0], Integer.valueOf(tempArray[1]));
                    }
                } else {
                    tempMap.put(tempArray[0], Integer.valueOf(tempArray[1]));
                }
            }
            Set<String> output = new HashSet<>();
            for (String s : tempMap.keySet()) {
                output.add(s + "." + tempMap.get(s));
            }
            return output;
        }
    }
}
