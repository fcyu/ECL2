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

    public DbTool(String db_name) {
        String id = "";
        String annotate;
        StringBuilder seq = new StringBuilder(99999);

        boolean new_pro = true;

        Pattern header_pattern = Pattern.compile(">([^\\s]*)(.*)");

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
                    annotate = head_matcher.group(2).trim();
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
}
