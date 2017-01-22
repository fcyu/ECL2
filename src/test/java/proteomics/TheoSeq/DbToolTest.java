package proteomics.TheoSeq;

import org.junit.BeforeClass;
import org.junit.Test;

import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.assertEquals;

public class DbToolTest {
    private static DbTool db_tool_obj;

    @BeforeClass
    public static void setUp() throws Exception {
        db_tool_obj = new DbTool(Thread.currentThread().getContextClassLoader().getResource("test.fasta").getPath());
    }

    @Test
    public void returnSeqMap() throws Exception {
        Map<String, String> pro_seq_map = db_tool_obj.getProSeqMap();
        Map<String, String> ground_truth = new HashMap<>();
        ground_truth.put("Pro1", "ASRIATAAAASKPSLNKF");
        ground_truth.put("Pro2", "STSVNPKLSKT");
        for (String k : pro_seq_map.keySet()) {
            assertEquals(pro_seq_map.get(k), ground_truth.get(k));
        }
    }

    @Test
    public void returnAnnotateMap() throws Exception {
        Map<String, String> pro_annotate_map = db_tool_obj.getProAnnotateMap();
        Map<String, String> ground_truth = new HashMap<>();
        ground_truth.put("Pro1", "test protein one");
        ground_truth.put("Pro2", "test protein two");
        ground_truth.put("Pro3", "test protein three");
        for (String k : pro_annotate_map.keySet()) {
            assertEquals(pro_annotate_map.get(k), ground_truth.get(k));
        }
    }
}