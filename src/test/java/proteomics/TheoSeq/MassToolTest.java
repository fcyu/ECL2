package proteomics.TheoSeq;

import org.junit.BeforeClass;
import org.junit.Test;
import proteomics.Types.AA;
import proteomics.Types.SparseBooleanVector;

import java.util.*;

import static org.junit.Assert.*;


public class MassToolTest {

    private static Map<Character, Double> fix_mod_map = new HashMap<>();

    @BeforeClass
    public static void setUp() throws Exception {
        fix_mod_map.put('G', 0d);
        fix_mod_map.put('A', 0d);
        fix_mod_map.put('S', 0d);
        fix_mod_map.put('P', 0d);
        fix_mod_map.put('V', 0d);
        fix_mod_map.put('T', 0d);
        fix_mod_map.put('C', 57.02146);
        fix_mod_map.put('I', 0d);
        fix_mod_map.put('L', 0d);
        fix_mod_map.put('N', 0d);
        fix_mod_map.put('D', 0d);
        fix_mod_map.put('Q', 0d);
        fix_mod_map.put('K', 0d);
        fix_mod_map.put('E', 0d);
        fix_mod_map.put('M', 0d);
        fix_mod_map.put('H', 0d);
        fix_mod_map.put('F', 0d);
        fix_mod_map.put('R', 0d);
        fix_mod_map.put('Y', 0d);
        fix_mod_map.put('W', 0d);
        fix_mod_map.put('U', 0d);
        fix_mod_map.put('O', 0d);
        fix_mod_map.put('n', 60d);
        fix_mod_map.put('c', 10d);
    }

    @Test
    public void calResidueMass() throws Exception {
        MassTool mass_tool_obj = new MassTool(1, fix_mod_map, "KR", "P", 1.0005, 0.6);
        assertEquals(2503.1357421875, mass_tool_obj.calResidueMass("nGASPVTCILNDQKEMHFRYWc"), 0.001);
    }

    @Test
    public void mzToBin() throws Exception {
        MassTool mass_tool_obj = new MassTool(1, fix_mod_map, "KR", "P", 1.0005, 0.6);
        assertEquals(11, mass_tool_obj.mzToBin(11), 1e-6);
        assertEquals(0, mass_tool_obj.mzToBin(0), 1e-6);
        assertEquals(-5, mass_tool_obj.mzToBin(-5), 1e-6);
    }

    @Test
    public void buildChainSet() throws Exception {
        // 1 missed-cleavage, N-term linkable
        MassTool mass_tool_obj = new MassTool(1, fix_mod_map, "KR", "P", 1.0005, 0.6);
        Set<String> result = mass_tool_obj.buildChainSet("MRGFASSASRIATAAAASKPSLNASTSVNPKLSKTMDYMRIFSVFVVTLWIIRVDARVFKTY", (short) 1);
        Set<String> ground_truth = new HashSet<>();
        ground_truth.add("nMRc");
        ground_truth.add("nMRGFASSASRc");
        ground_truth.add("nGFASSASRIATAAAASKPSLNASTSVNPKc");
        ground_truth.add("nIATAAAASKPSLNASTSVNPKc");
        ground_truth.add("nIATAAAASKPSLNASTSVNPKLSKc");
        ground_truth.add("nLSKTMDYMRc");
        ground_truth.add("nVFKTYc");
        assertEquals(ground_truth, result);

        // 2 missed-cleavage, N-term linkable
        mass_tool_obj = new MassTool(2, fix_mod_map, "KR", "P", 1.0005, 0.6);
        result = mass_tool_obj.buildChainSet("MRGFASSASRIATAAAASKPSLNASTSVNPKLSKTMDYMRIFSVFVVTLWIIRVDARVFKTY", (short) 1);
        ground_truth = new HashSet<>();
        ground_truth.add("nMRc");
        ground_truth.add("nIATAAAASKPSLNASTSVNPKc");
        ground_truth.add("nMRGFASSASRc");
        ground_truth.add("nGFASSASRIATAAAASKPSLNASTSVNPKc");
        ground_truth.add("nIATAAAASKPSLNASTSVNPKLSKc");
        ground_truth.add("nLSKTMDYMRc");
        ground_truth.add("nVFKTYc");
        ground_truth.add("nMRGFASSASRIATAAAASKPSLNASTSVNPKc");
        ground_truth.add("nGFASSASRIATAAAASKPSLNASTSVNPKLSKc");
        ground_truth.add("nIATAAAASKPSLNASTSVNPKLSKTMDYMRc");
        ground_truth.add("nLSKTMDYMRIFSVFVVTLWIIRc");
        ground_truth.add("nVDARVFKTYc");
        assertEquals(ground_truth, result);
    }

    @Test
    public void buildTheoVector() throws Exception { // todo: complete

    }

    @Test
    public void digestTrypsin() {
        // 0 missed cleavage
        MassTool mass_tool_obj = new MassTool(0, fix_mod_map, "KR", "P", 1.0005, 0.6);
        Map<Integer, List<int[]>> result = mass_tool_obj.digestTrypsin("FGTRHUYGKPHHYRPHGKHUUG");
        Map<Integer, List<int[]>> ground_truth = new HashMap<>();
        List<int[]> temp = new LinkedList<>();
        temp.add(new int[]{0, 4});
        temp.add(new int[]{4, 18});
        temp.add(new int[]{18, 22});
        ground_truth.put(0, temp);
        assertEquals(ground_truth.size(), result.size());
        for (int k : result.keySet()) {
            for (int i = 0; i < ground_truth.get(k).size(); ++i) {
                assertArrayEquals(ground_truth.get(k).get(i), result.get(k).get(i));
            }
        }

        result = mass_tool_obj.digestTrypsin("FGTRHUYGKPHHYRPHGKHUUR");
        ground_truth = new HashMap<>();
        temp = new LinkedList<>();
        temp.add(new int[]{0, 4});
        temp.add(new int[]{4, 18});
        temp.add(new int[]{18, 22});
        ground_truth.put(0, temp);
        assertEquals(ground_truth.size(), result.size());
        for (int k : result.keySet()) {
            for (int i = 0; i < ground_truth.get(k).size(); ++i) {
                assertArrayEquals(ground_truth.get(k).get(i), result.get(k).get(i));
            }
        }

        // 1 missed cleavage
        mass_tool_obj = new MassTool(1, fix_mod_map, "KR", "P", 1.0005, 0.6);
        result = mass_tool_obj.digestTrypsin("FGTRHUYGKPHHYRPHGKHUUG");
        ground_truth = new HashMap<>();
        temp = new LinkedList<>();
        temp.add(new int[]{0, 4});
        temp.add(new int[]{4, 18});
        temp.add(new int[]{18, 22});
        ground_truth.put(0, temp);
        temp = new LinkedList<>();
        temp.add(new int[]{0, 18});
        ground_truth.put(1, temp);
        assertEquals(ground_truth.size(), result.size());
        for (int k : result.keySet()) {
            for (int i = 0; i < ground_truth.get(k).size(); ++i) {
                assertArrayEquals(ground_truth.get(k).get(i), result.get(k).get(i));
            }
        }

        result = mass_tool_obj.digestTrypsin("FGTRHUYGKPHHYRPHGKHUUR");
        ground_truth = new HashMap<>();
        temp = new LinkedList<>();
        temp.add(new int[]{0, 4});
        temp.add(new int[]{4, 18});
        temp.add(new int[]{18, 22});
        ground_truth.put(0, temp);
        temp = new LinkedList<>();
        temp.add(new int[]{0, 18});
        ground_truth.put(1, temp);
        assertEquals(ground_truth.size(), result.size());
        for (int k : result.keySet()) {
            for (int i = 0; i < ground_truth.get(k).size(); ++i) {
                assertArrayEquals(ground_truth.get(k).get(i), result.get(k).get(i));
            }
        }
    }

    @Test
    public void seqToAAList() {
        MassTool mass_tool_obj = new MassTool(1, fix_mod_map, "KR", "P", 1.0005, 0.6);
        String seq = "nGHUKc";
        AA[] result = MassTool.seqToAAList(seq);
        AA[] ground_truth = new AA[]{new AA('n', 0), new AA('G', 0), new AA('H', 0), new AA('U', 0), new AA('K', 0), new AA('c', 0)};
        assertArrayEquals(ground_truth, result);

        seq = "nGH[3.02]UKc";
        result = MassTool.seqToAAList(seq);
        ground_truth = new AA[]{new AA('n', 0), new AA('G', 0), new AA('H', 3.02), new AA('U', 0), new AA('K', 0), new AA('c', 0)};
        assertArrayEquals(ground_truth, result);
    }
}