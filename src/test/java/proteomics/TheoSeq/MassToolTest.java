package proteomics.TheoSeq;

import org.junit.BeforeClass;
import org.junit.Test;
import proteomics.Types.AA;
import proteomics.Types.SparseBooleanVector;

import java.util.*;

import static org.junit.Assert.*;


public class MassToolTest {

    private static Map<Character, Float> fix_mod_map = new HashMap<>();

    @BeforeClass
    public static void setUp() throws Exception {
        fix_mod_map.put('G', 0f);
        fix_mod_map.put('A', 0f);
        fix_mod_map.put('S', 0f);
        fix_mod_map.put('P', 0f);
        fix_mod_map.put('V', 0f);
        fix_mod_map.put('T', 0f);
        fix_mod_map.put('C', 57.02146f);
        fix_mod_map.put('I', 0f);
        fix_mod_map.put('L', 0f);
        fix_mod_map.put('N', 0f);
        fix_mod_map.put('D', 0f);
        fix_mod_map.put('Q', 0f);
        fix_mod_map.put('K', 0f);
        fix_mod_map.put('E', 0f);
        fix_mod_map.put('M', 0f);
        fix_mod_map.put('H', 0f);
        fix_mod_map.put('F', 0f);
        fix_mod_map.put('R', 0f);
        fix_mod_map.put('Y', 0f);
        fix_mod_map.put('W', 0f);
        fix_mod_map.put('U', 0f);
        fix_mod_map.put('O', 0f);
        fix_mod_map.put('n', 60f);
        fix_mod_map.put('c', 10f);
    }

    @Test
    public void calResidueMass() throws Exception {
        MassTool mass_tool_obj = new MassTool(1, fix_mod_map, "KR", "P", 1.0005f, 0.6f);
        assertEquals(2503.1357421875, mass_tool_obj.calResidueMass(MassTool.seqToAAList("nGASPVTCILNDQKEMHFRYWc")), 0.001);
    }

    @Test
    public void mzToBin() throws Exception {
        MassTool mass_tool_obj = new MassTool(1, fix_mod_map, "KR", "P", 1.0005f, 0.6f);
        assertEquals(11, mass_tool_obj.mzToBin(11), 1e-6);
        assertEquals(0, mass_tool_obj.mzToBin(0), 1e-6);
        assertEquals(-5, mass_tool_obj.mzToBin(-5), 1e-6);
    }

    @Test
    public void buildChainSet() throws Exception {
        // 1 missed-cleavage, N-term linkable
        MassTool mass_tool_obj = new MassTool(1, fix_mod_map, "KR", "P", 1.0005f, 0.6f);
        Set<String> result = mass_tool_obj.buildChainSet("MRGFASSASRIATAAAASKPSLNASTSVNPKLSKTMDYMRIFSVFVVTLWIIRVDARVFKTY");
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
        mass_tool_obj = new MassTool(2, fix_mod_map, "KR", "P", 1.0005f, 0.6f);
        result = mass_tool_obj.buildChainSet("MRGFASSASRIATAAAASKPSLNASTSVNPKLSKTMDYMRIFSVFVVTLWIIRVDARVFKTY");
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
    public void buildChainIonArray() throws Exception {
        MassTool mass_tool_obj = new MassTool(1, fix_mod_map, "KR", "P", 1.0005f, 0.6f);
        float[][] result = mass_tool_obj.buildChainIonArray(mass_tool_obj.seqToAAList("nMKGRc"));
        float[][] ground_truth = new float[][]{
                {132.04776f + 60, 260.1427f + 60, 317.16415f + 60, 473.26526f + 60 + 10},
                {491.27582f + 60 + 10, 360.23535f + 10, 232.1404f + 10, 175.11893f + 10},
                {66.52752f + 30, 130.575f + 30, 159.08572f + 30, 237.13628f + 30 + 5},
                {246.14156f + 30 + 5, 180.6213f + 5, 116.57383f + 5, 88.063095f + 5},
                {44.687435f + 20, 87.38576f + 20, 106.392914f + 20, 158.42662f + 20 + 3.333f},
                {164.43013f + 20 + 3.333f, 120.74997f + 3.333f, 78.05165f + 3.333f, 59.044495f + 3.333f},
                {33.767395f + 15, 65.79114f + 15, 80.0465f + 15, 119.07178f + 15 + 2.5f},
                {123.57442f + 15 + 2.5f, 90.8143f + 2.5f, 58.79056f + 2.5f, 44.535194f + 2.5f},
                {27.215372f + 12, 52.834366f + 12, 64.23866f + 12, 95.458885f + 12 + 2},
                {99.06099f + 12 + 2, 72.85289f + 2, 47.2339f + 2, 35.829605f + 2},
                {22.847357f + 10, 44.196518f + 10, 53.700096f + 10, 79.71695f + 10 + 1.66666f},
                {82.718704f + 10 + 1.66666f, 60.878624f + 1.66666f, 39.529465f + 1.66666f, 30.025887f + 1.66666f}};
        for (int i = 0; i < result.length; ++i) {
            assertArrayEquals(ground_truth[i], result[i], 0.001f);
        }
    }

    @Test
    public void buildVector() throws Exception {
        MassTool mass_tool_obj = new MassTool(1, fix_mod_map, "KR", "P", 1.0005f, 0.6f);
        float[][] ion_matrix = new float[][]{ // nMKGRc
                {132.04776f, 260.1427f, 317.16415f, 473.26526f},
                {491.27582f, 360.23535f, 232.1404f, 175.11893f},
                {66.52752f, 130.575f, 159.08572f, 237.13628f},
                {246.14156f, 180.6213f, 116.57383f, 88.063095f},
                {44.687435f, 87.38576f, 106.392914f, 158.42662f},
                {164.43013f, 120.74997f, 78.05165f, 59.044495f},
                {33.767395f, 65.79114f, 80.0465f, 119.07178f},
                {123.57442f, 90.8143f, 58.79056f, 44.535194f},
                {27.215372f, 52.834366f, 64.23866f, 95.458885f},
                {99.06099f, 72.85289f, 47.2339f, 35.829605f},
                {22.847357f, 44.196518f, 53.700096f, 79.71695f},
                {82.718704f, 60.878624f, 39.529465f, 30.025887f}};

        SparseBooleanVector vector = mass_tool_obj.buildVector(ion_matrix, 1);
        SparseBooleanVector ground_truth = new SparseBooleanVector();
        ground_truth.put(132);
        ground_truth.put(260);
        ground_truth.put(360);
        ground_truth.put(232);
        ground_truth.put(473);
        ground_truth.put(491);
        ground_truth.put(317);
        ground_truth.put(175);
        assertEquals(ground_truth, vector);

        vector = mass_tool_obj.buildVector(ion_matrix, 2);
        assertEquals(ground_truth, vector);

        vector = mass_tool_obj.buildVector(ion_matrix, 3);
        ground_truth = new SparseBooleanVector();
        ground_truth.put(67);
        ground_truth.put(131);
        ground_truth.put(132);
        ground_truth.put(260);
        ground_truth.put(360);
        ground_truth.put(232);
        ground_truth.put(491);
        ground_truth.put(237);
        ground_truth.put(175);
        ground_truth.put(181);
        ground_truth.put(117);
        ground_truth.put(246);
        ground_truth.put(88);
        ground_truth.put(473);
        ground_truth.put(317);
        ground_truth.put(159);
        assertEquals(ground_truth, vector);

        vector = mass_tool_obj.buildVector(ion_matrix, 4);
        ground_truth = new SparseBooleanVector();
        ground_truth.put(67);
        ground_truth.put(131);
        ground_truth.put(132);
        ground_truth.put(260);
        ground_truth.put(164);
        ground_truth.put(360);
        ground_truth.put(232);
        ground_truth.put(106);
        ground_truth.put(491);
        ground_truth.put(237);
        ground_truth.put(45);
        ground_truth.put(78);
        ground_truth.put(175);
        ground_truth.put(181);
        ground_truth.put(117);
        ground_truth.put(246);
        ground_truth.put(87);
        ground_truth.put(88);
        ground_truth.put(473);
        ground_truth.put(121);
        ground_truth.put(59);
        ground_truth.put(317);
        ground_truth.put(158);
        ground_truth.put(159);
        assertEquals(ground_truth, vector);
    }

    @Test
    public void buildPseudoCLIonArray() throws Exception {
        MassTool mass_tool_obj = new MassTool(1, fix_mod_map, "KR", "P", 1.0005f, 0.6f);
        float[][] ion_matrix = new float[][]{ // nMKGRc
                {132.04776f, 260.1427f, 317.16415f, 473.26526f},
                {491.27582f, 360.23535f, 232.1404f, 175.11893f},
                {66.52752f, 130.575f, 159.08572f, 237.13628f},
                {246.14156f, 180.6213f, 116.57383f, 88.063095f},
                {44.687435f, 87.38576f, 106.392914f, 158.42662f},
                {164.43013f, 120.74997f, 78.05165f, 59.044495f},
                {33.767395f, 65.79114f, 80.0465f, 119.07178f},
                {123.57442f, 90.8143f, 58.79056f, 44.535194f},
                {27.215372f, 52.834366f, 64.23866f, 95.458885f},
                {99.06099f, 72.85289f, 47.2339f, 35.829605f},
                {22.847357f, 44.196518f, 53.700096f, 79.71695f},
                {82.718704f, 60.878624f, 39.529465f, 30.025887f}};

        float[][] result = mass_tool_obj.buildPseudoCLIonArray(ion_matrix, 2, 2, 10.5f, 6);
        float[][] ground_truth = new float[][]{
                {132.04776f, 270.6427f, 327.66415f, 483.76526f},
                {501.77582f, 370.73535f, 232.1404f, 175.11893f},
                {66.52752f, 130.575f + 5.25f, 159.08572f + 5.25f, 237.13628f + 5.25f},
                {246.14156f + 5.25f, 180.6213f + 5.25f, 116.57383f, 88.063095f},
                {0, 87.38576f + 3.5f, 106.392914f + 3.5f, 158.42662f + 3.5f},
                {164.43013f + 3.5f, 120.74997f + 3.5f, 0, 0},
                {0, 68.41614f, 82.6715f, 121.69678f},
                {126.19942f, 93.4393f, 0, 0},
                {0, 54.934364f, 66.33866f, 97.55888f},
                {101.16099f, 74.95289f, 0, 0},
                {0, 45.946518f, 55.450096f, 81.46695f},
                {84.468704f, 62.628624f, 0, 0}};
        for (int i = 0; i < result.length; ++i) {
            assertArrayEquals(ground_truth[i], result[i], 0);
        }

        result = mass_tool_obj.buildPseudoCLIonArray(ion_matrix, 2, 2, 10.5f, 6);
        ground_truth = new float[][]{
                {132.04776f, 270.6427f, 327.66415f, 483.76526f},
                {501.77582f, 370.73535f, 232.1404f, 175.11893f},
                {66.52752f, 135.825f, 164.33572f, 242.38628f},
                {251.39156f, 185.8713f, 116.57383f, 88.063095f},
                {0, 87.38576f + 3.5f, 106.392914f + 3.5f, 158.42662f + 3.5f},
                {164.43013f + 3.5f, 120.74997f + 3.5f, 0, 0},
                {0, 68.41614f, 82.6715f, 121.69678f},
                {126.19942f, 93.4393f, 0, 0},
                {0, 54.934364f, 66.33866f, 97.55888f},
                {101.16099f, 74.95289f, 0, 0},
                {0, 45.946518f, 55.450096f, 81.46695f},
                {84.468704f, 62.628624f, 0, 0}};
        for (int i = 0; i < result.length; ++i) {
            assertArrayEquals(ground_truth[i], result[i], 0);
        }
    }

    @Test
    public void digestTrypsin() {
        // 0 missed cleavage
        MassTool mass_tool_obj = new MassTool(0, fix_mod_map, "KR", "P", 1.0005f, 0.6f);
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
        mass_tool_obj = new MassTool(1, fix_mod_map, "KR", "P", 1.0005f, 0.6f);
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
        MassTool mass_tool_obj = new MassTool(1, fix_mod_map, "KR", "P", 1.0005f, 0.6f);
        String seq = "nGHUKc";
        AA[] result = MassTool.seqToAAList(seq);
        AA[] ground_truth = new AA[]{new AA('n', 0), new AA('G', 0), new AA('H', 0), new AA('U', 0), new AA('K', 0), new AA('c', 0)};
        assertArrayEquals(ground_truth, result);

        seq = "nGH[3.02]UKc";
        result = MassTool.seqToAAList(seq);
        ground_truth = new AA[]{new AA('n', 0), new AA('G', 0), new AA('H', 3.02f), new AA('U', 0), new AA('K', 0), new AA('c', 0)};
        assertArrayEquals(ground_truth, result);
    }

    @Test
    public void getModFreeSeq() {
        MassTool mass_tool_obj = new MassTool(1, fix_mod_map, "KR", "P", 1.0005f, 0.6f);
        String seq = "GH[3.02]UK";
        String result = mass_tool_obj.getModFreeSeq(seq);
        String ground_truth = "GHUK";
        assertEquals(ground_truth, result);
    }
}