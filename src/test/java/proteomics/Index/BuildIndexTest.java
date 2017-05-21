package proteomics.Index;

import org.junit.BeforeClass;
import org.junit.Test;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.ChainEntry;

import java.util.*;

import static org.junit.Assert.*;

public class BuildIndexTest {

    private static Map<String, String> parameter_map = new HashMap<>();

    @BeforeClass
    public static void setUp() throws Exception {
        // init parameter_map
        parameter_map.put("max_precursor_mass", "5000");
        parameter_map.put("min_chain_length", "5");
        parameter_map.put("max_chain_length", "50");
        parameter_map.put("db", Thread.currentThread().getContextClassLoader().getResource("test.fasta").getPath());
        parameter_map.put("decoy_pool_db", Thread.currentThread().getContextClassLoader().getResource("test.decoy.pool.fasta").getPath());
        parameter_map.put("missed_cleavage", "2");
        parameter_map.put("cl_mass", "138.0680796");
        parameter_map.put("mz_bin_size", "1.0005");
        parameter_map.put("mz_bin_offset", "0.4");
        parameter_map.put("var_mod1", "28.031300 Kn");
        parameter_map.put("var_mod2", "34.063117 Kn");
        parameter_map.put("var_mod3", "0.0 X");
        parameter_map.put("var_mod_max_num", "5");
        parameter_map.put("G", "0");
        parameter_map.put("A", "0");
        parameter_map.put("S", "0");
        parameter_map.put("P", "0");
        parameter_map.put("V", "0");
        parameter_map.put("T", "0");
        parameter_map.put("C", "57.02146");
        parameter_map.put("I", "0");
        parameter_map.put("L", "0");
        parameter_map.put("N", "0");
        parameter_map.put("D", "0");
        parameter_map.put("Q", "0");
        parameter_map.put("K", "0");
        parameter_map.put("E", "0");
        parameter_map.put("M", "0");
        parameter_map.put("H", "0");
        parameter_map.put("F", "0");
        parameter_map.put("R", "0");
        parameter_map.put("Y", "0");
        parameter_map.put("W", "0");
        parameter_map.put("U", "0");
        parameter_map.put("O", "0");
        parameter_map.put("n", "0");
        parameter_map.put("c", "0");
    }

    @Test
    public void getMassBinSeqMap() throws Exception {
        BuildIndex build_idx_obj = new BuildIndex(parameter_map);
        TreeMap<Integer, Set<String>> result = build_idx_obj.getMassBinSeqMap();
        TreeMap<Integer, Set<String>> ground_truth = new TreeMap<>();

        Set<String> temp;
        temp = new HashSet<>();
        temp.add("nIATAAAASKPSLNKc");
        temp.add("nNLSPKSAAAATAIKc");
        ground_truth.put(1341761, temp);

        temp = new HashSet<>();
        temp.add("n[28.03]IATAAAASKPSLNKc");
        temp.add("nIATAAAASKPSLNK[28.03]c");
        temp.add("n[28.03]NLSPKSAAAATAIKc");
        temp.add("nNLSPKSAAAATAIK[28.03]c");
        ground_truth.put(1369791, temp);

        temp = new HashSet<>();
        temp.add("n[34.06]IATAAAASKPSLNKc");
        temp.add("nIATAAAASKPSLNK[34.06]c");
        temp.add("n[34.06]NLSPKSAAAATAIKc");
        temp.add("nNLSPKSAAAATAIK[34.06]c");
        ground_truth.put(1375821, temp);

        temp = new HashSet<>();
        temp.add("n[28.03]NLSPKSAAAATAIK[28.03]c");
        temp.add("n[28.03]IATAAAASKPSLNK[28.03]c");
        ground_truth.put(1397821, temp);

        temp = new HashSet<>();
        temp.add("n[28.03]IATAAAASKPSLNK[34.06]c");
        temp.add("n[34.06]IATAAAASKPSLNK[28.03]c");
        temp.add("n[28.03]NLSPKSAAAATAIK[34.06]c");
        temp.add("n[34.06]NLSPKSAAAATAIK[28.03]c");
        ground_truth.put(1403851, temp);

        temp = new HashSet<>();
        temp.add("n[34.06]IATAAAASKPSLNK[34.06]c");
        temp.add("n[34.06]NLSPKSAAAATAIK[34.06]c");
        ground_truth.put(1409881, temp);


        temp = new HashSet<>();
        temp.add("nIATAAAASKPSLNKFc");
        temp.add("nKNLSPKSAAAATAIFc");
        ground_truth.put(1488830, temp);

        temp = new HashSet<>();
        temp.add("nIATAAAASK[28.03]PSLNKFc");
        temp.add("nIATAAAASKPSLNK[28.03]Fc");
        temp.add("n[28.03]IATAAAASKPSLNKFc");
        temp.add("nKNLSPK[28.03]SAAAATAIFc");
        temp.add("n[28.03]KNLSPKSAAAATAIFc");
        ground_truth.put(1516860, temp);

        temp = new HashSet<>();
        temp.add("nIATAAAASK[34.06]PSLNKFc");
        temp.add("n[34.06]IATAAAASKPSLNKFc");
        temp.add("nIATAAAASKPSLNK[34.06]Fc");
        temp.add("nKNLSPK[34.06]SAAAATAIFc");
        temp.add("n[34.06]KNLSPKSAAAATAIFc");
        ground_truth.put(1522890, temp);

        temp = new HashSet<>();
        temp.add("n[28.03]KNLSPK[28.03]SAAAATAIFc");
        temp.add("n[28.03]IATAAAASK[28.03]PSLNKFc");
        temp.add("n[28.03]IATAAAASKPSLNK[28.03]Fc");
        temp.add("n[28.03]K[28.03]NLSPKSAAAATAIFc");
        ground_truth.put(1544890, temp);

        temp = new HashSet<>();
        temp.add("n[28.03]KNLSPK[34.06]SAAAATAIFc");
        temp.add("n[28.03]IATAAAASK[34.06]PSLNKFc");
        temp.add("n[34.06]K[28.03]NLSPKSAAAATAIFc");
        temp.add("n[28.03]IATAAAASKPSLNK[34.06]Fc");
        temp.add("n[34.06]IATAAAASK[28.03]PSLNKFc");
        temp.add("n[34.06]IATAAAASKPSLNK[28.03]Fc");
        temp.add("n[34.06]KNLSPK[28.03]SAAAATAIFc");
        ground_truth.put(1550920, temp);

        temp = new HashSet<>();
        temp.add("n[34.06]KNLSPK[34.06]SAAAATAIFc");
        temp.add("n[34.06]K[34.06]NLSPKSAAAATAIFc");
        temp.add("n[34.06]IATAAAASK[34.06]PSLNKFc");
        temp.add("n[34.06]IATAAAASKPSLNK[34.06]Fc");
        ground_truth.put(1556950, temp);


        temp = new HashSet<>();
        temp.add("nASRIATAAAASKPSLNKc");
        temp.add("nNLSPKSAAAATAIRSAKc");
        ground_truth.put(1655931, temp);

        temp = new HashSet<>();
        temp.add("n[28.03]ASRIATAAAASKPSLNKc");
        temp.add("nASRIATAAAASKPSLNK[28.03]c");
        temp.add("nASRIATAAAASK[28.03]PSLNKc");
        temp.add("nNLSPKSAAAATAIRSAK[28.03]c");
        temp.add("n[28.03]NLSPKSAAAATAIRSAKc");
        temp.add("nNLSPK[28.03]SAAAATAIRSAKc");
        ground_truth.put(1683961, temp);

        temp = new HashSet<>();
        temp.add("n[34.06]ASRIATAAAASKPSLNKc");
        temp.add("nASRIATAAAASKPSLNK[34.06]c");
        temp.add("nASRIATAAAASK[34.06]PSLNKc");
        temp.add("nNLSPKSAAAATAIRSAK[34.06]c");
        temp.add("n[34.06]NLSPKSAAAATAIRSAKc");
        temp.add("nNLSPK[34.06]SAAAATAIRSAKc");
        ground_truth.put(1689991, temp);

        temp = new HashSet<>();
        temp.add("n[28.03]ASRIATAAAASKPSLNK[28.03]c");
        temp.add("nASRIATAAAASK[28.03]PSLNK[28.03]c");
        temp.add("nNLSPK[28.03]SAAAATAIRSAK[28.03]c");
        temp.add("n[28.03]NLSPKSAAAATAIRSAK[28.03]c");
        ground_truth.put(1711991, temp);

        temp = new HashSet<>();
        temp.add("nASRIATAAAASK[28.03]PSLNK[34.06]c");
        temp.add("n[34.06]ASRIATAAAASKPSLNK[28.03]c");
        temp.add("nASRIATAAAASK[34.06]PSLNK[28.03]c");
        temp.add("nNLSPK[28.03]SAAAATAIRSAK[34.06]c");
        temp.add("nNLSPK[34.06]SAAAATAIRSAK[28.03]c");
        temp.add("n[28.03]NLSPKSAAAATAIRSAK[34.06]c");
        temp.add("n[28.03]ASRIATAAAASKPSLNK[34.06]c");
        temp.add("n[34.06]NLSPKSAAAATAIRSAK[28.03]c");
        ground_truth.put(1718021, temp);

        temp = new HashSet<>();
        temp.add("n[34.06]NLSPKSAAAATAIRSAK[34.06]c");
        temp.add("nNLSPK[34.06]SAAAATAIRSAK[34.06]c");
        temp.add("n[34.06]ASRIATAAAASKPSLNK[34.06]c");
        temp.add("nASRIATAAAASK[34.06]PSLNK[34.06]c");
        ground_truth.put(1724051, temp);


        temp = new HashSet<>();
        temp.add("nASRIATAAAASKPSLNKFc");
        temp.add("nKNLSPKSAAAATAIRSAFc");
        ground_truth.put(1803000, temp);

        temp = new HashSet<>();
        temp.add("n[28.03]ASRIATAAAASKPSLNKFc");
        temp.add("nKNLSPK[28.03]SAAAATAIRSAFc");
        temp.add("nASRIATAAAASK[28.03]PSLNKFc");
        temp.add("nASRIATAAAASKPSLNK[28.03]Fc");
        temp.add("n[28.03]KNLSPKSAAAATAIRSAFc");
        ground_truth.put(1831030, temp);

        temp = new HashSet<>();
        temp.add("nKNLSPK[34.06]SAAAATAIRSAFc");
        temp.add("nASRIATAAAASKPSLNK[34.06]Fc");
        temp.add("nASRIATAAAASK[34.06]PSLNKFc");
        temp.add("n[34.06]KNLSPKSAAAATAIRSAFc");
        temp.add("n[34.06]ASRIATAAAASKPSLNKFc");
        ground_truth.put(1837060, temp);

        temp = new HashSet<>();
        temp.add("n[28.03]K[28.03]NLSPKSAAAATAIRSAFc");
        temp.add("n[28.03]ASRIATAAAASKPSLNK[28.03]Fc");
        temp.add("n[28.03]ASRIATAAAASK[28.03]PSLNKFc");
        temp.add("n[28.03]KNLSPK[28.03]SAAAATAIRSAFc");
        temp.add("nASRIATAAAASK[28.03]PSLNK[28.03]Fc");
        ground_truth.put(1859060, temp);

        temp = new HashSet<>();
        temp.add("n[34.06]ASRIATAAAASK[28.03]PSLNKFc");
        temp.add("n[34.06]KNLSPK[28.03]SAAAATAIRSAFc");
        temp.add("n[28.03]ASRIATAAAASKPSLNK[34.06]Fc");
        temp.add("nASRIATAAAASK[34.06]PSLNK[28.03]Fc");
        temp.add("n[28.03]KNLSPK[34.06]SAAAATAIRSAFc");
        temp.add("nASRIATAAAASK[28.03]PSLNK[34.06]Fc");
        temp.add("n[28.03]ASRIATAAAASK[34.06]PSLNKFc");
        temp.add("n[34.06]K[28.03]NLSPKSAAAATAIRSAFc");
        temp.add("n[34.06]ASRIATAAAASKPSLNK[28.03]Fc");
        ground_truth.put(1865090, temp);

        temp = new HashSet<>();
        temp.add("n[34.06]K[34.06]NLSPKSAAAATAIRSAFc");
        temp.add("n[34.06]ASRIATAAAASK[34.06]PSLNKFc");
        temp.add("nASRIATAAAASK[34.06]PSLNK[34.06]Fc");
        temp.add("n[34.06]ASRIATAAAASKPSLNK[34.06]Fc");
        temp.add("n[34.06]KNLSPK[34.06]SAAAATAIRSAFc");
        ground_truth.put(1871120, temp);

        assertEquals(ground_truth, result);
    }

    @Test
    public void getSeqEntryMap() throws Exception {
        BuildIndex build_idx_obj = new BuildIndex(parameter_map);
        MassTool mass_tool_obj = build_idx_obj.returnMassTool();
        Map<Character, Float> mass_table = mass_tool_obj.getMassTable();
        Map<String, ChainEntry> result = build_idx_obj.getSeqEntryMap();
        Map<String, ChainEntry> group_truth = new HashMap<>();

        Set<Short> link_site_set;

        String seq = "nASRIATAAAASKPSLNKc";
        link_site_set = new HashSet<>();
        link_site_set.add((short) 0);
        link_site_set.add((short) 12);
        ChainEntry chain_entry = new ChainEntry(seq, MassTool.seqToAAList(seq), mass_tool_obj.calResidueMass(MassTool.seqToAAList(seq)) + MassTool.H2O, link_site_set, true, false);
        group_truth.put(seq, chain_entry);

        seq = "nIATAAAASKPSLNK[28.03]c";
        link_site_set = new HashSet<>();
        link_site_set.add((short) 9);
        chain_entry = new ChainEntry(seq, MassTool.seqToAAList(seq), mass_tool_obj.calResidueMass(MassTool.seqToAAList(seq)) + MassTool.H2O, link_site_set, false, false);
        group_truth.put(seq, chain_entry);

        seq = "nNLSPK[34.06]SAAAATAIRSAK[34.06]c";
        link_site_set = new HashSet<>();
        link_site_set.add((short) 0);
        chain_entry = new ChainEntry(seq, MassTool.seqToAAList(seq), mass_tool_obj.calResidueMass(MassTool.seqToAAList(seq)) + MassTool.H2O, link_site_set, true, false);
        group_truth.put(seq, chain_entry);

        seq = "n[28.03]KNLSPK[28.03]SAAAATAIRSAFc";
        link_site_set = new HashSet<>();
        link_site_set.add((short) 1);
        chain_entry = new ChainEntry(seq, MassTool.seqToAAList(seq), mass_tool_obj.calResidueMass(MassTool.seqToAAList(seq)) + MassTool.H2O, link_site_set, true, true);
        group_truth.put(seq, chain_entry);

        for (String temp : group_truth.keySet()) {
            assertEquals(group_truth.get(temp), result.get(temp));
        }
    }
}