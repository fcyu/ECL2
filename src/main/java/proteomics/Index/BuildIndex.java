package proteomics.Index;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.ECL2;
import proteomics.TheoSeq.DbTool;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.AA;
import proteomics.Types.ChainEntry;
import proteomics.Types.VarModParam;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class BuildIndex {

    private static final Logger logger = LoggerFactory.getLogger(BuildIndex.class);
    private static final float ms1_bin_size = 0.001f;
    private static final Pattern varModParamPattern = Pattern.compile("([0-9.-]+)\\s+([A-Znc]+)");
    private static final int globalVarModMaxNum = 5; // Do not change this value. Otherwise, change generateLocalIdxModMassMap accordingly.
    private static final float varModMassResolution = 0.01f;

    public final float linker_mass;

    private final MassTool mass_tool_obj;
    private final Map<String, String> pro_annotate_map;
    private Map<Character, Float> fix_mod_map = new HashMap<>();
    private TreeMap<Integer, Set<String>> bin_seq_map = new TreeMap<>();
    private Map<Integer, Long> bin_candidate_num_map = new HashMap<>();
    private Map<String, ChainEntry> seq_entry_map = new HashMap<>();
    private Map<String, boolean[]> seq_term_map = new HashMap<>();
    private TreeMap<Float, Set<String>> uniprot_decoy_mass_seq_map = new TreeMap<>();
    private Set<String> for_check_duplicate = new HashSet<>();

    public BuildIndex(Map<String, String> parameter_map) {
        // initialize parameters
        int min_chain_length = Integer.valueOf(parameter_map.get("min_chain_length")) + 2; // n and c are counted in the sequence
        int max_chain_length = Integer.valueOf(parameter_map.get("max_chain_length")) + 2; // n and c are counted in the sequence
        String db_path = parameter_map.get("db");
        String uniprot_decoy_db = parameter_map.get("decoy_pool_db");
        int missed_cleavage = Integer.valueOf(parameter_map.get("missed_cleavage"));
        float mz_bin_size = Float.valueOf(parameter_map.get("mz_bin_size"));
        float one_minus_bin_offset = 1 - Float.valueOf(parameter_map.get("mz_bin_offset"));
        float max_precursor_mass = Float.valueOf(parameter_map.get("max_precursor_mass"));

        // Read fix modification
        fix_mod_map.put('G', Float.valueOf(parameter_map.get("G")));
        fix_mod_map.put('A', Float.valueOf(parameter_map.get("A")));
        fix_mod_map.put('S', Float.valueOf(parameter_map.get("S")));
        fix_mod_map.put('P', Float.valueOf(parameter_map.get("P")));
        fix_mod_map.put('V', Float.valueOf(parameter_map.get("V")));
        fix_mod_map.put('T', Float.valueOf(parameter_map.get("T")));
        fix_mod_map.put('C', Float.valueOf(parameter_map.get("C")));
        fix_mod_map.put('I', Float.valueOf(parameter_map.get("I")));
        fix_mod_map.put('L', Float.valueOf(parameter_map.get("L")));
        fix_mod_map.put('N', Float.valueOf(parameter_map.get("N")));
        fix_mod_map.put('D', Float.valueOf(parameter_map.get("D")));
        fix_mod_map.put('Q', Float.valueOf(parameter_map.get("Q")));
        fix_mod_map.put('K', Float.valueOf(parameter_map.get("K")));
        fix_mod_map.put('E', Float.valueOf(parameter_map.get("E")));
        fix_mod_map.put('M', Float.valueOf(parameter_map.get("M")));
        fix_mod_map.put('H', Float.valueOf(parameter_map.get("H")));
        fix_mod_map.put('F', Float.valueOf(parameter_map.get("F")));
        fix_mod_map.put('R', Float.valueOf(parameter_map.get("R")));
        fix_mod_map.put('Y', Float.valueOf(parameter_map.get("Y")));
        fix_mod_map.put('W', Float.valueOf(parameter_map.get("W")));
        fix_mod_map.put('U', Float.valueOf(parameter_map.get("U")));
        fix_mod_map.put('O', Float.valueOf(parameter_map.get("O")));
        fix_mod_map.put('n', Float.valueOf(parameter_map.get("n")));
        fix_mod_map.put('c', Float.valueOf(parameter_map.get("c")));

        if (Math.abs(fix_mod_map.get('K') - fix_mod_map.get('n')) > 1e-6) { // todo improve
            logger.error("K and N-term have different fixed modification. Exit.");
            System.exit(1);
        }

        linker_mass = Float.valueOf(parameter_map.get("cl_mass")) - 2 * fix_mod_map.get('K'); // todo improve

        // read protein database
        DbTool db_tool_obj = new DbTool(db_path);
        Map<String, String> pro_seq_map = db_tool_obj.getProSeqMap();
        pro_annotate_map = db_tool_obj.getProAnnotateMap();

        // define a new MassTool object
        mass_tool_obj = new MassTool(missed_cleavage, fix_mod_map, "KR", "P", mz_bin_size, one_minus_bin_offset);
        Map<String, Float> mass_table = mass_tool_obj.getMassTable();

        // generate seq_pro_map
        Map<String, Set<String>> seq_pro_map = buildSeqProMap(pro_seq_map, min_chain_length, max_chain_length);

        // read var mods
        Set<VarModParam> varModParams = new HashSet<>();
        for (String k : parameter_map.keySet()) {
            if (k.contentEquals("var_mod1")) {
                varModParams.addAll(getVarModParams(parameter_map.get(k)));
            } else if (k.contentEquals("var_mod2")) {
                varModParams.addAll(getVarModParams(parameter_map.get(k)));
            } else if (k.contentEquals("var_mod3")) {
                varModParams.addAll(getVarModParams(parameter_map.get(k)));
            } else if (k.contentEquals("var_mod4")) {
                varModParams.addAll(getVarModParams(parameter_map.get(k)));
            } else if (k.contentEquals("var_mod5")) {
                varModParams.addAll(getVarModParams(parameter_map.get(k)));
            } else if (k.contentEquals("var_mod6")) {
                varModParams.addAll(getVarModParams(parameter_map.get(k)));
            } else if (k.contentEquals("var_mod7")) {
                varModParams.addAll(getVarModParams(parameter_map.get(k)));
            } else if (k.contentEquals("var_mod8")) {
                varModParams.addAll(getVarModParams(parameter_map.get(k)));
            } else if (k.contentEquals("var_mod9")) {
                varModParams.addAll(getVarModParams(parameter_map.get(k)));
            }
        }


        int varModMaxNum = Math.min(globalVarModMaxNum, Integer.valueOf(parameter_map.get("var_mod_max_num")));

        // generate all peptide entries
        for (String seq : seq_pro_map.keySet()) {
            boolean proteinNTerm = seq_term_map.get(seq)[0];
            boolean proteinCTerm = seq_term_map.get(seq)[1];
            Set<String> proIdSet = seq_pro_map.get(seq);

            // mod free
            Set<Short> linkSiteSet = getLinkSiteSet(seq, proteinNTerm);
            if (!linkSiteSet.isEmpty()) {
                float totalMass = mass_tool_obj.calResidueMass(aaList) + mass_table.get("H2O");
                AA[] aaArray = MassTool.seqToAAList(seq);
                if (totalMass < max_precursor_mass - linker_mass) {
                    int bin = massToBin(totalMass);
                    if (bin_seq_map.containsKey(bin)) {
                        bin_seq_map.get(bin).add(seq);
                    } else {
                        Set<String> temp = new HashSet<>();
                        temp.add(seq);
                        bin_seq_map.put(bin, temp);
                    }
                    float[][] chainIonArray = mass_tool_obj.buildChainIonArray(aaList);
                    ChainEntry chainEntry = new ChainEntry(seq, totalMass, proIdSet, linkSiteSet, chainIonArray, proteinNTerm, proteinCTerm);
                    seq_entry_map.put(seq, chainEntry);
                }
            }

            // mod containing
            Set<String> varSeqSet = generateModSeq(seq, linkSiteSet, varModParams, varModMaxNum);
            for (String varSeq : varSeqSet) {
                linkSiteSet = getLinkSiteSet(varSeq, proteinNTerm);
                if (!linkSiteSet.isEmpty()) {
                    AA[] aaArray = MassTool.seqToAAList(varSeq);
                    float totalMass = mass_tool_obj.calResidueMass(aaArray) + mass_table.get("H2O");
                    if (totalMass < max_precursor_mass - linker_mass) {
                        int bin = massToBin(totalMass);
                        if (bin_seq_map.containsKey(bin)) {
                            bin_seq_map.get(bin).add(varSeq);
                        } else {
                            Set<String> temp = new HashSet<>();
                            temp.add(varSeq);
                            bin_seq_map.put(bin, temp);
                        }
                        float[][] chainIonArray = mass_tool_obj.buildChainIonArray(aaArray);
                        ChainEntry chainEntry = new ChainEntry(varSeq, totalMass, proIdSet, linkSiteSet, chainIonArray, proteinNTerm, proteinCTerm);
                        seq_entry_map.put(varSeq, chainEntry);
                    }
                }
            }
        }

        // summarize candidate numbers in each bin
        for (int bin_index : bin_seq_map.keySet()) {
            long candidate_num = 0;
            for (String seq : bin_seq_map.get(bin_index)) {
                candidate_num += seq_entry_map.get(seq).link_site_set.size();
            }
            bin_candidate_num_map.put(bin_index, candidate_num);
        }

        // read and index large decoy uniprot database
        if (ECL2.cal_evalue) {
            db_tool_obj = new DbTool(uniprot_decoy_db);
            Map<String, String> uniprot_decoy_pro_seq_map = db_tool_obj.getProSeqMap();
            for (String pro_seq : uniprot_decoy_pro_seq_map.values()) {
                Set<String> seq_set = mass_tool_obj.buildChainSet(pro_seq);
                for (String mod_free_seq : seq_set) {
                    if ((mod_free_seq.length() >= min_chain_length) && (mod_free_seq.length() <= max_chain_length)) {
                        if (!mod_free_seq.contains("B") && !mod_free_seq.contains("J") && !mod_free_seq.contains("X") && !mod_free_seq.contains("Z")) {
                            String temp_str = mod_free_seq.replace("L", "!").replace("I", "!").replace("K", "#").replace("Q", "#");
                            if (!for_check_duplicate.contains(temp_str)) {
                                for_check_duplicate.add(temp_str);
                                int temp_int = mod_free_seq.indexOf('K'); // only use the first K
                                if ((temp_int != -1) && (temp_int != mod_free_seq.length() - 2)) { // only consider chains whose link site is in the middle.
                                    // consider PTM free
                                    float mass = mass_tool_obj.calResidueMass(MassTool.seqToAAList(mod_free_seq)) + mass_table.get("H2O");
                                    if (mass < max_precursor_mass - linker_mass) {
                                        if (uniprot_decoy_mass_seq_map.containsKey(mass)) {
                                            if (uniprot_decoy_mass_seq_map.get(mass).size() < ECL2.random_score_point_t) { // limit the size of set for the sake of memory.
                                                uniprot_decoy_mass_seq_map.get(mass).add(mod_free_seq + "-" + temp_int);
                                            }
                                        } else {
                                            Set<String> temp = new HashSet<>();
                                            temp.add(mod_free_seq + "-" + temp_int);
                                            uniprot_decoy_mass_seq_map.put(mass, temp);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    public MassTool returnMassTool() {
        return mass_tool_obj;
    }

    public Map<String, String> getProAnnotateMap() {
        return pro_annotate_map;
    }

    public Map<Character, Float> getFixModMap() {
        return fix_mod_map;
    }

    public TreeMap<Integer, Set<String>> getMassBinSeqMap() {
        return bin_seq_map;
    }

    public Map<String, ChainEntry> getSeqEntryMap() {
        return seq_entry_map;
    }

    public float binToLeftMass(int bin_idx) {
        return bin_idx * ms1_bin_size;
    }

    public float binToRightMass(int bin_idx) {
        return (bin_idx + 1) * ms1_bin_size;
    }

    public int massToBin(float mass) {
        return (int) Math.floor(mass / ms1_bin_size);
    }

    public TreeMap<Float, Set<String>> getUniprotDecoyMassSeqMap() {
        return uniprot_decoy_mass_seq_map;
    }

    public Map<Integer, Long> getBinCandidateNumMap() {
        return bin_candidate_num_map;
    }

    private Map<String, Set<String>> buildSeqProMap(Map<String, String> pro_seq_map, int min_chain_length, int max_chain_length) {
        Map<String, Set<String>> seq_pro_map = new HashMap<>();
        for (String pro_id : pro_seq_map.keySet()) {
            String pro_seq = pro_seq_map.get(pro_id);
            Set<String> seq_set = mass_tool_obj.buildChainSet(pro_seq);
            for (String target_seq : seq_set) {
                if ((target_seq.length() < min_chain_length) || (target_seq.length() > max_chain_length) || target_seq.contains("B") || target_seq.contains("J") || target_seq.contains("X") || target_seq.contains("Z")) {
                    continue;
                }

                // Add the sequence to the check set for decoy duplicate check
                for_check_duplicate.add(target_seq.replace("L", "!").replace("I", "!").replace("K", "#").replace("Q", "#")); // "L" and "I"; "K"and "Q" have the same mass

                boolean n_term = false;
                boolean c_term = false;
                if (pro_seq.startsWith(getNCFreeSeq(target_seq))) {
                    n_term = true;
                }
                if (pro_seq.endsWith(getNCFreeSeq(target_seq))) {
                    c_term = true;
                }
                seq_term_map.put(target_seq, new boolean[]{n_term, c_term});

                if (seq_pro_map.containsKey(target_seq)) {
                    seq_pro_map.get(target_seq).add(pro_id);
                } else {
                    Set<String> pro_list = new HashSet<>();
                    pro_list.add(pro_id);
                    seq_pro_map.put(target_seq, pro_list);
                }
            }
        }

        Map<String, Set<String>> temp_map = new HashMap<>();
        for (String target_seq : seq_pro_map.keySet()) {
            String decoy_seq = "n" + (new StringBuilder(target_seq.substring(1, target_seq.length() - 2))).reverse().toString() + target_seq.substring(target_seq.length() - 2);

            // Check duplicate
            String new_decoy_seq = decoy_seq.replace("L", "!").replace("I", "!").replace("K", "#").replace("Q", "#");
            if (for_check_duplicate.contains(new_decoy_seq)) {
                // the decoy sequence is the same as the target sequence
                continue;
            }

            for_check_duplicate.add(new_decoy_seq);

            Set<String> target_id_set = seq_pro_map.get(target_seq);
            String target_pro_seq = pro_seq_map.get(target_id_set.iterator().next());
            boolean n_term = false;
            boolean c_term = false;
            if (target_pro_seq.startsWith(getNCFreeSeq(target_seq))) {
                n_term = true;
            }
            if (target_pro_seq.endsWith(getNCFreeSeq(target_seq))) {
                c_term = true;
            }
            seq_term_map.put(decoy_seq, new boolean[]{n_term, c_term});

            Set<String> decoy_id_set = new HashSet<>();
            for (String target_id : target_id_set)
            {
                decoy_id_set.add("DECOY_" + target_id);
            }
            temp_map.put(decoy_seq, decoy_id_set);
        }

        seq_pro_map.putAll(temp_map);

        return seq_pro_map;
    }

    private Set<String> generateModSeq(String seq, Set<Short> modFreeListSites, Set<VarModParam> varModParams, int varModMaxNum) { // todo: check
        // get all locations' var lists
        Map<Integer, List<Float>> idxModMassMap = new HashMap<>();
        for (int i = 0; i < seq.length(); ++i) {
            char aa = seq.charAt(i);
            for (VarModParam varModParam : varModParams) {
                if (varModParam.aa == aa) {
                    if (idxModMassMap.containsKey(i)) {
                        idxModMassMap.get(i).add(varModParam.modMass);
                    } else {
                        List<Float> temp = new LinkedList<>();
                        temp.add(varModParam.modMass);
                        idxModMassMap.put(i, temp);
                    }
                }
            }
        }

        // generate var containing sequences
        Set<String> varSeqSet = new HashSet<>();
        Integer[] allIdxArray = idxModMassMap.keySet().toArray(new Integer[idxModMassMap.size()]);
        Arrays.sort(allIdxArray);
        for (int i = 1; i <= Math.min(varModMaxNum, idxModMassMap.size()); ++i) {
            List<int[]> idxArrays = selectIdxes(allIdxArray, i);
            Set<String> varSetSubSet = new HashSet<>();
            for (int[] idxArray : idxArrays) {
                if (checkIdxArray(idxArray, modFreeListSites)) {
                    varSetSubSet.addAll(generateModSeqSub(seq, idxArray, idxModMassMap));
                }
            }
            if (!varSetSubSet.isEmpty()) {
                varSeqSet.addAll(checkKCTermMod(varSetSubSet)); // eliminate those sequence that the middle amino acids having the same mod mass and the n-term and the first amino acid or the c-term and the last amino acid have the same mod mass.
            }
        }

        return varSeqSet;
    }

    private List<int[]> selectIdxes(Integer[] allIdxArray, int num) {
        List<int[]> outputList = new LinkedList<>();
        Iterator<int[]> tempIterator = CombinatoricsUtils.combinationsIterator(allIdxArray.length, num);
        while (tempIterator.hasNext()) {
            int[] idxIdxArray = tempIterator.next();
            int[] idxArray = new int[idxIdxArray.length];
            for (int i = 0; i < idxIdxArray.length; ++i) {
                idxArray[i] = allIdxArray[idxIdxArray[i]];
            }
            Arrays.sort(idxArray);
            outputList.add(idxArray);
        }

        return outputList;
    }

    private boolean checkIdxArray(int[] idxArray, Set<Short> modFreeLinkSites) {
        for (int modFreeLinkSite : modFreeLinkSites) {
            if (Arrays.binarySearch(idxArray, modFreeLinkSite) < 0) {
                return true;
            }
        }
        return false;
    }

    private Set<String> generateModSeqSub(String seq, int[] idxArray, Map<Integer, List<Float>> idxModMassMap) {
        List<Map<Integer, Float>> localIdxModMassMaps = generateLocalIdxModMassMap(idxArray, idxModMassMap);

        Set<String> outputSet = new HashSet<>();
        for (Map<Integer, Float> localIdxModMassMap : localIdxModMassMaps) {
            StringBuilder sb = new StringBuilder(seq.length() * 10);
            for (int i = 0; i < seq.length(); ++i) {
                sb.append(seq.charAt(i));
                if (localIdxModMassMap.containsKey(i)) {
                    sb.append(String.format("[%.2f]", localIdxModMassMap.get(i)));
                }
            }
            outputSet.add(sb.toString());
        }

        return outputSet;
    }

    private List<Map<Integer, Float>> generateLocalIdxModMassMap(int[] idxArray, Map<Integer, List<Float>> idxModMassMap) {
        List<Map<Integer, Float>> outputList = new LinkedList<>();
        if (idxArray.length == globalVarModMaxNum) {
            for (int i0 = 0; i0 < idxModMassMap.get(idxArray[0]).size(); ++i0) {
                for (int i1 = 0; i1 < idxModMassMap.get(idxArray[1]).size(); ++i1) {
                    for (int i2 = 0; i2 < idxModMassMap.get(idxArray[2]).size(); ++i2) {
                        for (int i3 = 0; i3 < idxModMassMap.get(idxArray[3]).size(); ++i3) {
                            for (int i4 = 0; i4 < idxModMassMap.get(idxArray[4]).size(); ++i4) {
                                Map<Integer, Float> localIdxModMassMap = new HashMap<>();
                                localIdxModMassMap.put(idxArray[0], idxModMassMap.get(idxArray[0]).get(i0));
                                localIdxModMassMap.put(idxArray[1], idxModMassMap.get(idxArray[1]).get(i1));
                                localIdxModMassMap.put(idxArray[2], idxModMassMap.get(idxArray[2]).get(i2));
                                localIdxModMassMap.put(idxArray[3], idxModMassMap.get(idxArray[3]).get(i3));
                                localIdxModMassMap.put(idxArray[4], idxModMassMap.get(idxArray[4]).get(i4));
                                outputList.add(localIdxModMassMap);
                            }
                        }
                    }
                }
            }
        } else if (idxArray.length == 4) {
            for (int i0 = 0; i0 < idxModMassMap.get(idxArray[0]).size(); ++i0) {
                for (int i1 = 0; i1 < idxModMassMap.get(idxArray[1]).size(); ++i1) {
                    for (int i2 = 0; i2 < idxModMassMap.get(idxArray[2]).size(); ++i2) {
                        for (int i3 = 0; i3 < idxModMassMap.get(idxArray[3]).size(); ++i3) {
                            Map<Integer, Float> localIdxModMassMap = new HashMap<>();
                            localIdxModMassMap.put(idxArray[0], idxModMassMap.get(idxArray[0]).get(i0));
                            localIdxModMassMap.put(idxArray[1], idxModMassMap.get(idxArray[1]).get(i1));
                            localIdxModMassMap.put(idxArray[2], idxModMassMap.get(idxArray[2]).get(i2));
                            localIdxModMassMap.put(idxArray[3], idxModMassMap.get(idxArray[3]).get(i3));
                            outputList.add(localIdxModMassMap);
                        }
                    }
                }
            }
        } else if (idxArray.length == 3) {
            for (int i0 = 0; i0 < idxModMassMap.get(idxArray[0]).size(); ++i0) {
                for (int i1 = 0; i1 < idxModMassMap.get(idxArray[1]).size(); ++i1) {
                    for (int i2 = 0; i2 < idxModMassMap.get(idxArray[2]).size(); ++i2) {
                        Map<Integer, Float> localIdxModMassMap = new HashMap<>();
                        localIdxModMassMap.put(idxArray[0], idxModMassMap.get(idxArray[0]).get(i0));
                        localIdxModMassMap.put(idxArray[1], idxModMassMap.get(idxArray[1]).get(i1));
                        localIdxModMassMap.put(idxArray[2], idxModMassMap.get(idxArray[2]).get(i2));
                        outputList.add(localIdxModMassMap);
                    }
                }
            }
        } else if (idxArray.length == 2) {
            for (int i0 = 0; i0 < idxModMassMap.get(idxArray[0]).size(); ++i0) {
                for (int i1 = 0; i1 < idxModMassMap.get(idxArray[1]).size(); ++i1) {
                    Map<Integer, Float> localIdxModMassMap = new HashMap<>();
                    localIdxModMassMap.put(idxArray[0], idxModMassMap.get(idxArray[0]).get(i0));
                    localIdxModMassMap.put(idxArray[1], idxModMassMap.get(idxArray[1]).get(i1));
                    outputList.add(localIdxModMassMap);
                }
            }
        } else if (idxArray.length == 1) {
            for (int i0 = 0; i0 < idxModMassMap.get(idxArray[0]).size(); ++i0) {
                Map<Integer, Float> localIdxModMassMap = new HashMap<>();
                localIdxModMassMap.put(idxArray[0], idxModMassMap.get(idxArray[0]).get(i0));
                outputList.add(localIdxModMassMap);
            }
        }

        return outputList;
    }

    private Set<String> checkKCTermMod(Set<String> varSeqSet) { // eliminate those sequence that the middle amino acids having the same mod mass and the n-term and the first amino acid or the c-term and the last amino acid have the same mod mass. todo: check
        String[] varSeqArray = varSeqSet.toArray(new String[varSeqSet.size()]);
        Arrays.sort(varSeqArray); // Make sure that nK[].... is before n[]K..., so that n[]K... will be kept.
        int seqLength = MassTool.seqToAAList(varSeqArray[0]).length;
        AA[][] aaArrays = new AA[varSeqArray.length][seqLength];
        for (int i = 0; i < varSeqArray.length; ++i) {
            aaArrays[i] = MassTool.seqToAAList(varSeqArray[i]);
        }

        if (aaArrays.length > 1) {
            Set<String> keptSeqSet = new HashSet<>();
            for (int i = 0; i < aaArrays.length - 1; ++i) {
                boolean keep = true;
                for (int j = i + 1; j < aaArrays.length; ++j) {
                    if ((Math.abs(aaArrays[i][0].delta_mass - aaArrays[j][1].delta_mass) < varModMassResolution) && (Math.abs(aaArrays[i][1].delta_mass - aaArrays[j][0].delta_mass) < varModMassResolution) && (Math.abs(aaArrays[i][seqLength - 2].delta_mass - aaArrays[j][seqLength - 1].delta_mass) < varModMassResolution) && (Math.abs(aaArrays[i][seqLength - 1].delta_mass - aaArrays[j][seqLength - 2].delta_mass) < varModMassResolution)) {
                        keep = false;
                        for (int k = 2; k < seqLength - 2; ++k) {
                            if (Math.abs(aaArrays[i][k].delta_mass - aaArrays[j][k].delta_mass) > varModMassResolution) {
                                keep = true;
                                break;
                            }
                        }
                        if (!keep) {
                            break;
                        }
                    }
                }
                if (keep) {
                    keptSeqSet.add(varSeqArray[i]);
                }
            }
            keptSeqSet.add(varSeqArray[varSeqArray.length - 1]);
            return keptSeqSet;
        } else {
            return varSeqSet;
        }
    }

    private Set<Short> getLinkSiteSet(String seq, boolean n_term) {
        AA[] aa_list = MassTool.seqToAAList(seq);
        Set<Short> output = new HashSet<>();
        for (int i = 1; i < aa_list.length - 2; ++i) {
            if (aa_list[i].aa == 'K' && (Math.abs(aa_list[i].delta_mass) < varModMassResolution)) {
                output.add((short) i);
            }
        }
        if (n_term && !output.contains((short) 1) && (Math.abs(aa_list[0].delta_mass) < varModMassResolution)) {
            output.add((short) 0);
        }
        return output;
    }

    private String getNCFreeSeq(String seq) {
        return seq.substring(1, seq.length() - 1);
    }

    private Set<VarModParam> getVarModParams(String v) {
        Set<VarModParam> varModParams = new HashSet<>();

        Matcher varModMatcher = varModParamPattern.matcher(v);
        if (varModMatcher.matches()) {
            float modMass = Float.valueOf(varModMatcher.group(1));
            String aas = varModMatcher.group(2);
            if (Math.abs(modMass) < varModMassResolution) {
                return  varModParams;
            }
            for (int i = 0; i < aas.length(); ++i) {
                varModParams.add(new VarModParam(modMass, aas.charAt(i)));
            }
        } else {
            logger.error("Cannot parse variable modification parameter from {}.", v);
            System.exit(1);
        }

        return varModParams;
    }
}
