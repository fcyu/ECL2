package proteomics.Index;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.TheoSeq.DbTool;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.*;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class BuildIndex {

    private static final Logger logger = LoggerFactory.getLogger(BuildIndex.class);
    private static final Pattern varModParamPattern = Pattern.compile("([0-9.-]+)\\s+([A-Znc]+)\\s+([01])");
    private static final int globalVarModMaxNum = 5; // Do not change this value. Otherwise, change generateLocalIdxModMassMap accordingly.
    private static final float varModMassResolution = 0.01f;

    public final float linker_mass;

    private final MassTool mass_tool_obj;
    private final Map<String, String> pro_annotate_map;
    private Map<Character, Float> fix_mod_map = new HashMap<>(25, 1);
    private TreeMap<Integer, Set<String>> bin_seq_map = new TreeMap<>();
    private Map<Integer, Long> bin_candidate_num_map = new HashMap<>();
    private Map<String, ChainEntry> seq_entry_map = new HashMap<>();
    private Map<String, boolean[]> seq_term_map = new HashMap<>();
    private Set<String> for_check_duplicate = new HashSet<>();
    private Map<String, Set<String>> seqProMap;
    private final float ms1_bin_size ;

    public BuildIndex(Map<String, String> parameter_map) {
        // initialize parameters
        int min_chain_length = Integer.valueOf(parameter_map.get("min_chain_length")) + 2; // n and c are counted in the sequence
        int max_chain_length = Integer.valueOf(parameter_map.get("max_chain_length")) + 2; // n and c are counted in the sequence
        String db_path = parameter_map.get("db");
        int missed_cleavage = Integer.valueOf(parameter_map.get("missed_cleavage"));
        float mz_bin_size = Float.valueOf(parameter_map.get("mz_bin_size"));
        float one_minus_bin_offset = 1 - Float.valueOf(parameter_map.get("mz_bin_offset"));
        float max_precursor_mass = Float.valueOf(parameter_map.get("max_precursor_mass"));

        if (parameter_map.containsKey("ms1_bin_size")) {
            ms1_bin_size = Float.valueOf(parameter_map.get("ms1_bin_size"));
        } else {
            ms1_bin_size = 0.001f;
        }

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

        // generate seq_pro_map
        seqProMap = buildSeqProMap(pro_seq_map, min_chain_length, max_chain_length);

        // read var mods
        Set<VarModParam> varModParamSet = new HashSet<>(30, 1);
        Set<BinaryModParam> binaryModParamSet = new HashSet<>(10, 1);
        for (String k : parameter_map.keySet()) {
            if (k.contentEquals("var_mod1")) {
                getVarModParams(parameter_map.get(k), varModParamSet, binaryModParamSet);
            } else if (k.contentEquals("var_mod2")) {
                getVarModParams(parameter_map.get(k), varModParamSet, binaryModParamSet);
            } else if (k.contentEquals("var_mod3")) {
                getVarModParams(parameter_map.get(k), varModParamSet, binaryModParamSet);
            } else if (k.contentEquals("var_mod4")) {
                getVarModParams(parameter_map.get(k), varModParamSet, binaryModParamSet);
            } else if (k.contentEquals("var_mod5")) {
                getVarModParams(parameter_map.get(k), varModParamSet, binaryModParamSet);
            } else if (k.contentEquals("var_mod6")) {
                getVarModParams(parameter_map.get(k), varModParamSet, binaryModParamSet);
            } else if (k.contentEquals("var_mod7")) {
                getVarModParams(parameter_map.get(k), varModParamSet, binaryModParamSet);
            } else if (k.contentEquals("var_mod8")) {
                getVarModParams(parameter_map.get(k), varModParamSet, binaryModParamSet);
            } else if (k.contentEquals("var_mod9")) {
                getVarModParams(parameter_map.get(k), varModParamSet, binaryModParamSet);
            }
        }


        int varModMaxNum = Math.min(globalVarModMaxNum, Integer.valueOf(parameter_map.get("var_mod_max_num")));

        // generate all peptide entries
        for (String seq : seqProMap.keySet()) {
            boolean proteinNTerm = seq_term_map.get(seq)[0];
            boolean proteinCTerm = seq_term_map.get(seq)[1];

            // mod free
            Set<Short> linkSiteSet = getLinkSiteSet(seq, proteinNTerm, proteinCTerm);
            if (!linkSiteSet.isEmpty()) {
                float totalMass = (float) (mass_tool_obj.calResidueMass(seq) + MassTool.H2O);
                if (totalMass < max_precursor_mass - linker_mass) {
                    int bin = massToBin(totalMass);
                    if (bin_seq_map.containsKey(bin)) {
                        bin_seq_map.get(bin).add(seq);
                    } else {
                        Set<String> temp = new HashSet<>();
                        temp.add(seq);
                        bin_seq_map.put(bin, temp);
                    }
                    ChainEntry chainEntry = new ChainEntry(seq, totalMass, linkSiteSet, proteinNTerm, proteinCTerm, "0".hashCode());
                    seq_entry_map.put(seq, chainEntry);
                }
            }

            // mod containing
            Set<VarSequence> varSeqSet = generateModSeq(seq, linkSiteSet, varModParamSet, binaryModParamSet, varModMaxNum);
            for (VarSequence varSeq : varSeqSet) {
                linkSiteSet = new HashSet<>(5, 1);
                linkSiteSet.add(varSeq.linkSite);
                if (!linkSiteSet.isEmpty()) {
                    float totalMass = (float) (mass_tool_obj.calResidueMass(varSeq.seq) + MassTool.H2O);
                    if (totalMass < max_precursor_mass - linker_mass) {
                        int bin = massToBin(totalMass);
                        if (bin_seq_map.containsKey(bin)) {
                            bin_seq_map.get(bin).add(varSeq.seq);
                        } else {
                            Set<String> temp = new HashSet<>();
                            temp.add(varSeq.seq);
                            bin_seq_map.put(bin, temp);
                        }
                        ChainEntry chainEntry = new ChainEntry(varSeq.seq, totalMass, linkSiteSet, proteinNTerm, proteinCTerm, varSeq.binaryModType);
                        if (seq_entry_map.containsKey(varSeq.seq)) {
                            // Binary mod has the higher priority
                            if (seq_entry_map.get(varSeq.seq).binaryModType == "0".hashCode()) {
                                seq_entry_map.put(varSeq.seq, chainEntry);
                            }
                        } else {
                            seq_entry_map.put(varSeq.seq, chainEntry);
                        }
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
    }

    public Map<String, Set<String>> getSeqProMap() {
        return seqProMap;
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

    public Map<Integer, Long> getBinCandidateNumMap() {
        return bin_candidate_num_map;
    }

    private Map<String, Set<String>> buildSeqProMap(Map<String, String> pro_seq_map, int min_chain_length, int max_chain_length) {
        Map<String, Set<String>> seq_pro_map = new HashMap<>(pro_seq_map.size() * 150, 1);
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
                if (pro_seq.startsWith(target_seq.substring(1, target_seq.length() - 1))) {
                    n_term = true;
                }
                if (pro_seq.endsWith(target_seq.substring(1, target_seq.length() - 1))) {
                    c_term = true;
                }
                seq_term_map.put(target_seq, new boolean[]{n_term, c_term});

                if (seq_pro_map.containsKey(target_seq)) {
                    seq_pro_map.get(target_seq).add(pro_id);
                } else {
                    Set<String> pro_list = new HashSet<>(5, 1);
                    pro_list.add(pro_id);
                    seq_pro_map.put(target_seq, pro_list);
                }
            }
        }

        // generate decoy seq
        for (String pro_id : pro_seq_map.keySet()) {
            String pro_seq = pro_seq_map.get(pro_id);
            String decoy_pro_seq = (new StringBuilder(pro_seq)).reverse().toString();
            Set<String> decoy_seq_set = mass_tool_obj.buildChainSet(decoy_pro_seq);
            for (String decoy_seq : decoy_seq_set) {
                if ((decoy_seq.length() < min_chain_length) || (decoy_seq.length() > max_chain_length) || decoy_seq.contains("B") || decoy_seq.contains("J") || decoy_seq.contains("X") || decoy_seq.contains("Z")) {
                    continue;
                }

                // Check duplicate
                String new_decoy_seq = decoy_seq.replace("L", "!").replace("I", "!").replace("K", "#").replace("Q", "#");
                if (for_check_duplicate.contains(new_decoy_seq)) {
                    // the decoy sequence is the same as the target sequence
                    continue;
                }

                for_check_duplicate.add(new_decoy_seq);

                boolean n_term = false;
                boolean c_term = false;
                if (decoy_pro_seq.startsWith(decoy_seq.substring(1, decoy_seq.length() - 1))) {
                    n_term = true;
                }
                if (decoy_pro_seq.endsWith(decoy_seq.substring(1, decoy_seq.length() - 1))) {
                    c_term = true;
                }
                seq_term_map.put(decoy_seq, new boolean[]{n_term, c_term});

                if (seq_pro_map.containsKey(decoy_seq)) {
                    seq_pro_map.get(decoy_seq).add("DECOY_" + pro_id);
                } else {
                    Set<String> pro_list = new HashSet<>(5, 1);
                    pro_list.add("DECOY_" + pro_id);
                    seq_pro_map.put(decoy_seq, pro_list);
                }
            }
        }

        return seq_pro_map;
    }

    private Set<VarSequence> generateModSeq(String seq, Set<Short> modFreeListSites, Set<VarModParam> varModParamSet, Set<BinaryModParam> binaryModParamSet, int varModMaxNum) { // todo: check
        Set<VarSequence> varSeqSet = new HashSet<>();
        for (short linkSite : modFreeListSites) {
            // has binary mod
            for (BinaryModParam binaryModParam : binaryModParamSet) {
                // get all locations having binary mod
                Map<Integer, List<Float>> idxBinaryModMassMap = new HashMap<>(seq.length(), 1);
                for (int i = 0; i < seq.length(); ++i) {
                    if (i != linkSite) {
                        String aa = seq.substring(i, i + 1);
                        if (binaryModParam.aas.contains(aa)) {
                            List<Float> tempList = new LinkedList<>();
                            tempList.add(binaryModParam.modMass);
                            idxBinaryModMassMap.put(i, tempList);
                        }
                    }
                }
                if (!idxBinaryModMassMap.isEmpty()) {
                    // generate a sequence only containing the binary mod
                    StringBuilder sb = new StringBuilder(seq.length() * 10);
                    for (int i = 0; i < seq.length(); ++i) {
                        sb.append(seq.substring(i, i + 1));
                        if (idxBinaryModMassMap.containsKey(i)) {
                            sb.append(String.format("[%.2f]", idxBinaryModMassMap.get(i).get(0)));
                        }
                    }
                    varSeqSet.add(new VarSequence(sb.toString(), linkSite, binaryModParam.hashCode()));

                    if (idxBinaryModMassMap.size() < varModMaxNum) {
                        // generate sequences containing the binary mod and additional var mod
                        // get all locations having var mods
                        Map<Integer, List<Float>> idxVarModMassMap = new HashMap<>(seq.length(), 1);
                        for (int i = 0; i < seq.length(); ++i) {
                            if (i != linkSite) {
                                if (!idxBinaryModMassMap.containsKey(i)) {
                                    char aa = seq.charAt(i);
                                    for (VarModParam varModParam : varModParamSet) {
                                        if (varModParam.aa == aa) {
                                            if (idxVarModMassMap.containsKey(i)) {
                                                idxVarModMassMap.get(i).add(varModParam.modMass);
                                            } else {
                                                List<Float> temp = new LinkedList<>();
                                                temp.add(varModParam.modMass);
                                                idxVarModMassMap.put(i, temp);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        if (!idxVarModMassMap.isEmpty()) {
                            // generate var containing sequences
                            Map<Integer, List<Float>> idxBinaryVarModMassMap = new HashMap<>(seq.length(), 1);
                            idxBinaryVarModMassMap.putAll(idxBinaryModMassMap);
                            idxBinaryVarModMassMap.putAll(idxVarModMassMap);
                            Integer[] allIdxArray = idxVarModMassMap.keySet().toArray(new Integer[idxVarModMassMap.size()]);
                            Arrays.sort(allIdxArray);
                            for (int i = 1; i <= Math.min(varModMaxNum - idxBinaryModMassMap.size(), idxVarModMassMap.size()); ++i) {
                                List<int[]> idxCombinationList = generateIdxCombinations(allIdxArray, i);
                                Set<VarSequence> varSetSubSet = new HashSet<>(seq.length(), 1);
                                for (int[] idxCombination : idxCombinationList) {
                                    int[] allIdxCombination = new int[idxCombination.length + idxBinaryModMassMap.size()];
                                    int j = 0;
                                    for (int idx : idxBinaryModMassMap.keySet()) {
                                        allIdxCombination[j] = idx;
                                        ++j;
                                    }
                                    for (int k = 0; k < idxCombination.length; ++k) {
                                        allIdxCombination[j + k] = idxCombination[k];
                                    }
                                    Arrays.sort(allIdxCombination);
                                    varSetSubSet.addAll(generateModSeqSub(seq, allIdxCombination, idxBinaryVarModMassMap, linkSite, binaryModParam.hashCode()));
                                }
                                if (!varSetSubSet.isEmpty()) {
                                    varSeqSet.addAll(varSetSubSet); // eliminate those sequence that the middle amino acids having the same mod mass and the n-term and the first amino acid or the c-term and the last amino acid have the same mod mass.
                                }
                            }
                        }
                    }
                }
            }

            // does not have binary mod
            // get all locations' var lists
            Map<Integer, List<Float>> idxVarModMassMap = new HashMap<>(seq.length(), 1);
            for (int i = 0; i < seq.length(); ++i) {
                if (i != linkSite) {
                    char aa = seq.charAt(i);
                    for (VarModParam varModParam : varModParamSet) {
                        if (varModParam.aa == aa) {
                            if (idxVarModMassMap.containsKey(i)) {
                                idxVarModMassMap.get(i).add(varModParam.modMass);
                            } else {
                                List<Float> temp = new LinkedList<>();
                                temp.add(varModParam.modMass);
                                idxVarModMassMap.put(i, temp);
                            }
                        }
                    }
                }
            }
            if (!idxVarModMassMap.isEmpty()) {
                // generate var containing sequences
                Integer[] allIdxArray = idxVarModMassMap.keySet().toArray(new Integer[idxVarModMassMap.size()]);
                Arrays.sort(allIdxArray);
                for (int i = 1; i <= Math.min(varModMaxNum, idxVarModMassMap.size()); ++i) {
                    List<int[]> idxCombinationList = generateIdxCombinations(allIdxArray, i);
                    for (int[] idxCombination : idxCombinationList) {
                        varSeqSet.addAll(generateModSeqSub(seq, idxCombination, idxVarModMassMap, linkSite, "0".hashCode()));
                    }
                }
            }
        }

        return varSeqSet;
    }

    private List<int[]> generateIdxCombinations(Integer[] allIdxArray, int num) {
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

    private Set<VarSequence> generateModSeqSub(String seq, int[] idxCombination, Map<Integer, List<Float>> idxModMassMap, short linkSite, int binaryModType) {
        List<Map<Integer, Float>> localIdxModMassMaps = generateLocalIdxModMassMap(idxCombination, idxModMassMap);

        Set<VarSequence> outputSet = new HashSet<>();
        for (Map<Integer, Float> localIdxModMassMap : localIdxModMassMaps) {
            StringBuilder sb = new StringBuilder(seq.length() * 10);
            for (int i = 0; i < seq.length(); ++i) {
                sb.append(seq.charAt(i));
                if (localIdxModMassMap.containsKey(i)) {
                    sb.append(String.format("[%.2f]", localIdxModMassMap.get(i)));
                }
            }
            outputSet.add(new VarSequence(sb.toString(), linkSite, binaryModType));
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
                                Map<Integer, Float> localIdxModMassMap = new HashMap<>(6, 1);
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
                            Map<Integer, Float> localIdxModMassMap = new HashMap<>(5, 1);
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
                        Map<Integer, Float> localIdxModMassMap = new HashMap<>(4, 1);
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
                    Map<Integer, Float> localIdxModMassMap = new HashMap<>(3, 1);
                    localIdxModMassMap.put(idxArray[0], idxModMassMap.get(idxArray[0]).get(i0));
                    localIdxModMassMap.put(idxArray[1], idxModMassMap.get(idxArray[1]).get(i1));
                    outputList.add(localIdxModMassMap);
                }
            }
        } else if (idxArray.length == 1) {
            for (int i0 = 0; i0 < idxModMassMap.get(idxArray[0]).size(); ++i0) {
                Map<Integer, Float> localIdxModMassMap = new HashMap<>(2, 1);
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

    private Set<Short> getLinkSiteSet(String seq, boolean n_term, boolean c_term) {
        AA[] aa_list = MassTool.seqToAAList(seq);
        Set<Short> output = new HashSet<>(5, 1);
        for (int i = 1; i < aa_list.length - 2; ++i) {
            if (aa_list[i].aa == 'K' && (Math.abs(aa_list[i].delta_mass) < varModMassResolution)) {
                output.add((short) i);
            }
        }
        if (n_term && !output.contains((short) 1) && (Math.abs(aa_list[0].delta_mass) < varModMassResolution)) {
            output.add((short) 0);
        }
        if (c_term && aa_list[aa_list.length - 2].aa == 'K' && (Math.abs(aa_list[aa_list.length - 2].delta_mass) < varModMassResolution)) {
            output.add((short) (aa_list.length - 2));
        }
        return output;
    }

    private void getVarModParams(String v, Set<VarModParam> varModParamSet, Set<BinaryModParam> binaryModParamSet) {
        Matcher varModMatcher = varModParamPattern.matcher(v);
        if (varModMatcher.matches()) {
            float modMass = Float.valueOf(varModMatcher.group(1));
            String aas = varModMatcher.group(2);
            boolean isBinary = varModMatcher.group(3).contentEquals("1");
            if (Math.abs(modMass) > varModMassResolution) {
                if (isBinary) {
                    binaryModParamSet.add(new BinaryModParam(modMass, aas));
                } else {
                    for (int i = 0; i < aas.length(); ++i) {
                        if (!fix_mod_map.containsKey(aas.charAt(i))) {
                            varModParamSet.add(new VarModParam(modMass, aas.charAt(i)));
                        }
                    }
                }
            }
        } else {
            logger.error("Cannot parse variable modification parameter from {}.", v);
            System.exit(1);
        }
    }
}
