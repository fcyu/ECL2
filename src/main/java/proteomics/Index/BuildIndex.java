package proteomics.Index;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import ProteomicsLibrary.*;
import ProteomicsLibrary.Types.*;
import proteomics.Types.*;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class BuildIndex {

    private static final Logger logger = LoggerFactory.getLogger(BuildIndex.class);
    private static final Pattern varModParamPattern = Pattern.compile("([0-9.-]+)\\s+([A-Znc]+)\\s+([01])");
    private static final int globalVarModMaxNum = 5; // Do not change this value. Otherwise, change generateLocalIdxModMassMap accordingly.
    private static final double varModMassResolution = 0.01;

    public final double linker_mass;
    public final short linker_type;

    private final MassTool mass_tool_obj;
    private final Map<String, String> pro_annotate_map;
    private Map<Character, Double> fix_mod_map = new HashMap<>(25, 1);
    private TreeMap<Integer, Set<String>> bin_seq_map = new TreeMap<>();
    private Map<String, ChainEntry> seq_entry_map = new HashMap<>();
    private Map<String, Set<String>> seqProMap;
    private final double ms1_bin_size;
    private final double inverseMs1BinSize;

    public BuildIndex(Map<String, String> parameter_map) throws Exception {
        // initialize parameters
        int min_chain_length = Integer.valueOf(parameter_map.get("min_chain_length")) + 2; // n and c are counted in the sequence
        int max_chain_length = Integer.valueOf(parameter_map.get("max_chain_length")) + 2; // n and c are counted in the sequence
        String db_path = parameter_map.get("db");
        int missed_cleavage = Integer.valueOf(parameter_map.get("missed_cleavage"));
        double mz_bin_size = Double.valueOf(parameter_map.get("mz_bin_size"));
        double one_minus_bin_offset = 1 - Double.valueOf(parameter_map.get("mz_bin_offset"));
        ms1_bin_size = Double.valueOf(parameter_map.get("ms1_bin_size"));
        inverseMs1BinSize = 1 / ms1_bin_size;

        // Read fix modification
        fix_mod_map.put('G', Double.valueOf(parameter_map.get("G")));
        fix_mod_map.put('A', Double.valueOf(parameter_map.get("A")));
        fix_mod_map.put('S', Double.valueOf(parameter_map.get("S")));
        fix_mod_map.put('P', Double.valueOf(parameter_map.get("P")));
        fix_mod_map.put('V', Double.valueOf(parameter_map.get("V")));
        fix_mod_map.put('T', Double.valueOf(parameter_map.get("T")));
        fix_mod_map.put('C', Double.valueOf(parameter_map.get("C")));
        fix_mod_map.put('I', Double.valueOf(parameter_map.get("I")));
        fix_mod_map.put('L', Double.valueOf(parameter_map.get("L")));
        fix_mod_map.put('N', Double.valueOf(parameter_map.get("N")));
        fix_mod_map.put('D', Double.valueOf(parameter_map.get("D")));
        fix_mod_map.put('Q', Double.valueOf(parameter_map.get("Q")));
        fix_mod_map.put('K', Double.valueOf(parameter_map.get("K")));
        fix_mod_map.put('E', Double.valueOf(parameter_map.get("E")));
        fix_mod_map.put('M', Double.valueOf(parameter_map.get("M")));
        fix_mod_map.put('H', Double.valueOf(parameter_map.get("H")));
        fix_mod_map.put('F', Double.valueOf(parameter_map.get("F")));
        fix_mod_map.put('R', Double.valueOf(parameter_map.get("R")));
        fix_mod_map.put('Y', Double.valueOf(parameter_map.get("Y")));
        fix_mod_map.put('W', Double.valueOf(parameter_map.get("W")));
        fix_mod_map.put('U', Double.valueOf(parameter_map.get("U")));
        fix_mod_map.put('O', Double.valueOf(parameter_map.get("O")));
        fix_mod_map.put('n', Double.valueOf(parameter_map.get("n")));
        fix_mod_map.put('c', Double.valueOf(parameter_map.get("c")));

        linker_type = Short.valueOf(parameter_map.get("cl_type"));

        // check the cross-linker and the fix modification
        if (linker_type == 1) {
            if (Math.abs(fix_mod_map.get('K') - fix_mod_map.get('n')) > 1e-6) {
                linker_mass = 0;
                throw new Exception("The link sites have different fix modifications.");
            } else {
                linker_mass = Double.valueOf(parameter_map.get("cl_mass")) - fix_mod_map.get('K');
            }
        } else if (linker_type == 2) {
            linker_mass = Double.valueOf(parameter_map.get("cl_mass"));
        } else {
            linker_mass = 0;
            throw new Exception("The cross-linker type cannot be recognized.");
        }

        // read protein database
        Map<String, String> pro_seq_map;
        DbTool db_tool_obj = new DbTool(db_path, parameter_map.get("database_type"));
        if (parameter_map.get("append_contaminants").contentEquals("1")) {
            DbTool contaminantsDb = new DbTool(null, "contaminants");
            pro_seq_map = contaminantsDb.getProteinSequenceMap();
            pro_seq_map.putAll(db_tool_obj.getProteinSequenceMap()); // using the target sequence to replace contaminant sequence if there is conflict.
            pro_annotate_map = contaminantsDb.getProteinAnnotateMap();
            pro_annotate_map.putAll(db_tool_obj.getProteinAnnotateMap()); // using the target sequence to replace contaminant sequence if there is conflict.
        } else {
            pro_seq_map = db_tool_obj.getProteinSequenceMap();
            pro_annotate_map = db_tool_obj.getProteinAnnotateMap();
        }

        // define a new MassTool object
        mass_tool_obj = new MassTool(missed_cleavage, fix_mod_map, "KR", "P", true, mz_bin_size * 0.5, one_minus_bin_offset, "N14", "[]");

        // generate seq_pro_map
        Map<String, boolean[]> seq_term_map = new HashMap<>();
        seqProMap = buildSeqProMap(pro_seq_map, seq_term_map, min_chain_length, max_chain_length);

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
            Set<Short> linkSiteSet = getLinkSiteSet(seq, proteinNTerm, proteinCTerm, linker_type);
            if (!linkSiteSet.isEmpty()) {
                double totalMass = (mass_tool_obj.calResidueMass(seq) + mass_tool_obj.H2O);
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

            // mod containing
            Set<VarSequence> varSeqSet = generateModSeq(seq, linkSiteSet, varModParamSet, binaryModParamSet, varModMaxNum);
            for (VarSequence varSeq : varSeqSet) {
                linkSiteSet = new HashSet<>(5, 1);
                linkSiteSet.add(varSeq.linkSite);
                if (!linkSiteSet.isEmpty()) {
                    double totalMass = (mass_tool_obj.calResidueMass(varSeq.seq) + mass_tool_obj.H2O);
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

    public Map<String, Set<String>> getSeqProMap() {
        return seqProMap;
    }

    public MassTool returnMassTool() {
        return mass_tool_obj;
    }

    public Map<String, String> getProAnnotateMap() {
        return pro_annotate_map;
    }

    public Map<Character, Double> getFixModMap() {
        return fix_mod_map;
    }

    public TreeMap<Integer, Set<String>> getMassBinSeqMap() {
        return bin_seq_map;
    }

    public Map<String, ChainEntry> getSeqEntryMap() {
        return seq_entry_map;
    }

    public double binToLeftMass(int bin_idx) {
        return bin_idx * ms1_bin_size;
    }

    public double binToRightMass(int bin_idx) {
        return (bin_idx + 1) * ms1_bin_size;
    }

    public int massToBin(double mass) {
        return (int) Math.floor(mass * inverseMs1BinSize);
    }

    private Map<String, Set<String>> buildSeqProMap(Map<String, String> pro_seq_map, Map<String, boolean[]> seq_term_map, int min_chain_length, int max_chain_length) {
        Map<String, Set<String>> seq_pro_map = new HashMap<>(pro_seq_map.size() * 150, 1);
        Set<String> for_check_duplicate = new HashSet<>();
        for (String pro_id : pro_seq_map.keySet()) {
            String pro_seq = pro_seq_map.get(pro_id);
            Set<String> seq_set = mass_tool_obj.buildChainSet(pro_seq, linker_type);
            for (String target_seq : seq_set) {
                if ((target_seq.length() >= min_chain_length) && (target_seq.length() <= max_chain_length) && !MassTool.containsNonAAAndNC(target_seq)) {
                    if (!for_check_duplicate.contains(target_seq.replace('L', 'I'))) {
                        // Add the sequence to the check set for duplicate check
                        for_check_duplicate.add(target_seq.replace('L', 'I')); // "L" and "I" have the same mass

                        boolean n_term = false;
                        boolean c_term = false;
                        if (pro_seq.startsWith(target_seq.substring(1, target_seq.length() - 1))) {
                            n_term = true;
                        }
                        if (pro_seq.startsWith("M") && pro_seq.substring(1).startsWith(target_seq.substring(1, target_seq.length() - 1))) { // consider the first "M" being cut situation.
                            n_term = true;
                        }
                        if (pro_seq.endsWith(target_seq.substring(1, target_seq.length() - 1))) {
                            c_term = true;
                        }
                        seq_term_map.put(target_seq, new boolean[]{n_term, c_term});

                        Set<String> pro_set = new HashSet<>(5, 1);
                        pro_set.add(pro_id);
                        seq_pro_map.put(target_seq, pro_set);
                    } else if (seq_pro_map.containsKey(target_seq)) {
                        // considering the case that the sequence has multiple proteins. In the above if clock, such a protein wasn't recorded.
                        seq_pro_map.get(target_seq).add(pro_id);
                    }
                }
            }
        }

        // generate decoy seq
        for (String pro_id : pro_seq_map.keySet()) {
            String pro_seq = pro_seq_map.get(pro_id);
            String decoy_pro_seq = (new StringBuilder(pro_seq)).reverse().toString();
            Set<String> decoy_seq_set = mass_tool_obj.buildChainSet(decoy_pro_seq, linker_type);
            if (pro_seq.startsWith("M")) {
                decoy_seq_set.addAll(mass_tool_obj.buildChainSet((new StringBuilder(pro_seq.substring(1))).reverse().toString(), linker_type));
            }
            for (String decoy_seq : decoy_seq_set) {
                if ((decoy_seq.length() >= min_chain_length) && (decoy_seq.length() <= max_chain_length) && !MassTool.containsNonAAAndNC(decoy_seq)) {
                    if (!for_check_duplicate.contains(decoy_seq.replace('L', 'I'))) {
                        for_check_duplicate.add(decoy_seq.replace('L', 'I'));

                        boolean n_term = false;
                        boolean c_term = false;
                        if (decoy_pro_seq.startsWith(decoy_seq.substring(1, decoy_seq.length() - 1))) { // here, we don't consider the first "M" being cut situation.
                            n_term = true;
                        }
                        if (decoy_pro_seq.endsWith(decoy_seq.substring(1, decoy_seq.length() - 1))) {
                            c_term = true;
                        }
                        seq_term_map.put(decoy_seq, new boolean[]{n_term, c_term});

                        Set<String> pro_set = new HashSet<>(5, 1);
                        pro_set.add("DECOY_" + pro_id);
                        seq_pro_map.put(decoy_seq, pro_set);
                    }
                    if (seq_pro_map.containsKey(decoy_seq)) {
                        seq_pro_map.get(decoy_seq).add("DECOY_" + pro_id);
                    }
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
                Map<Integer, List<Double>> idxBinaryModMassMap = new HashMap<>(seq.length(), 1);
                for (int i = 0; i < seq.length(); ++i) {
                    if (i != linkSite) {
                        String aa = seq.substring(i, i + 1);
                        if (binaryModParam.aas.contains(aa)) {
                            List<Double> tempList = new LinkedList<>();
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
                            sb.append(String.format(Locale.US, "[%.3f]", idxBinaryModMassMap.get(i).get(0)));
                        }
                    }
                    varSeqSet.add(new VarSequence(sb.toString(), linkSite, binaryModParam.hashCode()));

                    if (idxBinaryModMassMap.size() < varModMaxNum) {
                        // generate sequences containing the binary mod and additional var mod
                        // get all locations having var mods
                        Map<Integer, List<Double>> idxVarModMassMap = new HashMap<>(seq.length(), 1);
                        for (int i = 0; i < seq.length(); ++i) {
                            if (i != linkSite) {
                                if (!idxBinaryModMassMap.containsKey(i)) {
                                    char aa = seq.charAt(i);
                                    for (VarModParam varModParam : varModParamSet) {
                                        if (varModParam.aa == aa) {
                                            if (idxVarModMassMap.containsKey(i)) {
                                                idxVarModMassMap.get(i).add(varModParam.modMass);
                                            } else {
                                                List<Double> temp = new LinkedList<>();
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
                            Map<Integer, List<Double>> idxBinaryVarModMassMap = new HashMap<>(seq.length(), 1);
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
            Map<Integer, List<Double>> idxVarModMassMap = new HashMap<>(seq.length(), 1);
            for (int i = 0; i < seq.length(); ++i) {
                if (i != linkSite) {
                    char aa = seq.charAt(i);
                    for (VarModParam varModParam : varModParamSet) {
                        if (varModParam.aa == aa) {
                            if (idxVarModMassMap.containsKey(i)) {
                                idxVarModMassMap.get(i).add(varModParam.modMass);
                            } else {
                                List<Double> temp = new LinkedList<>();
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

    private Set<VarSequence> generateModSeqSub(String seq, int[] idxCombination, Map<Integer, List<Double>> idxModMassMap, short linkSite, int binaryModType) {
        List<Map<Integer, Double>> localIdxModMassMaps = generateLocalIdxModMassMap(idxCombination, idxModMassMap);

        Set<VarSequence> outputSet = new HashSet<>();
        for (Map<Integer, Double> localIdxModMassMap : localIdxModMassMaps) {
            StringBuilder sb = new StringBuilder(seq.length() * 10);
            for (int i = 0; i < seq.length(); ++i) {
                sb.append(seq.charAt(i));
                if (localIdxModMassMap.containsKey(i)) {
                    sb.append(String.format(Locale.US, "[%.3f]", localIdxModMassMap.get(i)));
                }
            }
            outputSet.add(new VarSequence(sb.toString(), linkSite, binaryModType));
        }

        return outputSet;
    }

    private List<Map<Integer, Double>> generateLocalIdxModMassMap(int[] idxArray, Map<Integer, List<Double>> idxModMassMap) {
        List<Map<Integer, Double>> outputList = new LinkedList<>();
        if (idxArray.length == globalVarModMaxNum) {
            for (int i0 = 0; i0 < idxModMassMap.get(idxArray[0]).size(); ++i0) {
                for (int i1 = 0; i1 < idxModMassMap.get(idxArray[1]).size(); ++i1) {
                    for (int i2 = 0; i2 < idxModMassMap.get(idxArray[2]).size(); ++i2) {
                        for (int i3 = 0; i3 < idxModMassMap.get(idxArray[3]).size(); ++i3) {
                            for (int i4 = 0; i4 < idxModMassMap.get(idxArray[4]).size(); ++i4) {
                                Map<Integer, Double> localIdxModMassMap = new HashMap<>(6, 1);
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
                            Map<Integer, Double> localIdxModMassMap = new HashMap<>(5, 1);
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
                        Map<Integer, Double> localIdxModMassMap = new HashMap<>(4, 1);
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
                    Map<Integer, Double> localIdxModMassMap = new HashMap<>(3, 1);
                    localIdxModMassMap.put(idxArray[0], idxModMassMap.get(idxArray[0]).get(i0));
                    localIdxModMassMap.put(idxArray[1], idxModMassMap.get(idxArray[1]).get(i1));
                    outputList.add(localIdxModMassMap);
                }
            }
        } else if (idxArray.length == 1) {
            for (int i0 = 0; i0 < idxModMassMap.get(idxArray[0]).size(); ++i0) {
                Map<Integer, Double> localIdxModMassMap = new HashMap<>(2, 1);
                localIdxModMassMap.put(idxArray[0], idxModMassMap.get(idxArray[0]).get(i0));
                outputList.add(localIdxModMassMap);
            }
        }

        return outputList;
    }

    private Set<String> checkKCTermMod(Set<String> varSeqSet) { // eliminate those sequence that the middle amino acids having the same mod mass and the n-term and the first amino acid or the c-term and the last amino acid have the same mod mass. todo: check
        String[] varSeqArray = varSeqSet.toArray(new String[varSeqSet.size()]);
        Arrays.sort(varSeqArray); // Make sure that nK[].... is before n[]K..., so that n[]K... will be kept.
        int seqLength = MassTool.seqToAAList(varSeqArray[0], "[]").length;
        AA[][] aaArrays = new AA[varSeqArray.length][seqLength];
        for (int i = 0; i < varSeqArray.length; ++i) {
            aaArrays[i] = MassTool.seqToAAList(varSeqArray[i], "[]");
        }

        if (aaArrays.length > 1) {
            Set<String> keptSeqSet = new HashSet<>();
            for (int i = 0; i < aaArrays.length - 1; ++i) {
                boolean keep = true;
                for (int j = i + 1; j < aaArrays.length; ++j) {
                    if ((Math.abs(aaArrays[i][0].ptmDeltaMass - aaArrays[j][1].ptmDeltaMass) < varModMassResolution) && (Math.abs(aaArrays[i][1].ptmDeltaMass - aaArrays[j][0].ptmDeltaMass) < varModMassResolution) && (Math.abs(aaArrays[i][seqLength - 2].ptmDeltaMass - aaArrays[j][seqLength - 1].ptmDeltaMass) < varModMassResolution) && (Math.abs(aaArrays[i][seqLength - 1].ptmDeltaMass - aaArrays[j][seqLength - 2].ptmDeltaMass) < varModMassResolution)) {
                        keep = false;
                        for (int k = 2; k < seqLength - 2; ++k) {
                            if (Math.abs(aaArrays[i][k].ptmDeltaMass - aaArrays[j][k].ptmDeltaMass) > varModMassResolution) {
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

    private Set<Short> getLinkSiteSet(String seq, boolean n_term, boolean c_term, short linker_type) {
        AA[] aa_list = MassTool.seqToAAList(seq, "[]");
        Set<Short> output = new HashSet<>(5, 1);
        for (int i = 1; i < aa_list.length - 2; ++i) {
            if (linker_type == 1 && aa_list[i].aa == 'K' && (Math.abs(aa_list[i].ptmDeltaMass) < varModMassResolution)) {
                output.add((short) i);
            } else if (linker_type == 2 && aa_list[i].aa == 'C' && (Math.abs(aa_list[i].ptmDeltaMass) < varModMassResolution)) {
                output.add((short) i);
            }
        }
        if (linker_type == 1 && n_term && !output.contains((short) 1) && (Math.abs(aa_list[0].ptmDeltaMass) < varModMassResolution)) {
            output.add((short) 0);
        }
        if (linker_type == 1 && c_term && aa_list[aa_list.length - 2].aa == 'K' && (Math.abs(aa_list[aa_list.length - 2].ptmDeltaMass) < varModMassResolution)) {
            output.add((short) (aa_list.length - 2));
        }
        return output;
    }

    private void getVarModParams(String v, Set<VarModParam> varModParamSet, Set<BinaryModParam> binaryModParamSet) throws Exception {
        Matcher varModMatcher = varModParamPattern.matcher(v);
        if (varModMatcher.matches()) {
            double modMass = Double.valueOf(varModMatcher.group(1));
            String aas = varModMatcher.group(2);
            boolean isBinary = varModMatcher.group(3).contentEquals("1");
            if (Math.abs(modMass) > varModMassResolution) {
                if (isBinary) {
                    binaryModParamSet.add(new BinaryModParam(modMass, aas));
                } else {
                    for (int i = 0; i < aas.length(); ++i) {
                        if (Math.abs(fix_mod_map.get(aas.charAt(i))) <= varModMassResolution) {
                            varModParamSet.add(new VarModParam(modMass, aas.charAt(i)));
                        }
                    }
                }
            }
        } else {
            throw new Exception(String.format(Locale.US, "Cannot parse variable modification parameter from %s.", v));
        }
    }
}
