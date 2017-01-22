package proteomics.Search;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.ECL2;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.*;
import proteomics.Index.BuildIndex;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class Search {

    private static final Logger logger = LoggerFactory.getLogger(Search.class);

    public static final double single_chain_t = 0;

    private final int max_common_ion_charge;
    public final float ms1_tolerance;
    public final int ms1_tolerance_unit;
    private final Map<String, ChainEntry> chain_entry_map;
    private final Map<String, Float> fix_mod_map;
    private final MassTool mass_tool_obj;
    private final Map<String, Float> mass_table;
    private final TreeMap<Integer, Set<String>> bin_seq_map;
    private final BuildIndex build_index_obj;
    private int[] C13_correction_range;
    public final boolean consider_two_identical_chains;
    private Map<Integer, Long> bin_candidate_num_map;

    /////////////////////////////////////////public methods////////////////////////////////////////////////////////////
    public Search(BuildIndex build_index_obj, Map<String, String> parameter_map) throws Exception {
        this.build_index_obj = build_index_obj;
        chain_entry_map = build_index_obj.getSeqEntryMap();
        fix_mod_map = build_index_obj.getFixModMap();
        mass_tool_obj = build_index_obj.returnMassTool();
        max_common_ion_charge = Integer.valueOf(parameter_map.get("max_common_ion_charge"));
        ms1_tolerance_unit = Integer.valueOf(parameter_map.get("ms1_tolerance_unit"));
        ms1_tolerance = Float.valueOf(parameter_map.get("ms1_tolerance"));
        mass_table = mass_tool_obj.getMassTable();
        bin_seq_map = build_index_obj.getMassBinSeqMap();
        consider_two_identical_chains = parameter_map.get("consider_two_identical_chains").contentEquals("1");
        bin_candidate_num_map = build_index_obj.getBinCandidateNumMap();

        String[] temp = parameter_map.get("C13_correction_range").split(",");
        C13_correction_range = new int[temp.length];
        for (int i = 0; i < temp.length; ++i) {
            C13_correction_range[i] = Integer.valueOf(temp[i].trim());
        }
        Arrays.sort(C13_correction_range);
    }

    HashMap<Integer, ResultEntry> doSearch(NavigableMap<Integer, Set<SpectrumEntry>> bin_spectra_map) throws Exception {
        // get the necessary max ms1 mass
        int max_spectrum_mass_bin = bin_spectra_map.lastKey();
        float max_spectrum_mass = build_index_obj.binToRightMass(max_spectrum_mass_bin);
        int max_chain_bin_idx = Math.min(build_index_obj.massToBin(max_spectrum_mass + C13_correction_range[C13_correction_range.length - 1] * mass_table.get("C13_DIFF") - build_index_obj.linker_mass) + 1 - bin_seq_map.firstKey(), bin_seq_map.lastKey());
        int min_chain_bin_idx = bin_seq_map.firstKey();
        float min_additional_mass = build_index_obj.binToLeftMass(min_chain_bin_idx) + build_index_obj.linker_mass; // this is the minimum additional mass added to a peptide chain.

        if (max_chain_bin_idx < min_chain_bin_idx) {
            return new HashMap<>();
        }

        Map<Integer, List<DebugEntry>> debug_scan_entry_map = new HashMap<>();
        Map<Integer, Map<String, Double>> dev_scan_chain_score_map = new HashMap<>();

        Map<SpectrumEntry, TreeMap<Integer, ChainResultEntry>> scan_bin_chain_map = new HashMap<>();
        for (int bin_idx : bin_seq_map.keySet()) {
            if (bin_idx > max_chain_bin_idx) {
                break;
            }

            float bin_left_mass = build_index_obj.binToLeftMass(bin_idx - 1);
            float left_mass = bin_left_mass + min_additional_mass + C13_correction_range[0] * mass_table.get("C13_DIFF");
            if (ms1_tolerance_unit == 0) {
                left_mass -= ms1_tolerance;
            } else {
                left_mass -= (left_mass + build_index_obj.binToRightMass(max_chain_bin_idx) + build_index_obj.linker_mass) * ms1_tolerance * 1e-6f;
            }

            NavigableMap<Integer, Set<SpectrumEntry>> sub_spectra_map = bin_spectra_map.subMap(Math.max(build_index_obj.massToBin(left_mass), bin_spectra_map.firstKey()), true, max_spectrum_mass_bin, true);

            if (sub_spectra_map.isEmpty()) {
                continue;
            }

            Set<String> seq_set = bin_seq_map.get(bin_idx);
            for (String seq : seq_set) {
                ChainEntry chain_entry = chain_entry_map.get(seq);
                linearScan(sub_spectra_map, chain_entry, bin_idx, scan_bin_chain_map, debug_scan_entry_map, dev_scan_chain_score_map);
            }
        }

        // write debug file
        if (ECL2.debug) {
            for (int scan_num : debug_scan_entry_map.keySet()) {
                try (BufferedWriter writer = new BufferedWriter(new FileWriter(scan_num + ".csv"))) {
                    List<DebugEntry> temp = debug_scan_entry_map.get(scan_num);
                    Collections.sort(temp, Comparator.reverseOrder());
                    writer.write("chain,link_site,mass,score\n");
                    for (DebugEntry t : temp) {
                        writer.write(String.format("%s,%d,%f,%f\n", t.chain, t.link_site, t.mass, t.score));
                    }
                } catch (IOException ex) {
                    logger.error(ex.getMessage());
                    ex.printStackTrace();
                    System.exit(1);
                }
            }
        }

        HashMap<Integer, ResultEntry> scan_result_map = new HashMap<>();
        for (SpectrumEntry spectrum_entry : scan_bin_chain_map.keySet()) {
            float max_mass = spectrum_entry.mass_without_linker_mass / 2;
            float left_tol = ms1_tolerance;
            float right_tol = ms1_tolerance;
            if (ms1_tolerance_unit == 1) {
                left_tol = spectrum_entry.precursor_mass - (spectrum_entry.precursor_mass / (1 + ms1_tolerance * 1e-6f));
                right_tol = (spectrum_entry.precursor_mass / (1 - ms1_tolerance * 1e-6f)) - spectrum_entry.precursor_mass;
            }
            max_mass += right_tol;

            int max_v = build_index_obj.massToBin(max_mass) + 1;

            long candidate_num = 0;

            TreeMap<Integer, ChainResultEntry> bin_chain_map = scan_bin_chain_map.get(spectrum_entry);
            for (int idx_1 : bin_chain_map.keySet()) {
                if (idx_1 > max_v) {
                    break;
                }

                long bin1_candidate_num = bin_candidate_num_map.get(idx_1);

                float left_mass_2;
                float right_mass_2;
                int left_idx_2;
                int right_idx_2;
                NavigableMap<Integer, ChainResultEntry> sub_map = new TreeMap<>();

                // consider C13 correction from -2 to +1
                for (int i : C13_correction_range) {
                    left_mass_2 = spectrum_entry.mass_without_linker_mass + i * mass_table.get("C13_DIFF") - build_index_obj.binToRightMass(idx_1) - left_tol;
                    right_mass_2 = spectrum_entry.mass_without_linker_mass + i * mass_table.get("C13_DIFF") - build_index_obj.binToLeftMass(idx_1) + right_tol;
                    left_idx_2 = build_index_obj.massToBin(left_mass_2);
                    right_idx_2 = build_index_obj.massToBin(right_mass_2) + 1;
                    sub_map.putAll(bin_chain_map.subMap(left_idx_2, true, right_idx_2, true));
                }

                if (!sub_map.isEmpty()) {
                    ChainResultEntry chain_score_entry_1 = bin_chain_map.get(idx_1);
                    for (int idx_2 : sub_map.keySet()) {
                        if (idx_1 == idx_2) {
                            candidate_num += bin1_candidate_num * (bin1_candidate_num + 1) / 2;
                        } else {
                            candidate_num += bin1_candidate_num * bin_candidate_num_map.get(idx_2);
                        }

                        ChainResultEntry chain_score_entry_2 = bin_chain_map.get(idx_2);
                        double score = 0;
                        if (chain_score_entry_1.getNormalizedSeq().contentEquals(chain_score_entry_2.getNormalizedSeq())) {
                            if (consider_two_identical_chains) {
                                score = chain_score_entry_1.getScore() + chain_score_entry_2.getScore();
                            }
                        } else {
                            score = chain_score_entry_1.getScore() + chain_score_entry_2.getScore();
                        }

                        double second_score = 0;
                        double temp_1 = -1;
                        if (chain_score_entry_1.getSecondSeq() != null) {
                            if (chain_score_entry_1.getNormalizedSecondSeq().contentEquals(chain_score_entry_2.getNormalizedSeq())) {
                                if (consider_two_identical_chains) {
                                    temp_1 = chain_score_entry_1.getSecondScore() + chain_score_entry_2.getScore();
                                }
                            } else {
                                temp_1 = chain_score_entry_1.getSecondScore() + chain_score_entry_2.getScore();
                            }
                        }
                        double temp_2 = -1;
                        if (chain_score_entry_2.getSecondSeq() != null) {
                            if (chain_score_entry_1.getNormalizedSeq().contentEquals(chain_score_entry_2.getNormalizedSecondSeq())) {
                                if (consider_two_identical_chains) {
                                    temp_2 = chain_score_entry_1.getScore() + chain_score_entry_2.getSecondScore();
                                }
                            } else {
                                temp_2 = chain_score_entry_1.getScore() + chain_score_entry_2.getSecondScore();
                            }
                        }

                        if (temp_1 > 0) {
                            if (temp_1 >= temp_2) {
                                second_score = temp_1;
                            }
                        } else if (temp_2 > 0) {
                            if (temp_2 > temp_1) {
                                second_score = temp_2;
                            }
                        }

                        if (scan_result_map.containsKey(spectrum_entry.scan_num)) {
                            ResultEntry result_entry = scan_result_map.get(spectrum_entry.scan_num);
                            if (ECL2.cal_evalue && (result_entry.getScoreCount() < ECL2.score_point_t)) {
                                for (double s1 : chain_score_entry_1.getScoreList()) {
                                    for (double s2 : chain_score_entry_2.getScoreList()) {
                                        result_entry.addToScoreHistogram(s1 + s2);
                                    }
                                }
                            }
                            if (score > result_entry.getScore()) {
                                result_entry.setSecondScore(Math.max(result_entry.getScore(), second_score));
                                result_entry.setScore(score);
                                result_entry.setChain1(chain_score_entry_1.getSeq());
                                result_entry.setChain2(chain_score_entry_2.getSeq());
                                result_entry.setLinkSite1(chain_score_entry_1.getLinkSite());
                                result_entry.setLinkSite2(chain_score_entry_2.getLinkSite());
                            } else if (Math.max(score, second_score) > result_entry.getSecondScore()) {
                                result_entry.setSecondScore(Math.max(score, second_score));
                            }
                        } else {
                            ResultEntry result_entry = new ResultEntry(spectrum_entry.spectrum_id, spectrum_entry.precursor_mz, spectrum_entry.precursor_mass, spectrum_entry.rt, spectrum_entry.precursor_charge);
                            if (ECL2.cal_evalue && (result_entry.getScoreCount() < ECL2.score_point_t)) {
                                for (double s1 : chain_score_entry_1.getScoreList()) {
                                    for (double s2 : chain_score_entry_2.getScoreList()) {
                                        result_entry.addToScoreHistogram(s1 + s2);
                                    }
                                }
                            }
                            result_entry.setScore(score);
                            result_entry.setSecondScore(second_score);
                            result_entry.setChain1(chain_score_entry_1.getSeq());
                            result_entry.setChain2(chain_score_entry_2.getSeq());
                            result_entry.setLinkSite1(chain_score_entry_1.getLinkSite());
                            result_entry.setLinkSite2(chain_score_entry_2.getLinkSite());
                            scan_result_map.put(spectrum_entry.scan_num, result_entry);
                        }
                    }
                }
            }

            // set candidate number
            if (scan_result_map.containsKey(spectrum_entry.scan_num)) {
                scan_result_map.get(spectrum_entry.scan_num).setCandidateNum(candidate_num);
            }

            if (ECL2.dev) {
                if (scan_result_map.containsKey(spectrum_entry.scan_num)) {
                    ResultEntry result_entry = scan_result_map.get(spectrum_entry.scan_num);
                    Map<String, Double> chain_score_map = dev_scan_chain_score_map.get(spectrum_entry.scan_num);
                    Double chain_score_1 = chain_score_map.get(result_entry.getChain1() + "-" + result_entry.getLinkSite1());
                    Double chain_score_2 = chain_score_map.get(result_entry.getChain2() + "-" + result_entry.getLinkSite2());
                    List<Double> scores = new LinkedList<>(chain_score_map.values());
                    Collections.sort(scores, Comparator.reverseOrder());
                    int chain_rank_1 = scores.indexOf(chain_score_1) + 1;
                    int chain_rank_2 = scores.indexOf(chain_score_2) + 1;
                    result_entry.setChainDetails(chain_score_1, chain_rank_1, chain_score_2, chain_rank_2);
                }
            }
        }
        return scan_result_map;
    }

    public FinalResultEntry convertResultEntry(int scanNum, ResultEntry result_entry) {
        int rank = 1;
        String chain_seq_1 = result_entry.getChain1();
        String chain_seq_2 = result_entry.getChain2();
        ChainEntry chain_entry_1 = chain_entry_map.get(chain_seq_1);
        ChainEntry chain_entry_2 = chain_entry_map.get(chain_seq_2);

        float spectrum_mass = result_entry.spectrum_mass;
        float theo_mass = chain_entry_1.chain_mass + chain_entry_2.chain_mass + build_index_obj.linker_mass;

        int C13_Diff_num = getC13Num(spectrum_mass, theo_mass);
        spectrum_mass += C13_Diff_num * mass_table.get("C13_DIFF");
        float ppm = (spectrum_mass - theo_mass) * 1e6f / theo_mass;

        String pro_1 = "";
        for (String temp : chain_entry_1.pro_id) {
            pro_1 += temp + ";";
        }
        String pro_2 = "";
        for (String temp : chain_entry_2.pro_id) {
            pro_2 += temp + ";";
        }

        int hit_type; // 0 = T-T; 1 = D-D; 2 = T-D;
        if (pro_1.startsWith("DECOY") && pro_2.startsWith("DECOY")) {
            hit_type = 1;
        } else if (!pro_1.startsWith("DECOY") && !pro_2.startsWith("DECOY")) {
            hit_type = 0;
        } else {
            hit_type = 2;
        }

        String cl_type = "intra_protein";
        boolean keep = false;
        for (String temp_1 : chain_entry_1.pro_id) {
            if (temp_1.startsWith("DECOY_")) {
                temp_1 = temp_1.substring(6);
            }
            for (String temp_2 : chain_entry_2.pro_id) {
                if (temp_2.startsWith("DECOY_")) {
                    temp_2 = temp_2.substring(6);
                }
                if (temp_1.contentEquals(temp_2)) {
                    keep = true;
                    break;
                }
            }
        }

        if (!keep) {
            cl_type = "inter_protein";
        }

        double delta_c = 1 - (result_entry.getSecondScore() / result_entry.getScore());

        // add fix modification to the sequences.
        String final_seq_1 = addFixMod(chain_seq_1, result_entry.getLinkSite1());
        String final_seq_2 = addFixMod(chain_seq_2, result_entry.getLinkSite2());

        return new FinalResultEntry(scanNum, result_entry.spectrum_id, rank, result_entry.charge, result_entry.spectrum_mz, result_entry.spectrum_mass, theo_mass, result_entry.rt, ppm, result_entry.getScore(), delta_c, final_seq_1, result_entry.getLinkSite1(), pro_1, final_seq_2, result_entry.getLinkSite2(), pro_2, cl_type, hit_type, C13_Diff_num, result_entry.getEValue(), result_entry.getScoreCount(), result_entry.getRSquare(), result_entry.getSlope(), result_entry.getIntercept(), result_entry.getStartIdx(), result_entry.getEndIdx(), result_entry.getChainScore1(), result_entry.getChainRank1(), result_entry.getChainScore2(), result_entry.getChainRank2(), result_entry.getCandidateNum());
    }

    private void linearScan(NavigableMap<Integer, Set<SpectrumEntry>> sub_spectra_map, ChainEntry chain_entry, int chain_bin_idx, Map<SpectrumEntry, TreeMap<Integer, ChainResultEntry>> scan_bin_chain_map, Map<Integer, List<DebugEntry>> debug_scan_entry_map, Map<Integer, Map<String, Double>> dev_scan_chain_score_map) {
        for (int spectra_bin_idx : sub_spectra_map.keySet()) {
            Map<Float, SparseBooleanVector> massVectorMap = new HashMap<>();
            for (int link_site_1 : chain_entry.link_site_set) {
                for (SpectrumEntry spectrum_entry : sub_spectra_map.get(spectra_bin_idx)) {
                    SparseBooleanVector theo_mz;
                    if (massVectorMap.containsKey(spectrum_entry.precursor_mass)) {
                        theo_mz = massVectorMap.get(spectrum_entry.precursor_mass);
                    } else {
                        int precursor_charge = spectrum_entry.precursor_charge;
                        theo_mz = mass_tool_obj.buildVector(mass_tool_obj.buildPseudoCLIonArray(chain_entry.chain_ion_array, link_site_1, max_common_ion_charge, spectrum_entry.precursor_mass - chain_entry.chain_mass, precursor_charge), precursor_charge);
                    }

                    // Calculate dot produce
                    double dot_product = theo_mz.dot(spectrum_entry.pl_map_xcorr) * 0.25;

                    if (ECL2.debug) {
                        if (debug_scan_entry_map.containsKey(spectrum_entry.scan_num)) {
                            debug_scan_entry_map.get(spectrum_entry.scan_num).add(new DebugEntry(chain_entry.seq, link_site_1, chain_entry.chain_mass, dot_product));
                        } else {
                            List<DebugEntry> temp = new LinkedList<>();
                            temp.add(new DebugEntry(chain_entry.seq, link_site_1, chain_entry.chain_mass, dot_product));
                            debug_scan_entry_map.put(spectrum_entry.scan_num, temp);
                        }
                    }

                    if (ECL2.dev && (dot_product > single_chain_t)) {
                        if (dev_scan_chain_score_map.containsKey(spectrum_entry.scan_num)) {
                            dev_scan_chain_score_map.get(spectrum_entry.scan_num).put(chain_entry.seq + "-" + link_site_1, dot_product);
                        } else {
                            Map<String, Double> temp_map = new HashMap<>();
                            temp_map.put(chain_entry.seq + "-" + link_site_1, dot_product);
                            dev_scan_chain_score_map.put(spectrum_entry.scan_num, temp_map);
                        }
                    }

                    // Record result
                    if (dot_product > single_chain_t) {
                        if (scan_bin_chain_map.containsKey(spectrum_entry)) {
                            Map<Integer, ChainResultEntry> temp_map = scan_bin_chain_map.get(spectrum_entry);
                            if (temp_map.containsKey(chain_bin_idx)) {
                                ChainResultEntry chain_result_entry = temp_map.get(chain_bin_idx);
                                chain_result_entry.addToScoreList(dot_product);
                                if (dot_product > chain_result_entry.getScore()) {
                                    chain_result_entry.setSecondScore(chain_result_entry.getScore());
                                    chain_result_entry.setSecondSeq(chain_result_entry.getSeq());
                                    chain_result_entry.setScore(dot_product);
                                    chain_result_entry.setSeq(chain_entry.seq);
                                    chain_result_entry.setLinkSite(link_site_1);
                                } else if (dot_product > chain_result_entry.getSecondScore()) {
                                    chain_result_entry.setSecondScore(dot_product);
                                    chain_result_entry.setSecondSeq(chain_entry.seq);
                                }
                            } else {
                                ChainResultEntry chain_result_entry = new ChainResultEntry();
                                chain_result_entry.addToScoreList(dot_product);
                                chain_result_entry.setSeq(chain_entry.seq);
                                chain_result_entry.setLinkSite(link_site_1);
                                chain_result_entry.setScore(dot_product);
                                chain_result_entry.setSecondScore(0);
                                temp_map.put(chain_bin_idx, chain_result_entry);
                            }
                        } else {
                            TreeMap<Integer, ChainResultEntry> temp_map = new TreeMap<>();
                            ChainResultEntry chain_result_entry = new ChainResultEntry();
                            chain_result_entry.addToScoreList(dot_product);
                            chain_result_entry.setSeq(chain_entry.seq);
                            chain_result_entry.setLinkSite(link_site_1);
                            chain_result_entry.setScore(dot_product);
                            chain_result_entry.setSecondScore(0);
                            temp_map.put(chain_bin_idx, chain_result_entry);
                            scan_bin_chain_map.put(spectrum_entry, temp_map);
                        }
                    }
                }
            }
        }
    }

    private int getC13Num(float exp_mass, float theo_mass) {
        float min_diff = 10;
        int num = 0;

        for (int i : C13_correction_range) {
            float temp = Math.abs(exp_mass + i * mass_table.get("C13_DIFF") - theo_mass);
            if (temp < min_diff) {
                min_diff = temp;
                num = i;
            }
        }

        return num;
    }

    private String addFixMod(String seq, int linkSite) {
        AA[] aaList = mass_tool_obj.seqToAAList(seq);
        StringBuilder sb = new StringBuilder(seq.length() * 3);
        for (int i = 0; i < aaList.length; ++i) {
            AA aa = aaList[i];
            if (i == linkSite) { // priority order: linkSite > fixMod > varMod
                sb.append(aa.aa);
            } else if (Math.abs(fix_mod_map.get(aa.aa)) > 1e-6) {
                sb.append(String.format("%s[%.2f]", aa.aa, fix_mod_map.get(aa.aa)));
            } else {
                if (Math.abs(aa.delta_mass) > 1e-6) {
                    sb.append(String.format("%s[%.2f]", aa.aa, aa.delta_mass));
                } else {
                    sb.append(aa.aa);
                }
            }
        }
        return sb.toString();
    }
}