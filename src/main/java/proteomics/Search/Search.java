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

    public final float ms1_tolerance;
    public final int ms1_tolerance_unit;
    private final Map<String, ChainEntry> chain_entry_map;
    private final Map<Character, Float> fix_mod_map;
    private final MassTool mass_tool_obj;
    private final TreeMap<Integer, Set<String>> bin_seq_map;
    private final BuildIndex build_index_obj;
    private final int[] C13_correction_range;
    private final Map<Integer, Long> bin_candidate_num_map;
    private final float single_chain_t;
    private final boolean cal_evalue;

    /////////////////////////////////////////public methods////////////////////////////////////////////////////////////
    public Search(BuildIndex build_index_obj, Map<String, String> parameter_map) {
        this.build_index_obj = build_index_obj;
        chain_entry_map = build_index_obj.getSeqEntryMap();
        fix_mod_map = build_index_obj.getFixModMap();
        mass_tool_obj = build_index_obj.returnMassTool();
        ms1_tolerance_unit = Integer.valueOf(parameter_map.get("ms1_tolerance_unit"));
        ms1_tolerance = Float.valueOf(parameter_map.get("ms1_tolerance"));
        bin_seq_map = build_index_obj.getMassBinSeqMap();
        bin_candidate_num_map = build_index_obj.getBinCandidateNumMap();

        if (parameter_map.containsKey("single_chain_t")) {
            single_chain_t = Float.valueOf(parameter_map.get("single_chain_t"));
        } else {
            single_chain_t = 0.1f;
        }

        if (parameter_map.containsKey("cal_evalue") && parameter_map.get("cal_evalue").trim().contentEquals("0")) {
            cal_evalue = false;
        } else {
            cal_evalue = true;
        }

        String[] temp = parameter_map.get("C13_correction_range").split(",");
        C13_correction_range = new int[temp.length];
        for (int i = 0; i < temp.length; ++i) {
            C13_correction_range[i] = Integer.valueOf(temp[i].trim());
        }
        Arrays.sort(C13_correction_range);
    }

    ResultEntry doSearch(SpectrumEntry spectrumEntry, SparseVector xcorrPL, int specMaxBinIdx) throws IOException {
        int max_chain_bin_idx = Math.min(build_index_obj.massToBin(spectrumEntry.mass_without_linker_mass + C13_correction_range[C13_correction_range.length - 1] * 1.00335483f) + 1 - bin_seq_map.firstKey(), bin_seq_map.lastKey());
        int min_chain_bin_idx = bin_seq_map.firstKey();

        if (max_chain_bin_idx < min_chain_bin_idx) {
            return null;
        }

        // set MS1 tolerance for further usage.
        float leftMs1Tol = ms1_tolerance;
        float rightMs1Tol = ms1_tolerance;
        if (ms1_tolerance_unit == 1) {
            leftMs1Tol = spectrumEntry.precursor_mass - (spectrumEntry.precursor_mass / (1 + ms1_tolerance * 1e-6f));
            rightMs1Tol = (spectrumEntry.precursor_mass / (1 - ms1_tolerance * 1e-6f)) - spectrumEntry.precursor_mass;
        }

        // for debug
        List<DebugEntry> debugEntryList = new LinkedList<>();
        Map<String, Double> devChainScoreMap = new HashMap<>();

        // start search
        TreeMap<Integer, ChainResultEntry> binChainMap = new TreeMap<>();
        Set<Integer> checkedBinSet = new HashSet<>(bin_seq_map.size() + 1, 1);
        long candidate_num = 0;
        ResultEntry resultEntry = new ResultEntry(spectrumEntry.spectrum_id, spectrumEntry.precursor_mz, spectrumEntry.precursor_mass, spectrumEntry.rt, spectrumEntry.precursor_charge, cal_evalue, binChainMap);
        for (int idx_1 : bin_seq_map.keySet()) {
            long bin1_candidate_num = bin_candidate_num_map.get(idx_1);

            float left_mass_2;
            float right_mass_2;
            int left_idx_2;
            int right_idx_2;
            NavigableMap<Integer, Set<String>> sub_map = new TreeMap<>();

            // consider C13 correction
            for (int i : C13_correction_range) {
                left_mass_2 = spectrumEntry.mass_without_linker_mass + i * 1.00335483f - build_index_obj.binToRightMass(idx_1) - leftMs1Tol;
                right_mass_2 = spectrumEntry.mass_without_linker_mass + i * 1.00335483f - build_index_obj.binToLeftMass(idx_1) + rightMs1Tol;
                left_idx_2 = build_index_obj.massToBin(left_mass_2);
                right_idx_2 = build_index_obj.massToBin(right_mass_2) + 1;
                sub_map.putAll(bin_seq_map.subMap(left_idx_2, true, right_idx_2, true));
            }

            if (!sub_map.isEmpty()) {
                if (checkedBinSet.contains(idx_1)) {
                    // The first bin has reach the middle point of the whole range. All pairs have been checked.
                    break;
                }

                // this bin hasn't been visited. Linear scan first.
                for (String seq : bin_seq_map.get(idx_1)) {
                    ChainEntry chainEntry = chain_entry_map.get(seq);
                    linearScan(spectrumEntry, xcorrPL, specMaxBinIdx, chainEntry, idx_1, binChainMap, debugEntryList, devChainScoreMap);
                }
                checkedBinSet.add(idx_1);

                if (binChainMap.containsKey(idx_1)) { // There may be no chain with chain score > single_chain_t
                    ChainResultEntry chain_score_entry_1 = binChainMap.get(idx_1);
                    for (int idx_2 : sub_map.keySet()) {
                        if (!checkedBinSet.contains(idx_2)) {
                            // this bin hasn't been visited. Linear scan first.
                            for (String seq : bin_seq_map.get(idx_2)) {
                                ChainEntry chainEntry = chain_entry_map.get(seq);
                                linearScan(spectrumEntry, xcorrPL, specMaxBinIdx, chainEntry, idx_2, binChainMap, debugEntryList, devChainScoreMap);
                            }
                            checkedBinSet.add(idx_2);
                        }

                        if (binChainMap.containsKey(idx_2)) { // There may be no chain with chain score > single_chain_t
                            ChainResultEntry chain_score_entry_2 = binChainMap.get(idx_2);

                            if (idx_1 == idx_2) {
                                candidate_num += bin1_candidate_num * (bin1_candidate_num + 1) / 2;
                            } else {
                                candidate_num += bin1_candidate_num * bin_candidate_num_map.get(idx_2);
                            }

                            // only two sequences with the same binary mod type can be linked.
                            if (chain_entry_map.get(chain_score_entry_1.getSeq()).binaryModType == chain_entry_map.get(chain_score_entry_2.getSeq()).binaryModType) {
                                double score;
                                if (chain_score_entry_1.getPtmFreeSeq().contentEquals(chain_score_entry_2.getPtmFreeSeq())) {
                                    score = (chain_score_entry_1.getScore() + chain_score_entry_2.getScore()) / 2;
                                } else {
                                    score = chain_score_entry_1.getScore() + chain_score_entry_2.getScore();
                                }

                                // calculate second score
                                double second_score = 0;
                                double temp_1 = -1;
                                if (chain_score_entry_1.getSecondSeq() != null) {
                                    if (chain_score_entry_1.getSecondPtmFreeSeq().contentEquals(chain_score_entry_2.getPtmFreeSeq())) {
                                        temp_1 = (chain_score_entry_1.getSecondScore() + chain_score_entry_2.getScore()) / 2;
                                    } else {
                                        temp_1 = chain_score_entry_1.getSecondScore() + chain_score_entry_2.getScore();
                                    }
                                }
                                double temp_2 = -1;
                                if (chain_score_entry_2.getSecondSeq() != null) {
                                    if (chain_score_entry_1.getPtmFreeSeq().contentEquals(chain_score_entry_2.getSecondPtmFreeSeq())) {
                                        temp_2 = (chain_score_entry_1.getScore() + chain_score_entry_2.getSecondScore()) / 2;
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

                                if (cal_evalue && (resultEntry.getScoreCount() < ECL2.score_point_t)) {
                                    for (double s1 : chain_score_entry_1.getScoreList()) {
                                        for (double s2 : chain_score_entry_2.getScoreList()) {
                                            resultEntry.addToScoreHistogram(s1 + s2);
                                        }
                                    }
                                }
                                if (score > resultEntry.getScore()) {
                                    resultEntry.setSecondScore(Math.max(resultEntry.getScore(), second_score));
                                    resultEntry.setScore(score);
                                    resultEntry.setChain1(chain_score_entry_1.getSeq());
                                    resultEntry.setChain2(chain_score_entry_2.getSeq());
                                    resultEntry.setLinkSite1(chain_score_entry_1.getLinkSite());
                                    resultEntry.setLinkSite2(chain_score_entry_2.getLinkSite());
                                } else if (Math.max(score, second_score) > resultEntry.getSecondScore()) {
                                    resultEntry.setSecondScore(Math.max(score, second_score));
                                }
                            }
                        }
                    }
                }
            }
        }

        // set candidate number
        resultEntry.setCandidateNum(candidate_num);

        // write debug file
        if (ECL2.debug) {
            BufferedWriter writer = new BufferedWriter(new FileWriter(spectrumEntry.spectrum_id + ".csv"));
            debugEntryList.sort(Comparator.reverseOrder());
            writer.write("chain,link_site,mass,score\n");
            for (DebugEntry t : debugEntryList) {
                writer.write(String.format(Locale.US, "%s,%d,%f,%f\n", addFixMod(t.chain, t.link_site), t.link_site, t.mass, t.score));
            }
            writer.close();
        }

        if (ECL2.dev && (candidate_num > 0)) {
            // get ranks of each chain
            Double chain_score_1 = devChainScoreMap.get(resultEntry.getChain1() + "-" + resultEntry.getLinkSite1());
            Double chain_score_2 = devChainScoreMap.get(resultEntry.getChain2() + "-" + resultEntry.getLinkSite2());
            List<Double> scores = new LinkedList<>(devChainScoreMap.values());
            scores.sort(Comparator.reverseOrder());
            int chain_rank_1 = scores.indexOf(chain_score_1) + 1;
            int chain_rank_2 = scores.indexOf(chain_score_2) + 1;
            resultEntry.setChainDetails(chain_score_1, chain_rank_1, chain_score_2, chain_rank_2);
        }

        if ((resultEntry.getChain1() != null) && (resultEntry.getChain2() != null)) {
            return resultEntry;
        } else {
            return null;
        }
    }

    public FinalResultEntry convertResultEntry(int scanNum, ResultEntry result_entry, Map<String, Set<String>> seqProMap, SpectrumEntry spectrumEntry) {
        int rank = 1;
        String chain_seq_1 = result_entry.getChain1();
        String chain_seq_2 = result_entry.getChain2();
        ChainEntry chain_entry_1 = chain_entry_map.get(chain_seq_1);
        ChainEntry chain_entry_2 = chain_entry_map.get(chain_seq_2);

        float spectrum_mass = result_entry.spectrum_mass;
        float theo_mass = chain_entry_1.chain_mass + chain_entry_2.chain_mass + build_index_obj.linker_mass;

        int C13_Diff_num = getC13Num(spectrum_mass, theo_mass);
        spectrum_mass += C13_Diff_num * 1.00335483f;
        float ppm = (spectrum_mass - theo_mass) * 1e6f / theo_mass;

        Set<String> pro1Set = new TreeSet<>();
        boolean isDecoy1 = false;
        for (String temp : seqProMap.get(chain_entry_1.seq.replaceAll("[^A-Znc]", ""))) {
            pro1Set.add(temp);
            if (temp.startsWith("DECOY")) { // there is no overlapped peptide between target and decoy.
                isDecoy1 = true;
            }
        }
        Set<String> pro2Set = new TreeSet<>();
        boolean isDecoy2 = false;
        for (String temp : seqProMap.get(chain_entry_2.seq.replaceAll("[^A-Znc]", ""))) {
            pro2Set.add(temp);
            if (temp.startsWith("DECOY")) { // there is no overlapped peptide between target and decoy.
                isDecoy2 = true;
            }
        }

        int hit_type; // 0 = T-T; 1 = D-D; 2 = T-D;
        if (isDecoy1 && isDecoy2) {
            hit_type = 1;
        } else if (!isDecoy1 && !isDecoy2) {
            hit_type = 0;
        } else {
            hit_type = 2;
        }

        String cl_type = "intra_protein";
        boolean keep = false;
        for (String temp_1 : seqProMap.get(chain_entry_1.seq.replaceAll("[^A-Znc]", ""))) {
            if (temp_1.startsWith("DECOY_")) {
                temp_1 = temp_1.substring(6);
            }
            for (String temp_2 : seqProMap.get(chain_entry_2.seq.replaceAll("[^A-Znc]", ""))) {
                if (temp_2.startsWith("DECOY_")) {
                    temp_2 = temp_2.substring(6);
                }
                if (temp_1.contentEquals(temp_2)) { // bias to intra protein cross-linking
                    keep = true;
                    break;
                }
            }
            if (keep) {
                break;
            }
        }

        if (!keep) {
            cl_type = "inter_protein";
        }

        double delta_c = 1 - (result_entry.getSecondScore() / result_entry.getScore());

        // add fix modification to the sequences.
        String final_seq_1 = addFixMod(chain_seq_1, result_entry.getLinkSite1());
        String final_seq_2 = addFixMod(chain_seq_2, result_entry.getLinkSite2());

        return new FinalResultEntry(scanNum, result_entry.spectrum_id, rank, result_entry.charge, result_entry.spectrum_mz, result_entry.spectrum_mass, theo_mass, result_entry.rt, ppm, result_entry.getScore(), delta_c, final_seq_1, result_entry.getLinkSite1(), String.join(";", pro1Set), final_seq_2, result_entry.getLinkSite2(), String.join(";", pro2Set), cl_type, hit_type, C13_Diff_num, result_entry.getEValue(), result_entry.getScoreCount(), result_entry.getRSquare(), result_entry.getSlope(), result_entry.getIntercept(), result_entry.getStartIdx(), result_entry.getEndIdx(), result_entry.getChainScore1(), result_entry.getChainRank1(), result_entry.getChainScore2(), result_entry.getChainRank2(), result_entry.getCandidateNum(), cal_evalue, spectrumEntry.mgfTitle);
    }

    private void linearScan(SpectrumEntry spectrumEntry, SparseVector xcorrPL, int specMaxBinIdx, ChainEntry chainEntry, int binInx, TreeMap<Integer, ChainResultEntry> binChainMap, List<DebugEntry> debugEntryList, Map<String, Double> devChainScoreMap) {
        for (short link_site_1 : chainEntry.link_site_set) {
            int precursor_charge = spectrumEntry.precursor_charge;
            SparseBooleanVector theo_mz = mass_tool_obj.buildTheoVector(chainEntry.seq, link_site_1, spectrumEntry.precursor_mass - chainEntry.chain_mass, precursor_charge, specMaxBinIdx);

            // Calculate dot produce
            double dot_product = theo_mz.dot(xcorrPL) * 0.005;

            if (ECL2.debug) {
                debugEntryList.add(new DebugEntry(chainEntry.seq, link_site_1, chainEntry.chain_mass, dot_product));
            }

            if (ECL2.dev && (dot_product > single_chain_t)) {
                devChainScoreMap.put(chainEntry.seq + "-" + link_site_1, dot_product);
            }

            // Record result
            if (dot_product > single_chain_t) {
                if (binChainMap.containsKey(binInx)) {
                    ChainResultEntry chain_result_entry = binChainMap.get(binInx);
                    chain_result_entry.addToScoreList(dot_product, cal_evalue);
                    if (dot_product > chain_result_entry.getScore()) {
                        chain_result_entry.setSecondScore(chain_result_entry.getScore());
                        chain_result_entry.setSecondSeq(chain_result_entry.getSeq());
                        chain_result_entry.setScore(dot_product);
                        chain_result_entry.setSeq(chainEntry.seq);
                        chain_result_entry.setLinkSite(link_site_1);
                    } else if (dot_product > chain_result_entry.getSecondScore()) {
                        chain_result_entry.setSecondScore(dot_product);
                        chain_result_entry.setSecondSeq(chainEntry.seq);
                    }
                } else {
                    ChainResultEntry chain_result_entry = new ChainResultEntry();
                    chain_result_entry.addToScoreList(dot_product, cal_evalue);
                    chain_result_entry.setSeq(chainEntry.seq);
                    chain_result_entry.setLinkSite(link_site_1);
                    chain_result_entry.setScore(dot_product);
                    chain_result_entry.setSecondScore(0);
                    binChainMap.put(binInx, chain_result_entry);
                }
            }
        }
    }

    private int getC13Num(float exp_mass, float theo_mass) {
        float min_diff = 10;
        int num = 0;

        for (int i : C13_correction_range) {
            float temp = Math.abs(exp_mass + i * 1.00335483f - theo_mass);
            if (temp < min_diff) {
                min_diff = temp;
                num = i;
            }
        }

        return num;
    }

    private String addFixMod(String seq, int linkSite) {
        AA[] aaList = MassTool.seqToAAList(seq);
        StringBuilder sb = new StringBuilder(seq.length() * 3);
        for (int i = 0; i < aaList.length; ++i) {
            AA aa = aaList[i];
            if (i == linkSite) { // priority order: linkSite > fixMod > varMod
                sb.append(aa.aa);
            } else if (Math.abs(fix_mod_map.get(aa.aa)) > 1e-6) {
                sb.append(String.format(Locale.US, "%c[%.3f]", aa.aa, fix_mod_map.get(aa.aa)));
            } else {
                if (Math.abs(aa.delta_mass) > 1e-6) {
                    sb.append(String.format(Locale.US, "%c[%.3f]", aa.aa, aa.delta_mass));
                } else {
                    sb.append(aa.aa);
                }
            }
        }
        return sb.toString();
    }
}