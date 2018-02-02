package proteomics.Search;

import proteomics.ECL2;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.*;
import proteomics.Index.BuildIndex;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class Search {

    public final double ms1_tolerance;
    public final int ms1_tolerance_unit;
    private final Map<String, ChainEntry> chain_entry_map;
    private final MassTool mass_tool_obj;
    private final TreeMap<Integer, Set<String>> bin_seq_map;
    private final BuildIndex build_index_obj;
    private final int[] C13_correction_range;
    final double single_chain_t;
    private final boolean cal_evalue;

    /////////////////////////////////////////public methods////////////////////////////////////////////////////////////
    public Search(BuildIndex build_index_obj, Map<String, String> parameter_map) {
        this.build_index_obj = build_index_obj;
        chain_entry_map = build_index_obj.getSeqEntryMap();
        mass_tool_obj = build_index_obj.returnMassTool();
        ms1_tolerance_unit = Integer.valueOf(parameter_map.get("ms1_tolerance_unit"));
        ms1_tolerance = Double.valueOf(parameter_map.get("ms1_tolerance"));
        bin_seq_map = build_index_obj.getMassBinSeqMap();
        single_chain_t = Double.valueOf(parameter_map.get("single_chain_t"));
        if (parameter_map.get("cal_evalue").contentEquals("0")) {
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

    ResultEntry doSearch(SparseVector xcorrPL, TreeMap<Integer, List<Double>> binScoresMap, double massWithoutLinker, int precursorCharge, double precursorMass, String scanId) throws IOException {
        int max_chain_bin_idx = Math.min(build_index_obj.massToBin(massWithoutLinker + C13_correction_range[C13_correction_range.length - 1] * MassTool.C13_DIFF) + 1 - bin_seq_map.firstKey(), bin_seq_map.lastKey());
        int min_chain_bin_idx = bin_seq_map.firstKey();

        if (max_chain_bin_idx < min_chain_bin_idx) {
            return null;
        }

        // set MS1 tolerance for further usage.
        double leftMs1Tol = ms1_tolerance;
        double rightMs1Tol = ms1_tolerance;
        if (ms1_tolerance_unit == 1) {
            leftMs1Tol = precursorMass - (precursorMass / (1 + ms1_tolerance * 1e-6));
            rightMs1Tol = (precursorMass / (1 - ms1_tolerance * 1e-6)) - precursorMass;
        }

        // for debug
        List<DebugEntry> debugEntryList = new LinkedList<>();
        Map<String, Double> devChainScoreMap = new HashMap<>();

        // start search
        TreeMap<Integer, ChainResultEntry> binChainMap = new TreeMap<>();
        Set<Integer> checkedBinSet = new HashSet<>(bin_seq_map.size() + 1, 1);
        long candidate_num = 0;
        ResultEntry resultEntry = new ResultEntry(scanId, cal_evalue, binChainMap);
        int stopIdx = (int) Math.ceil(build_index_obj.massToBin((massWithoutLinker + MassTool.C13_DIFF) * 0.5 + rightMs1Tol));
        for (int idx_1 : bin_seq_map.keySet()) {
            if (idx_1 > stopIdx) {
                break;
            }

            double left_mass_2;
            double right_mass_2;
            int left_idx_2;
            int right_idx_2;
            NavigableMap<Integer, Set<String>> sub_map = new TreeMap<>();

            // consider C13 correction
            for (int i : C13_correction_range) {
                left_mass_2 = massWithoutLinker + i * MassTool.C13_DIFF - build_index_obj.binToRightMass(idx_1) - leftMs1Tol;
                right_mass_2 = massWithoutLinker + i * MassTool.C13_DIFF - build_index_obj.binToLeftMass(idx_1) + rightMs1Tol;
                left_idx_2 = build_index_obj.massToBin(left_mass_2);
                right_idx_2 = build_index_obj.massToBin(right_mass_2) + 1;
                sub_map.putAll(bin_seq_map.subMap(left_idx_2, true, right_idx_2, true));
            }

            if (!sub_map.isEmpty()) {
                if (!checkedBinSet.contains(idx_1)) {
                    // This bin hasn't been visited. Linear scan first.
                    for (String seq : bin_seq_map.get(idx_1)) {
                        ChainEntry chainEntry = chain_entry_map.get(seq);
                        linearScan(xcorrPL, chainEntry, idx_1, binChainMap, binScoresMap, debugEntryList, devChainScoreMap, precursorCharge, precursorMass);
                    }
                    checkedBinSet.add(idx_1);
                }

                if (binChainMap.containsKey(idx_1)) { // There may be no chain with chain score > single_chain_t
                    ChainResultEntry chain_score_entry_1 = binChainMap.get(idx_1);
                    for (int idx_2 : sub_map.keySet()) {
                        if (!checkedBinSet.contains(idx_2)) {
                            // this bin hasn't been visited. Linear scan first.
                            for (String seq : bin_seq_map.get(idx_2)) {
                                ChainEntry chainEntry = chain_entry_map.get(seq);
                                linearScan(xcorrPL, chainEntry, idx_2, binChainMap, binScoresMap, debugEntryList, devChainScoreMap, precursorCharge, precursorMass);
                            }
                            checkedBinSet.add(idx_2);
                        }

                        if (binChainMap.containsKey(idx_2)) { // There may be no chain with chain score > single_chain_t
                            ChainResultEntry chain_score_entry_2 = binChainMap.get(idx_2);

                            // only two sequences with the same binary mod type can be linked.
                            if (chain_entry_map.get(chain_score_entry_1.getSeq()).binaryModType == chain_entry_map.get(chain_score_entry_2.getSeq()).binaryModType) {
                                double score = chain_score_entry_1.getScore() + chain_score_entry_2.getScore();
                                ++candidate_num;

                                // calculate second score
                                double second_score = 0;
                                double temp_1 = -1;
                                if (chain_score_entry_1.getSecondSeq() != null) {
                                    temp_1 = chain_score_entry_1.getSecondScore() + chain_score_entry_2.getScore();
                                }
                                double temp_2 = -1;
                                if (chain_score_entry_2.getSecondSeq() != null) {
                                    temp_2 = chain_score_entry_1.getScore() + chain_score_entry_2.getSecondScore();
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
                                    for (double s1 : binScoresMap.get(idx_1)) {
                                        for (double s2 : binScoresMap.get(idx_2)) {
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
            BufferedWriter writer = new BufferedWriter(new FileWriter(scanId + ".csv"));
            debugEntryList.sort(Comparator.reverseOrder());
            writer.write("chain,link_site,mass,score\n");
            for (DebugEntry t : debugEntryList) {
                writer.write(String.format(Locale.US, "%s,%d,%f,%f\n", t.chain, t.link_site, t.mass, t.score));
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

    private void linearScan(SparseVector xcorrPL, ChainEntry chainEntry, int binInx, TreeMap<Integer, ChainResultEntry> binChainMap, TreeMap<Integer, List<Double>> binScoresMap, List<DebugEntry> debugEntryList, Map<String, Double> devChainScoreMap, int precursorCharge, double precursorMass) {
        for (short link_site_1 : chainEntry.link_site_set) {
            // generate theoretical fragment ion bins and calculate XCorr.
            double dot_product = mass_tool_obj.generateTheoFragmentAndCalXCorr(chainEntry.seq, link_site_1, precursorMass - chainEntry.chain_mass, precursorCharge, xcorrPL);

            if (ECL2.debug) {
                debugEntryList.add(new DebugEntry(chainEntry.seq, link_site_1, chainEntry.chain_mass, dot_product));
            }

            if (ECL2.dev && (dot_product > single_chain_t)) {
                devChainScoreMap.put(chainEntry.seq + "-" + link_site_1, dot_product);
            }

            // Record result
            if (dot_product > single_chain_t) {
                if (cal_evalue) {
                    if (binScoresMap.containsKey(binInx)) {
                        binScoresMap.get(binInx).add(dot_product);
                    } else {
                        List<Double> tempList = new ArrayList<>();
                        tempList.add(dot_product);
                        binScoresMap.put(binInx, tempList);
                    }
                }
                if (binChainMap.containsKey(binInx)) {
                    ChainResultEntry chain_result_entry = binChainMap.get(binInx);
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
                    chain_result_entry.setSeq(chainEntry.seq);
                    chain_result_entry.setLinkSite(link_site_1);
                    chain_result_entry.setScore(dot_product);
                    chain_result_entry.setSecondScore(0);
                    binChainMap.put(binInx, chain_result_entry);
                }
            }
        }
    }
}