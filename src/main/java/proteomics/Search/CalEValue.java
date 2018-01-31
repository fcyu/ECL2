package proteomics.Search;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.ECL2;
import proteomics.Index.BuildIndex;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.*;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class CalEValue {

    private static final Logger logger = LoggerFactory.getLogger(CalEValue.class);
    private static final float maxTolerance = 20;
    private static final float toleranceStep = 1;

    static void calEValue(int scan_num, ResultEntry result_entry, BuildIndex buildIndexObj, TreeMap<Integer, List<Double>> binScoresMap, float linker_mass, float originalTolerance, SparseVector xcorrPL, double singleChainT) throws IOException {
        int gap_num = ECL2.score_point_t - result_entry.getScoreCount();
        float tolerance = originalTolerance;
        float massWithoutLinker = result_entry.spectrum_mass - linker_mass;
        int maxBinIdx = buildIndexObj.massToBin(massWithoutLinker * 0.5f);
        while (gap_num > 0 && tolerance <= maxTolerance) {
            gap_num = generateRandomRandomScores(gap_num, tolerance, toleranceStep, binScoresMap, result_entry.spectrum_mass, massWithoutLinker, result_entry.charge, xcorrPL, buildIndexObj, buildIndexObj.getMassBinSeqMap(), buildIndexObj.getSeqEntryMap(), buildIndexObj.returnMassTool(), result_entry, maxBinIdx, singleChainT);
            tolerance += toleranceStep;
        }

        if (gap_num > 0) {
            logger.debug("Scan {}: Estimated e-value may not be accurate ({} data points).", scan_num, ECL2.score_point_t - gap_num);
        }

        int[] score_histogram = result_entry.getScoreHistogram();
        double inverseHistogramBinSize = ResultEntry.getInverseHistogramBinSize();

        int min_nonzero_idx = 0;
        for (int i = 0; i < score_histogram.length; ++i) {
            if (score_histogram[i] > 0) {
                min_nonzero_idx = i;
                break;
            }
        }
        if (min_nonzero_idx == score_histogram.length - 1) {
            logger.debug("Fail to estimate e-value for Scan {} (an empty score histogram).", scan_num);
            result_entry.setEValue(9999);
            return;
        }

        int max_nonzero_idx = 0;
        for (int i = score_histogram.length - 1; i >= 0; --i) {
            if (score_histogram[i] > 0) {
                max_nonzero_idx = i;
                break;
            }
        }
        if (max_nonzero_idx - min_nonzero_idx < 5) {
            logger.debug("Fail to estimate e-value for Scan {} (doesn't have a useful score histogram).", scan_num);
            result_entry.setEValue(9999);
            if (ECL2.debug) {
                BufferedWriter writer = new BufferedWriter(new FileWriter(scan_num + ".evalue.csv"));
                writer.write("histogram\n");
                for (int i = 0; i < max_nonzero_idx; ++i) {
                    writer.write(String.format(Locale.US, "%d\n", score_histogram[i]));
                }
                writer.close();
            }
            return;
        }

        // generate survival array
        int[] survival_count_array = new int[max_nonzero_idx + 1];
        survival_count_array[max_nonzero_idx] = score_histogram[max_nonzero_idx];
        for (int i = max_nonzero_idx - 1; i >= 0; --i) {
            survival_count_array[i] = survival_count_array[i + 1] + score_histogram[i];
        }

        // find the next max nonzero idx. from this down to the changing point are from null.
        int null_end_idx = max_nonzero_idx;
        for (int i = max_nonzero_idx; i > min_nonzero_idx + 4; --i) {
            if (score_histogram[i] == 0 && score_histogram[i - 1] == 0 && score_histogram[i - 2] != 0) {
                null_end_idx = i - 2;
                break;
            }
        }

        // log transform. eliminate the non-null points
        double[] ln_survival_count_array = new double[null_end_idx + 1];
        for (int i = 0; i <= null_end_idx; ++i) {
            ln_survival_count_array[i] = Math.log(survival_count_array[i]);
        }

        // find the start point
        int start_idx;
        if (null_end_idx > 3 * inverseHistogramBinSize) {
            start_idx = Math.max(min_nonzero_idx, (int) (0.75 * null_end_idx));
        } else {
            start_idx = Math.max(min_nonzero_idx, (int) (0.5 * null_end_idx));
        }

        // linear regression
        double max_r_square = 0;
        double optimal_slope = 1;
        double optimal_intercept = 0;
        int optimal_start_idx = 0;
        while (start_idx >= min_nonzero_idx) {
            double x_sum = 0;
            double y_sum = 0;
            double x_mean = 0;
            double y_mean = 0;
            double xx_mean = 0;
            double yy_mean = 0;
            double xy_mean = 0;
            int point_num = 0;
            double r_square;
            double slope;
            double intercept;

            for (int i = start_idx; i <= null_end_idx; i++) {
                if (survival_count_array[i] > 0) {
                    x_sum += i;
                    y_sum += ln_survival_count_array[i];
                    ++point_num;
                }
            }

            if (point_num > 0) {
                x_mean = x_sum / point_num;
                y_mean = y_sum / point_num;
            }

            for (int i = start_idx; i <= null_end_idx; i++) {
                if (survival_count_array[i] > 0) {
                    double d_x = i - x_mean;
                    double d_y = ln_survival_count_array[i] - y_mean;
                    xx_mean += d_x * d_x;
                    yy_mean += d_y * d_y;
                    xy_mean += d_x * d_y;
                }
            }

            if (xx_mean > 0) {
                slope = xy_mean / xx_mean;
            } else {
                slope = 0;
            }

            intercept = y_mean - slope * x_mean;

            double ssr = 0.0;
            for (int i = start_idx; i <= null_end_idx; i++) {
                if (survival_count_array[i] > 0) {
                    double fit = slope * i + intercept;
                    ssr += (fit - y_mean) * (fit - y_mean);
                }
            }

            r_square = ssr / yy_mean;
            if (slope < 0) {
                optimal_slope = slope;
                optimal_intercept = intercept;
                max_r_square = r_square;
                optimal_start_idx = start_idx;
                break;
            }
            --start_idx;
        }

        if (optimal_slope >= 0) {
            result_entry.setEValue(9999);
            logger.debug("Estimating E-value failed. Scan: {}, mass: {}, slope: {}, intercept: {}, R square: {}, point num: {}.",scan_num, result_entry.spectrum_mass, optimal_slope, optimal_intercept, max_r_square, result_entry.getScoreCount());
        } else {
            result_entry.setEValue(Math.exp((optimal_slope * Math.round(result_entry.getScore() * inverseHistogramBinSize) + optimal_intercept) + Math.log((double) result_entry.getCandidateNum() / (double) result_entry.getScoreCount()))); // double point precision limitation.
            result_entry.setEValueDetails((float) max_r_square, (float) optimal_slope, (float) optimal_intercept, optimal_start_idx, null_end_idx);
        }

        if (ECL2.debug) {
            BufferedWriter writer = new BufferedWriter(new FileWriter(scan_num + ".evalue.csv"));
            writer.write(String.format(Locale.US, "histogram,survival,ln(survival),slope=%.4f,intercept=%.4f,rsquare=%.4f,start=%d,end=%d\n", optimal_slope, optimal_intercept, max_r_square, optimal_start_idx, null_end_idx));
            for (int i = 0; i <= max_nonzero_idx; ++i) {
                if (i < ln_survival_count_array.length) {
                    writer.write(String.format(Locale.US, "%d,%d,%.4f\n", score_histogram[i], survival_count_array[i], ln_survival_count_array[i]));
                } else {
                    writer.write(String.format(Locale.US, "%d,%d\n", score_histogram[i], survival_count_array[i]));
                }
            }
            writer.close();
        }
    }

    private static int generateRandomRandomScores(int gap_num, float tolerance, float toleranceStep, TreeMap<Integer, List<Double>> binScoresMap, float precursorMass, float massWithoutLinker, int precursorCharge, SparseVector xcorrPL, BuildIndex buildIndex, TreeMap<Integer, Set<String>> binSequencesMap, Map<String, ChainEntry> seqEntryMap, MassTool massTool, ResultEntry resultEntry, int maxBinIdx, double singleChainT) {
        for (int binIdx1 : binSequencesMap.keySet()) {
            if (binIdx1 < maxBinIdx) {
                int leftBinIdx1 = buildIndex.massToBin(massWithoutLinker - tolerance - toleranceStep) - binIdx1;
                int rightBinIdx1 = buildIndex.massToBin(massWithoutLinker - tolerance) - binIdx1;
                int leftBinIdx2 = buildIndex.massToBin(massWithoutLinker + tolerance) - binIdx1;
                int rightBinIdx2 = buildIndex.massToBin(massWithoutLinker + tolerance + toleranceStep) - binIdx1;
                TreeMap<Integer, Set<String>> subMap = new TreeMap<>();
                subMap.putAll(binSequencesMap.subMap(leftBinIdx1, true, rightBinIdx1, false));
                subMap.putAll(binSequencesMap.subMap(leftBinIdx2, false, rightBinIdx2, true));
                if (!subMap.isEmpty()) {
                    for (double score1 : subFunction(binIdx1, binScoresMap, binSequencesMap, seqEntryMap, massTool, precursorMass, precursorCharge, xcorrPL, singleChainT)) {
                        for (int binIdx2 : subMap.keySet()) {
                            for (double score2 : subFunction(binIdx2, binScoresMap, binSequencesMap, seqEntryMap, massTool, precursorMass, precursorCharge, xcorrPL, singleChainT)) {
                                resultEntry.addToScoreHistogram(score1 + score2);
                                --gap_num;
                                if (gap_num <= 0) {
                                    return gap_num;
                                }
                            }
                        }
                    }
                }
            }
        }
        return gap_num;
    }

    private static List<Double> subFunction(int binIdx, TreeMap<Integer, List<Double>> binScoresMap, TreeMap<Integer, Set<String>> binSequencesMap, Map<String, ChainEntry> seqEntryMap, MassTool massTool, float precursorMass, int precursorCharge, SparseVector xcorrPL, double singleChainT) {
        List<Double> scoreList = new ArrayList<>();
        if (binScoresMap.containsKey(binIdx)) {
            scoreList = binScoresMap.get(binIdx);
        } else {
            for (String seq : binSequencesMap.get(binIdx)) {
                ChainEntry chainEntry = seqEntryMap.get(seq);
                for (short linkSite : chainEntry.link_site_set) {
                    double xcorr = massTool.generateTheoFragmentAndCalXCorr(seq, linkSite, precursorMass - chainEntry.chain_mass, precursorCharge, xcorrPL);
                    if (xcorr > singleChainT) {
                        scoreList.add(xcorr);
                    }
                }
            }
            binScoresMap.put(binIdx, scoreList);
        }
        return scoreList;
    }
}
