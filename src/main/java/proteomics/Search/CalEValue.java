package proteomics.Search;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.ECL2;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.ResultEntry;
import proteomics.Types.SparseBooleanVector;
import proteomics.Types.SparseVector;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.NavigableMap;
import java.util.Set;
import java.util.TreeMap;

public class CalEValue {

    private static final Logger logger = LoggerFactory.getLogger(CalEValue.class);

    private ResultEntry result_entry;
    private TreeMap<Float, Set<String>> uniprot_decoy_mass_seq_map;
    private float linker_mass;
    private MassTool mass_tool_obj;
    private int max_common_ion_charge;
    private SparseVector pl_map_xcorr;
    private boolean consider_two_identical_chains;
    private float e_value_precursor_mass_tol;

    CalEValue(int scan_num, ResultEntry result_entry, SparseVector pl_map_xcorr, TreeMap<Float, Set<String>> uniprot_decoy_mass_seq_map, MassTool mass_tool_obj, float linker_mass, int max_common_ion_charge, boolean consider_two_identical_chains, float e_value_precursor_mass_tol) {
        this.result_entry = result_entry;
        this.uniprot_decoy_mass_seq_map = uniprot_decoy_mass_seq_map;
        this.linker_mass = linker_mass;
        this.mass_tool_obj = mass_tool_obj;
        this.max_common_ion_charge = max_common_ion_charge;
        this.pl_map_xcorr = pl_map_xcorr;
        this.consider_two_identical_chains = consider_two_identical_chains;
        this.e_value_precursor_mass_tol = e_value_precursor_mass_tol;

        int gap_num = ECL2.score_point_t - result_entry.getScoreCount();
        if (gap_num >= 0) {
            generateRandomRandomScores(gap_num);
        }

        int[] score_histogram = result_entry.getScoreHistogram();
        double histogram_bin_size = result_entry.getHistogramBinSize();

        // find the max non-zero idx. It may be an outlier.
        int max_nonzero_idx = 0;
        for (int i = 0; i < score_histogram.length; ++i) {
            if (score_histogram[i] > 0) {
                max_nonzero_idx = i;
            }
        }

        // generate survival array
        int[] survival_count_array = new int[max_nonzero_idx + 1];
        survival_count_array[max_nonzero_idx] = score_histogram[max_nonzero_idx];
        for (int i = max_nonzero_idx - 1; i >= 0; --i) {
            survival_count_array[i] = survival_count_array[i + 1] + score_histogram[i];
        }

        // find the next max nonzero idx. from this down to the changing point are from null.
        int top_count = (int) (survival_count_array[0] * 0.01); // these are from non-null.
        int null_end_idx = 0;
        for (int i = max_nonzero_idx - 1; i >= 0; --i) {
            if (survival_count_array[i] > top_count) {
                null_end_idx = i;
                break;
            }
        }

        // log transform. eliminate the non-null points
        double[] ln_survival_count_array = new double[null_end_idx + 1];
        for (int i = 0; i <= null_end_idx; ++i) {
            ln_survival_count_array[i] = Math.log(survival_count_array[i]);
        }

        // try and find the best linear regression result
        int start_idx = (int) (null_end_idx * 0.5);
        double max_r_square = 0;
        double optimal_slope = 1;
        double optimal_intercept = 0;
        int optimal_start_idx = 0;
        while (start_idx >= 0) {
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
            if (r_square > max_r_square) {
                optimal_slope = slope;
                optimal_intercept = intercept;
                max_r_square = r_square;
                optimal_start_idx = start_idx;
            }
            --start_idx;
        }

        if (optimal_slope >= 0) {
            result_entry.setEValue(9999);
            logger.debug("Estimating E-value failed. Scan: {}, mass: {}, slope: {}, intercept: {}, R square: {}, point num: {}.",scan_num, result_entry.spectrum_mass, optimal_slope, optimal_intercept, max_r_square, result_entry.getScoreCount());
        } else {
            double p_value = Math.exp(optimal_slope * Math.round(result_entry.getScore() / histogram_bin_size) + optimal_intercept) / result_entry.getScoreCount();
            result_entry.setEValue(p_value * result_entry.getCandidateNum());
            result_entry.setEValueDetails((float) max_r_square, (float) optimal_slope, (float) optimal_intercept, optimal_start_idx, null_end_idx);
        }

        if (ECL2.debug) {
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(scan_num + ".evalue.csv"))) {
                writer.write(String.format("histogram,survival,ln(survival),slope=%.4f,intercept=%.4f,rsquare=%.4f,start=%d,end=%d\n", optimal_slope, optimal_intercept, max_r_square, optimal_start_idx, null_end_idx));
                for (int i = 0; i < max_nonzero_idx; ++i) {
                    if (i < ln_survival_count_array.length) {
                        writer.write(String.format("%d,%d,%.4f\n", score_histogram[i], survival_count_array[i], ln_survival_count_array[i]));
                    } else {
                        writer.write(String.format("%d,%d\n", score_histogram[i], survival_count_array[i]));
                    }
                }
            } catch (IOException ex) {
                logger.error(ex.getMessage());
                ex.printStackTrace();
                System.exit(1);
            }
        }
    }

    private int generateRandomRandomScores(int gap_num) {
        float precursor_mass = result_entry.spectrum_mass;
        float max_mass = (precursor_mass - linker_mass) / 2;
        for (float mass_1 : uniprot_decoy_mass_seq_map.keySet()) {
            if (mass_1 <= max_mass) {
                float left_mass = precursor_mass - linker_mass - mass_1 - e_value_precursor_mass_tol;
                float right_mass = precursor_mass - linker_mass - mass_1 + e_value_precursor_mass_tol;
                NavigableMap<Float, Set<String>> sub_map = uniprot_decoy_mass_seq_map.subMap(left_mass, true, right_mass, true);
                if (!sub_map.isEmpty()) {
                    for (String seq_1_link_site : uniprot_decoy_mass_seq_map.get(mass_1)) {
                        String[] temp = seq_1_link_site.split("-");
                        String seq_1 = temp[0];
                        SparseBooleanVector theo_mz_1 = mass_tool_obj.buildVector(mass_tool_obj.buildPseudoCLIonArray(mass_tool_obj.buildChainIonArray(mass_tool_obj.seqToAAList(seq_1)), link_site, max_common_ion_charge, precursor_mass - mass_1, result_entry.charge), result_entry.charge);
                        short link_site = Short.valueOf(temp[1]);
                        double score1 = theo_mz_1.dot(pl_map_xcorr) * 0.25;
                        if (score1 > Search.single_chain_t) {
                            gap_num = generateMoreScoresSub(sub_map, seq_1, score1, precursor_mass, gap_num);
                            if (gap_num < 0) {
                                return gap_num;
                            }
                        }
                    }
                }
            }
        }
        return gap_num;
    }

    private int generateMoreScoresSub(NavigableMap<Float, Set<String>> sub_map, String seq_1, double score1, float precursor_mass, int gap_num) {
        for (float mass_2 : sub_map.keySet()) {
            for (String seq_2_link_site : sub_map.get(mass_2)) {
                String[] temp = seq_2_link_site.split("-");
                String seq_2 = temp[0];
                short link_site = Short.valueOf(temp[1]);
                if (!seq_1.contentEquals(seq_2) || consider_two_identical_chains) {
                    SparseBooleanVector theo_mz_2 = mass_tool_obj.buildVector(mass_tool_obj.buildPseudoCLIonArray(mass_tool_obj.buildChainIonArray(mass_tool_obj.seqToAAList(seq_2)), link_site, max_common_ion_charge, precursor_mass - mass_2, result_entry.charge), result_entry.charge);
                    double score2 = theo_mz_2.dot(pl_map_xcorr) * 0.25;
                    if (score2 > Search.single_chain_t) {
                        double score = score1 + score2;
                        result_entry.addToScoreHistogram(score);
                        --gap_num;
                        if (gap_num < 0) {
                            return gap_num;
                        }
                    }
                }
            }
        }
        return gap_num;
    }
}
