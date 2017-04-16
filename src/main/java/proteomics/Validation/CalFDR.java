package proteomics.Validation;


import proteomics.ECL2;
import proteomics.Types.FinalResultEntry;

import java.util.List;

public class CalFDR {

    private final float precision;

    private double min_score = 0;
    private double max_score = 0;
    private float[] qvalue_array = null;
    private List<FinalResultEntry> results;

    public CalFDR(List<FinalResultEntry> results) {
        this.results = results;
        if (ECL2.cal_evalue) {
            precision = 0.1f;
        } else {
            precision = 0.001f;
        }

        // find the max score
        for (FinalResultEntry entry : results) {
            if (ECL2.cal_evalue) {
                if (entry.negative_log10_evalue > max_score) {
                    max_score = entry.negative_log10_evalue;
                }
                if (entry.negative_log10_evalue < min_score) {
                    min_score = entry.negative_log10_evalue;
                }
            } else {
                if (entry.score > max_score) {
                    max_score = entry.score;
                }
                if (entry.score < min_score) {
                    min_score = entry.score;
                }
            }
        }

        final int array_length = 1 + (int) Math.ceil((max_score - min_score) / precision);
        float[] decoy_count_vector = new float[array_length];
        float[] target_count_vector = new float[array_length];
        float[] fuse_count_vector = new float[array_length];
        float[] fdr_array = new float[array_length];
        qvalue_array = new float[array_length];

        for (FinalResultEntry re : results) {
            if (re.hit_type == 1) {
                int idx;
                if (ECL2.cal_evalue) {
                    idx = (int) Math.floor((re.negative_log10_evalue - min_score) / precision);
                } else {
                    idx = (int) Math.floor((re.score - min_score) / precision);
                }
                ++decoy_count_vector[idx];
            } else if (re.hit_type == 0) {
                int idx;
                if (ECL2.cal_evalue) {
                    idx = (int) Math.floor((re.negative_log10_evalue - min_score) / precision);
                } else {
                    idx = (int) Math.floor((re.score - min_score) / precision);
                }
                ++target_count_vector[idx];
            } else {
                int idx;
                if (ECL2.cal_evalue) {
                    idx = (int) Math.floor((re.negative_log10_evalue - min_score) / precision);
                } else {
                    idx = (int) Math.floor((re.score - min_score) / precision);
                }
                ++fuse_count_vector[idx];
            }
        }

        // Calculate FDR
        for (int idx_1 = 0; idx_1 < array_length - 1; ++idx_1) {
            int decoy_count = 0;
            int fuse_count = 0;
            int target_count = 0;
            for (int idx_2 = idx_1 + 1; idx_2 < array_length; ++idx_2) {
                decoy_count += decoy_count_vector[idx_2];
                target_count += target_count_vector[idx_2];
                fuse_count += fuse_count_vector[idx_2];
            }

            float fdr;
            if (target_count == 0) {
                fdr = 0;
            } else if (fuse_count < decoy_count) {
                fdr = 0;
            } else {
                fdr = (float) (fuse_count - decoy_count) / (float) target_count;
            }

            fdr = Math.min(fdr, 1); // Adjust those fdrs that are larger than 1
            fdr_array[idx_1] = fdr;
        }

        // Convert FDR to qvalue
        float last_q_value = fdr_array[0];
        qvalue_array[0] = last_q_value;

        for (int idx_1 = 1; idx_1 < array_length; ++idx_1) {
            float q_value = fdr_array[idx_1];
            if (q_value > last_q_value) {
                qvalue_array[idx_1] = last_q_value;
            } else {
                qvalue_array[idx_1] = q_value;
                last_q_value = q_value;
            }
        }
    }

    public List<FinalResultEntry> includeStats() {
        for (FinalResultEntry re : results) {
            if (re.hit_type == 0) {
                if (ECL2.cal_evalue) {
                    int idx = (int) Math.floor((re.negative_log10_evalue - min_score) / precision);
                    re.qvalue = qvalue_array[idx];
                } else {
                    int idx = (int) Math.floor((re.score - min_score) / precision);
                    re.qvalue = qvalue_array[idx];
                }
            }
        }

        return results;
    }
}