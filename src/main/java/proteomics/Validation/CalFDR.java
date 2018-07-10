package proteomics.Validation;


import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.sql.*;
import java.util.*;

public class CalFDR {

    private static final Logger logger = LoggerFactory.getLogger(CalFDR.class);

    public static Map<String, Entry> calFDR(String sqlPath, boolean cal_evalue, String clType) throws SQLException {
        double min_score = 999;
        double max_score = -999;
        double inversePrecision;
        if (cal_evalue) {
            inversePrecision = 10;
        } else {
            inversePrecision = 1000;
        }


        String scoreType = "score";
        if (cal_evalue) {
            scoreType = "eValue";
        }

        Map<String, Entry> scanIdEntryMap = new HashMap<>();
        Connection sqlConnection = DriverManager.getConnection(sqlPath);
        Statement sqlStatement = sqlConnection.createStatement();
        ResultSet sqlResultSet = sqlStatement.executeQuery(String.format(Locale.US, "SELECT scanId, %s, hitType FROM spectraTable WHERE clType='%s'", scoreType, clType));
        while (sqlResultSet.next()) {
            int hitType = sqlResultSet.getInt("hitType");
            if (!sqlResultSet.wasNull()) {
                String scanId = sqlResultSet.getString("scanId");
                double score = sqlResultSet.getDouble(scoreType);
                if (cal_evalue) {
                    score = -1 * Math.log10(score);
                }
                scanIdEntryMap.put(scanId, new Entry(scanId, score, hitType));
                if (score > max_score) {
                    max_score = score;
                }
                if (score < min_score) {
                    min_score = score;
                }
            }
        }
        sqlResultSet.close();
        sqlStatement.close();
        sqlConnection.close();

        if (scanIdEntryMap.size() < 500) {
            logger.warn("For {}, There are only {} hits in total, which is not enough to estimate an accurate FDR/q-value.", clType, scanIdEntryMap.size());
        }

        final int array_length = 1 + (int) Math.ceil((max_score - min_score) * inversePrecision);
        if (array_length <= 0) { // there is no hit
            return scanIdEntryMap;
        }
        double[] decoy_count_vector = new double[array_length];
        double[] target_count_vector = new double[array_length];
        double[] fuse_count_vector = new double[array_length];
        double[] fdr_array = new double[array_length];
        double[] qvalue_array = new double[array_length];

        for (Entry re : scanIdEntryMap.values()) {
            if (re.hitType == 1) {
                int idx = (int) Math.round((re.score - min_score) * inversePrecision);
                ++decoy_count_vector[idx];
            } else if (re.hitType == 0) {
                int idx = (int) Math.round((re.score - min_score) * inversePrecision);
                ++target_count_vector[idx];
            } else {
                int idx = (int) Math.round((re.score - min_score) * inversePrecision);
                ++fuse_count_vector[idx];
            }
        }

        // Calculate FDR
        boolean printWarning = false;
        for (int idx_1 = 0; idx_1 < array_length - 1; ++idx_1) {
            int decoy_count = 0;
            int fuse_count = 0;
            int target_count = 0;
            for (int idx_2 = idx_1 + 1; idx_2 < array_length; ++idx_2) {
                decoy_count += decoy_count_vector[idx_2];
                target_count += target_count_vector[idx_2];
                fuse_count += fuse_count_vector[idx_2];
            }

            double fdr;
            if (target_count == 0) {
                fdr = 0;
            } else if (fuse_count < decoy_count) {
                printWarning = true;
                fdr = 0;
            } else {
                fdr = (double) (fuse_count - decoy_count) / (double) target_count;
            }

            fdr = Math.min(fdr, 1); // Adjust those fdrs that are larger than 1
            fdr_array[idx_1] = fdr;
        }

        if (printWarning) {
            logger.warn("For {}, there are less T-D hits than D-D hits. The estimated FDR/q-value is not accurate.", clType);
        }

        // Convert FDR to qvalue
        double last_q_value = fdr_array[0];
        qvalue_array[0] = last_q_value;

        for (int idx_1 = 1; idx_1 < array_length; ++idx_1) {
            double q_value = fdr_array[idx_1];
            if (q_value > last_q_value) {
                qvalue_array[idx_1] = last_q_value;
            } else {
                qvalue_array[idx_1] = q_value;
                last_q_value = q_value;
            }
        }

        for (Entry entry : scanIdEntryMap.values()) {
            entry.qValue = qvalue_array[(int) Math.round((entry.score - min_score) * inversePrecision)];
        }

        return scanIdEntryMap;
    }

    public static class Entry implements Comparable<Entry> {

        public final String scanId;
        public final double score;
        public final int hitType;

        public double qValue;

        public Entry(String scanId, double score, int hitType) {
            this.scanId = scanId;
            this.score = score;
            this.hitType = hitType;
        }

        public boolean equals(Object other) {
            if (other instanceof Entry) {
                return ((Entry) other).scanId.contentEquals(scanId) && ((Entry) other).score == score && ((Entry) other).hitType == hitType;
            } else {
                return false;
            }
        }

        public int compareTo(Entry other) {
            if (score > other.score) {
                return 1;
            } else if (score < other.score) {
                return -1;
            } else {
                if (hitType == 0 && other.hitType != 0) {
                    return 1;
                } else if (hitType != 0 && other.hitType == 0) {
                    return -1;
                } else {
                    return 0;
                }
            }
        }
    }
}