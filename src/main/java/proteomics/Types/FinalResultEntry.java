package proteomics.Types;


public class FinalResultEntry implements Comparable<FinalResultEntry> {

    public final int scan_num;
    public final String spectrum_id;
    public final int rank;
    public final int charge;
    public final float spectrum_mz;
    public final float spectrum_mass;
    public final float peptide_mass;
    public final float rt;
    public final float ppm;
    public final double score;
    public final double delta_c;
    public final String seq_1;
    public final int link_site_1;
    public final String pro_id_1;
    public final String seq_2;
    public final int link_site_2;
    public final String pro_id_2;
    public final String cl_type;
    public final int hit_type; // 0 = T-T; 1 = D-D; 2 = T-D;
    public final int C13_correction;
    public final double e_value;
    public final double negative_log10_evalue;
    public float qvalue = -1;
    public final String mgfTitle;

    private final String toString;
    private final int hashCode;

    public final long candidate_num;

    public int point_count;
    public float r_square;
    public float slope;
    public float intercept;
    public int start_idx;
    public int end_idx;

    public double chain_score_1;
    public int chain_rank_1;
    public double chain_score_2;
    public int chain_rank_2;

    private final boolean cal_evalue;

    public FinalResultEntry(int scan_num, String spectrum_id, int rank, int charge, float spectrum_mz, float spectrum_mass, float peptide_mass, float rt, float ppm, double score, double delta_c, String seq_1, int link_site_1, String pro_id_1, String seq_2, int link_site_2, String pro_id_2, String cl_type, int hit_type, int C13_correction, double e_value, int point_count, float r_square, float slope, float intercept, int start_idx, int end_idx, double chain_score_1, int chain_rank_1, double chain_score_2, int chain_rank_2, long candidate_num, boolean cal_evalue, String mgfTitle) {
        this.scan_num = scan_num;
        this.spectrum_id = spectrum_id;
        this.rank = rank;
        this.charge = charge;
        this.spectrum_mz = spectrum_mz;
        this.spectrum_mass = spectrum_mass;
        this.peptide_mass = peptide_mass;
        this.rt = rt;
        this.ppm = ppm;
        this.score = score;
        this.delta_c = delta_c;
        if (seq_1.compareTo(seq_2) < 0) {
            this.seq_1 = seq_1;
            this.link_site_1 = link_site_1;
            this.pro_id_1 = pro_id_1;
            this.seq_2 = seq_2;
            this.link_site_2 = link_site_2;
            this.pro_id_2 = pro_id_2;
        } else {
            this.seq_1 = seq_2;
            this.link_site_1 = link_site_2;
            this.pro_id_1 = pro_id_2;
            this.seq_2 = seq_1;
            this.link_site_2 = link_site_1;
            this.pro_id_2 = pro_id_1;
        }
        this.cl_type = cl_type;
        this.hit_type = hit_type;
        this.C13_correction = C13_correction;
        this.e_value = e_value;
        negative_log10_evalue = -1 * Math.log10(e_value);

        this.point_count = point_count;
        this.r_square = r_square;
        this.slope = slope;
        this.intercept = intercept;
        this.start_idx = start_idx;
        this.end_idx = end_idx;

        this.chain_score_1 = chain_score_1;
        this.chain_rank_1 = chain_rank_1;
        this.chain_score_2 = chain_score_2;
        this.chain_rank_2 = chain_rank_2;

        this.candidate_num = candidate_num;

        toString = scan_num + "-" + seq_1 + "-" + link_site_1 + "-" + seq_2 + link_site_2;
        hashCode = toString.hashCode();

        this.cal_evalue = cal_evalue;

        this.mgfTitle = mgfTitle;
    }

    public int compareTo(FinalResultEntry other) {
        if (cal_evalue) {
            if (this.negative_log10_evalue > other.negative_log10_evalue) {
                return 1;
            } else if (this.negative_log10_evalue < other.negative_log10_evalue) {
                return -1;
            } else {
                return 0;
            }
        } else {
            if (this.score > other.score) {
                return 1;
            } else if (this.score < other.score) {
                return -1;
            } else {
                return 0;
            }
        }
    }

    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        if (other instanceof FinalResultEntry) {
            return ((FinalResultEntry) other).hashCode == hashCode;
        } else {
            return false;
        }
    }
}
