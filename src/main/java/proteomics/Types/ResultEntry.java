package proteomics.Types;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.TreeMap;


public class ResultEntry{

    private static final Logger logger = LoggerFactory.getLogger(ResultEntry.class);
    public static final double inverseHistogramBinSize = 20;
    private static final double max_score = 20;

    public final String spectrum_id;
    public final float spectrum_mass;
    public final float spectrum_mz;
    public final float rt;
    public final int charge;

    private final TreeMap<Integer, ChainResultEntry> binChainMap;

    private String chain_seq_1;
    private String chain_seq_2;
    private double score;
    private double second_score;
    private int link_site_1;
    private int link_site_2;

    private int[] score_histogram;
    private long candidate_num;
    private int score_count = 0;
    private double e_value = 9999;
    private float r_square;
    private float slope;
    private float intercept;
    private int start_idx;
    private int end_idx;

    private double chain_score_1;
    private int chain_rank_1;
    private double chain_score_2;
    private int chain_rank_2;

    public ResultEntry(String spectrum_id, float spectrum_mz, float spectrum_mass, float rt, int charge, boolean cal_evalue, TreeMap<Integer, ChainResultEntry> binChainMap) {
        if (cal_evalue) {
            score_histogram = new int[(int) Math.round(max_score * inverseHistogramBinSize) + 1]; // start from zero score.
        }
        this.spectrum_id = spectrum_id;
        this.spectrum_mz = spectrum_mz;
        this.spectrum_mass = spectrum_mass;
        this.rt = rt;
        this.charge = charge;
        this.binChainMap = binChainMap;
    }

    public void setChain1(String chain_seq_1) {
        this.chain_seq_1 = chain_seq_1;
    }

    public void setChain2(String chain_seq_2) {
        this.chain_seq_2 = chain_seq_2;
    }

    public void setScore(double score) {
        this.score = score;
    }

    public void setSecondScore(double second_score) {
        this.second_score = second_score;
    }

    public void setLinkSite1(int link_site_1) {
        this.link_site_1 = link_site_1;
    }

    public void setLinkSite2(int link_site_2) {
        this.link_site_2 = link_site_2;
    }

    public void setEValue(double e_value) {
        this.e_value = e_value;
    }

    public void setCandidateNum(long candidate_num) {
        this.candidate_num = candidate_num;
    }

    public void addToScoreHistogram(double score) {
        try {
            ++score_histogram[(int) Math.round(score * inverseHistogramBinSize)];
            ++score_count;
        } catch (ArrayIndexOutOfBoundsException ex) {
            logger.warn("Score {} is out of the range [0, {}].", score, max_score);
        }
    }

    public void setEValueDetails(float r_square, float slope, float intercept, int start_idx, int end_idx) {
        this.r_square = r_square;
        this.slope = slope;
        this.intercept = intercept;
        this.start_idx = start_idx;
        this.end_idx = end_idx;
    }

    public void setChainDetails(double chain_score_1, int chain_rank_1, double chain_score_2, int chain_rank_2) {
        this.chain_score_1 = chain_score_1;
        this.chain_rank_1 = chain_rank_1;
        this.chain_score_2 = chain_score_2;
        this.chain_rank_2 = chain_rank_2;
    }

    public double getScore() {
        return score;
    }

    public double getSecondScore() {
        return second_score;
    }

    public String getChain1() {
        return chain_seq_1;
    }

    public String getChain2() {
        return chain_seq_2;
    }

    public int getLinkSite1() {
        return link_site_1;
    }

    public int getLinkSite2() {
        return link_site_2;
    }

    public int[] getScoreHistogram() {
        return score_histogram;
    }

    public static double getInverseHistogramBinSize() {
        return inverseHistogramBinSize;
    }

    public double getEValue() {
        return e_value;
    }

    public int getScoreCount() {
        return score_count;
    }

    public float getRSquare() {
        return r_square;
    }

    public float getSlope() {
        return slope;
    }

    public float getIntercept() {
        return intercept;
    }

    public int getStartIdx() {
        return start_idx;
    }

    public int getEndIdx() {
        return end_idx;
    }

    public double getChainScore1() {
        return chain_score_1;
    }

    public int getChainRank1() {
        return chain_rank_1;
    }

    public double getChainScore2() {
        return chain_score_2;
    }

    public int getChainRank2() {
        return chain_rank_2;
    }

    public long getCandidateNum() {
        return candidate_num;
    }

    public TreeMap<Integer, ChainResultEntry> getBinChainMap() {
        return binChainMap;
    }
}
