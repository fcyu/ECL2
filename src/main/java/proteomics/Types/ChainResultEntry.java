package proteomics.Types;

import proteomics.ECL2;

import java.util.LinkedList;
import java.util.List;

public class ChainResultEntry implements Comparable<ChainResultEntry>{

    private String seq;
    private String normalized_seq;
    private String second_seq;
    private String normalized_second_seq;
    private int link_site;
    private double score;
    private double second_score;
    private List<Double> score_list = new LinkedList<>();

    public ChainResultEntry() {}

    public void setSeq(String seq) {
        this.seq = seq;
        normalized_seq = seq.replaceAll("I", "L").replaceAll("K", "Q");
    }

    public void setSecondSeq(String second_seq) {
        this.second_seq = second_seq;
        normalized_second_seq = second_seq.replaceAll("I", "L").replaceAll("K", "Q");
    }

    public void addToScoreList(double score) {
        if (ECL2.cal_evalue && (score_list.size() < ECL2.chain_score_point_t)) {
            score_list.add(score);
        }
    }

    public void setLinkSite(int link_site) {
        this.link_site = link_site;
    }

    public void setScore(double score) {
        this.score = score;
    }

    public void setSecondScore(double second_score) {
        this.second_score = second_score;
    }

    public String getSeq() {
        return seq;
    }

    public String getNormalizedSeq() {
        return normalized_seq;
    }

    public String getSecondSeq() {
        return second_seq;
    }

    public String getNormalizedSecondSeq() {
        return normalized_second_seq;
    }

    public int getLinkSite() {
        return link_site;
    }

    public double getScore() {
        return score;
    }

    public double getSecondScore() {
        return second_score;
    }

    public int compareTo(ChainResultEntry other) {
        if (other.score > score) {
            return -1;
        } else if (other.score < score) {
            return 1;
        } else {
            return 0;
        }
    }

    public List<Double> getScoreList() {
        return score_list;
    }
}
