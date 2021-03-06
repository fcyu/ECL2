package proteomics.Types;

import ProteomicsLibrary.DbTool;

public class ChainResultEntry implements Comparable<ChainResultEntry>{

    private String seq;
    private String ptmFreeSeq;
    private String second_seq;
    private String secondPtmFreeSeq;
    private int link_site;
    private int secondLinkSite;
    private double score;
    private double second_score;

    public ChainResultEntry() {}

    public void setSeq(String seq) {
        this.seq = seq;
        ptmFreeSeq = DbTool.getPtmFreePeptide(seq);
    }

    public void setSecondSeq(String second_seq) {
        this.second_seq = second_seq;
        secondPtmFreeSeq = DbTool.getPtmFreePeptide(second_seq);
    }

    public void setLinkSite(int link_site) {
        this.link_site = link_site;
    }

    public void setSecondLinkSite(int secondLinkSite) {
        this.secondLinkSite = secondLinkSite;
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

    public String getPtmFreeSeq() {
        return ptmFreeSeq;
    }

    public String getSecondSeq() {
        return second_seq;
    }

    public String getSecondPtmFreeSeq() {
        return secondPtmFreeSeq;
    }

    public int getLinkSite() {
        return link_site;
    }

    public int getSecondLinkSite() {
        return secondLinkSite;
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
}
