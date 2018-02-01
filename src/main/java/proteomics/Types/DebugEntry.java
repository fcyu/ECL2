package proteomics.Types;

public class DebugEntry implements Comparable<DebugEntry> {

    public final String chain;
    public final int link_site;
    public final double mass;
    public final double score;

    public DebugEntry(String chain, int link_site, double mass, double score) {
        this.chain = chain;
        this.link_site = link_site;
        this.mass = mass;
        this.score = score;
    }

    public int compareTo(DebugEntry other) {
        if (score > other.score) {
            return 1;
        } else if (score < other.score) {
            return -1;
        } else {
            return 0;
        }
    }
}
