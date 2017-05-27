package proteomics.Types;


public class VarSequence {

    public final String seq;
    public final short linkSite;
    public final int binaryModType;
    private final String toString;

    public VarSequence(String seq, short linkSite, int binaryModType) {
        this.seq = seq;
        this.linkSite = linkSite;
        this.binaryModType = binaryModType;
        toString = seq + "-" + linkSite + "-" + binaryModType;
    }

    public boolean equals(Object other) {
        if (other instanceof VarSequence) {
            VarSequence temp = (VarSequence) other;
            return temp.seq.contentEquals(seq) && (temp.linkSite == linkSite) && (temp.binaryModType == binaryModType);
        } else {
            return false;
        }
    }

    public String toString() {
        return toString;
    }

    public int hashCode() {
        return toString.hashCode();
    }
}
