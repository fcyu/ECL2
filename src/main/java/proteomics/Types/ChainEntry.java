package proteomics.Types;

import java.util.*;

public class ChainEntry {

    public final String seq;
    public final double chain_mass;
    public final Set<Short> link_site_set;
    public final boolean n_term;
    public final boolean c_term;
    public final int binaryModType;

    public ChainEntry(String seq, double chain_mass, Set<Short> link_site_set, boolean n_term, boolean c_term, int binaryModType) {
        this.seq = seq;
        this.chain_mass = chain_mass;
        this.link_site_set = link_site_set;
        this.n_term = n_term;
        this.c_term = c_term;
        this.binaryModType = binaryModType;
    }

    @Override
    public boolean equals(Object other) {
        if (other instanceof ChainEntry) {
            ChainEntry temp = (ChainEntry) other;
            return seq.contentEquals(temp.seq);
        } else {
            return false;
        }
    }
}