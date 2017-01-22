package proteomics.Types;

import java.util.*;

public class ChainEntry {

    public final String seq;
    public final float chain_mass;
    public final Set<String> pro_id;
    public final Set<Integer> link_site_set;
    public final float[][] chain_ion_array;
    public final boolean n_term;
    public final boolean c_term;
    public final boolean is_decoy;

    public ChainEntry(String seq, float chain_mass, Set<String> pro_id, Set<Integer> link_site_set, float[][] chain_ion_array, boolean n_term, boolean c_term) {
        this.seq = seq;
        this.chain_mass = chain_mass;
        this.pro_id = pro_id;
        this.link_site_set = link_site_set;
        this.chain_ion_array = chain_ion_array;
        this.n_term = n_term;
        this.c_term = c_term;
        is_decoy = pro_id.iterator().next().startsWith("DECOY");
    }

    @Override
    public boolean equals(Object other) {
        if (other instanceof ChainEntry) {
            ChainEntry temp = (ChainEntry) other;
            if (seq.contentEquals(temp.seq) && (Math.abs(chain_mass - temp.chain_mass) <= 0.001) && pro_id.equals(temp.pro_id) && link_site_set.equals(temp.link_site_set) && (n_term == temp.n_term) && (c_term == temp.c_term) && (is_decoy == temp.is_decoy)) {
                for (int i = 0; i < chain_ion_array.length; ++i) {
                    if (!Arrays.equals(chain_ion_array[i], temp.chain_ion_array[i])) {
                        return false;
                    }
                }
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    }
}