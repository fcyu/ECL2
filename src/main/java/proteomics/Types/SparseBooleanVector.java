package proteomics.Types;

import ProteomicsLibrary.Types.SparseVector;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;

import java.util.Set;

public class SparseBooleanVector {

    private Multiset<Integer> sparse_vector = HashMultiset.create();

    SparseBooleanVector(Multiset<Integer> sparse_vector) {
        this.sparse_vector = sparse_vector;
    }

    public SparseBooleanVector() {}

    public void put(int idx) {
        sparse_vector.add(idx);
    }

    public double dot(SparseVector other) {
        double output = 0;
        for (Integer i : sparse_vector) {
            output += other.get(i);
        }
        return output;
    }

    public Set<Integer> getIdxSet() {
        return sparse_vector.elementSet();
    }

    public boolean equals(Object other) {
        if (other instanceof SparseBooleanVector) {
            SparseBooleanVector temp = (SparseBooleanVector) other;
            if (temp.sparse_vector.size() == sparse_vector.size()) {
                for (int v : sparse_vector) {
                    if (!temp.sparse_vector.contains(v)) {
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
