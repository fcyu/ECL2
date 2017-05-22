package proteomics.Types;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class SparseBooleanVector {

    private Set<Integer> sparse_vector = new HashSet<>();

    SparseBooleanVector(Set<Integer> sparse_vector) {
        this.sparse_vector = sparse_vector;
    }

    public SparseBooleanVector() {}

    public void put(int idx) {
        sparse_vector.add(idx);
    }

    double norm2square() {
        return sparse_vector.size();
    }

    public double dot(SparseVector other) { // caution: it will change the original Sparse Boolean Vector.
        double output = 0;
        sparse_vector.retainAll(other.getIdxSet());
        for (int i : sparse_vector) {
            output += other.get(i);
        }
        return output;
    }

    double dot(SparseBooleanVector other) {
        Set<Integer> intersectedKeys = new HashSet<>(sparse_vector);
        intersectedKeys.retainAll(other.sparse_vector);
        return intersectedKeys.size();
    }

    public int size() {
        return sparse_vector.size();
    }

    public boolean contains(int v) {
        return sparse_vector.contains(v);
    }

    public Set<Integer> getIdxSet() {
        return sparse_vector;
    }

    public boolean equals(Object other) {
        if (other instanceof SparseBooleanVector) {
            SparseBooleanVector temp = (SparseBooleanVector) other;
            if (temp.size() == sparse_vector.size()) {
                for (int v : sparse_vector) {
                    if (!temp.contains(v)) {
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
