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

    public double dot(SparseVector other) {
        double output = 0;
        Map<Integer, Float> other_vector = other.getVectorMap();
        for (int i : sparse_vector) {
            if (other_vector.containsKey(i)) {
                output += other_vector.get(i);
            }
        }
        return output;
    }

    double dot(SparseBooleanVector other) {
        double dot_product = 0;
        for (int k : other.sparse_vector) {
            if (sparse_vector.contains(k)) {
                ++dot_product;
            }
        }
        return dot_product;
    }

    public int size() {
        return sparse_vector.size();
    }

    public boolean contains(int v) {
        return sparse_vector.contains(v);
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
