package proteomics.Types;

import java.util.*;

public class SparseVector {

    private Map<Integer, Float> sparse_vector = new HashMap<>();

    SparseVector(Map<Integer, Float> sparse_vector) {
        for (int i : sparse_vector.keySet()) {
            this.sparse_vector.put(i, sparse_vector.get(i));
        }
    }

    public SparseVector() {}

    public void add(int i, float v) {
        if (Math.abs(v) > 1e-6) {
            if (sparse_vector.containsKey(i)) {
                sparse_vector.put(i, sparse_vector.get(i) + v);
            } else {
                sparse_vector.put(i, v);
            }
        }
    }

    public void put(int i, float v) {
        if (Math.abs(v) > 1e-6) {
            sparse_vector.put(i, v);
        }
    }

    public float get(int i) {
        if (sparse_vector.containsKey(i)) {
            return sparse_vector.get(i);
        } else {
            return 0;
        }
    }

    public Set<Integer> getIdxSet() {
        return sparse_vector.keySet();
    }

    float getMaxValue() {
        List<Float> intensity_list = new LinkedList<>(sparse_vector.values());
        intensity_list.sort(Collections.reverseOrder());
        return intensity_list.get(0);
    }

    float getMinValue() {
        List<Float> intensity_list = new LinkedList<>(sparse_vector.values());
        Collections.sort(intensity_list);
        return intensity_list.get(0);
    }

    double norm2square() {
        float output = 0;
        for (float v : sparse_vector.values()) {
            output += v * v;
        }
        return output;
    }

    double dot(SparseVector other) {
        double output = 0;
        Map<Integer, Float> other_vector = other.sparse_vector;
        Set<Integer> intersectedKeys = new HashSet<>(sparse_vector.keySet());
        intersectedKeys.retainAll(other_vector.keySet());
        for (int i : intersectedKeys) {
            output += sparse_vector.get(i) * other_vector.get(i);
        }
        return output;
    }

    Map<Integer, Float> getVectorMap() {
        return sparse_vector;
    }

    Set<Integer> getNonzeroIdx() {
        return sparse_vector.keySet();
    }

    boolean isNonzero(int i) {
        return get(i) != 0;
    }

    public int getNonzeroNum() {
        return sparse_vector.size();
    }
}
