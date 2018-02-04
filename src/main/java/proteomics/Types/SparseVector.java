package proteomics.Types;

import java.util.*;

public class SparseVector {

    private Map<Integer, Double> sparse_vector = new HashMap<>();

    SparseVector(Map<Integer, Double> sparse_vector) {
        for (int i : sparse_vector.keySet()) {
            this.sparse_vector.put(i, sparse_vector.get(i));
        }
    }

    public SparseVector() {}

    public void add(int i, double v) {
        if (Math.abs(v) > 1e-6) {
            if (sparse_vector.containsKey(i)) {
                sparse_vector.put(i, sparse_vector.get(i) + v);
            } else {
                sparse_vector.put(i, v);
            }
        }
    }

    public void put(int i, double v) {
        if (Math.abs(v) > 1e-6) {
            sparse_vector.put(i, v);
        }
    }

    public double get(int i) {
        if (sparse_vector.containsKey(i)) {
            return sparse_vector.get(i);
        } else {
            return 0;
        }
    }

    public Set<Integer> getIdxSet() {
        return sparse_vector.keySet();
    }

    double getMaxValue() {
        List<Double> intensity_list = new ArrayList<>(sparse_vector.values());
        intensity_list.sort(Collections.reverseOrder());
        return intensity_list.get(0);
    }

    double getMinValue() {
        List<Double> intensity_list = new ArrayList<>(sparse_vector.values());
        Collections.sort(intensity_list);
        return intensity_list.get(0);
    }

    double norm2square() {
        double output = 0;
        for (double v : sparse_vector.values()) {
            output += v * v;
        }
        return output;
    }

    double dot(SparseVector other) {
        double output = 0;
        Map<Integer, Double> other_vector = other.sparse_vector;
        Set<Integer> intersectedKeys = new HashSet<>(sparse_vector.keySet());
        intersectedKeys.retainAll(other_vector.keySet());
        for (int i : intersectedKeys) {
            output += sparse_vector.get(i) * other_vector.get(i);
        }
        return output;
    }

    Map<Integer, Double> getVectorMap() {
        return sparse_vector;
    }

    Set<Integer> getNonzeroIdx() {
        return sparse_vector.keySet();
    }

    public int getNonzeroNum() {
        return sparse_vector.size();
    }

    public int getMaxIdx() {
        int maxIdx = 0;
        for (int idx : sparse_vector.keySet()) {
            if (idx > maxIdx) {
                maxIdx = idx;
            }
        }
        return(maxIdx);
    }
}
