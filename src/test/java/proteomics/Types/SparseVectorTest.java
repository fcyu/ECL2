package proteomics.Types;

import org.junit.Before;
import org.junit.Test;

import ProteomicsLibrary.Types.SparseVector;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import static org.junit.Assert.*;

public class SparseVectorTest {

    private SparseVector vector;

    @Before
    public void setUp() throws Exception {
        Map<Integer, Double> temp = new HashMap<>();
        temp.put(1, 0.5);
        temp.put(3, 1.2);
        temp.put(11, 3.2);
        vector = new SparseVector(temp);
    }

    @Test
    public void add() throws Exception {
        vector.add(1, 2.1);
        vector.add(2, 3);
        assertEquals(2.6, vector.get(1), 1e-4);
        assertEquals(3, vector.get(2), 1e-4);
    }

    @Test
    public void put() throws Exception {
        vector.put(1, 2.1);
        vector.put(2, 3);
        assertEquals(2.1, vector.get(1), 1e-4);
        assertEquals(3, vector.get(2), 1e-4);
    }

    @Test
    public void get() throws Exception {
        assertEquals(0.5, vector.get(1), 1e-4);
        assertEquals(1.2, vector.get(3), 1e-4);
        assertEquals(3.2, vector.get(11), 1e-4);
    }

    @Test
    public void idxSet() throws Exception {
        Set<Integer> result = vector.getNonzeroIdx();
        Set<Integer> ground_truth = new HashSet<>();
        ground_truth.add(1);
        ground_truth.add(3);
        ground_truth.add(11);
        assertEquals(ground_truth, result);
    }

    @Test
    public void getMaxValue() throws Exception {
        assertEquals(3.2, vector.getMaxValue(), 1e-4);
    }

    @Test
    public void getMinValue() throws Exception {
        assertEquals(0.5, vector.getMinValue(), 1e-4);
    }

    @Test
    public void norm2square() throws Exception {
        assertEquals(11.93, vector.norm2square(), 1e-4);
    }

    @Test
    public void dot() throws Exception {
        SparseVector other = new SparseVector();
        other.put(1, 3);
        other.put(2, 5);
        other.put(11, 2);
        other.put(30, 4);
        assertEquals(7.9, vector.dot(other), 1e-4);
    }

    @Test
    public void getVectorMap() throws Exception {
        Map<Integer, Double> ground_truth = new HashMap<>();
        ground_truth.put(1, 0.5);
        ground_truth.put(3, 1.2);
        ground_truth.put(11, 3.2);
        assertEquals(ground_truth, vector.getVectorMap());
    }

    @Test
    public void getNonzeroIdx() throws Exception {
        Set<Integer> ground_truth = new HashSet<>();
        ground_truth.add(1);
        ground_truth.add(3);
        ground_truth.add(11);
        assertEquals(ground_truth, vector.getNonzeroIdx());
    }
}