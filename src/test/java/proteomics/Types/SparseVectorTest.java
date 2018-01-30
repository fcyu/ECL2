package proteomics.Types;

import org.junit.Before;
import org.junit.Test;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import static org.junit.Assert.*;

public class SparseVectorTest {

    private SparseVector vector;

    @Before
    public void setUp() throws Exception {
        Map<Integer, Float> temp = new HashMap<>();
        temp.put(1, 0.5f);
        temp.put(3, 1.2f);
        temp.put(11, 3.2f);
        vector = new SparseVector(temp);
    }

    @Test
    public void add() throws Exception {
        vector.add(1, 2.1f);
        vector.add(2, 3f);
        assertEquals(2.6f, vector.get(1), 0);
        assertEquals(3f, vector.get(2), 0);
    }

    @Test
    public void put() throws Exception {
        vector.put(1, 2.1f);
        vector.put(2, 3f);
        assertEquals(2.1f, vector.get(1), 0);
        assertEquals(3f, vector.get(2), 0);
    }

    @Test
    public void get() throws Exception {
        assertEquals(0.5f, vector.get(1), 0);
        assertEquals(1.2f, vector.get(3), 0);
        assertEquals(3.2f, vector.get(11), 0);
    }

    @Test
    public void idxSet() throws Exception {
        Set<Integer> result = vector.getIdxSet();
        Set<Integer> ground_truth = new HashSet<>();
        ground_truth.add(1);
        ground_truth.add(3);
        ground_truth.add(11);
        assertEquals(ground_truth, result);
    }

    @Test
    public void getMaxValue() throws Exception {
        assertEquals(3.2f, vector.getMaxValue(), 0);
    }

    @Test
    public void getMinValue() throws Exception {
        assertEquals(0.5f, vector.getMinValue(), 0);
    }

    @Test
    public void norm2square() throws Exception {
        assertEquals(11.93f, vector.norm2square(), 0);
    }

    @Test
    public void dot() throws Exception {
        SparseVector other = new SparseVector();
        other.put(1, 3);
        other.put(2, 5);
        other.put(11, 2);
        other.put(30, 4);
        assertEquals(7.9, vector.dot(other), 1e-6);
    }

    @Test
    public void getVectorMap() throws Exception {
        Map<Integer, Float> ground_truth = new HashMap<>();
        ground_truth.put(1, 0.5f);
        ground_truth.put(3, 1.2f);
        ground_truth.put(11, 3.2f);
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

    @Test
    public void getNonzeroNum() throws Exception {
        assertEquals(3, vector.getNonzeroNum());
    }
}