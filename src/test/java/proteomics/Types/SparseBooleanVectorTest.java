package proteomics.Types;

import org.junit.Before;
import org.junit.Test;

import java.util.HashSet;
import java.util.Set;

import static org.junit.Assert.*;

public class SparseBooleanVectorTest {

    private SparseBooleanVector vector;

    @Before
    public void setUp() throws Exception {
        Set<Integer> temp = new HashSet<>();
        temp.add(1);
        temp.add(3);
        temp.add(23);
        vector = new SparseBooleanVector(temp);
    }

    @Test
    public void put() throws Exception {
        Set<Integer> temp = new HashSet<>();
        temp.add(1);
        temp.add(3);
        temp.add(23);
        temp.add(2);
        SparseBooleanVector ground_truth = new SparseBooleanVector(temp);

        vector.put(2);
        assertEquals(ground_truth, vector);
    }

    @Test
    public void norm2square() throws Exception {
        assertEquals(3, vector.norm2square(), 0);
    }

    @Test
    public void dot() throws Exception {
        SparseVector other = new SparseVector();
        other.put(1, 1.5f);
        other.put(2, 3);
        assertEquals(1.5f, vector.dot(other), 0);
    }

    @Test
    public void dot1() throws Exception {
        SparseBooleanVector other = new SparseBooleanVector();
        other.put(1);
        other.put(2);
        assertEquals(1, vector.dot(other), 0);
    }

    @Test
    public void size() throws Exception {
        assertEquals(3, vector.size());
    }

    @Test
    public void contains() throws Exception {
        assertTrue(vector.contains(1));
        assertFalse(vector.contains(2));
    }

    @Test
    public void equals() throws Exception {
        SparseBooleanVector other_1 = new SparseBooleanVector();
        other_1.put(1);
        other_1.put(3);
        other_1.put(23);
        assertTrue(vector.equals(other_1));

        SparseBooleanVector other_2 = new SparseBooleanVector();
        other_1.put(1);
        assertFalse(vector.equals(other_2));
    }
}