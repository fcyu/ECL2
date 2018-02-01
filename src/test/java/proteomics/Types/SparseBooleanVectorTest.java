package proteomics.Types;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;
import org.junit.Before;
import org.junit.Test;

import java.util.HashSet;
import java.util.Set;

import static org.junit.Assert.*;

public class SparseBooleanVectorTest {

    private SparseBooleanVector vector;

    @Before
    public void setUp() throws Exception {
        Multiset<Integer> temp = HashMultiset.create();
        temp.add(1);
        temp.add(3);
        temp.add(23);
        vector = new SparseBooleanVector(temp);
    }

    @Test
    public void put() throws Exception {
        Multiset<Integer> temp = HashMultiset.create();
        temp.add(1);
        temp.add(3);
        temp.add(23);
        temp.add(2);
        SparseBooleanVector ground_truth = new SparseBooleanVector(temp);

        vector.put(2);
        assertEquals(ground_truth, vector);
    }

    @Test
    public void dot() throws Exception {
        SparseVector other = new SparseVector();
        other.put(1, 1.5);
        other.put(2, 3);
        assertEquals(1.5, vector.dot(other), 0);
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