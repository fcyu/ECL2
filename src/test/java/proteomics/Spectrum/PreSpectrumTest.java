package proteomics.Spectrum;

import org.junit.BeforeClass;
import org.junit.Test;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.SparseVector;

import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

public class PreSpectrumTest {

    private static Map<Character, Float> fix_mod_map = new HashMap<>();
    private static MassTool mass_tool_obj;
    private static PreSpectrum pre_spectrum_obj;

    @BeforeClass
    public static void setUp() throws Exception {
        float mz_bin_size = 0.02f;
        float one_minus_bin_offset = 1;

        fix_mod_map.put('G', 0f);
        fix_mod_map.put('A', 0f);
        fix_mod_map.put('S', 0f);
        fix_mod_map.put('P', 0f);
        fix_mod_map.put('V', 0f);
        fix_mod_map.put('T', 0f);
        fix_mod_map.put('C', 57.02146f);
        fix_mod_map.put('I', 0f);
        fix_mod_map.put('L', 0f);
        fix_mod_map.put('N', 0f);
        fix_mod_map.put('D', 0f);
        fix_mod_map.put('Q', 0f);
        fix_mod_map.put('K', 0f);
        fix_mod_map.put('E', 0f);
        fix_mod_map.put('M', 0f);
        fix_mod_map.put('H', 0f);
        fix_mod_map.put('F', 0f);
        fix_mod_map.put('R', 0f);
        fix_mod_map.put('Y', 0f);
        fix_mod_map.put('W', 0f);
        fix_mod_map.put('U', 0f);
        fix_mod_map.put('O', 0f);
        fix_mod_map.put('n', 0f);

        mass_tool_obj = new MassTool(0, fix_mod_map, "KR", "P", mz_bin_size, one_minus_bin_offset);
        pre_spectrum_obj = new PreSpectrum(mass_tool_obj);
    }

    @Test
    public void digitizeSpec() throws Exception {
        TreeMap<Float, Float> peaks_map = new TreeMap<>();
        peaks_map.put(113.0716f, 11347.8955f);
        peaks_map.put(114.0665f, 2344.7898f);
        peaks_map.put(141.0654f, 1897.8801f);

        float[] results = pre_spectrum_obj.digitizeSpec(peaks_map, 1000f);

        float[] ground_truth = new float[mass_tool_obj.mzToBin(1000f) + 1];
        ground_truth[5654] = 11347.8955f;
        ground_truth[5704] = 2344.7898f;
        ground_truth[7054] = 1897.8801f;
        assertArrayEquals(ground_truth, results, 0);
    }

    @Test
    public void normalizeSpec() {
        Map<Double, Double> peaks_map = new TreeMap<>();
        peaks_map.put(113.0716, 11347.8955);
        peaks_map.put(114.0665, 2344.7898);
        peaks_map.put(116.0709, 19463.0645);
        peaks_map.put(127.0508, 1947.5690);
        peaks_map.put(158.0924, 52436.9258);
        peaks_map.put(158.1044, 2455.5698);
        peaks_map.put(167.0381, 1459.7645);
        peaks_map.put(167.0818, 1813.0936);
        peaks_map.put(175.1191, 83178.7266);
        peaks_map.put(175.1348, 3220.7153);
        peaks_map.put(176.1223, 3157.8721);
        peaks_map.put(182.0925, 2318.5164);
        peaks_map.put(184.1080, 2771.4836);
        peaks_map.put(185.0920, 44454.5898);
        peaks_map.put(249.0988, 3125.5686);
        peaks_map.put(251.1497, 2475.2751);
        peaks_map.put(255.1467, 19911.5371);
        peaks_map.put(329.1923, 2416.5129);
        peaks_map.put(331.1594, 2136.6970);
        peaks_map.put(340.8202, 4608.2510);
        peaks_map.put(343.1721, 4315.1021);
        peaks_map.put(360.1987, 63162.6563);

        TreeMap<Float, Float> result = pre_spectrum_obj.normalizeSpec(peaks_map, 400);

        TreeMap<Float, Float> ground_truth = new TreeMap<>();
        ground_truth.put(113.0716f, 0.76357f);
        ground_truth.put(114.0665f, 0.3471f);
        ground_truth.put(116.0709f, 1.00f);
        ground_truth.put(127.0508f, 0.3163f);
        ground_truth.put(158.0924f, 1.00f);
        ground_truth.put(158.1044f, 0.2164f);
        ground_truth.put(167.0381f, 0.1325f);
        ground_truth.put(167.0818f, 0.1476f);
        ground_truth.put(175.1191f, 1.00f);
        ground_truth.put(175.1348f, 0.1968f);
        ground_truth.put(176.1223f, 0.1948f);
        ground_truth.put(182.0925f, 0.1670f);
        ground_truth.put(184.1080f, 0.1825f);
        ground_truth.put(185.0920f, 0.7311f);
        ground_truth.put(249.0988f, 0.3962f);
        ground_truth.put(251.1497f, 0.3526f);
        ground_truth.put(255.1467f, 1.00f);
        ground_truth.put(329.1923f, 0.7241f);
        ground_truth.put(331.1594f, 0.6809f);
        ground_truth.put(340.8202f, 1.00f);
        ground_truth.put(343.1721f, 0.9677f);
        ground_truth.put(360.1987f, 1.00f);

        for (float mz : result.keySet()) {
            assertEquals(ground_truth.get(mz), result.get(mz), 0.001);
        }
    }

    @Test
    public void prepareXcorr() { // todo: complete [lowest priority]
        TreeMap<Float, Float> peaks_map = new TreeMap<>();
        peaks_map.put(184.1080f, 1f);

        SparseVector result = pre_spectrum_obj.prepareXcorr(peaks_map, 200f);

        SparseVector ground_truth = new SparseVector();
    }
}