package proteomics.Spectrum;

import proteomics.TheoSeq.MassTool;
import proteomics.Types.SparseVector;

import java.util.*;

public class PreSpectrum {

    private static final float DEFAULT_INTENSITY = 1; // DO NOT change. Otherwise, change the whole project accordingly.
    private static final int XCORR_OFFSET = 75;

    private final MassTool mass_tool_obj;

    public PreSpectrum(MassTool mass_tool_obj) {
        this.mass_tool_obj = mass_tool_obj;
    }

    TreeMap<Float, Float> preSpectrum (Map<Double, Double> peaks_map, float precursor_mass) {
        return normalizeSpec(peaks_map, precursor_mass);
    }

    // prepare for XCorr
    public SparseVector prepareXcorr(TreeMap<Float, Float> pl_map, float precursor_mass) {
        float[] pl_array = digitizeSpec(pl_map, precursor_mass);

        SparseVector xcorr_pl = new SparseVector();
        float my_sum = 0;
        int offset_range = 2 * XCORR_OFFSET + 1;
        for (int i = 0; i < XCORR_OFFSET; ++i) {
            my_sum += pl_array[i];
        }

        float factor = 1 / (float) (offset_range - 1); // caution: 1/150 rather than 1/151
        for (int i = XCORR_OFFSET; i < pl_array.length + XCORR_OFFSET; ++i) {
            if (i < pl_array.length) {
                my_sum += pl_array[i];
            }
            if (i >= offset_range) {
                my_sum -= pl_array[i - offset_range];
            }
            float temp = (pl_array[i - XCORR_OFFSET] - (my_sum - pl_array[i - XCORR_OFFSET]) * factor);
            if (Math.abs(temp) > 1e-6) {
                xcorr_pl.put(i - XCORR_OFFSET, temp);
            }
        }

        return xcorr_pl;
    }

    TreeMap<Float, Float> normalizeSpec(Map<Double, Double> pl_map, float precursor_mass) {
        // sqrt the intensity and find the highest intensity.
        TreeMap<Float, Float> sqrt_pl_map = new TreeMap<>();
        for (double mz : pl_map.keySet()) {
            if ((pl_map.get(mz) > 1e-6) && (mz < precursor_mass)) {
                float sqrt_intensity = (float) Math.sqrt(pl_map.get(mz));
                sqrt_pl_map.put((float) mz, sqrt_intensity);
            }
        }

        // divide the spectrum into 10 windows and normalize each windows to DEFAULT_INTENSITY
        TreeMap<Float, Float> windowed_pl_map = new TreeMap<>();
        float min_mz = sqrt_pl_map.firstKey();
        float max_mz = sqrt_pl_map.lastKey();
        float window_size = (max_mz - min_mz) / 10 + 1;
        for (int i = 0; i < 10; ++i) {
            // find the max intensity in each window
            float left_mz = Math.min(min_mz + i * window_size, max_mz);
            float right_mz = Math.min(left_mz + window_size, max_mz);
            NavigableMap<Float, Float> sub_map;
            if (right_mz < max_mz) {
                sub_map = sqrt_pl_map.subMap(left_mz, true, right_mz, false);
            } else {
                sub_map = sqrt_pl_map.subMap(left_mz, true, right_mz, true);
            }
            if (!sub_map.isEmpty()) {
                Float[] intensity_array = sub_map.values().toArray(new Float[sub_map.size()]);
                Arrays.sort(intensity_array);
                float temp_1 = DEFAULT_INTENSITY / intensity_array[intensity_array.length - 1];
                float temp_2 = (float) 0.05 * intensity_array[intensity_array.length - 1];
                for (float mz : sub_map.keySet()) {
                    if (sub_map.get(mz) > temp_2) {
                        windowed_pl_map.put(mz, sub_map.get(mz) * temp_1);
                    }
                }
            }
        }

        return windowed_pl_map;
    }

    float[] digitizeSpec(TreeMap<Float, Float> pl, float precursor_mass) {
        float[] digitized_pl = new float[mass_tool_obj.mzToBin(precursor_mass) + 1];
        for (float mz : pl.keySet()) {
            int idx = mass_tool_obj.mzToBin(mz);
            digitized_pl[idx] = Math.max(pl.get(mz), digitized_pl[idx]);
        }

        return digitized_pl;
    }
}
