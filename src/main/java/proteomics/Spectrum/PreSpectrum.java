package proteomics.Spectrum;

import proteomics.ECL2;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.SparseVector;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class PreSpectrum {

    private static final int XCORR_OFFSET = 75;

    private final MassTool mass_tool_obj;
    private final boolean flankingPeaks;

    public PreSpectrum(MassTool mass_tool_obj, boolean flankingPeaks) {
        this.mass_tool_obj = mass_tool_obj;
        this.flankingPeaks =flankingPeaks;
    }

    public SparseVector preSpectrum (Map<Double, Double> peaks_map, double precursor_mass, String scanId) throws IOException {
        // sqrt the intensity
        Map<Double, Double> sqrt_pl_map = new HashMap<>(peaks_map.size() + 1, 1);
        for (double mz : peaks_map.keySet()) {
            if ((peaks_map.get(mz) > 1e-6) && (mz < precursor_mass)) {
                sqrt_pl_map.put(mz, Math.sqrt(peaks_map.get(mz)));
            }
        }

        // digitize the spectrum
        double[] pl_array = digitizeSpec(sqrt_pl_map);

        // normalize the spectrum
        double[] normalizedSpectrum =  normalizeSpec(pl_array);

        if (ECL2.debug) {
            BufferedWriter writer = new BufferedWriter(new FileWriter(scanId + ".normalized.spectrum.csv"));
            writer.write("bin_idx,intensity\n");
            for (int i = 0; i < normalizedSpectrum.length; ++i) {
                writer.write(i + "," + normalizedSpectrum[i] + "\n");
            }
            writer.close();
        }

        return prepareXcorr(normalizedSpectrum);
    }

    // prepare for XCorr
    private SparseVector prepareXcorr(double[] pl_array) {
        SparseVector xcorr_pl = new SparseVector();
        int offset_range = 2 * XCORR_OFFSET + 1;
        double factor = 1 / (double) (offset_range - 1); // caution: 1/150 rather than 1/151
        double my_sum = 0;
        for (int i = 0; i < XCORR_OFFSET; ++i) {
            my_sum += pl_array[i];
        }

        double[] tempArray = new double[pl_array.length];
        for (int i = XCORR_OFFSET; i < pl_array.length + XCORR_OFFSET; ++i) {
            if (i < pl_array.length) {
                my_sum += pl_array[i];
            }
            if (i >= offset_range) {
                my_sum -= pl_array[i - offset_range];
            }
            tempArray[i - XCORR_OFFSET] = (my_sum - pl_array[i - XCORR_OFFSET]) * factor;
        }

        for (int i = 1; i < pl_array.length; ++i) {
            double temp = pl_array[i] - tempArray[i];
            if (flankingPeaks) {
                temp += (pl_array[i - 1] - tempArray[i - 1]) * 0.5;
                if (i + 1 < pl_array.length) {
                    temp += (pl_array[i + 1] - tempArray[i + 1]) * 0.5;
                }
            }
            if (Math.abs(temp) > 1e-6) {
                xcorr_pl.put(i, temp);
            }
        }

        return xcorr_pl;
    }

    private double[] normalizeSpec(double[] plArray) {
        double maxIntensity = 0;
        for (double intensity : plArray) {
            if (intensity > maxIntensity) {
                maxIntensity = intensity;
            }
        }

        double[] normalizedPlArray = new double[plArray.length];
        int windowSize = (plArray.length / 10) + 1;
        for (int i = 0; i < 10; ++i) {
            // find the max intensity in each window
            double maxWindowIntensity = 0;
            for (int j = 0; j < windowSize; ++j) {
                int idx = i * windowSize + j;
                if (idx < plArray.length) {
                    if (plArray[idx] > maxWindowIntensity) {
                        maxWindowIntensity = plArray[idx];
                    }
                }
            }

            if (maxWindowIntensity > 0) {
                double temp1 = 50 / maxWindowIntensity;
                double temp2 = 0.05 * maxIntensity; // caution: Xolik does not have this
                for (int j = 0; j < windowSize; ++j) {
                    int idx = i * windowSize + j;
                    if (idx < plArray.length) {
                        if (plArray[idx] > temp2) {
                            normalizedPlArray[idx] = plArray[idx] * temp1;
                        }
                    }
                }
            }
        }

        return normalizedPlArray;
    }

    private double[] digitizeSpec(Map<Double, Double> pl) {
        Double[] mzArray = pl.keySet().toArray(new Double[pl.size()]);
        Arrays.sort(mzArray);
        double[] digitized_pl = new double[mass_tool_obj.mzToBin(mzArray[mzArray.length - 1]) + 1];
        for (double mz : pl.keySet()) {
            int idx = mass_tool_obj.mzToBin(mz);
            digitized_pl[idx] = Math.max(pl.get(mz), digitized_pl[idx]);
        }

        return digitized_pl;
    }
}
