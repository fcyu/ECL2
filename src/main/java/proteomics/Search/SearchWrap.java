package proteomics.Search;

import proteomics.ECL2;
import proteomics.Index.BuildIndex;
import proteomics.Spectrum.PreSpectrum;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.*;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.concurrent.Callable;

public class SearchWrap implements Callable<FinalResultEntry> {

    private final Search search_obj;
    private final SpectrumEntry spectrumEntry;
    private final BuildIndex build_index_obj;
    private final PreSpectrum preSpectrumObj;
    private final Map<String, Set<String>> seqProMap;
    private final boolean cal_evalue;
    private final double delta_c_t;

    public SearchWrap(Search search_obj, SpectrumEntry spectrumEntry, BuildIndex build_index_obj, MassTool mass_tool_obj, Map<String, Set<String>> seqProMap, boolean cal_evalue, float delta_c_t, boolean flankingPeaks) {
        this.search_obj = search_obj;
        this.spectrumEntry = spectrumEntry;
        this.build_index_obj = build_index_obj;
        preSpectrumObj = new PreSpectrum(mass_tool_obj, flankingPeaks);
        this.seqProMap = seqProMap;
        this.cal_evalue = cal_evalue;
        this.delta_c_t = delta_c_t;
    }

    @Override
    public FinalResultEntry call() throws IOException {
        SparseVector xcorrPL = preSpectrumObj.preSpectrum(spectrumEntry.originalPlMap, spectrumEntry.precursor_mass, spectrumEntry.scan_num);
        if (ECL2.debug) {
            BufferedWriter writer = new BufferedWriter(new FileWriter(spectrumEntry.scan_num + ".xcorr.spectrum.csv"));
            writer.write("bin_idx,intensity\n");
            for (int idx : xcorrPL.getIdxSet()) {
                writer.write(idx + "," + xcorrPL.get(idx) + "\n");
            }
            writer.close();
        }
        TreeMap<Integer, List<Double>> binScoresMap = new TreeMap<>();
        ResultEntry resultEntry =  search_obj.doSearch(spectrumEntry, xcorrPL, binScoresMap);
        if (resultEntry != null) {
            if (1 - (resultEntry.getSecondScore() / resultEntry.getScore()) >= delta_c_t) {
                if (cal_evalue) {
                    double originalTolerance;
                    if (search_obj.ms1_tolerance_unit == 1) {
                        originalTolerance = resultEntry.spectrum_mass * search_obj.ms1_tolerance * 1e-6f;
                    } else {
                        originalTolerance = search_obj.ms1_tolerance;
                    }
                    CalEValue.calEValue(spectrumEntry.scan_num, resultEntry, build_index_obj, binScoresMap, build_index_obj.linker_mass, originalTolerance, xcorrPL, search_obj.single_chain_t);
                    if (resultEntry.getEValue() != 9999) {
                        return search_obj.convertResultEntry(spectrumEntry.scan_num, resultEntry, seqProMap, spectrumEntry);
                    } else {
                        return null;
                    }
                } else {
                    return search_obj.convertResultEntry(spectrumEntry.scan_num, resultEntry, seqProMap, spectrumEntry);
                }
            } else {
                return null;
            }
        } else {
            return null;
        }
    }
}
