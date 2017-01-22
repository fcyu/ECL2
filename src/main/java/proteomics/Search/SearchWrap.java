package proteomics.Search;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.ECL2;
import proteomics.Index.BuildIndex;
import proteomics.Spectrum.PreSpectrum;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.FinalResultEntry;
import proteomics.Types.ResultEntry;
import proteomics.Types.SparseVector;
import proteomics.Types.SpectrumEntry;

import java.util.*;
import java.util.concurrent.Callable;

public class SearchWrap implements Callable<List<FinalResultEntry>> {

    private static final Logger logger = LoggerFactory.getLogger(SearchWrap.class);

    private final Search search_obj;
    private final Map<Integer, SpectrumEntry> numSpectrumMap;
    private final Integer[] scanNumArray;
    private final BuildIndex build_index_obj;
    private final MassTool mass_tool_obj;
    private final int max_common_ion_charge;
    private final PreSpectrum preSpectrumObj;

    public SearchWrap(Search search_obj, Map<Integer, SpectrumEntry> numSpectrumMap, Integer[] scanNumArray, BuildIndex build_index_obj, MassTool mass_tool_obj, int max_common_ion_charge) {
        this.search_obj = search_obj;
        this.numSpectrumMap = numSpectrumMap;
        this.scanNumArray = scanNumArray;
        this.build_index_obj = build_index_obj;
        this.mass_tool_obj = mass_tool_obj;
        this.max_common_ion_charge = max_common_ion_charge;
        preSpectrumObj = new PreSpectrum(mass_tool_obj);
    }

    @Override
    public List<FinalResultEntry> call() {
        List<FinalResultEntry> outputList = new LinkedList<>();
        for (int scanNum : scanNumArray) {
            SpectrumEntry spectrumEntry = numSpectrumMap.get(scanNum);
            SparseVector xcorrPL = preSpectrumObj.prepareXcorr(spectrumEntry.originalPlMap, spectrumEntry.precursor_mass);
            ResultEntry resultEntry =  search_obj.doSearch(spectrumEntry, xcorrPL);
            if (resultEntry != null) {
                if (1 - (resultEntry.getSecondScore() / resultEntry.getScore()) > ECL2.delta_c_t) {
                    if (ECL2.cal_evalue) {
                        float e_value_precursor_mass_tol;
                        if (search_obj.ms1_tolerance_unit == 1) {
                            e_value_precursor_mass_tol = resultEntry.spectrum_mass * search_obj.ms1_tolerance * 1e-6f;
                        } else {
                            e_value_precursor_mass_tol = search_obj.ms1_tolerance;
                        }
                        new CalEValue(scanNum, resultEntry, xcorrPL, build_index_obj.getUniprotDecoyMassSeqMap(), mass_tool_obj, build_index_obj.linker_mass, max_common_ion_charge, search_obj.consider_two_identical_chains, e_value_precursor_mass_tol);
                        if (resultEntry.getEValue() != 9999) {
                            outputList.add(search_obj.convertResultEntry(scanNum, resultEntry));
                        }
                    }
                }
            }
        }
        return outputList;
    }
}
