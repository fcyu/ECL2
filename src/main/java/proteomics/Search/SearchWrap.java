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

import java.util.concurrent.Callable;

public class SearchWrap implements Callable<FinalResultEntry> {

    private static final Logger logger = LoggerFactory.getLogger(SearchWrap.class);

    private final Search search_obj;
    private final SpectrumEntry spectrumEntry;
    private final BuildIndex build_index_obj;
    private final MassTool mass_tool_obj;
    private final int max_common_ion_charge;
    private final PreSpectrum preSpectrumObj;

    public SearchWrap(Search search_obj, SpectrumEntry spectrumEntry, BuildIndex build_index_obj, MassTool mass_tool_obj, int max_common_ion_charge) {
        this.search_obj = search_obj;
        this.spectrumEntry = spectrumEntry;
        this.build_index_obj = build_index_obj;
        this.mass_tool_obj = mass_tool_obj;
        this.max_common_ion_charge = max_common_ion_charge;
        preSpectrumObj = new PreSpectrum(mass_tool_obj);
    }

    @Override
    public FinalResultEntry call() {
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
                    new CalEValue(spectrumEntry.scan_num, resultEntry, xcorrPL, build_index_obj.getUniprotDecoyMassSeqMap(), mass_tool_obj, build_index_obj.linker_mass, max_common_ion_charge, search_obj.consider_two_identical_chains, e_value_precursor_mass_tol);
                    if (resultEntry.getEValue() != 9999) {
                        return search_obj.convertResultEntry(spectrumEntry.scan_num, resultEntry);
                    } else {
                        return null;
                    }
                } else {
                    return search_obj.convertResultEntry(spectrumEntry.scan_num, resultEntry);
                }
            } else {
                return null;
            }
        } else {
            return null;
        }
    }
}
