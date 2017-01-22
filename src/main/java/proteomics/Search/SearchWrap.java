package proteomics.Search;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.ECL2;
import proteomics.Index.BuildIndex;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.FinalResultEntry;
import proteomics.Types.ResultEntry;
import proteomics.Types.SpectrumEntry;

import java.util.*;
import java.util.concurrent.Callable;

public class SearchWrap implements Callable<List<FinalResultEntry>> {

    private static final Logger logger = LoggerFactory.getLogger(SearchWrap.class);

    private final Search search_obj;
    private final NavigableMap<Integer, Set<SpectrumEntry>> bin_spectra_map;
    private final Map<Integer, SpectrumEntry> num_spectrum_map;
    private final BuildIndex build_index_obj;
    private final MassTool mass_tool_obj;
    private final int max_common_ion_charge;

    public SearchWrap(Search search_obj, NavigableMap<Integer, Set<SpectrumEntry>> bin_spectra_map, Map<Integer, SpectrumEntry> num_spectrum_map, BuildIndex build_index_obj, MassTool mass_tool_obj, int max_common_ion_charge) {
        this.search_obj = search_obj;
        this.bin_spectra_map = bin_spectra_map;
        this.num_spectrum_map = num_spectrum_map;
        this.build_index_obj = build_index_obj;
        this.mass_tool_obj = mass_tool_obj;
        this.max_common_ion_charge = max_common_ion_charge;
    }

    @Override
    public List<FinalResultEntry> call() {
        try {
            if (!bin_spectra_map.isEmpty()) {
                List<FinalResultEntry> outputList = new LinkedList<>();
                Map<Integer, ResultEntry> search_result =  search_obj.doSearch(bin_spectra_map);
                for (int scan_num : search_result.keySet()) {
                    ResultEntry result_entry = search_result.get(scan_num);
                    if (1 - (result_entry.getSecondScore() / result_entry.getScore()) > ECL2.delta_c_t) {
                        if (ECL2.cal_evalue) {
                            float e_value_precursor_mass_tol;
                            if (search_obj.ms1_tolerance_unit == 1) {
                                e_value_precursor_mass_tol = result_entry.spectrum_mass * search_obj.ms1_tolerance * 1e-6f;
                            } else {
                                e_value_precursor_mass_tol = search_obj.ms1_tolerance;
                            }
                            new CalEValue(scan_num, result_entry, num_spectrum_map.get(scan_num).pl_map_xcorr, build_index_obj.getUniprotDecoyMassSeqMap(), mass_tool_obj, build_index_obj.linker_mass, max_common_ion_charge, search_obj.consider_two_identical_chains, e_value_precursor_mass_tol);
                            if (result_entry.getEValue() != 9999) {
                                outputList.add(search_obj.convertResultEntry(scan_num, result_entry));
                            }
                        }
                    }
                }

                System.gc();
                return outputList;
            }
        } catch (Exception ex) {
            logger.error(ex.getMessage());
            ex.printStackTrace();
            System.exit(1);
        }
        return null;
    }
}
