package proteomics.Search;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.ECL2;
import proteomics.Index.BuildIndex;
import proteomics.Spectrum.PreSpectrum;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.*;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;

public class SearchWrap implements Callable<FinalResultEntry> {

    private static final Logger logger = LoggerFactory.getLogger(SearchWrap.class);

    private final Search search_obj;
    private final SpectrumEntry spectrumEntry;
    private final BuildIndex build_index_obj;
    private final MassTool mass_tool_obj;
    private final int max_common_ion_charge;
    private final PreSpectrum preSpectrumObj;
    private final Map<String, Set<String>> seqProMap;

    public SearchWrap(Search search_obj, SpectrumEntry spectrumEntry, BuildIndex build_index_obj, MassTool mass_tool_obj, int max_common_ion_charge, Map<String, Set<String>> seqProMap) {
        this.search_obj = search_obj;
        this.spectrumEntry = spectrumEntry;
        this.build_index_obj = build_index_obj;
        this.mass_tool_obj = mass_tool_obj;
        this.max_common_ion_charge = max_common_ion_charge;
        preSpectrumObj = new PreSpectrum(mass_tool_obj);
        this.seqProMap = seqProMap;
    }

    @Override
    public FinalResultEntry call() {
        SparseVector xcorrPL = preSpectrumObj.preSpectrum(spectrumEntry.originalPlMap, spectrumEntry.precursor_mass, spectrumEntry.scan_num);
        if (ECL2.debug) {
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(spectrumEntry.scan_num + ".xcorr.spectrum.csv"))) {
                writer.write("bin_idx,intensity\n");
                for (int idx : xcorrPL.getIdxSet()) {
                    writer.write(idx + "," + xcorrPL.get(idx) + "\n");
                }
            } catch (IOException ex) {
                ex.printStackTrace();
                logger.error(ex.getMessage());
                System.exit(1);
            }
        }
        ResultEntry resultEntry =  search_obj.doSearch(spectrumEntry, xcorrPL);
        if (resultEntry != null) {
            if (ECL2.debug) {
                SparseBooleanVector chainVector1 = mass_tool_obj.buildTheoVector(resultEntry.getChain1(), (short) resultEntry.getLinkSite1(), spectrumEntry.precursor_mass - (float) (mass_tool_obj.calResidueMass(resultEntry.getChain1()) + MassTool.H2O), spectrumEntry.precursor_charge, max_common_ion_charge, xcorrPL.getMaxIdx());
                SparseBooleanVector chainVector2 = mass_tool_obj.buildTheoVector(resultEntry.getChain2(), (short) resultEntry.getLinkSite2(), spectrumEntry.precursor_mass - (float) (mass_tool_obj.calResidueMass(resultEntry.getChain2()) + MassTool.H2O), spectrumEntry.precursor_charge, max_common_ion_charge, xcorrPL.getMaxIdx());
                try (BufferedWriter writer = new BufferedWriter(new FileWriter(spectrumEntry.scan_num + ".chain.spectrum.csv"))) {
                    writer.write(resultEntry.getChain1() + " bin idx\n");
                    for (int idx : chainVector1.getIdxSet()) {
                        writer.write(idx + "\n");
                    }
                    writer.write(resultEntry.getChain2() + " bin idx\n");
                    for (int idx : chainVector2.getIdxSet()) {
                        writer.write(idx + "\n");
                    }
                } catch (IOException ex) {
                    ex.printStackTrace();
                    logger.error(ex.getMessage());
                    System.exit(1);
                }
            }
            if (1 - (resultEntry.getSecondScore() / resultEntry.getScore()) > ECL2.delta_c_t) {
                if (ECL2.cal_evalue) {
                    float e_value_precursor_mass_tol;
                    if (search_obj.ms1_tolerance_unit == 1) {
                        e_value_precursor_mass_tol = resultEntry.spectrum_mass * search_obj.ms1_tolerance * 1e-6f;
                    } else {
                        e_value_precursor_mass_tol = search_obj.ms1_tolerance;
                    }
                    new CalEValue(spectrumEntry.scan_num, resultEntry, xcorrPL, build_index_obj.getUniprotDecoyMassSeqMap(), mass_tool_obj, build_index_obj.linker_mass, max_common_ion_charge, e_value_precursor_mass_tol);
                    if (resultEntry.getEValue() != 9999) {
                        return search_obj.convertResultEntry(spectrumEntry.scan_num, resultEntry, seqProMap);
                    } else {
                        return null;
                    }
                } else {
                    return search_obj.convertResultEntry(spectrumEntry.scan_num, resultEntry, seqProMap);
                }
            } else {
                return null;
            }
        } else {
            return null;
        }
    }
}
