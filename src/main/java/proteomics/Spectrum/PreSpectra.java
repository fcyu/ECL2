package proteomics.Spectrum;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Index.BuildIndex;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.SparseVector;
import proteomics.Types.SpectrumEntry;
import uk.ac.ebi.pride.tools.jmzreader.*;
import uk.ac.ebi.pride.tools.jmzreader.model.*;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLParsingException;
import uk.ac.ebi.pride.tools.mzxml_parser.mzxml.model.Scan;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.*;

public class PreSpectra {

    private static final Logger logger = LoggerFactory.getLogger(PreSpectra.class);

    private TreeMap<Integer, Set<SpectrumEntry>> bin_spectra_map = new TreeMap<>();
    private Map<Integer, SpectrumEntry> num_spectrum_map = new HashMap<>();

    public PreSpectra(MzXMLFile spectra_parser, BuildIndex build_index_obj, Map<String, String> parameter_map, MassTool mass_tool_obj, FilterSpectra.MassScan[] scan_num_array, int idx_start, int idx_end) {
        int min_peak_num = Integer.valueOf(parameter_map.get("min_peak_num"));
        Map<String, Float> mass_table = mass_tool_obj.getMassTable();

        PreSpectrum pre_spectrum_obj = new PreSpectrum(mass_tool_obj);

        PrintStream original_stream = System.out;
        PrintStream null_stream = new PrintStream(new OutputStream() {
            @Override
            public void write(int b) throws IOException {}
        });
        System.setOut(null_stream);

        int idx = idx_start;
        try {
            while ((idx < idx_end) && (idx < scan_num_array.length)) {
                Spectrum spectrum = spectra_parser.getSpectrumById(scan_num_array[idx].scan_id);
                int precursor_charge = spectrum.getPrecursorCharge();
                float precursor_mz = spectrum.getPrecursorMZ().floatValue();
                float precursor_mass = precursor_mz * precursor_charge - precursor_charge * mass_table.get("PROTON");
                Map<Double, Double> raw_mz_intensity_map = spectrum.getPeakList();
                SparseVector pl_map = pre_spectrum_obj.preSpectrum(raw_mz_intensity_map, precursor_mass);
                if (pl_map.getNonzeroNum() <= min_peak_num) {
                    ++idx;
                    continue;
                }

                int scan_num = Integer.valueOf(scan_num_array[idx].scan_id);

                Scan scan = spectra_parser.getScanByNum((long) scan_num);
                float rt = scan.getRetentionTime().getSeconds();

                SpectrumEntry spectrum_entry = new SpectrumEntry(scan_num, spectrum.getId(), precursor_mz, precursor_mass, precursor_charge, rt, pl_map, build_index_obj.linker_mass);
                num_spectrum_map.put(scan_num, spectrum_entry);
                int bin_idx = build_index_obj.massToBin(precursor_mass);

                if (bin_spectra_map.containsKey(bin_idx)) {
                    Set<SpectrumEntry> spectrum_set = bin_spectra_map.get(bin_idx);
                    spectrum_set.add(spectrum_entry);
                } else {
                    Set<SpectrumEntry> scan_num_set = new HashSet<>();
                    scan_num_set.add(spectrum_entry);
                    bin_spectra_map.put(bin_idx, scan_num_set);
                }

                ++idx;
            }
        } catch (JMzReaderException |MzXMLParsingException ex) {
            logger.error(ex.getMessage());
            ex.printStackTrace();
            System.exit(1);
        }

        System.setOut(original_stream);
    }

    public TreeMap<Integer, Set<SpectrumEntry>> returnMassNumMap() {
        return bin_spectra_map;
    }

    public Map<Integer, SpectrumEntry> getNumSpectrumMap() {
        return num_spectrum_map;
    }
}
