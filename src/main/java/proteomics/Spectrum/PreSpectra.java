package proteomics.Spectrum;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.ECL2;
import proteomics.Index.BuildIndex;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.SpectrumEntry;
import uk.ac.ebi.pride.tools.jmzreader.model.*;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLParsingException;
import uk.ac.ebi.pride.tools.mzxml_parser.mzxml.model.Scan;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class PreSpectra {

    private static final Logger logger = LoggerFactory.getLogger(PreSpectra.class);
    private static final Pattern pattern = Pattern.compile("^[0-9]+$");

    private Map<Integer, SpectrumEntry> num_spectrum_map = new HashMap<>();
    private Set<Integer> debug_scan_num_set = new HashSet<>();

    public PreSpectra(MzXMLFile spectra_parser, BuildIndex build_index_obj, Map<String, String> parameter_map, MassTool mass_tool_obj) {
        int min_ms1_charge = Integer.valueOf(parameter_map.get("min_ms1_charge"));
        int max_ms1_charge = Integer.valueOf(parameter_map.get("max_ms1_charge"));
        float min_precursor_mass =  Float.valueOf(parameter_map.get("min_precursor_mass"));
        float max_precursor_mass = Float.valueOf(parameter_map.get("max_precursor_mass"));
        int min_peak_num = Integer.valueOf(parameter_map.get("min_peak_num"));

        //  In DEBUG mode, filter out unlisted scan num
        if (ECL2.debug) {
            for (String k : parameter_map.keySet()) {
                Matcher matcher = pattern.matcher(k);
                if (matcher.find()) {
                    debug_scan_num_set.add(Integer.valueOf(k));
                }
            }
        }

        PreSpectrum pre_spectrum_obj = new PreSpectrum(mass_tool_obj);

        PrintStream original_stream = System.out;
        PrintStream null_stream = new PrintStream(new OutputStream() {
            @Override
            public void write(int b) throws IOException {}
        });
        System.setOut(null_stream);

        Iterator<Spectrum> spectrumIterator = spectra_parser.getSpectrumIterator();
        try {
            while (spectrumIterator.hasNext()) {
                Spectrum spectrum = spectrumIterator.next();

                if (ECL2.debug && !debug_scan_num_set.contains(Integer.valueOf(spectrum.getId()))) {
                    continue;
                }

                if (spectrum.getMsLevel() != 2) {
                    continue;
                }

                if (spectrum.getPrecursorCharge() == null) {
                    logger.warn("Scan {} doesn't have charge information. Skip.", spectrum.getId());
                    continue;
                }
                int precursor_charge = spectrum.getPrecursorCharge();

                if ((precursor_charge < min_ms1_charge) || (precursor_charge > max_ms1_charge)) {
                    continue;
                }

                float precursor_mz = spectrum.getPrecursorMZ().floatValue();
                float precursor_mass = precursor_mz * precursor_charge - precursor_charge * 1.00727646688f;

                if ((precursor_mass > max_precursor_mass) || (precursor_mass < min_precursor_mass)) {
                    continue;
                }

                Map<Double, Double> raw_mz_intensity_map = spectrum.getPeakList();

                if (ECL2.debug) {
                    try (BufferedWriter writer = new BufferedWriter(new FileWriter(Integer.valueOf(spectrum.getId()) + ".raw.spectrum.csv"))) {
                        writer.write("mz,intensity\n");
                        for (double mz : raw_mz_intensity_map.keySet()) {
                            if (Math.abs(raw_mz_intensity_map.get(mz)) > 1e-6) {
                                writer.write(mz + "," + raw_mz_intensity_map.get(mz) + "\n");
                            }
                        }
                    } catch (IOException ex) {
                        ex.printStackTrace();
                        logger.error(ex.getMessage());
                        System.exit(1);
                    }
                }

                if (raw_mz_intensity_map.size() < min_peak_num) {
                    logger.debug("Scan {} doesn't contain enough peak number ({}). Skip.", spectrum.getId(), min_peak_num);
                    continue;
                }

                TreeMap<Float, Float> originalPlMap = pre_spectrum_obj.preSpectrum(raw_mz_intensity_map, precursor_mass);

                if (ECL2.debug) {
                    try (BufferedWriter writer = new BufferedWriter(new FileWriter(Integer.valueOf(spectrum.getId()) + ".normalized.spectrum.csv"))) {
                        writer.write("mz,intensity\n");
                        for (float mz : originalPlMap.keySet()) {
                            writer.write(mz + "," + originalPlMap.get(mz) + "\n");
                        }
                    } catch (IOException ex) {
                        ex.printStackTrace();
                        logger.error(ex.getMessage());
                        System.exit(1);
                    }
                }

                if (originalPlMap.size() <= min_peak_num) {
                    continue;
                }

                int scan_num = Integer.valueOf(spectrum.getId());

                Scan scan = spectra_parser.getScanByNum((long) scan_num);
                float rt = scan.getRetentionTime().getSeconds();

                SpectrumEntry spectrum_entry = new SpectrumEntry(scan_num, spectrum.getId(), precursor_mz, precursor_mass, precursor_charge, rt, originalPlMap, build_index_obj.linker_mass);
                num_spectrum_map.put(scan_num, spectrum_entry);
            }
        } catch (MzXMLParsingException ex) {
            logger.error(ex.getMessage());
            ex.printStackTrace();
            System.exit(1);
        }

        System.setOut(original_stream);
    }

    public Map<Integer, SpectrumEntry> getNumSpectrumMap() {
        return num_spectrum_map;
    }
}
