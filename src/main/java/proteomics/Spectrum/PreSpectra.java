package proteomics.Spectrum;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.ECL2;
import proteomics.Index.BuildIndex;
import proteomics.Types.SpectrumEntry;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.model.*;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;
import uk.ac.ebi.pride.tools.mzxml_parser.mzxml.model.Scan;
import uk.ac.ebi.pride.tools.mgf_parser.model.Ms2Query;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class PreSpectra {

    private static final Logger logger = LoggerFactory.getLogger(PreSpectra.class);
    private static final Pattern pattern = Pattern.compile("^[0-9]+$");
    private static final Pattern scanNumPattern1 = Pattern.compile("Scan:([0-9]+) ", Pattern.CASE_INSENSITIVE);
    private static final Pattern scanNumPattern2 = Pattern.compile("scan=([0-9]+)", Pattern.CASE_INSENSITIVE);
    private static final Pattern scanNumPattern3 = Pattern.compile("^[^.]+\\.([0-9]+)\\.[0-9]+\\.[0-9]");

    private Map<Integer, SpectrumEntry> num_spectrum_map = new HashMap<>();
    private Set<Integer> debug_scan_num_set = new HashSet<>();

    public PreSpectra(JMzReader spectra_parser, BuildIndex build_index_obj, Map<String, String> parameter_map, String ext) {
        int min_ms1_charge = Integer.valueOf(parameter_map.get("min_ms1_charge"));
        int max_ms1_charge = Integer.valueOf(parameter_map.get("max_ms1_charge"));
        float min_precursor_mass =  Float.valueOf(parameter_map.get("min_precursor_mass"));
        float max_precursor_mass = Float.valueOf(parameter_map.get("max_precursor_mass"));

        //  In DEBUG mode, filter out unlisted scan num
        if (ECL2.debug) {
            for (String k : parameter_map.keySet()) {
                Matcher matcher = pattern.matcher(k);
                if (matcher.find()) {
                    debug_scan_num_set.add(Integer.valueOf(k));
                }
            }
        }

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

                double precursor_mz = spectrum.getPrecursorMZ();
                float precursor_mass = (float) (precursor_mz * precursor_charge - precursor_charge * 1.00727646688);

                if ((precursor_mass > max_precursor_mass) || (precursor_mass < min_precursor_mass)) {
                    continue;
                }

                Map<Double, Double> raw_mz_intensity_map = spectrum.getPeakList();

                int peakCount = 0;
                for (double intensity : raw_mz_intensity_map.values()) {
                    if (intensity > 1e-6) {
                        ++peakCount;
                    }
                }
                if (peakCount < 10) {
                    continue;
                }

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

                int scan_num = -1;
                String mgfTitle = "";
                float rt = -1;
                try {
                    if (ext.toLowerCase().contentEquals("mgf")) {
                        mgfTitle = ((Ms2Query) spectrum).getTitle();
                        Matcher matcher1 = scanNumPattern1.matcher(mgfTitle);
                        Matcher matcher2 = scanNumPattern2.matcher(mgfTitle);
                        Matcher matcher3 = scanNumPattern3.matcher(mgfTitle);
                        if (matcher1.find()) {
                            scan_num = Integer.valueOf(matcher1.group(1));
                        } else if (matcher2.find()) {
                            scan_num = Integer.valueOf(matcher2.group(1));
                        } else if (matcher3.find()) {
                            scan_num = Integer.valueOf(matcher3.group(1));
                        } else {
                            logger.error("Cannot get scan number from the MGF title {}. Please report your MGF title to the author.", mgfTitle);
                            System.exit(1);
                        }
                    } else {
                        scan_num = Integer.valueOf(spectrum.getId());
                        Scan scan = ((MzXMLFile) spectra_parser).getScanByNum((long) scan_num);
                        rt = scan.getRetentionTime().getSeconds();
                    }
                } catch (Exception ex) {
                    ex.printStackTrace();
                    logger.error(ex.getMessage());
                    System.exit(1);
                }

                SpectrumEntry spectrum_entry = new SpectrumEntry(scan_num, spectrum.getId(), (float) precursor_mz, precursor_mass, precursor_charge, rt, raw_mz_intensity_map, build_index_obj.linker_mass, mgfTitle);
                num_spectrum_map.put(scan_num, spectrum_entry);
            }
        } catch (Exception ex) {
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
