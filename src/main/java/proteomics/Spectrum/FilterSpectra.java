package proteomics.Spectrum;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.ECL2;
import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class FilterSpectra {

    private static final Logger logger = LoggerFactory.getLogger(FilterSpectra.class);
    private static final Pattern pattern = Pattern.compile("^[0-9]+$");

    private final MassScan[] scan_num_array;
    private Set<Integer> debug_scan_num_set = new HashSet<>();

    public FilterSpectra(MzXMLFile spectra_reader, Map<String, String> parameter_map, Map<String, Float> mass_table) {
        int min_ms1_charge = Integer.valueOf(parameter_map.get("min_ms1_charge"));
        int max_ms1_charge = Integer.valueOf(parameter_map.get("max_ms1_charge"));
        float min_precursor_mass =  Float.valueOf(parameter_map.get("min_precursor_mass"));
        float max_precursor_mass = Float.valueOf(parameter_map.get("max_precursor_mass"));
        int min_peak_num = Integer.valueOf(parameter_map.get("min_peak_num"));
        List<MassScan> temp_mass_scan = new LinkedList<>();
        Iterator<Spectrum> spectra_iterator = spectra_reader.getSpectrumIterator();

        //  In DEBUG mode, filter out unlisted scan num
        if (ECL2.debug) {
            for (String k : parameter_map.keySet()) {
                Matcher matcher = pattern.matcher(k);
                if (matcher.find()) {
                    debug_scan_num_set.add(Integer.valueOf(k));
                }
            }
        }

        while (spectra_iterator.hasNext()) {
            Spectrum spectrum = spectra_iterator.next();

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
            float precursor_mass = precursor_mz * precursor_charge - precursor_charge * mass_table.get("PROTON");
            if ((precursor_mass > max_precursor_mass) || (precursor_mass < min_precursor_mass)) {
                continue;
            }

            if (spectrum.getPeakList().size() < min_peak_num) {
                logger.debug("Scan {} doesn't contain enough peak number ({}). Skip.", spectrum.getId(), min_peak_num);
                continue;
            }

            temp_mass_scan.add(new MassScan(precursor_mass, spectrum.getId()));
        }

        Collections.sort(temp_mass_scan);
        scan_num_array = temp_mass_scan.toArray(new MassScan[temp_mass_scan.size()]);
    }

    public MassScan[] getScanNumArray() {
        return scan_num_array;
    }

    public class MassScan implements Comparable<MassScan> {

        final float mass;
        final String scan_id;

        MassScan(float mass, String scan_id) {
            this.mass = mass;
            this.scan_id = scan_id;
        }

        public int compareTo(MassScan other) {
            if (other.mass > mass) {
                return -1;
            } else if (other.mass < mass) {
                return 1;
            } else {
                return 0;
            }
        }
    }
}