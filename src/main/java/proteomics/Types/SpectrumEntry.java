package proteomics.Types;

import java.util.Map;

public class SpectrumEntry {
    public final int scan_num;
    public final String spectrum_id;
    public final float precursor_mz;
    public final float precursor_mass;
    public final float rt;
    public final float mass_without_linker_mass;
    public final int precursor_charge;
    public final String mgfTitle;
    public Map<Double, Double> originalPlMap;

    public SpectrumEntry(int scan_num, String spectrum_id, float precursor_mz, float precursor_mass, int precursor_charge, float rt, Map<Double, Double> originalPlMap, float linker_mass, String mgfTitle) {
        this.scan_num = scan_num;
        this.spectrum_id = spectrum_id;
        this.precursor_mz = precursor_mz;
        this.precursor_mass = precursor_mass;
        mass_without_linker_mass = precursor_mass - linker_mass;
        this.precursor_charge = precursor_charge;
        this.rt = rt;
        this.mgfTitle = mgfTitle;
        this.originalPlMap = originalPlMap;
    }

    public int hashCode() {
        return scan_num;
    }

    public boolean equals(Object other) {
        if (other instanceof SpectrumEntry) {
            SpectrumEntry temp = (SpectrumEntry) other;
            return (temp.scan_num == scan_num);
        } else {
            return false;
        }
    }
}