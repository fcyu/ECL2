package proteomics.Types;

public class BinaryVarMod {

    public final String site;
    public final float mod_mass;
    private final String to_string;

    public BinaryVarMod(String site, float mod_mass) {
        this.site = site;
        this.mod_mass = mod_mass;
        to_string = String.format("%.2f@%s", mod_mass, site);
    }

    public boolean equals(Object other) {
        if (other instanceof BinaryVarMod) {
            BinaryVarMod temp = (BinaryVarMod) other;
            return ((temp.site.contentEquals(site)) && (temp.mod_mass == mod_mass));
        } else {
            return false;
        }
    }

    public String toString() {
        return to_string;
    }

    public int hashCode() {
        return to_string.hashCode();
    }
}
