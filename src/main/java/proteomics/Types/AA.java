package proteomics.Types;

import java.util.Locale;

public class AA {

    public final char aa;
    public final float delta_mass;
    private final String toString;
    private final int hashCode;

    public AA(char aa, float delta_mass) {
        this.aa = aa;
        this.delta_mass = delta_mass;
        if (Math.abs(delta_mass) > 1e-6) {
            toString = String.format(Locale.US, "%c[%.3f]", aa, delta_mass);
        } else {
            toString = String.valueOf(aa);
        }
        hashCode = toString.hashCode();
    }

    public String toString() {
        return toString;
    }

    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        if (other instanceof AA) {
            AA temp = (AA) other;
            return temp.toString().contentEquals(this.toString());
        } else {
            return false;
        }
    }
}
