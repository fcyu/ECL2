package proteomics.Types;

import java.util.Locale;

public class AA {

    public final char aa;
    public final float delta_mass;
    private final int hashCode;

    public AA(char aa, float delta_mass) {
        this.aa = aa;
        this.delta_mass = delta_mass;
        String toString;
        if (Math.abs(delta_mass) > 1e-6) {
            toString = String.format(Locale.US, "%c[%.3f]", aa, delta_mass);
        } else {
            toString = String.valueOf(aa);
        }
        hashCode = toString.hashCode();
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
