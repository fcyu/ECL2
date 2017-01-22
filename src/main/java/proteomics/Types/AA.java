package proteomics.Types;

public class AA {

    public final String aa;
    public final float delta_mass;

    public AA(String aa, float delta_mass) {
        this.aa = aa;
        this.delta_mass = delta_mass;
    }

    public String toString() {
        if (Math.abs(delta_mass) > 1e-6) {
            return String.format("%s[%.2f]", aa, delta_mass);
        } else {
            return aa;
        }
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
