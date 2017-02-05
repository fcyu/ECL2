package proteomics.Types;

public class AA {

    public final char aa;
    public final float delta_mass;

    public AA(char aa, float delta_mass) {
        this.aa = aa;
        this.delta_mass = delta_mass;
    }

    public String toString() {
        if (Math.abs(delta_mass) > 1e-6) {
            return String.format("%c[%.2f]", aa, delta_mass);
        } else {
            return String.valueOf(aa);
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
