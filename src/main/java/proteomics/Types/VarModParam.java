package proteomics.Types;


public class VarModParam {

    public final float modMass;
    public final char aa;
    private final String toString;
    private final int hashCode;

    public VarModParam(float modMass, char aa) {
        this.modMass = modMass;
        this.aa = aa;
        toString = modMass + "@" + aa;
        hashCode = toString.hashCode();
    }

    public String toString() {
        return toString;
    }

    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        if (other instanceof VarModParam) {
            VarModParam temp = (VarModParam) other;
            return (Math.abs(temp.modMass - modMass) <= 0.01) && (temp.aa == aa);
        } else {
            return false;
        }
    }
}
