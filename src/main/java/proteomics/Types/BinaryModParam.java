package proteomics.Types;


public class BinaryModParam {
    public final float modMass;
    public final String aas;
    public final String toString;

    public BinaryModParam(float modMass, String aas) {
        this.modMass = modMass;
        this.aas = aas;
        toString = modMass + "@" + aas + "(binary)";
    }

    public String toString() {
        return toString;
    }

    public int hashCode() {
        return toString.hashCode();
    }

    public boolean equals(Object other) {
        if (other instanceof BinaryModParam) {
            BinaryModParam temp = (BinaryModParam) other;
            return (Math.abs(temp.modMass - modMass) <= 0.01) && (temp.aas.contentEquals(aas));
        } else {
            return false;
        }
    }
}
