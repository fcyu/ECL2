package proteomics.Types;


public class BinaryModParam {
    public final float modMass;
    public final String aas;

    private final int hashCode;

    public BinaryModParam(float modMass, String aas) {
        this.modMass = modMass;
        this.aas = aas;
        String toString = modMass + "@" + aas + "(binary)";
        hashCode = toString.hashCode();
    }

    public int hashCode() {
        return  hashCode;
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
