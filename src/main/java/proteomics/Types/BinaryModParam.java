/*
 * Copyright 2016-2019 The Hong Kong University of Science and Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package proteomics.Types;


public class BinaryModParam {
    public final double modMass;
    public final String aas;

    private final int hashCode;

    public BinaryModParam(double modMass, String aas) {
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
