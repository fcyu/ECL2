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


public class VarSequence {

    public final String seq;
    public final short linkSite;
    public final int binaryModType;
    private final int hashCode;

    public VarSequence(String seq, short linkSite, int binaryModType) {
        this.seq = seq;
        this.linkSite = linkSite;
        this.binaryModType = binaryModType;
        String toString = seq + "-" + linkSite + "-" + binaryModType;
        hashCode = toString.hashCode();
    }

    public boolean equals(Object other) {
        if (other instanceof VarSequence) {
            VarSequence temp = (VarSequence) other;
            return temp.seq.contentEquals(seq) && (temp.linkSite == linkSite) && (temp.binaryModType == binaryModType);
        } else {
            return false;
        }
    }

    public int hashCode() {
        return hashCode;
    }
}
