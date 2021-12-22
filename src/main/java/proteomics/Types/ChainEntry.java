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

import java.util.*;

public class ChainEntry {

    public final String seq;
    public final double chain_mass;
    public final Set<Short> link_site_set;
    public final boolean n_term;
    public final boolean c_term;
    public final int binaryModType;

    public ChainEntry(String seq, double chain_mass, Set<Short> link_site_set, boolean n_term, boolean c_term, int binaryModType) {
        this.seq = seq;
        this.chain_mass = chain_mass;
        this.link_site_set = link_site_set;
        this.n_term = n_term;
        this.c_term = c_term;
        this.binaryModType = binaryModType;
    }

    @Override
    public boolean equals(Object other) {
        if (other instanceof ChainEntry) {
            ChainEntry temp = (ChainEntry) other;
            return seq.contentEquals(temp.seq);
        } else {
            return false;
        }
    }
}