package proteomics.TheoSeq;

import proteomics.Types.AA;
import proteomics.Types.SparseBooleanVector;

import java.util.*;
import java.util.regex.*;

public class MassTool {

    private static final int max_charge = 6;

    private static final Pattern mod_aa_pattern = Pattern.compile("([A-Znc])(\\[([0-9\\.\\-]+)\\])?");

    private final Map<String, Float> mass_table = new HashMap<>();
    private final int missed_cleavage;
    private final String cut_site;
    private final String protect_site;
    private float mz_bin_size = 1.0005f;
    private float one_minus_bin_offset = 0.4f;

    public MassTool(int missed_cleavage, Map<String, Float> fix_mod_map, String cut_site, String protect_site, float mz_bin_size, float one_minus_bin_offset) {
        this.mz_bin_size = mz_bin_size;
        this.missed_cleavage = missed_cleavage;
        this.cut_site = cut_site;
        this.protect_site = protect_site;
        this.one_minus_bin_offset = one_minus_bin_offset;
        mass_table.put("G", 57.021464f + fix_mod_map.get("G"));
        mass_table.put("A", 71.037114f + fix_mod_map.get("A"));
        mass_table.put("S", 87.032028f + fix_mod_map.get("S"));
        mass_table.put("P", 97.052764f + fix_mod_map.get("P"));
        mass_table.put("V", 99.068414f + fix_mod_map.get("V"));
        mass_table.put("T", 101.047678f + fix_mod_map.get("I"));
        mass_table.put("C", 103.009184f + fix_mod_map.get("C"));
        mass_table.put("I", 113.084064f + fix_mod_map.get("I"));
        mass_table.put("L", 113.084064f + fix_mod_map.get("L"));
        mass_table.put("N", 114.042927f + fix_mod_map.get("N"));
        mass_table.put("D", 115.026943f + fix_mod_map.get("D"));
        mass_table.put("Q", 128.058578f + fix_mod_map.get("Q"));
        mass_table.put("K", 128.094963f + fix_mod_map.get("K"));
        mass_table.put("E", 129.042593f + fix_mod_map.get("E"));
        mass_table.put("M", 131.040485f + fix_mod_map.get("M"));
        mass_table.put("H", 137.058912f + fix_mod_map.get("H"));
        mass_table.put("F", 147.068414f + fix_mod_map.get("F"));
        mass_table.put("R", 156.101111f + fix_mod_map.get("R"));
        mass_table.put("Y", 163.063329f + fix_mod_map.get("Y"));
        mass_table.put("W", 186.079313f + fix_mod_map.get("W"));
        mass_table.put("U", 150.953636f + fix_mod_map.get("U"));
        mass_table.put("O", 237.147727f + fix_mod_map.get("O"));
        mass_table.put("n", fix_mod_map.get("n"));
        mass_table.put("c", fix_mod_map.get("c"));
        mass_table.put("C13_DIFF", 1.00335483f);
        mass_table.put("H2O", 18.010564684f);
        mass_table.put("NH3", 17.026549106f);
        mass_table.put("PROTON", 1.00727646688f);
        mass_table.put("Hatom", 1.007825032f);
        mass_table.put("Natom", 14.00307401f);
        mass_table.put("Oatom", 15.99491462f);
        mass_table.put("Patom", 30.97376151f);
        mass_table.put("Satom", 31.97207069f);
    }

    public int mzToBin(float mz) {
        return (int) Math.floor(mz / mz_bin_size + one_minus_bin_offset);
    }

    public float calResidueMass(AA[] aa_array) { // n and c are also AA.
        float total_mass = 0;
        for (AA aa : aa_array) {
            total_mass += mass_table.get(aa.aa) + aa.delta_mass;
        }

        return total_mass;
    }

    public Set<String> buildChainSet(String pro_seq) {
        Map<Integer, List<int[]>> digest_range_map = digestTrypsin(pro_seq);
        Set<String> chain_seq_set = new HashSet<>();

        for (int i = 0; i <= missed_cleavage; ++i) {
            for (int[] digest_range_1 : digest_range_map.get(i)) {
                String sub_string = pro_seq.substring(digest_range_1[0], digest_range_1[1]);
                if (sub_string.substring(0, sub_string.length() - 1).contains("K")) {
                    // If there is a K in middle, this peptide is a chain.
                    chain_seq_set.add("n" + sub_string + "c");
                }
            }
            if (digest_range_map.get(i).size() > 0) {
                int[] digest_range = digest_range_map.get(i).get(0);
                String sub_string = pro_seq.substring(digest_range[0], digest_range[1]);
                chain_seq_set.add("n" + sub_string + "c");
            }
        }
        return chain_seq_set;
    }

    public float[][] buildChainIonArray(AA[] aa_array) {
        // [NOTE] The b/y-ions charge 0
        float[][] chain_ion_array = new float[2 * max_charge][aa_array.length - 2];
        AA n = aa_array[0];
        AA c = aa_array[aa_array.length - 1];
        float b_ion_mass = mass_table.get(n.aa) + n.delta_mass;
        float y_ion_mass = calResidueMass(aa_array) + mass_table.get("H2O");

        for (int charge = 1; charge <= max_charge; ++charge) {
            float b_ion_mass_charge = (b_ion_mass + charge * mass_table.get("PROTON")) / charge;
            float y_ion_mass_charge = (y_ion_mass + charge * mass_table.get("PROTON")) / charge;
            int charge_idx = charge - 1;

            for (int idx = 1; idx < aa_array.length - 1; ++idx) {
                // y-ion
                chain_ion_array[2 * charge_idx + 1][idx - 1] = y_ion_mass_charge;

                AA aa = aa_array[idx];

                // b-ion
                if (idx == aa_array.length - 2) {
                    b_ion_mass_charge += (mass_table.get(aa.aa) + aa.delta_mass + mass_table.get(c.aa) + c.delta_mass) / charge;
                } else {
                    b_ion_mass_charge += (mass_table.get(aa.aa) + aa.delta_mass) / charge;
                }
                chain_ion_array[2 * charge_idx][idx - 1] = b_ion_mass_charge;

                // Calculate next y-ion
                if (idx == 1) {
                    y_ion_mass_charge -= (mass_table.get(aa.aa) + aa.delta_mass + mass_table.get(n.aa) + n.delta_mass) / charge;
                } else {
                    y_ion_mass_charge -= (mass_table.get(aa.aa) + aa.delta_mass) / charge;
                }
            }
        }

        return chain_ion_array;
    }

    public Map<String, Float> getMassTable() {
        return mass_table;
    }

    public SparseBooleanVector buildVector(float[][] ion_matrix, int precursor_charge) {
        SparseBooleanVector output_vector = new SparseBooleanVector();

        int max_row = Math.min(max_charge, precursor_charge - 1) * 2;
        if (precursor_charge == 1) {
            max_row = 2;
        }

        for (int i = 0; i < max_row; ++i) {
            for (int j = 0; j < ion_matrix[0].length; ++j) {
                if (ion_matrix[i][j] > 1e-6) {
                    output_vector.put(mzToBin(ion_matrix[i][j]));
                }
            }
        }

        return output_vector;
    }

    public float[][] buildPseudoCLIonArray(float[][] seq_ion, int link_site, int max_common_ion_charge, float additional_mass, int precursor_charge) {
        int charge = Math.min(precursor_charge - 1, max_charge);
        float[][] pseudo_cl_ion_array = new float[charge * 2][seq_ion[0].length];

        if (link_site > 0) {
            link_site -= 1;
        }

        for (int i = 0; i < charge; ++i) {
            // common ion
            if (i < max_common_ion_charge) {
                System.arraycopy(seq_ion[i * 2], 0, pseudo_cl_ion_array[i * 2], 0, link_site);
                System.arraycopy(seq_ion[i * 2 + 1], link_site + 1, pseudo_cl_ion_array[i * 2 + 1], link_site + 1, seq_ion[0].length - link_site - 1);
            }

            // xlink ion
            float addition_mz = additional_mass / (i + 1);
            for (int j = 0; j < seq_ion[0].length; ++j) {
                if (j < link_site) {
                    pseudo_cl_ion_array[i * 2 + 1][j] = seq_ion[i * 2 + 1][j] + addition_mz;
                } else if (j == link_site) {
                    pseudo_cl_ion_array[i * 2][j] = seq_ion[i * 2][j] + addition_mz;
                    pseudo_cl_ion_array[i * 2 + 1][j] = seq_ion[i * 2 + 1][j] + addition_mz;
                } else {
                    pseudo_cl_ion_array[i * 2][j] = seq_ion[i * 2][j] + addition_mz;
                }
            }
        }

        return pseudo_cl_ion_array;
    }

    public AA[] seqToAAList(String seq) {
        Matcher matcher = mod_aa_pattern.matcher(seq);
        List<AA> temp = new LinkedList<>();
        while (matcher.find()) {
            String aa = matcher.group(1);
            float delta_mass = 0;
            if (matcher.group(3) != null) {
                delta_mass = Float.valueOf(matcher.group(3));
            }
            temp.add(new AA(aa, delta_mass));
        }
        return temp.toArray(new AA[temp.size()]);
    }

    public String getModFreeSeq(String mod_seq) {
        Matcher matcher = mod_aa_pattern.matcher(mod_seq);
        StringBuilder sb = new StringBuilder(mod_seq.length());
        while (matcher.find()) {
            sb.append(matcher.group(1));
        }
        return sb.toString();
    }

    Map<Integer, List<int[]>> digestTrypsin(String pro_seq) {
        // Cut a protein
        List<Integer> cut_point_list = new LinkedList<>();
        int length = pro_seq.length();
        Pattern cut_pattern = Pattern.compile("[" + cut_site + "](?![" + protect_site + "])");
        int idx_start = 0;
        Matcher match_obj = cut_pattern.matcher(pro_seq);
        cut_point_list.add(0);
        while (idx_start < length) {
            if (match_obj.find()) {
                int cut_point = match_obj.end();
                cut_point_list.add(cut_point);
                idx_start = cut_point;
            } else {
                cut_point_list.add(length);
                break;
            }
        }

        Collections.sort(cut_point_list);

        // Deal with missed cleavage
        Map<Integer, List<int[]>> digest_range_map = new HashMap<>();
        for (int time = 0; time <= missed_cleavage; ++time) {
            List<int[]> temp = new LinkedList<>();
            int left_point;
            int right_point;
            for (int i = 0; i + 1 + time < cut_point_list.size(); ++i) {
                left_point = cut_point_list.get(i);
                right_point = cut_point_list.get(i + 1 + time);
                temp.add(new int[]{left_point, right_point});
            }
            digest_range_map.put(time, temp);
        }

        return digest_range_map;
    }
}
