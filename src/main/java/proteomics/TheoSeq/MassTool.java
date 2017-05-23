package proteomics.TheoSeq;

import proteomics.Types.AA;
import proteomics.Types.SparseBooleanVector;

import java.util.*;
import java.util.regex.*;

public class MassTool {

    private static final int max_charge = 6;
    private static final Pattern mod_aa_pattern = Pattern.compile("([A-Znc])(\\[([0-9\\.\\-]+)\\])?");

    public static final double H2O = 18.010564684;

    private final Map<Character, Float> mass_table = new HashMap<>();
    private final int missed_cleavage;
    private final String cut_site;
    private final String protect_site;
    private double mz_bin_size = 1.0005;
    private double one_minus_bin_offset = 0.4;

    public MassTool(int missed_cleavage, Map<Character, Float> fix_mod_map, String cut_site, String protect_site, float mz_bin_size, float one_minus_bin_offset) {
        this.mz_bin_size = mz_bin_size;
        this.missed_cleavage = missed_cleavage;
        this.cut_site = cut_site;
        this.protect_site = protect_site;
        this.one_minus_bin_offset = one_minus_bin_offset;
        mass_table.put('G', 57.021464f + fix_mod_map.get('G'));
        mass_table.put('A', 71.037114f + fix_mod_map.get('A'));
        mass_table.put('S', 87.032028f + fix_mod_map.get('S'));
        mass_table.put('P', 97.052764f + fix_mod_map.get('P'));
        mass_table.put('V', 99.068414f + fix_mod_map.get('V'));
        mass_table.put('T', 101.047678f + fix_mod_map.get('I'));
        mass_table.put('C', 103.009184f + fix_mod_map.get('C'));
        mass_table.put('I', 113.084064f + fix_mod_map.get('I'));
        mass_table.put('L', 113.084064f + fix_mod_map.get('L'));
        mass_table.put('N', 114.042927f + fix_mod_map.get('N'));
        mass_table.put('D', 115.026943f + fix_mod_map.get('D'));
        mass_table.put('Q', 128.058578f + fix_mod_map.get('Q'));
        mass_table.put('K', 128.094963f + fix_mod_map.get('K'));
        mass_table.put('E', 129.042593f + fix_mod_map.get('E'));
        mass_table.put('M', 131.040485f + fix_mod_map.get('M'));
        mass_table.put('H', 137.058912f + fix_mod_map.get('H'));
        mass_table.put('F', 147.068414f + fix_mod_map.get('F'));
        mass_table.put('R', 156.101111f + fix_mod_map.get('R'));
        mass_table.put('Y', 163.063329f + fix_mod_map.get('Y'));
        mass_table.put('W', 186.079313f + fix_mod_map.get('W'));
        mass_table.put('U', 150.953636f + fix_mod_map.get('U'));
        mass_table.put('O', 132.08988f + fix_mod_map.get('O'));
        mass_table.put('n', fix_mod_map.get('n'));
        mass_table.put('c', fix_mod_map.get('c'));
    }

    public int mzToBin(double mz) {
        return (int) (mz / mz_bin_size + one_minus_bin_offset);
    }

    public float calResidueMass(String seq) { // n and c are also AA.
        double total_mass = 0;
        Matcher matcher = mod_aa_pattern.matcher(seq);
        while (matcher.find()) {
            char aa = matcher.group(1).charAt(0);
            float delta_mass = 0;
            if (matcher.group(3) != null) {
                delta_mass = Float.valueOf(matcher.group(3));
            }
            total_mass += mass_table.get(aa) + delta_mass;
        }

        return (float) total_mass;
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
                if (digest_range_1[1] == pro_seq.length()) {
                    // This is the end of the protein. No digestion site, so "K" in any position including C-term can be linked.
                    if (sub_string.contains("K")) {
                        chain_seq_set.add("n" + sub_string + "c");
                    }
                }
            }
            // Add N-term peptide
            if (digest_range_map.get(i).size() > 0) {
                int[] digest_range = digest_range_map.get(i).get(0);
                String sub_string = pro_seq.substring(digest_range[0], digest_range[1]);
                chain_seq_set.add("n" + sub_string + "c");
            }
        }
        return chain_seq_set;
    }

    public Map<Character, Float> getMassTable() {
        return mass_table;
    }

    public SparseBooleanVector buildTheoVector(String seq, short linkSite, float additional_mass, int precursor_charge, int max_common_ion_charge, int specMaxBinIdx) {
        linkSite = (short) Math.max(1, linkSite);

        int localMaxCharge = Math.min(max_charge, Math.max(precursor_charge - 1, 1));
        int localMaxCommonIonCharge = Math.min(max_common_ion_charge, localMaxCharge);

        SparseBooleanVector outputVector = new SparseBooleanVector();

        Matcher matcher = mod_aa_pattern.matcher(seq);
        List<AA> temp = new LinkedList<>();
        while (matcher.find()) {
            char aa = matcher.group(1).charAt(0);
            float delta_mass = 0;
            if (matcher.group(3) != null) {
                delta_mass = Float.valueOf(matcher.group(3));
            }
            temp.add(new AA(aa, delta_mass));
        }
        AA[] aaArray = temp.toArray(new AA[temp.size()]);

        // traverse the sequence to get b-ion
        double bIonMass = aaArray[0].delta_mass; // add N-term modification
        for (int i = 1; i < aaArray.length - 2; ++i) {
            bIonMass += mass_table.get(aaArray[i].aa) + aaArray[i].delta_mass;
            if (i < linkSite) {
                for (int charge = 1; charge <= localMaxCommonIonCharge; ++charge) {
                    int idx = mzToBin(bIonMass / charge + 1.00727646688);
                    if (idx <= specMaxBinIdx) {
                        outputVector.put(idx);
                    }
                }
            } else {
                for (int charge = 1; charge <= localMaxCharge; ++charge) {
                    int idx = mzToBin((bIonMass + additional_mass) / charge + 1.00727646688);
                    if (idx <= specMaxBinIdx) {
                        outputVector.put(idx);
                    }
                }
            }
        }
        // calculate the last b-ion with C-term modification
        bIonMass +=  mass_table.get(aaArray[aaArray.length - 2].aa) + aaArray[aaArray.length - 2].delta_mass + mass_table.get(aaArray[aaArray.length - 1].aa) + aaArray[aaArray.length - 1].delta_mass;
        for (int charge = 1; charge <= localMaxCharge; ++charge) {
            int idx = mzToBin((bIonMass + additional_mass) / charge + 1.00727646688);
            if (idx <= specMaxBinIdx) {
                outputVector.put(idx); // for the fragment containing all amino acids, the additional mass is always included.
            }
        }

        // traverse the sequence with reversed order to get y-ion
        // the whole sequence
        double yIonMass = bIonMass + H2O;
        for (int charge = 1; charge <= localMaxCharge; ++charge) {
            int idx = mzToBin((yIonMass + additional_mass) / charge + 1.00727646688);
            if (idx <= specMaxBinIdx) {
                outputVector.put(idx); // for the fragment containing all amino acids, the additional mass is always included.
            }
        }
        // delete the first amino acid and N-term modification
        yIonMass -= mass_table.get(aaArray[0].aa) + aaArray[0].delta_mass + mass_table.get(aaArray[1].aa) + aaArray[1].delta_mass;
        if (1 >= linkSite) {
            for (int charge = 1; charge <= localMaxCommonIonCharge; ++charge) {
                int idx = mzToBin(yIonMass / charge + 1.00727646688);
                if (idx <= specMaxBinIdx) {
                    outputVector.put(idx);
                }
            }
        } else {
            for (int charge = 1; charge <= localMaxCharge; ++charge) {
                int idx = mzToBin((yIonMass + additional_mass) / charge + 1.00727646688);
                if (idx <= specMaxBinIdx) {
                    outputVector.put(idx);
                }
            }
        }
        // rest of the sequence
        for (int i = 2; i < aaArray.length - 2; ++i) {
            yIonMass -= mass_table.get(aaArray[i].aa) + aaArray[i].delta_mass;
            if (i >= linkSite) { // caution: here, it is different from b-ion
                for (int charge = 1; charge <= localMaxCommonIonCharge; ++charge) {
                    int idx = mzToBin(yIonMass / charge + 1.00727646688);
                    if (idx <= specMaxBinIdx) {
                        outputVector.put(idx);
                    }
                }
            } else {
                for (int charge = 1; charge <= localMaxCharge; ++charge) {
                    int idx = mzToBin((yIonMass + additional_mass) / charge + 1.00727646688);
                    if (idx <= specMaxBinIdx) {
                        outputVector.put(idx);
                    }
                }
            }
        }

        return outputVector;
    }

    public static AA[] seqToAAList(String seq) {
        Matcher matcher = mod_aa_pattern.matcher(seq);
        List<AA> temp = new LinkedList<>();
        while (matcher.find()) {
            char aa = matcher.group(1).charAt(0);
            float delta_mass = 0;
            if (matcher.group(3) != null) {
                delta_mass = Float.valueOf(matcher.group(3));
            }
            temp.add(new AA(aa, delta_mass));
        }
        return temp.toArray(new AA[temp.size()]);
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
