package proteomics.TheoSeq;

import proteomics.Types.AA;
import proteomics.Types.SparseVector;

import java.util.*;
import java.util.regex.*;

public class MassTool {

    private static final int max_charge = 6;
    private static final Pattern mod_aa_pattern = Pattern.compile("([A-Znc])(\\[([0-9\\.\\-]+)\\])?");
    public static final double PROTON = 1.00727646688;
    public static final double C13_DIFF = 1.00335483;

    public static double H2O;
    public Map<String, Double> elementTable = new HashMap<>();

    private final Map<Character, Double> mass_table = new HashMap<>(25, 1);
    private final int missed_cleavage;
    private final String cut_site;
    private final String protect_site;
    private final double inverseMzBinSize;
    private final double one_minus_bin_offset;

    public MassTool(int missed_cleavage, Map<Character, Double> fix_mod_map, String cut_site, String protect_site, double mz_bin_size, double one_minus_bin_offset) {
        inverseMzBinSize = 1 / mz_bin_size;
        this.missed_cleavage = missed_cleavage;
        this.cut_site = cut_site;
        this.protect_site = protect_site;
        this.one_minus_bin_offset = one_minus_bin_offset;

        elementTable.put("-", 0d);
        elementTable.put("H", 1.0078246);
        elementTable.put("He", 3.01603);
        elementTable.put("Li", 6.015121);
        elementTable.put("Be", 9.012182);
        elementTable.put("B", 10.012937);
        elementTable.put("C", 12.0000000);
        elementTable.put("N", 14.0030732);
        elementTable.put("O", 15.9949141);
        elementTable.put("F", 18.9984032);
        elementTable.put("Ne", 19.992435);
        elementTable.put("Na", 22.989767);
        elementTable.put("Mg", 23.985042);
        elementTable.put("Al", 26.981539);
        elementTable.put("Si", 27.976927);
        elementTable.put("P", 30.973762);
        elementTable.put("S", 31.972070);
        elementTable.put("Cl", 34.9688531);
        elementTable.put("Ar", 35.967545);
        elementTable.put("K", 38.963707);
        elementTable.put("Ca", 39.962591);
        elementTable.put("Sc", 44.955910);
        elementTable.put("Ti", 45.952629);
        elementTable.put("V", 49.947161);
        elementTable.put("Cr", 49.946046);
        elementTable.put("Mn", 54.938047);
        elementTable.put("Fe", 53.939612);
        elementTable.put("Co", 58.933198);
        elementTable.put("Ni", 57.935346);
        elementTable.put("Cu", 62.939598);
        elementTable.put("Zn", 63.929145);
        elementTable.put("Ga", 68.925580);
        elementTable.put("Ge", 69.924250);
        elementTable.put("As", 74.921594);
        elementTable.put("Se", 73.922475);
        elementTable.put("Br", 78.918336);
        elementTable.put("Kr", 77.914);
        elementTable.put("Rb", 84.911794);
        elementTable.put("Sr", 83.913430);
        elementTable.put("Y", 88.905849);
        elementTable.put("Zr", 89.904703);
        elementTable.put("Nb", 92.906377);
        elementTable.put("Mo", 91.906808);
        elementTable.put("Tc", 98.0);
        elementTable.put("Ru", 95.907599);
        elementTable.put("Rh", 102.905500);
        elementTable.put("Pd", 101.905634);
        elementTable.put("Ag", 106.905092);
        elementTable.put("Cd", 105.906461);
        elementTable.put("In", 112.904061);
        elementTable.put("Sn", 111.904826);
        elementTable.put("Sb", 120.903821);
        elementTable.put("Te", 119.904048);
        elementTable.put("I", 126.904473);
        elementTable.put("Xe", 123.905894);
        elementTable.put("Cs", 132.905429);
        elementTable.put("Ba", 129.906282);
        elementTable.put("La", 137.90711);
        elementTable.put("Ce", 135.907140);
        elementTable.put("Pr", 140.907647);
        elementTable.put("Nd", 141.907719);
        elementTable.put("Pm", 145.0);
        elementTable.put("Sm", 143.911998);
        elementTable.put("Eu", 150.919847);
        elementTable.put("Gd", 151.919786);
        elementTable.put("Tb", 158.925342);
        elementTable.put("Dy", 155.925277);
        elementTable.put("Ho", 164.930319);
        elementTable.put("Er", 161.928775);
        elementTable.put("Tm", 168.934212);
        elementTable.put("Yb", 167.933894);
        elementTable.put("Lu", 174.940770);
        elementTable.put("Hf", 173.940044);
        elementTable.put("Ta", 179.947462);
        elementTable.put("W", 179.946701);
        elementTable.put("Re", 184.952951);
        elementTable.put("Os", 183.952488);
        elementTable.put("Ir", 190.960584);
        elementTable.put("Pt", 189.959917);
        elementTable.put("Au", 196.966543);
        elementTable.put("Hg", 195.965807);
        elementTable.put("Tl", 202.972320);
        elementTable.put("Pb", 203.973020);
        elementTable.put("Bi", 208.980374);
        elementTable.put("Po", 209.0);
        elementTable.put("At", 210.0);
        elementTable.put("Rn", 222.0);
        elementTable.put("Fr", 223.0);
        elementTable.put("Ra", 226.025);
        // elementTable.put("Ac", 227.028); // conflict with Unimod bricks
        elementTable.put("Th", 232.038054);
        elementTable.put("Pa", 231.0359);
        elementTable.put("U", 234.040946);
        elementTable.put("Np", 237.048);
        elementTable.put("Pu", 244.0);
        elementTable.put("Am", 243.0);
        elementTable.put("Cm", 247.0);
        elementTable.put("Bk", 247.0);
        elementTable.put("Cf", 251.0);
        elementTable.put("Es", 252.0);
        elementTable.put("Fm", 257.0);
        elementTable.put("Md", 258.0);
        elementTable.put("No", 259.0);
        elementTable.put("Lr", 260.0);
        elementTable.put("13C", 13.0033554);
        elementTable.put("15N", 15.0001088);
        elementTable.put("18O", 17.9991616);
        elementTable.put("2H", 2.0141021);
        elementTable.put("dHex", elementTable.get("C") * 6 + elementTable.get("O") * 4 + elementTable.get("H") * 10);
        elementTable.put("Hep", elementTable.get("C") * 7 + elementTable.get("O") * 6 + elementTable.get("H") * 12);
        elementTable.put("Hex", elementTable.get("C") * 6 + elementTable.get("O") * 5 + elementTable.get("H") * 10);
        elementTable.put("HexA", elementTable.get("C") * 6 + elementTable.get("O") * 6 + elementTable.get("H") * 8);
        elementTable.put("HexN", elementTable.get("C") * 6 + elementTable.get("O") * 4 + elementTable.get("H") * 11 + elementTable.get("N"));
        elementTable.put("HexNAc", elementTable.get("C") * 8 + elementTable.get("O") * 5 + + elementTable.get("N") + elementTable.get("H") * 13);
        elementTable.put("Kdn", elementTable.get("C") * 9 + elementTable.get("H") * 14 + elementTable.get("O") * 8);
        elementTable.put("Kdo", elementTable.get("C") * 8 + elementTable.get("H") * 12 + elementTable.get("O") * 7);
        elementTable.put("NeuAc", elementTable.get("C") * 11 + elementTable.get("H") * 17 + elementTable.get("O") * 8 + elementTable.get("N"));
        elementTable.put("NeuGc", elementTable.get("C") * 11 + elementTable.get("H") * 17 + elementTable.get("O") * 9 + elementTable.get("N"));
        elementTable.put("Pent", elementTable.get("C") * 5 + elementTable.get("O") * 4 + elementTable.get("H") * 8);
        elementTable.put("Phos", elementTable.get("O") * 3 + elementTable.get("H") + elementTable.get("P"));
        elementTable.put("Sulf", elementTable.get("S") + elementTable.get("O") * 3);
        elementTable.put("Water", elementTable.get("H") * 2 + elementTable.get("O"));
        elementTable.put("Me", elementTable.get("C") + elementTable.get("H") * 2);
        elementTable.put("Ac", elementTable.get("C") * 2 + elementTable.get("H") * 2 + elementTable.get("O")); // Caution! This is not Actinium

        mass_table.put('G', (elementTable.get("C") * 2 + elementTable.get("H") * 3 + elementTable.get("N") + elementTable.get("O") + fix_mod_map.get('G')));
        mass_table.put('A', (elementTable.get("C") * 3 + elementTable.get("H") * 5 + elementTable.get("N") + elementTable.get("O") + fix_mod_map.get('A')));
        mass_table.put('S', (elementTable.get("C") * 3 + elementTable.get("H") * 5 + elementTable.get("N") + elementTable.get("O") * 2 + fix_mod_map.get('S')));
        mass_table.put('P', (elementTable.get("C") * 5 + elementTable.get("H") * 7 + elementTable.get("N") + elementTable.get("O") + fix_mod_map.get('P')));
        mass_table.put('V', (elementTable.get("C") * 5 + elementTable.get("H") * 9 + elementTable.get("N") + elementTable.get("O") + fix_mod_map.get('V')));
        mass_table.put('T', (elementTable.get("C") * 4 + elementTable.get("H") * 7 + elementTable.get("N") + elementTable.get("O") * 2 + fix_mod_map.get('I')));
        mass_table.put('C', (elementTable.get("C") * 3 + elementTable.get("H") * 5 + elementTable.get("N") + elementTable.get("O") + elementTable.get("S") + fix_mod_map.get('C')));
        mass_table.put('I', (elementTable.get("C") * 6 + elementTable.get("H") * 11 + elementTable.get("N") + elementTable.get("O") + fix_mod_map.get('I')));
        mass_table.put('L', (elementTable.get("C") * 6 + elementTable.get("H") * 11 + elementTable.get("N") + elementTable.get("O") + fix_mod_map.get('L')));
        mass_table.put('N', (elementTable.get("C") * 4 + elementTable.get("H") * 6 + elementTable.get("N") * 2 + elementTable.get("O") * 2 + fix_mod_map.get('N')));
        mass_table.put('D', (elementTable.get("C") * 4 + elementTable.get("H") * 5 + elementTable.get("N") + elementTable.get("O") * 3 + fix_mod_map.get('D')));
        mass_table.put('Q', (elementTable.get("C") * 5 + elementTable.get("H") * 8 + elementTable.get("N") * 2 + elementTable.get("O") * 2 + fix_mod_map.get('Q')));
        mass_table.put('K', (elementTable.get("C") * 6 + elementTable.get("H") * 12 + elementTable.get("N") * 2 + elementTable.get("O") + fix_mod_map.get('K')));
        mass_table.put('E', (elementTable.get("C") * 5 + elementTable.get("H") * 7 + elementTable.get("N") + elementTable.get("O") * 3 + fix_mod_map.get('E')));
        mass_table.put('M', (elementTable.get("C") * 5 + elementTable.get("H") * 9 + elementTable.get("N") + elementTable.get("O") + elementTable.get("S") + fix_mod_map.get('M')));
        mass_table.put('H', (elementTable.get("C") * 6 + elementTable.get("H") * 7 + elementTable.get("N") * 3 + elementTable.get("O") + fix_mod_map.get('H')));
        mass_table.put('F', (elementTable.get("C") * 9 + elementTable.get("H") * 9 + elementTable.get("N") + elementTable.get("O") + fix_mod_map.get('F')));
        mass_table.put('R', (elementTable.get("C") * 6 + elementTable.get("H") * 12 + elementTable.get("N") * 4 + elementTable.get("O") + fix_mod_map.get('R')));
        mass_table.put('Y', (elementTable.get("C") * 9 + elementTable.get("H") * 9 + elementTable.get("N") + elementTable.get("O") * 2 + fix_mod_map.get('Y')));
        mass_table.put('W', (elementTable.get("C") * 11 + elementTable.get("H") * 10 + elementTable.get("N") * 2 + elementTable.get("O") + fix_mod_map.get('W')));
        mass_table.put('U', (elementTable.get("C") * 3 + elementTable.get("H") * 7 + elementTable.get("N") + elementTable.get("O") * 2 + elementTable.get("Se") + fix_mod_map.get('U')));
        mass_table.put('O', (elementTable.get("C") * 12 + elementTable.get("H") * 21 + elementTable.get("N") * 3 + elementTable.get("O") * 3 + fix_mod_map.get('O')));
        mass_table.put('n', fix_mod_map.get('n'));
        mass_table.put('c', fix_mod_map.get('c'));
        H2O = elementTable.get("H") * 2 + elementTable.get("O");
    }

    public int mzToBin(double mz) {
        return (int) (mz * inverseMzBinSize + one_minus_bin_offset);
    }

    public double calResidueMass(String seq) { // n and c are also AA.
        double total_mass = 0;
        Matcher matcher = mod_aa_pattern.matcher(seq);
        while (matcher.find()) {
            char aa = matcher.group(1).charAt(0);
            double delta_mass = 0;
            if (matcher.group(3) != null) {
                delta_mass = Double.valueOf(matcher.group(3));
            }
            total_mass += mass_table.get(aa) + delta_mass;
        }

        return total_mass;
    }

    public Set<String> buildChainSet(String pro_seq, short linker_type) {
        Map<Integer, List<int[]>> digest_range_map = digestTrypsin(pro_seq);
        Set<String> chain_seq_set = new HashSet<>();

        for (int i = 0; i <= missed_cleavage; ++i) {
            for (int[] digest_range_1 : digest_range_map.get(i)) {
                String sub_string = pro_seq.substring(digest_range_1[0], digest_range_1[1]);
                if (linker_type == 1 && sub_string.substring(0, sub_string.length() - 1).contains("K")) {
                    chain_seq_set.add("n" + sub_string + "c");
                } else if (linker_type == 2 && sub_string.substring(0, sub_string.length() - 1).contains("C")) {
                    chain_seq_set.add("n" + sub_string + "c");
                }

                if (digest_range_1[1] == pro_seq.length()) {
                    // This is the end of the protein. No digestion site, so the link-sites in any position including C-term can be linked.
                    if (linker_type == 1 && sub_string.contains("K")) {
                        chain_seq_set.add("n" + sub_string + "c");
                    } else if (linker_type == 2 && sub_string.contains("C")) {
                        chain_seq_set.add("n" + sub_string + "c");
                    }
                }
            }
            if (linker_type == 1) {
                // Add N-term peptide
                if (digest_range_map.get(i).size() > 0) {
                    int[] digest_range = digest_range_map.get(i).get(0);
                    String sub_string = pro_seq.substring(digest_range[0], digest_range[1]);
                    chain_seq_set.add("n" + sub_string + "c");
                }
            }
        }
        return chain_seq_set;
    }

    public Map<Character, Double> getmass_table() {
        return mass_table;
    }

    public double generateTheoFragmentAndCalXCorr(String seq, short linkSite, double additional_mass, int precursor_charge, SparseVector xcorrPL) {
        linkSite = (short) Math.max(1, linkSite);

        int localMaxCharge = Math.min(max_charge, Math.max(precursor_charge - 1, 1));
        double[] inverseChargeArray = new double[localMaxCharge];
        for (int charge = 1; charge <= localMaxCharge; ++charge) {
            inverseChargeArray[charge - 1] = (double) 1 / (double) charge;
        }

        AA[] aaArray = seqToAAList(seq);

        double xcorr = 0;

        // traverse the sequence to get b-ion
        double bIonMass = mass_table.get(aaArray[0].aa) + aaArray[0].delta_mass; // add N-term modification
        for (int i = 1; i < aaArray.length - 2; ++i) {
            bIonMass += mass_table.get(aaArray[i].aa) + aaArray[i].delta_mass;
            if (i < linkSite) {
                for (double inverseCharge : inverseChargeArray) {
                    xcorr += xcorrPL.get(mzToBin(bIonMass * inverseCharge + PROTON));
                }
            } else {
                for (double inverseCharge : inverseChargeArray) {
                    xcorr += xcorrPL.get(mzToBin((bIonMass + additional_mass) * inverseCharge + PROTON));
                }
            }
        }
        // calculate the last b-ion with C-term modification
        bIonMass +=  mass_table.get(aaArray[aaArray.length - 2].aa) + aaArray[aaArray.length - 2].delta_mass + mass_table.get(aaArray[aaArray.length - 1].aa) + aaArray[aaArray.length - 1].delta_mass;
        for (double inverseCharge : inverseChargeArray) {
            xcorr += xcorrPL.get(mzToBin((bIonMass + additional_mass) * inverseCharge + PROTON)); // for the fragment containing all amino acids, the additional mass is always included.
        }

        // traverse the sequence with reversed order to get y-ion
        // the whole sequence
        double yIonMass = bIonMass + H2O;
        for (double inverseCharge : inverseChargeArray) {
            xcorr += xcorrPL.get(mzToBin((yIonMass + additional_mass) * inverseCharge + PROTON)); // for the fragment containing all amino acids, the additional mass is always included.
        }
        // delete the first amino acid and N-term modification
        yIonMass -= mass_table.get(aaArray[0].aa) + aaArray[0].delta_mass + mass_table.get(aaArray[1].aa) + aaArray[1].delta_mass;
        if (1 >= linkSite) {
            for (double inverseCharge : inverseChargeArray) {
                xcorr += xcorrPL.get(mzToBin(yIonMass * inverseCharge + PROTON));
            }
        } else {
            for (double inverseCharge : inverseChargeArray) {
                xcorr += xcorrPL.get(mzToBin((yIonMass + additional_mass) * inverseCharge + PROTON));
            }
        }
        // rest of the sequence
        for (int i = 2; i < aaArray.length - 2; ++i) {
            yIonMass -= mass_table.get(aaArray[i].aa) + aaArray[i].delta_mass;
            if (i >= linkSite) { // caution: here, it is different from b-ion
                for (double inverseCharge : inverseChargeArray) {
                    xcorr += xcorrPL.get(mzToBin(yIonMass * inverseCharge + PROTON));
                }
            } else {
                for (double inverseCharge : inverseChargeArray) {
                    xcorr += xcorrPL.get(mzToBin((yIonMass + additional_mass) * inverseCharge + PROTON));
                }
            }
        }

        return xcorr * 0.005;
    }

    public static AA[] seqToAAList(String seq) {
        Matcher matcher = mod_aa_pattern.matcher(seq);
        List<AA> temp = new LinkedList<>();
        while (matcher.find()) {
            char aa = matcher.group(1).charAt(0);
            double delta_mass = 0;
            if (matcher.group(3) != null) {
                delta_mass = Double.valueOf(matcher.group(3));
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
        Map<Integer, List<int[]>> digest_range_map = new HashMap<>(5, 1);
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
