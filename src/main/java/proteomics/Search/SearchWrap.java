package proteomics.Search;

import ProteomicsLibrary.DbTool;
import proteomics.ECL2;
import proteomics.Index.BuildIndex;
import ProteomicsLibrary.PrepareSpectrum;
import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.Types.*;
import proteomics.Types.*;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.JMzReaderException;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.locks.ReentrantLock;

import static proteomics.ECL2.isotopeCorrectionArray;

public class SearchWrap implements Callable<SearchWrap.Entry> {

    private final Search search_obj;
    private final BuildIndex build_index_obj;
    private final PrepareSpectrum preSpectrumObj;
    private final boolean cal_evalue;
    private final double delta_c_t;
    private final JMzReader spectraParser;
    private final ReentrantLock lock;
    private final String scanId;
    private final int precursorCharge;
    private final double massWithoutLinker;
    private final double precursorMass;
    private final String sqlPath;
    private final boolean flankingPeaks;
    private final boolean linkSamePeptide;

    public SearchWrap(Search search_obj, BuildIndex build_index_obj, MassTool mass_tool_obj, boolean cal_evalue, double delta_c_t, boolean flankingPeaks, JMzReader spectraParser, ReentrantLock lock, String scanId, int precursorCharge, double massWithoutLinker, double precursorMass, String sqlPath, boolean linkSamePeptide) {
        this.search_obj = search_obj;
        this.build_index_obj = build_index_obj;
        preSpectrumObj = new PrepareSpectrum(mass_tool_obj);
        this.cal_evalue = cal_evalue;
        this.delta_c_t = delta_c_t;
        this.spectraParser = spectraParser;
        this.lock = lock;
        this.scanId = scanId;
        this.precursorCharge = precursorCharge;
        this.massWithoutLinker = massWithoutLinker;
        this.precursorMass = precursorMass;
        this.sqlPath = sqlPath;
        this.flankingPeaks = flankingPeaks;
        this.linkSamePeptide = linkSamePeptide;
    }

    @Override
    public Entry call() throws IOException, JMzReaderException, SQLException {
        Map<Double, Double> rawPLMap;
        try {
            lock.lock();
            // Reading peak list.
            rawPLMap = spectraParser.getSpectrumById(scanId).getPeakList();
        } finally {
            lock.unlock();
        }

        SparseVector xcorrPL = preSpectrumObj.preSpectrumCometStyle(rawPLMap, precursorMass, flankingPeaks);

        if (xcorrPL.isEmpty()) {
            return null;
        }

        if (ECL2.debug) {
            BufferedWriter writer = new BufferedWriter(new FileWriter(scanId + ".xcorr.spectrum.csv"));
            writer.write("bin_idx,intensity\n");
            for (int idx : xcorrPL.getNonzeroIdx()) {
                writer.write(idx + "," + xcorrPL.get(idx) + "\n");
            }
            writer.close();
        }
        TreeMap<Integer, List<Double>> binScoresMap = new TreeMap<>();
        ResultEntry resultEntry =  search_obj.doSearch(xcorrPL, binScoresMap, massWithoutLinker, precursorCharge, precursorMass, scanId, linkSamePeptide);
        if (resultEntry != null) {
            if (1 - (resultEntry.getSecondScore() / resultEntry.getScore()) >= delta_c_t) {
                if (cal_evalue) {
                    double originalTolerance;
                    if (search_obj.ms1_tolerance_unit == 1) {
                        originalTolerance = precursorMass * search_obj.ms1_tolerance * 1e-6;
                    } else {
                        originalTolerance = search_obj.ms1_tolerance;
                    }
                    CalEValue.calEValue(scanId, resultEntry, build_index_obj, binScoresMap, precursorCharge, massWithoutLinker, precursorMass, originalTolerance, xcorrPL, search_obj.single_chain_t);
                    if (resultEntry.getEValue() != 9999) {
                        return recordResult(resultEntry);
                    } else {
                        return null;
                    }
                } else {
                    return recordResult(resultEntry);
                }
            } else {
                return null;
            }
        } else {
            return null;
        }
    }

    private Entry recordResult(ResultEntry result_entry) throws SQLException {
        Map<String, ChainEntry> chainEntryMap = build_index_obj.getSeqEntryMap();
        Map<String, Set<String>> seqProMap = build_index_obj.getSeqProMap();

        int rank = 1;
        String chain_seq_1 = result_entry.getChain1();
        String chain_seq_2 = result_entry.getChain2();
        ChainEntry chain_entry_1 = chainEntryMap.get(chain_seq_1);
        ChainEntry chain_entry_2 = chainEntryMap.get(chain_seq_2);

        double theo_mass = chain_entry_1.chain_mass + chain_entry_2.chain_mass + build_index_obj.linker_mass;

        int C13_Diff_num = getC13Num(precursorMass, theo_mass);
        double precursorMassNew = precursorMass + C13_Diff_num * MassTool.C13_DIFF;
        double ppm = (precursorMassNew - theo_mass) * 1e6 / theo_mass;

        Set<String> pro1Set = new TreeSet<>();
        boolean isDecoy1 = false;
        for (String temp : seqProMap.get(DbTool.getPtmFreePeptide(chain_entry_1.seq))) {
            pro1Set.add(temp);
            if (temp.startsWith("DECOY")) { // there is no overlapped peptide between target and decoy.
                isDecoy1 = true;
            }
        }
        Set<String> pro2Set = new TreeSet<>();
        boolean isDecoy2 = false;
        for (String temp : seqProMap.get(DbTool.getPtmFreePeptide(chain_entry_2.seq))) {
            pro2Set.add(temp);
            if (temp.startsWith("DECOY")) { // there is no overlapped peptide between target and decoy.
                isDecoy2 = true;
            }
        }

        int hit_type; // 0 = T-T; 1 = D-D; 2 = T-D;
        if (isDecoy1 && isDecoy2) {
            hit_type = 1;
        } else if (!isDecoy1 && !isDecoy2) {
            hit_type = 0;
        } else {
            hit_type = 2;
        }

        String cl_type = "intra_protein";
        boolean keep = false;
        for (String temp_1 : seqProMap.get(DbTool.getPtmFreePeptide(chain_entry_1.seq))) {
            if (temp_1.startsWith("DECOY_")) {
                temp_1 = temp_1.substring(6);
            }
            for (String temp_2 : seqProMap.get(DbTool.getPtmFreePeptide(chain_entry_2.seq))) {
                if (temp_2.startsWith("DECOY_")) {
                    temp_2 = temp_2.substring(6);
                }
                if (temp_1.contentEquals(temp_2)) { // bias to intra protein cross-linking
                    keep = true;
                    break;
                }
            }
            if (keep) {
                break;
            }
        }

        if (!keep) {
            cl_type = "inter_protein";
        }

        double delta_c = 1 - (result_entry.getSecondScore() / result_entry.getScore());

        // add fix modification to the sequences.
        String final_seq_1 = addFixMod(chain_seq_1, result_entry.getLinkSite1());
        String final_seq_2 = addFixMod(chain_seq_2, result_entry.getLinkSite2());

        Connection localSqlConnection = DriverManager.getConnection(sqlPath);
        Statement localSqlStatement = localSqlConnection.createStatement();
        ResultSet localSqlResultSet = localSqlStatement.executeQuery(String.format(Locale.US, "SELECT scanNum, scanId, precursorCharge, precursorMz, precursorMass, rt, massWithoutLinker, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, score, hitType FROM spectraTable WHERE scanId='%s'", scanId));
        if (localSqlResultSet.next()) {
            int scanNum = localSqlResultSet.getInt("scanNum");
            int precursorCharge = localSqlResultSet.getInt("precursorCharge");
            double precursorMz = localSqlResultSet.getDouble("precursorMz");
            int rt = localSqlResultSet.getInt("rt");
            double massWithoutLinker = localSqlResultSet.getDouble("massWithoutLinker");
            String mgfTitle = localSqlResultSet.getString("mgfTitle");
            int isotopeCorrectionNum = localSqlResultSet.getInt("isotopeCorrectionNum");
            double ms1PearsonCorrelationCoefficient = localSqlResultSet.getDouble("ms1PearsonCorrelationCoefficient");
            // sqlStr = String.format(Locale.US, "REPLACE INTO spectraTable (scanNum, scanId, precursorCharge, precursorMz, precursorMass, rt, massWithoutLinker, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, theoMass, score, deltaC, rank, ppm, seq1, linkSite1, proId1, seq2, linkSite2, proId2, clType, hitType, eValue, candidateNum, pointCount, rSquare, slope, intercept, startIdx, endIdx, chainScore1, chainRank1, chainScore2, chainRank2) VALUES (%d, '%s', %d, %f, %f, %d, %f, '%s', %d, %f, %f, %f, %f, %d, %f, '%s', %d, '%s', '%s', %d, '%s', '%s', %d, %f, %d, %d, %f, %f, %f, %d, %d, %f, %d, %f, %d)", scanNum, scanId, precursorCharge, precursorMz, precursorMass, rt, massWithoutLinker, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, theo_mass, result_entry.getScore(), delta_c, rank, ppm, final_seq_1, result_entry.getLinkSite1(), String.join(";", pro1Set), final_seq_2, result_entry.getLinkSite2(), String.join(";", pro2Set), cl_type, hit_type, result_entry.getEValue(), result_entry.getCandidateNum(), result_entry.getScoreCount(), result_entry.getRSquare(), result_entry.getSlope(), result_entry.getIntercept(), result_entry.getStartIdx(), result_entry.getEndIdx(), result_entry.getChainScore1(), result_entry.getChainRank1(), result_entry.getChainScore2(), result_entry.getChainRank2());
            Entry entry = new Entry(scanNum, scanId, precursorCharge, precursorMz, precursorMass, rt, massWithoutLinker, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, theo_mass, result_entry.getScore(), delta_c, rank, ppm, final_seq_1, result_entry.getLinkSite1(), String.join(";", pro1Set), final_seq_2, result_entry.getLinkSite2(), String.join(";", pro2Set), cl_type, hit_type, result_entry.getEValue(), result_entry.getCandidateNum(), result_entry.getScoreCount(), result_entry.getRSquare(), result_entry.getSlope(), result_entry.getIntercept(), result_entry.getStartIdx(), result_entry.getEndIdx(), result_entry.getChainScore1(), result_entry.getChainRank1(), result_entry.getChainScore2(), result_entry.getChainRank2());

            localSqlResultSet.close();
            localSqlStatement.close();
            localSqlConnection.close();
            return entry;
        }  else {
            throw new NullPointerException(String.format(Locale.US, "There is no record %s in the spectraTable.", scanId));
        }
    }

    private int getC13Num(double exp_mass, double theo_mass) {
        double min_diff = 10;
        int num = 0;

        for (int i : isotopeCorrectionArray) {
            double temp = Math.abs(exp_mass + i * MassTool.C13_DIFF - theo_mass);
            if (temp < min_diff) {
                min_diff = temp;
                num = i;
            }
        }

        return num;
    }

    private String addFixMod(String seq, int linkSite) {
        Map<Character, Double> fix_mod_map = build_index_obj.getFixModMap();
        AA[] aaList = MassTool.seqToAAList(seq, "[]");
        StringBuilder sb = new StringBuilder(seq.length() * 3);
        for (int i = 0; i < aaList.length; ++i) {
            AA aa = aaList[i];
            if (i == linkSite) { // priority order: linkSite > fixMod > varMod
                sb.append(aa.aa);
            } else if (Math.abs(fix_mod_map.get(aa.aa)) > 1e-6) {
                sb.append(String.format(Locale.US, "%c[%.3f]", aa.aa, fix_mod_map.get(aa.aa)));
            } else {
                if (Math.abs(aa.ptmDeltaMass) > 1e-6) {
                    sb.append(String.format(Locale.US, "%c[%.3f]", aa.aa, aa.ptmDeltaMass));
                } else {
                    sb.append(aa.aa);
                }
            }
        }
        return sb.toString();
    }


    public class Entry {

        public final int scanNum;
        public final String scanId;
        public final int precursorCharge;
        public final double precursorMz;
        public final double precursorMass;
        public final int rt;
        public final double massWithoutLinker;
        public final String mgfTitle;
        public final int isotopeCorrectionNum;
        public final double ms1PearsonCorrelationCoefficient;
        public final double theoMass;
        public final double score;
        public final double delteC;
        public final int rank;
        public final double ppm;
        public final String seq1;
        public final int linkSite1;
        public final String pro1;
        public final String seq2;
        public final int linkSite2;
        public final String pro2;
        public final String clType;
        public final int hitType;
        public final double eValue;
        public final long candidateNum;
        public final int scoreCount;
        public final double rSquare;
        public final double slope;
        public final double intercept;
        public final int startIdx;
        public final int endIdx;
        public final double chainScore1;
        public final int chainRank1;
        public final double chainScore2;
        public final int chainRank2;

        public Entry(int scanNum, String scanId, int precursorCharge, double precursorMz, double precursorMass, int rt, double massWithoutLinker, String mgfTitle, int isotopeCorrectionNum, double ms1PearsonCorrelationCoefficient, double theoMass, double score, double delteC, int rank, double ppm, String seq1, int linkSite1, String pro1, String seq2, int linkSite2, String pro2, String clType, int hitType, double eValue, long candidateNum, int scoreCount, double rSquare, double slope, double intercept, int startIdx, int endIdx, double chainScore1, int chainRank1, double chainScore2, int chainRank2) {
            this.scanNum = scanNum;
            this.scanId = scanId;
            this.precursorCharge = precursorCharge;
            this.precursorMz = precursorMz;
            this.precursorMass = precursorMass;
            this.rt = rt;
            this.massWithoutLinker = massWithoutLinker;
            this.mgfTitle = mgfTitle;
            this.isotopeCorrectionNum = isotopeCorrectionNum;
            this.ms1PearsonCorrelationCoefficient = ms1PearsonCorrelationCoefficient;
            this.theoMass = theoMass;
            this.score = score;
            this.delteC = delteC;
            this.rank = rank;
            this.ppm = ppm;
            this.seq1 = seq1;
            this.linkSite1 = linkSite1;
            this.pro1 = pro1;
            this.seq2 = seq2;
            this.linkSite2 = linkSite2;
            this.pro2 = pro2;
            this.clType = clType;
            this.hitType = hitType;
            this.eValue = eValue;
            this.candidateNum = candidateNum;
            this.scoreCount = scoreCount;
            this.rSquare = rSquare;
            this.slope = slope;
            this.intercept = intercept;
            this.startIdx = startIdx;
            this.endIdx = endIdx;
            this.chainScore1 = chainScore1;
            this.chainRank1 = chainRank1;
            this.chainScore2 = chainScore2;
            this.chainRank2 = chainRank2;
        }
    }
}
