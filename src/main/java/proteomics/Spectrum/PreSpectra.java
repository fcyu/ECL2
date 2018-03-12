package proteomics.Spectrum;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.ECL2;
import proteomics.Index.BuildIndex;
import ProteomicsLibrary.MassTool;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.JMzReaderException;
import uk.ac.ebi.pride.tools.jmzreader.model.*;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLParsingException;
import uk.ac.ebi.pride.tools.mzxml_parser.mzxml.model.Scan;
import uk.ac.ebi.pride.tools.mgf_parser.model.Ms2Query;

import java.io.*;
import java.sql.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class PreSpectra {

    private static final Logger logger = LoggerFactory.getLogger(PreSpectra.class);
    private static final Pattern pattern = Pattern.compile("^[0-9]+$");
    private static final Pattern scanNumPattern1 = Pattern.compile("Scan:([0-9]+) ", Pattern.CASE_INSENSITIVE);
    private static final Pattern scanNumPattern2 = Pattern.compile("scan=([0-9]+)", Pattern.CASE_INSENSITIVE);
    private static final Pattern scanNumPattern3 = Pattern.compile("^[^.]+\\.([0-9]+)\\.[0-9]+\\.[0-9]");

    private Map<Integer, TreeMap<Integer, TreeSet<DevEntry>>> scanDevEntryMap = new HashMap<>();
    private int usefulSpectraNum = 0;

    public PreSpectra(JMzReader spectra_parser, double ms1Tolerance, double leftInverseMs1Tolerance, double rightInverseMs1Tolerance, int ms1ToleranceUnit, BuildIndex build_index_obj, Map<String, String> parameter_map, String ext, String sqlPath) throws MzXMLParsingException, IOException, SQLException, JMzReaderException {
        Set<Integer> debug_scan_num_set = new HashSet<>();
        //  In DEBUG mode, filter out unlisted scan num
        if (ECL2.debug) {
            for (String k : parameter_map.keySet()) {
                Matcher matcher = pattern.matcher(k);
                if (matcher.find()) {
                    debug_scan_num_set.add(Integer.valueOf(k));
                }
            }
        }

        IsotopeDistribution isotopeDistribution = new IsotopeDistribution(build_index_obj.returnMassTool().getElementTable(), 0, "N14");

        // prepare SQL database
        Connection sqlConnection = DriverManager.getConnection(sqlPath);
        Statement sqlStatement = sqlConnection.createStatement();
        sqlStatement.executeUpdate("PRAGMA journal_mode=WAL");
        sqlStatement.executeUpdate("DROP TABLE IF EXISTS spectraTable");
        sqlStatement.executeUpdate("CREATE TABLE spectraTable (scanNum INTEGER NOT NULL, scanId TEXT PRIMARY KEY, precursorCharge INTEGER NOT NULL, precursorMz REAL NOT NULL, precursorMass REAL NOT NULL, rt INTEGER NOT NULL, massWithoutLinker REAL NOT NULL, mgfTitle TEXT NOT NULL, isotopeCorrectionNum INTEGER NOT NULL, ms1PearsonCorrelationCoefficient REAL NOT NULL, theoMass REAL, score REAL, deltaC REAL, rank INTEGER, ppm REAL, seq1 TEXT, linkSite1 INTEGER, proId1 TEXT, seq2 TEXT, linkSite2 INTEGER, proId2 TEXT, clType TEXT, hitType INTEGER, eValue REAL, candidateNum INTEGER, pointCount INTEGER, rSquare REAL, slope REAL, intercept REAL, startIdx INTEGER, endIdx INTEGER, chainScore1 REAL, chainRank1 INTEGER, chainScore2 REAL, chainRank2 INTEGER)");
        sqlStatement.close();

        PreparedStatement sqlPrepareStatement = sqlConnection.prepareStatement("INSERT INTO spectraTable (scanNum, scanId, precursorCharge, precursorMz, precursorMass, rt, massWithoutLinker, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
        sqlConnection.setAutoCommit(false);

        Iterator<Spectrum> spectrumIterator = spectra_parser.getSpectrumIterator();
        String parentId = null;
        while (spectrumIterator.hasNext()) {
            try {
                Spectrum spectrum = spectrumIterator.next();

                if (ECL2.debug && !debug_scan_num_set.contains(Integer.valueOf(spectrum.getId()))) {
                    continue;
                }

                if (ext.toLowerCase().contentEquals("mzxml")) {
                    if (spectrum.getMsLevel() == 1) {
                        parentId = spectrum.getId();
                        continue;
                    }
                }

                if (spectrum.getPrecursorCharge() == null) {
                    logger.warn("Scan {} doesn't have charge information. Skip.", spectrum.getId());
                    continue;
                }
                int precursor_charge = spectrum.getPrecursorCharge();
                double precursor_mz = spectrum.getPrecursorMZ();
                double precursor_mass = (precursor_mz * precursor_charge - precursor_charge * MassTool.PROTON);

                Map<Double, Double> raw_mz_intensity_map = spectrum.getPeakList();

                if (raw_mz_intensity_map.size() < 5) {
                    continue;
                }

                if (ECL2.debug) {
                    BufferedWriter writer = new BufferedWriter(new FileWriter(Integer.valueOf(spectrum.getId()) + ".raw.spectrum.csv"));
                    writer.write("mz,intensity\n");
                    for (double mz : raw_mz_intensity_map.keySet()) {
                        if (Math.abs(raw_mz_intensity_map.get(mz)) > 1e-6) {
                            writer.write(mz + "," + raw_mz_intensity_map.get(mz) + "\n");
                        }
                    }
                    writer.close();
                }

                int scan_num = -1;
                String mgfTitle = "";
                int rt = -1;
                int isotopeCorrectionNum = 0;
                double pearsonCorrelationCoefficient = -1;
                TreeMap<Integer, TreeSet<DevEntry>> chargeDevEntryMap = new TreeMap<>();
                if (ext.toLowerCase().contentEquals("mgf")) {
                    mgfTitle = ((Ms2Query) spectrum).getTitle();
                    Matcher matcher1 = scanNumPattern1.matcher(mgfTitle);
                    Matcher matcher2 = scanNumPattern2.matcher(mgfTitle);
                    Matcher matcher3 = scanNumPattern3.matcher(mgfTitle);
                    if (matcher1.find()) {
                        scan_num = Integer.valueOf(matcher1.group(1));
                    } else if (matcher2.find()) {
                        scan_num = Integer.valueOf(matcher2.group(1));
                    } else if (matcher3.find()) {
                        scan_num = Integer.valueOf(matcher3.group(1));
                    } else {
                        throw new NullPointerException(String.format(Locale.US, "Cannot get scan number from the MGF title %s. Please report your MGF title to the author.", mgfTitle));
                    }
                } else {
                    scan_num = Integer.valueOf(spectrum.getId());
                    Scan scan = ((MzXMLFile) spectra_parser).getScanByNum((long) scan_num);
                    rt = scan.getRetentionTime().getSeconds();
                    if (parameter_map.get("C13_correction").contentEquals("1")) {
                        TreeMap<Double, Double> parentPeakList = new TreeMap<>(spectra_parser.getSpectrumById(parentId).getPeakList());
                        // We do not try to correct the precursor charge if there is one.
                        Entry entry = getIsotopeCorrectionNum(precursor_mz, precursor_charge, 1 / (double) precursor_charge, parentPeakList, ms1Tolerance, leftInverseMs1Tolerance, rightInverseMs1Tolerance, ms1ToleranceUnit, isotopeDistribution, chargeDevEntryMap);
                        if (entry.pearsonCorrelationCoefficient >= 0.7) { // If the Pearson correlation coefficient is smaller than 0.7, there is not enough evidence to change the original precursor mz.
                            isotopeCorrectionNum = entry.isotopeCorrectionNum;
                            pearsonCorrelationCoefficient = entry.pearsonCorrelationCoefficient;
                        }
                        precursor_mass += isotopeCorrectionNum * MassTool.C13_DIFF;
                    }
                }

                sqlPrepareStatement.setInt(1, scan_num);
                sqlPrepareStatement.setString(2, spectrum.getId());
                sqlPrepareStatement.setInt(3, precursor_charge);
                sqlPrepareStatement.setDouble(4, precursor_mz);
                sqlPrepareStatement.setDouble(5, precursor_mass);
                sqlPrepareStatement.setInt(6, rt);
                sqlPrepareStatement.setDouble(7, precursor_mass - build_index_obj.linker_mass);
                sqlPrepareStatement.setString(8, mgfTitle);
                sqlPrepareStatement.setInt(9, isotopeCorrectionNum);
                sqlPrepareStatement.setDouble(10, pearsonCorrelationCoefficient);
                sqlPrepareStatement.executeUpdate();
                ++usefulSpectraNum;

                if (ECL2.dev) {
                    scanDevEntryMap.put(scan_num, chargeDevEntryMap);
                }
            } catch (RuntimeException ex) {
                logger.error(ex.toString());
            }
        }
        sqlConnection.commit();
        sqlConnection.setAutoCommit(true);
        sqlPrepareStatement.close();
        sqlConnection.close();
        logger.info("Useful MS/MS spectra number: {}.", usefulSpectraNum);
    }

    public Map<Integer, TreeMap<Integer, TreeSet<DevEntry>>> getScanDevEntryMap() {
        return scanDevEntryMap;
    }

    public int getUsefulSpectraNum() {
        return usefulSpectraNum;
    }

    private Entry getIsotopeCorrectionNum(double precursorMz, int charge, double inverseCharge, TreeMap<Double, Double> parentPeakList, double ms1Tolerance, double leftInverseMs1Tolerance, double rightInverseMs1Tolerance, int ms1ToleranceUnit, IsotopeDistribution isotopeDistribution, TreeMap<Integer, TreeSet<DevEntry>> chargeDevEntryMap) {
        Entry entry = new Entry(0, 0);
        double leftTol = ms1Tolerance * 2;
        double rightTol = ms1Tolerance * 2;
        if (ms1ToleranceUnit == 1) {
            leftTol = (precursorMz - precursorMz * leftInverseMs1Tolerance) * 2;
            rightTol = (precursorMz * rightInverseMs1Tolerance - precursorMz) * 2;
        }
        for (int isotopeCorrectionNum : ECL2.isotopeCorrectionArray) {
            double[][] expMatrix = new double[3][2];
            for (int i = 0; i < 3; ++i) {
                expMatrix[i][0] = precursorMz + (isotopeCorrectionNum + i) * MassTool.C13_DIFF * inverseCharge;
                NavigableMap<Double, Double> subMap = parentPeakList.subMap(expMatrix[i][0] - leftTol, true, expMatrix[i][0] + rightTol, true);
                for (double intensity : subMap.values()) {
                    expMatrix[i][1] = Math.max(expMatrix[i][1], intensity);
                }
            }

            if (Math.abs(expMatrix[0][1]) > 1) { // In bottom-up proteomics, the precursor mass won't be so large that we cannot observe the monoisotopic peak.
                Map<String, Integer> elementMap = isotopeDistribution.getElementMapFromMonoMass((expMatrix[0][0] - MassTool.PROTON) * charge);
                List<IsotopeDistribution.Peak> theoIsotopeDistribution = isotopeDistribution.calculate(elementMap);
                double pearsonCorrelationCoefficient = scaleAndCalPearsonCorrelationCoefficient(expMatrix, theoIsotopeDistribution, charge, isotopeCorrectionNum);
                if (isotopeCorrectionNum == 0 && pearsonCorrelationCoefficient > 0.8) { // Unless there is a strong evidence, we prefer the original MZ.
                    entry.pearsonCorrelationCoefficient = pearsonCorrelationCoefficient;
                    entry.isotopeCorrectionNum = isotopeCorrectionNum;
                    break;
                } else if (pearsonCorrelationCoefficient > entry.pearsonCorrelationCoefficient) {
                    entry.pearsonCorrelationCoefficient = pearsonCorrelationCoefficient;
                    entry.isotopeCorrectionNum = isotopeCorrectionNum;
                }
                if (ECL2.dev) {
                    double[][] theoMatrix = new double[expMatrix.length][2];
                    for (int i = 0; i < expMatrix.length; ++i) {
                        IsotopeDistribution.Peak peak = theoIsotopeDistribution.get(i);
                        theoMatrix[i][0] = peak.mass * inverseCharge + MassTool.PROTON;
                        theoMatrix[i][1] = peak.realArea;
                    }
                    if (chargeDevEntryMap.containsKey(charge)) {
                        chargeDevEntryMap.get(charge).add(new DevEntry(isotopeCorrectionNum, pearsonCorrelationCoefficient, expMatrix, theoMatrix));
                    } else {
                        TreeSet<DevEntry> tempSet = new TreeSet<>();
                        tempSet.add(new DevEntry(isotopeCorrectionNum, pearsonCorrelationCoefficient, expMatrix, theoMatrix));
                        chargeDevEntryMap.put(charge, tempSet);
                    }
                }
            }
        }
        return entry;
    }

    private double scaleAndCalPearsonCorrelationCoefficient(double[][] expMatrix, List<IsotopeDistribution.Peak> theoIsotopeDistribution, int precursorCharge, int isotopeCorrection) {
        // get theo peaks.
        int peakNum = Math.min(expMatrix.length, theoIsotopeDistribution.size());
        double[][] theoMatrix = new double[peakNum][2];
        for (int i = 0; i < peakNum; ++i) {
            IsotopeDistribution.Peak peak = theoIsotopeDistribution.get(i);
            theoMatrix[i][0] = peak.mass / precursorCharge + MassTool.PROTON;
            theoMatrix[i][1] = peak.realArea;
        }

        // scale theo peaks
        double scale = expMatrix[-1 * isotopeCorrection][1] / theoMatrix[-1 * isotopeCorrection][1];
        if (Math.abs(scale) > 1e-6) {
            for (int i = 0; i < peakNum; ++i) {
                theoMatrix[i][1] *= scale;
            }
            // calculate Pearson correlation coefficient.
            double theoIntensityMean = 0;
            double expIntensityMean = 0;
            for (int i = 0; i < peakNum; ++i) {
                theoIntensityMean += theoMatrix[i][1];
                expIntensityMean += expMatrix[i][1];
            }
            theoIntensityMean /= peakNum;
            expIntensityMean /= peakNum;
            double temp1 = 0;
            double temp2 = 0;
            double temp3 = 0;
            for (int i = 0; i < peakNum; ++i) {
                temp1 += (theoMatrix[i][1] - theoIntensityMean) * (expMatrix[i][1] - expIntensityMean);
                temp2 += Math.pow(theoMatrix[i][1] - theoIntensityMean, 2);
                temp3 += Math.pow(expMatrix[i][1] - expIntensityMean, 2);
            }
            return (temp1 == 0 || temp2 == 0) ? 0 : temp1 / (Math.sqrt(temp2 * temp3));
        } else {
            return 0;
        }
    }


    private class Entry {

        double pearsonCorrelationCoefficient;
        int isotopeCorrectionNum;

        Entry(double pearsonCorrelationCoefficient, int isotopeCorrectionNum) {
            this.pearsonCorrelationCoefficient = pearsonCorrelationCoefficient;
            this.isotopeCorrectionNum = isotopeCorrectionNum;
        }
    }


    public static class DevEntry implements Comparable<DevEntry> {

        public final int isotopeCorrectionNum;
        public final double pearsonCorrelationCoefficient;
        public final double[][] expMatrix;
        public final double[][] theoMatrix;
        private final String toString;
        private final int hashCode;

        public DevEntry(int isotopeCorrectionNum, double pearsonCorrelationCoefficient, double[][] expMatrix, double[][] theoMatrix) {
            this.isotopeCorrectionNum = isotopeCorrectionNum;
            this.pearsonCorrelationCoefficient = pearsonCorrelationCoefficient;
            this.expMatrix = expMatrix;
            this.theoMatrix = theoMatrix;
            toString = isotopeCorrectionNum + "-" + pearsonCorrelationCoefficient;
            hashCode = toString.hashCode();
        }

        public String toString() {
            return toString;
        }

        public int hashCode() {
            return hashCode;
        }

        public boolean equals(Object other) {
            if (other instanceof DevEntry) {
                DevEntry temp = (DevEntry) other;
                return isotopeCorrectionNum == temp.isotopeCorrectionNum && pearsonCorrelationCoefficient == temp.pearsonCorrelationCoefficient;
            } else {
                return false;
            }
        }

        public int compareTo(DevEntry other) {
            if (isotopeCorrectionNum > other.isotopeCorrectionNum) {
                return 1;
            } else if (isotopeCorrectionNum < other.isotopeCorrectionNum) {
                return -1;
            } else {
                return 0;
            }
        }
    }
}
