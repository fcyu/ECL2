package proteomics.Spectrum;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.ECL2;
import proteomics.Index.BuildIndex;
import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.IsotopeDistribution;
import static ProteomicsLibrary.Utilities.*;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.model.*;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;
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

    private int usefulSpectraNum = 0;

    public PreSpectra(JMzReader spectra_parser, double ms1Tolerance, double leftInverseMs1Tolerance, double rightInverseMs1Tolerance, int ms1ToleranceUnit, BuildIndex build_index_obj, Map<String, String> parameter_map, String ext, String sqlPath) throws Exception {
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
                if (ext.toLowerCase().contentEquals("mgf")) {
                    mgfTitle = ((Ms2Query) spectrum).getTitle();
                    scan_num = getScanNum(mgfTitle);
                } else {
                    scan_num = Integer.valueOf(spectrum.getId());
                    Scan scan = ((MzXMLFile) spectra_parser).getScanByNum((long) scan_num);
                    rt = scan.getRetentionTime().getSeconds();
                    if (parameter_map.get("C13_correction").contentEquals("1")) {
                        TreeMap<Double, Double> parentPeakList = new TreeMap<>(spectra_parser.getSpectrumById(parentId).getPeakList());
                        // We do not try to correct the precursor charge if there is one.
                        IsotopeDistribution.Entry entry = isotopeDistribution.getIsotopeCorrectionNum(precursor_mz, ms1Tolerance, ms1ToleranceUnit, precursor_charge, parentPeakList);
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

    public int getUsefulSpectraNum() {
        return usefulSpectraNum;
    }
}
