package proteomics.Spectrum;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.ECL2;
import proteomics.Index.BuildIndex;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
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

    public PreSpectra(JMzReader spectra_parser, BuildIndex build_index_obj, Map<String, String> parameter_map, String ext, String sqlPath) throws MzXMLParsingException, IOException, SQLException {
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

        // prepare SQL database
        Connection sqlConnection = DriverManager.getConnection(sqlPath);
        Statement sqlStatement = sqlConnection.createStatement();
        sqlStatement.executeUpdate("DROP TABLE IF EXISTS spectraTable");
        sqlStatement.executeUpdate("CREATE TABLE spectraTable (scanNum INTEGER NOT NULL, scanId TEXT PRIMARY KEY, precursorCharge INTEGER NOT NULL, precursorMz REAL NOT NULL, precursorMass REAL NOT NULL, rt INTEGER NOT NULL, massWithoutLinker REAL NOT NULL, mgfTitle TEXT NOT NULL, isotopeCorrectionNum INTEGER NOT NULL, ms1PearsonCorrelationCoefficient REAL NOT NULL, theoMass REAL, score REAL, deltaC REAL, rank INTEGER, ppm REAL, seq1 TEXT, linkSite1 INTEGER, proId1 TEXT, seq2 TEXT, linkSite2 INTEGER, proId2 TEXT, clType TEXT, hitType INTEGER, eValue REAL, candidateNum INTEGER, pointCount INTEGER, rSquare REAL, slope REAL, intercept REAL, startIdx INTEGER, endIdx INTEGER, chainScore1 REAL, chainRank1 INTEGER, chainScore2 REAL, chainRank2 INTEGER)");
        sqlStatement.close();

        PreparedStatement sqlPrepareStatement = sqlConnection.prepareStatement("INSERT INTO spectraTable (scanNum, scanId, precursorCharge, precursorMz, precursorMass, rt, massWithoutLinker, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
        sqlConnection.setAutoCommit(false);
        int usefulSpectraNum = 0;

        Iterator<Spectrum> spectrumIterator = spectra_parser.getSpectrumIterator();
        while (spectrumIterator.hasNext()) {
            Spectrum spectrum = spectrumIterator.next();

            if (ECL2.debug && !debug_scan_num_set.contains(Integer.valueOf(spectrum.getId()))) {
                continue;
            }

            if (spectrum.getMsLevel() != 2) {
                continue;
            }

            if (spectrum.getPrecursorCharge() == null) {
                logger.warn("Scan {} doesn't have charge information. Skip.", spectrum.getId());
                continue;
            }
            int precursor_charge = spectrum.getPrecursorCharge();
            double precursor_mz = spectrum.getPrecursorMZ();
            double precursor_mass = (precursor_mz * precursor_charge - precursor_charge * 1.00727646688);

            Map<Double, Double> raw_mz_intensity_map = spectrum.getPeakList();

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
                    logger.error("Cannot get scan number from the MGF title {}. Please report your MGF title to the author.", mgfTitle);
                    System.exit(1);
                }
            } else {
                scan_num = Integer.valueOf(spectrum.getId());
                Scan scan = ((MzXMLFile) spectra_parser).getScanByNum((long) scan_num);
                rt = scan.getRetentionTime().getSeconds();
            }

            sqlPrepareStatement.setInt(1, scan_num);
            sqlPrepareStatement.setString(2, spectrum.getId());
            sqlPrepareStatement.setInt(3, precursor_charge);
            sqlPrepareStatement.setDouble(4, precursor_mz);
            sqlPrepareStatement.setDouble(5, precursor_mass);
            sqlPrepareStatement.setInt(6, rt);
            sqlPrepareStatement.setDouble(7, precursor_mass - build_index_obj.linker_mass);
            sqlPrepareStatement.setString(8, mgfTitle);
            sqlPrepareStatement.setInt(9, 0); // todo: isotopeCorrectionNum
            sqlPrepareStatement.setDouble(10, -1); // todo: ms1PearsonCorrelationCoefficient
            sqlPrepareStatement.executeUpdate();
            ++usefulSpectraNum;
        }
        sqlConnection.commit();
        sqlConnection.setAutoCommit(true);
        sqlPrepareStatement.close();
        sqlConnection.close();
        logger.info("Useful MS/MS spectra number: {}.", usefulSpectraNum);
    }
}
