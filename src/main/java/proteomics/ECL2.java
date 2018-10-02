package proteomics;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Index.BuildIndex;
import proteomics.Parameter.Parameter;
import proteomics.Search.SearchWrap;
import proteomics.Search.Search;
import proteomics.Spectrum.PreSpectra;
import ProteomicsLibrary.MassTool;
import proteomics.Validation.CalFDR;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;

import java.io.*;
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.sql.*;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.locks.ReentrantLock;

public class ECL2 {

    public static final int score_point_t = 15000;
    public static final int[] isotopeCorrectionArray = new int[]{-2, -1, 0}; // do not change it

    private static final Logger logger = LoggerFactory.getLogger(ECL2.class);
    public static final String version = "2.1.8";

    public static boolean debug;
    public static boolean dev;

    public static void main(String[] args) {
        long startTime = System.nanoTime();

        // Process inputs
        if (args.length != 2) {
            help();
        }

        // Set parameters
        String parameter_path = args[0].trim();
        String spectra_path = args[1].trim();

        logger.info("ECL2 Version {}.", version);
        logger.info("Author: Fengchao Yu. Email: fyuab@connect.ust.hk");

        String dbName = null;
        String hostName = "unknown-host";
        try {
            hostName = InetAddress.getLocalHost().getHostName();
            logger.info("Computer: {}.", hostName);
        } catch (UnknownHostException ex) {
            logger.warn("Cannot get the computer's name.");
        }

        try {
            dbName = String.format(Locale.US, "ECL2.%s.%s.temp.db", hostName, new SimpleDateFormat("yyyy-MM-dd-HH-mm-ss").format(Calendar.getInstance().getTime()));

            logger.info("Parameter file: {}.", parameter_path);
            logger.info("Spectra file: {}.", spectra_path);

            new ECL2(parameter_path, spectra_path, dbName);
        } catch (Exception ex) {
            ex.printStackTrace();
            logger.error(ex.toString());
        } finally {
            if (dbName != null) {
                (new File(dbName)).delete();
                (new File(dbName + "-wal")).delete();
                (new File(dbName + "-shm")).delete();
            }
        }

        double totalHour = (double) (System.nanoTime() - startTime) * 1e-9 / 3600;
        logger.info("Running time: {} hours.", totalHour);
        logger.info("Done!");
    }

    private ECL2(String parameter_path, String spectra_path, String dbName) throws Exception{
        // Get the parameter map
        Parameter parameter = new Parameter(parameter_path);
        Map<String, String> parameter_map = parameter.returnParameterMap();

        // print all parameters
        logger.info("Parameters:");
        for (String k : parameter_map.keySet()) {
            logger.info("{} = {}", k, parameter_map.get(k));
        }

        boolean cal_evalue = true;
        if (parameter_map.get("cal_evalue").contentEquals("0")) {
            cal_evalue = false;
        }

        double delta_c_t = Double.valueOf(parameter_map.get("delta_c_t"));

        boolean flankingPeaks = true;
        if (parameter_map.get("flanking_peaks").contentEquals("0")) {
            flankingPeaks = false;
        }

        debug = parameter_map.get("debug").contentEquals("1");
        dev = parameter_map.get("dev").contentEquals("1");

        if (debug) {
            logger.info("In DEBUG mode.");
        }

        if (dev) {
            logger.info("In DEV mode.");
        }

        // Prepare search
        logger.info("Protein database: {}.", parameter_map.get("db"));
        logger.info("Indexing protein database...");
        BuildIndex build_index_obj = new BuildIndex(parameter_map);
        MassTool mass_tool_obj = build_index_obj.returnMassTool();

        logger.info("Peptide chains number: {}.", build_index_obj.getSeqEntryMap().size());

        String sqlPath = "jdbc:sqlite:" + dbName;
        Class.forName("org.sqlite.JDBC").newInstance();

        logger.info("Reading spectra...");
        JMzReader spectra_parser;
        File spectra_file = new File(spectra_path);
        if ((!spectra_file.exists() || (spectra_file.isDirectory()))) {
            throw new FileNotFoundException("The spectra file not found.");
        }
        String[] temp = spectra_path.split("\\.");
        String ext = temp[temp.length - 1];
        if (ext.contentEquals("mzXML")) {
            spectra_parser = new MzXMLFile(spectra_file);
        } else if (ext.toLowerCase().contentEquals("mgf")) {
            spectra_parser = new MgfFile(spectra_file);
        } else {
            throw new Exception(String.format(Locale.US, "Unsupported data format %s. ECL2 only support mzXML and MGF.", ext));
        }

        double ms1Tolerance = Double.valueOf(parameter_map.get("ms1_tolerance"));
        double leftInverseMs1Tolerance = 1 / (1 + ms1Tolerance * 1e-6);
        double rightInverseMs1Tolerance = 1 / (1 - ms1Tolerance * 1e-6);
        int ms1ToleranceUnit = Integer.valueOf(parameter_map.get("ms1_tolerance_unit"));

        PreSpectra preSpectra = new PreSpectra(spectra_parser, ms1Tolerance, leftInverseMs1Tolerance, rightInverseMs1Tolerance, ms1ToleranceUnit, build_index_obj, parameter_map, ext, sqlPath);

        logger.info("Start searching...");

        int thread_num = Integer.valueOf(parameter_map.get("thread_num"));
        if (thread_num == 0) {
            thread_num = 3 + Runtime.getRuntime().availableProcessors();
        }
        if (debug) {
            thread_num = 1;
        }
        ExecutorService thread_pool = Executors.newFixedThreadPool(thread_num);
        Search search_obj = new Search(build_index_obj, parameter_map, ms1Tolerance, leftInverseMs1Tolerance, rightInverseMs1Tolerance, ms1ToleranceUnit);
        ArrayList<Future<SearchWrap.Entry>> taskList = new ArrayList<>(preSpectra.getUsefulSpectraNum() + 10);
        Connection sqlConnection = DriverManager.getConnection(sqlPath);
        Statement sqlStatement = sqlConnection.createStatement();
        ResultSet sqlResultSet = sqlStatement.executeQuery("SELECT scanId, precursorCharge, massWithoutLinker, precursorMass FROM spectraTable");
        ReentrantLock lock = new ReentrantLock();
        while (sqlResultSet.next()) {
            String scanId = sqlResultSet.getString("scanId");
            int precursorCharge = sqlResultSet.getInt("precursorCharge");
            double massWithoutLinker = sqlResultSet.getDouble("massWithoutLinker");
            double precursorMass = sqlResultSet.getDouble("precursorMass");
            taskList.add(thread_pool.submit(new SearchWrap(search_obj, build_index_obj, mass_tool_obj, cal_evalue, delta_c_t, flankingPeaks, spectra_parser, lock, scanId, precursorCharge, massWithoutLinker, precursorMass, sqlPath, parameter_map.get("link_same_peptide").contentEquals("1"))));
        }
        sqlResultSet.close();
        sqlStatement.close();

        // check progress every minute, record results,and delete finished tasks.
        PreparedStatement sqlPreparedStatement = sqlConnection.prepareStatement("REPLACE INTO spectraTable (scanNum, scanId, precursorCharge, precursorMz, precursorMass, rt, massWithoutLinker, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, theoMass, score, deltaC, rank, ppm, seq1, linkSite1, proId1, seq2, linkSite2, proId2, clType, hitType, eValue, candidateNum, pointCount, rSquare, slope, intercept, startIdx, endIdx, chainScore1, chainRank1, chainScore2, chainRank2) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
        sqlConnection.setAutoCommit(false);
        int lastProgress = 0;
        int resultCount = 0;
        int totalCount = taskList.size();
        int count = 0;
        while (count < totalCount) {
            // record search results and delete finished ones.
            List<Future<SearchWrap.Entry>> toBeDeleteTaskList = new ArrayList<>(totalCount - count);
            for (Future<SearchWrap.Entry> task : taskList) {
                if (task.isDone()) {
                    if (task.get() != null) {
                        SearchWrap.Entry entry = task.get();
                        sqlPreparedStatement.setInt(1, entry.scanNum);
                        sqlPreparedStatement.setString(2, entry.scanId);
                        sqlPreparedStatement.setInt(3, entry.precursorCharge);
                        sqlPreparedStatement.setDouble(4, entry.precursorMz);
                        sqlPreparedStatement.setDouble(5, entry.precursorMass);
                        sqlPreparedStatement.setInt(6, entry.rt);
                        sqlPreparedStatement.setDouble(7, entry.massWithoutLinker);
                        sqlPreparedStatement.setString(8, entry.mgfTitle);
                        sqlPreparedStatement.setInt(9, entry.isotopeCorrectionNum);
                        sqlPreparedStatement.setDouble(10, entry.ms1PearsonCorrelationCoefficient);
                        sqlPreparedStatement.setDouble(11, entry.theoMass);
                        sqlPreparedStatement.setDouble(12, entry.score);
                        sqlPreparedStatement.setDouble(13, entry.delteC);
                        sqlPreparedStatement.setInt(14, entry.rank);
                        sqlPreparedStatement.setDouble(15, entry.ppm);
                        sqlPreparedStatement.setString(16, entry.seq1);
                        sqlPreparedStatement.setInt(17, entry.linkSite1);
                        sqlPreparedStatement.setString(18, entry.pro1);
                        sqlPreparedStatement.setString(19, entry.seq2);
                        sqlPreparedStatement.setInt(20, entry.linkSite2);
                        sqlPreparedStatement.setString(21, entry.pro2);
                        sqlPreparedStatement.setString(22, entry.clType);
                        sqlPreparedStatement.setInt(23, entry.hitType);
                        sqlPreparedStatement.setDouble(24, entry.eValue);
                        sqlPreparedStatement.setLong(25, entry.candidateNum);
                        sqlPreparedStatement.setInt(26, entry.scoreCount);
                        sqlPreparedStatement.setDouble(27, entry.rSquare);
                        sqlPreparedStatement.setDouble(28, entry.slope);
                        sqlPreparedStatement.setDouble(29, entry.intercept);
                        sqlPreparedStatement.setInt(30, entry.startIdx);
                        sqlPreparedStatement.setInt(31, entry.endIdx);
                        sqlPreparedStatement.setDouble(32, entry.chainScore1);
                        sqlPreparedStatement.setInt(33, entry.chainRank1);
                        sqlPreparedStatement.setDouble(34, entry.chainScore2);
                        sqlPreparedStatement.setInt(35, entry.chainRank2);
                        sqlPreparedStatement.executeUpdate();
                        ++resultCount;
                    }
                    toBeDeleteTaskList.add(task);
                } else if (task.isCancelled()) {
                    toBeDeleteTaskList.add(task);
                }
            }
            count += toBeDeleteTaskList.size();
            taskList.removeAll(toBeDeleteTaskList);
            taskList.trimToSize();

            sqlConnection.commit();

            int progress = count * 20 / totalCount;
            if (progress != lastProgress) {
                logger.info("Searching {}%...", progress * 5);
                lastProgress = progress;
            }

            if (count == totalCount) {
                break;
            }

            if (debug) {
                Thread.sleep(1000);
            } else {
                Thread.sleep(6000);
            }
        }

        // shutdown threads
        thread_pool.shutdown();
        if (!thread_pool.awaitTermination(60, TimeUnit.SECONDS)) {
            thread_pool.shutdownNow();
            if (!thread_pool.awaitTermination(60, TimeUnit.SECONDS))
                System.err.println("Pool did not terminate");
        }

        sqlConnection.commit();
        sqlConnection.setAutoCommit(true);
        sqlConnection.close();
        if (lock.isLocked()) {
            lock.unlock();
        }

        if (resultCount == 0) {
            logger.warn("There is no useful results.");
        } else {
            // save result
            logger.info("Estimating q-value...");
            Map<String, CalFDR.Entry> inter_result = CalFDR.calFDR(sqlPath, cal_evalue, "inter_protein");
            Map<String, CalFDR.Entry> intra_result = CalFDR.calFDR(sqlPath, cal_evalue, "intra_protein");

            logger.info("Saving results...");
            Map<String, CalFDR.Entry> allResult = new HashMap<>();
            allResult.putAll(inter_result);
            allResult.putAll(intra_result);
            saveTargetResult(allResult, build_index_obj.getProAnnotateMap(), spectra_path, cal_evalue, sqlPath);
        }
    }

    private static void saveTargetResult(Map<String, CalFDR.Entry> result, Map<String, String> pro_annotate_map, String id_file_name, boolean cal_evalue, String sqlPath) throws IOException, SQLException, NullPointerException {
        BufferedWriter intraTargetWriter = new BufferedWriter(new FileWriter(id_file_name + ".intra.target.csv"));
        BufferedWriter intraDecoyWriter = new BufferedWriter(new FileWriter(id_file_name + ".intra.decoy.csv"));
        BufferedWriter interTargetWriter = new BufferedWriter(new FileWriter(id_file_name + ".inter.target.csv"));
        BufferedWriter interDecoyWriter = new BufferedWriter(new FileWriter(id_file_name + ".inter.decoy.csv"));

        Connection sqlConnection = DriverManager.getConnection(sqlPath);
        Statement sqlStatement = sqlConnection.createStatement();
        ResultSet sqlResultSet = null;

        if (dev) {
            intraTargetWriter.write("scan_num,spectrum_id,spectrum_mz,spectrum_mass,peptide_mass,rt,C13_correction,charge,score,delta_C,ppm,peptide,protein,protein_annotation_1,protein_annotation_2,e_value,q_value,mgf_title,,MS1_pearson_correlation_coefficient,candidate_num,point_num,r_square,slope,intercept,start_idx,end_idx,chain_score_1,chain_rank_1,chain_score_2,chain_rank_2\n");
            intraDecoyWriter.write("scan_num,spectrum_id,spectrum_mz,spectrum_mass,peptide_mass,rt,C13_correction,charge,score,delta_C,ppm,peptide,protein,protein_annotation_1,protein_annotation_2,e_value,mgf_title,,MS1_pearson_correlation_coefficient,candidate_num,point_num,r_square,slope,intercept,start_idx,end_idx,chain_score_1,chain_rank_1,chain_score_2,chain_rank_2\n");
            interTargetWriter.write("scan_num,spectrum_id,spectrum_mz,spectrum_mass,peptide_mass,rt,C13_correction,charge,score,delta_C,ppm,peptide,protein,protein_annotation_1,protein_annotation_2,e_value,q_value,mgf_title,,MS1_pearson_correlation_coefficient,candidate_num,point_num,r_square,slope,intercept,start_idx,end_idx,chain_score_1,chain_rank_1,chain_score_2,chain_rank_2\n");
            interDecoyWriter.write("scan_num,spectrum_id,spectrum_mz,spectrum_mass,peptide_mass,rt,C13_correction,charge,score,delta_C,ppm,peptide,protein,protein_annotation_1,protein_annotation_2,e_value,mgf_title,,MS1_pearson_correlation_coefficient,candidate_num,point_num,r_square,slope,intercept,start_idx,end_idx,chain_score_1,chain_rank_1,chain_score_2,chain_rank_2\n");
        } else {
            intraTargetWriter.write("scan_num,spectrum_id,spectrum_mz,spectrum_mass,peptide_mass,rt,C13_correction,charge,score,delta_C,ppm,peptide,protein,protein_annotation_1,protein_annotation_2,e_value,q_value,mgf_title,\n");
            intraDecoyWriter.write("scan_num,spectrum_id,spectrum_mz,spectrum_mass,peptide_mass,rt,C13_correction,charge,score,delta_C,ppm,peptide,protein,protein_annotation_1,protein_annotation_2,e_value,mgf_title,\n");
            interTargetWriter.write("scan_num,spectrum_id,spectrum_mz,spectrum_mass,peptide_mass,rt,C13_correction,charge,score,delta_C,ppm,peptide,protein,protein_annotation_1,protein_annotation_2,e_value,q_value,mgf_title,\n");
            interDecoyWriter.write("scan_num,spectrum_id,spectrum_mz,spectrum_mass,peptide_mass,rt,C13_correction,charge,score,delta_C,ppm,peptide,protein,protein_annotation_1,protein_annotation_2,e_value,mgf_title,\n");
        }
        List<CalFDR.Entry> entryList = new ArrayList<>(result.values());
        entryList.sort(Comparator.reverseOrder());
        for (CalFDR.Entry entry : entryList) {
            sqlResultSet = sqlStatement.executeQuery(String.format(Locale.US, "SELECT scanNum, precursorMz, precursorMass, rt, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, precursorCharge, theoMass, score, deltaC, ppm, seq1, linkSite1, seq2, linkSite2, proId1, proId2, eValue, mgfTitle, candidateNum, pointCount, rSquare, slope, intercept, startIdx, endIdx, chainScore1, chainRank1, chainScore2, chainRank2, hitType, clType FROM spectraTable WHERE scanId='%s'", entry.scanId));
            if (sqlResultSet.next()) {
                List<String> proAnnotationList1 = new ArrayList<>();
                for (String s : sqlResultSet.getString("proId1").split(";")) {
                    if (s.startsWith("DECOY_")) {
                        proAnnotationList1.add("DECOY");
                    } else {
                        proAnnotationList1.add(pro_annotate_map.get(s));
                    }
                }
                List<String> proAnnotationList2 = new ArrayList<>();
                for (String s : sqlResultSet.getString("proId2").split(";")) {
                    if (s.startsWith("DECOY_")) {
                        proAnnotationList2.add("DECOY");
                    } else {
                        proAnnotationList2.add(pro_annotate_map.get(s));
                    }
                }

                if (dev) {
                    if (sqlResultSet.getString("clType").contentEquals("intra_protein")) {
                        if (sqlResultSet.getInt("hitType") == 0) {
                            intraTargetWriter.write(sqlResultSet.getInt("scanNum") + "," + entry.scanId + "," + sqlResultSet.getDouble("precursorMz") + "," + sqlResultSet.getDouble("precursorMass") + "," + sqlResultSet.getDouble("theoMass") + "," + sqlResultSet.getInt("rt") + "," + sqlResultSet.getInt("isotopeCorrectionNum") + "," + sqlResultSet.getInt("precursorCharge") + "," + sqlResultSet.getDouble("score") + "," + sqlResultSet.getDouble("deltaC") + "," + sqlResultSet.getDouble("ppm") + "," + sqlResultSet.getString("seq1") + "-" + sqlResultSet.getInt("linkSite1") + "-" + sqlResultSet.getString("seq2") + "-" + sqlResultSet.getInt("linkSite2") + "," + sqlResultSet.getString("proId1") + "-" + sqlResultSet.getString("proId2") + ",\"" + String.join(";", proAnnotationList1).replaceAll("\"", "") + "\",\"" + String.join(";", proAnnotationList2).replaceAll("\"", "") + "\"," + (cal_evalue ? String.format(Locale.US, "%E", sqlResultSet.getDouble("eValue")) : "-") + "," + entry.qValue + ",\"" + sqlResultSet.getString("mgfTitle").replaceAll("\"", "") + "\",," + sqlResultSet.getDouble("ms1PearsonCorrelationCoefficient") + "," + sqlResultSet.getInt("candidateNum") + "," + sqlResultSet.getInt("pointCount") + "," + sqlResultSet.getDouble("rSquare") + "," + sqlResultSet.getDouble("slope") + "," + sqlResultSet.getDouble("intercept") + "," + sqlResultSet.getInt("startIdx") + "," + sqlResultSet.getInt("endIdx") + "," + sqlResultSet.getDouble("chainScore1") + "," + sqlResultSet.getInt("chainRank1") + "," + sqlResultSet.getDouble("chainScore2") + "," + sqlResultSet.getInt("chainRank2") + "\n");
                        } else {
                            intraDecoyWriter.write(sqlResultSet.getInt("scanNum") + "," + entry.scanId + "," + sqlResultSet.getDouble("precursorMz") + "," + sqlResultSet.getDouble("precursorMass") + "," + sqlResultSet.getDouble("theoMass") + "," + sqlResultSet.getInt("rt") + "," + sqlResultSet.getInt("isotopeCorrectionNum") + "," + sqlResultSet.getInt("precursorCharge") + "," + sqlResultSet.getDouble("score") + "," + sqlResultSet.getDouble("deltaC") + "," + sqlResultSet.getDouble("ppm") + "," + sqlResultSet.getString("seq1") + "-" + sqlResultSet.getInt("linkSite1") + "-" + sqlResultSet.getString("seq2") + "-" + sqlResultSet.getInt("linkSite2") + "," + sqlResultSet.getString("proId1") + "-" + sqlResultSet.getString("proId2") + ",\"" + String.join(";", proAnnotationList1).replaceAll("\"", "") + "\",\"" + String.join(";", proAnnotationList2).replaceAll("\"", "") + "\"," + (cal_evalue ? String.format(Locale.US, "%E", sqlResultSet.getDouble("eValue")) : "-") + ",\"" + sqlResultSet.getString("mgfTitle").replaceAll("\"", "") + "\",," + sqlResultSet.getDouble("ms1PearsonCorrelationCoefficient") + "," + sqlResultSet.getInt("candidateNum") + "," + sqlResultSet.getInt("pointCount") + "," + sqlResultSet.getDouble("rSquare") + "," + sqlResultSet.getDouble("slope") + "," + sqlResultSet.getDouble("intercept") + "," + sqlResultSet.getInt("startIdx") + "," + sqlResultSet.getInt("endIdx") + "," + sqlResultSet.getDouble("chainScore1") + "," + sqlResultSet.getInt("chainRank1") + "," + sqlResultSet.getDouble("chainScore2") + "," + sqlResultSet.getInt("chainRank2") + "\n");
                        }
                    } else {
                        if (sqlResultSet.getInt("hitType") == 0) {
                            interTargetWriter.write(sqlResultSet.getInt("scanNum") + "," + entry.scanId + "," + sqlResultSet.getDouble("precursorMz") + "," + sqlResultSet.getDouble("precursorMass") + "," + sqlResultSet.getDouble("theoMass") + "," + sqlResultSet.getInt("rt") + "," + sqlResultSet.getInt("isotopeCorrectionNum") + "," + sqlResultSet.getInt("precursorCharge") + "," + sqlResultSet.getDouble("score") + "," + sqlResultSet.getDouble("deltaC") + "," + sqlResultSet.getDouble("ppm") + "," + sqlResultSet.getString("seq1") + "-" + sqlResultSet.getInt("linkSite1") + "-" + sqlResultSet.getString("seq2") + "-" + sqlResultSet.getInt("linkSite2") + "," + sqlResultSet.getString("proId1") + "-" + sqlResultSet.getString("proId2") + ",\"" + String.join(";", proAnnotationList1).replaceAll("\"", "") + "\",\"" + String.join(";", proAnnotationList2).replaceAll("\"", "") + "\"," + (cal_evalue ? String.format(Locale.US, "%E", sqlResultSet.getDouble("eValue")) : "-") + "," + entry.qValue + ",\"" + sqlResultSet.getString("mgfTitle").replaceAll("\"", "") + "\",," + sqlResultSet.getDouble("ms1PearsonCorrelationCoefficient") + "," + sqlResultSet.getInt("candidateNum") + "," + sqlResultSet.getInt("pointCount") + "," + sqlResultSet.getDouble("rSquare") + "," + sqlResultSet.getDouble("slope") + "," + sqlResultSet.getDouble("intercept") + "," + sqlResultSet.getInt("startIdx") + "," + sqlResultSet.getInt("endIdx") + "," + sqlResultSet.getDouble("chainScore1") + "," + sqlResultSet.getInt("chainRank1") + "," + sqlResultSet.getDouble("chainScore2") + "," + sqlResultSet.getInt("chainRank2") + "\n");
                        } else {
                            interDecoyWriter.write(sqlResultSet.getInt("scanNum") + "," + entry.scanId + "," + sqlResultSet.getDouble("precursorMz") + "," + sqlResultSet.getDouble("precursorMass") + "," + sqlResultSet.getDouble("theoMass") + "," + sqlResultSet.getInt("rt") + "," + sqlResultSet.getInt("isotopeCorrectionNum") + "," + sqlResultSet.getInt("precursorCharge") + "," + sqlResultSet.getDouble("score") + "," + sqlResultSet.getDouble("deltaC") + "," + sqlResultSet.getDouble("ppm") + "," + sqlResultSet.getString("seq1") + "-" + sqlResultSet.getInt("linkSite1") + "-" + sqlResultSet.getString("seq2") + "-" + sqlResultSet.getInt("linkSite2") + "," + sqlResultSet.getString("proId1") + "-" + sqlResultSet.getString("proId2") + ",\"" + String.join(";", proAnnotationList1).replaceAll("\"", "") + "\",\"" + String.join(";", proAnnotationList2).replaceAll("\"", "") + "\"," + (cal_evalue ? String.format(Locale.US, "%E", sqlResultSet.getDouble("eValue")) : "-") + ",\"" + sqlResultSet.getString("mgfTitle").replaceAll("\"", "") + "\",," + sqlResultSet.getDouble("ms1PearsonCorrelationCoefficient") + "," + sqlResultSet.getInt("candidateNum") + "," + sqlResultSet.getInt("pointCount") + "," + sqlResultSet.getDouble("rSquare") + "," + sqlResultSet.getDouble("slope") + "," + sqlResultSet.getDouble("intercept") + "," + sqlResultSet.getInt("startIdx") + "," + sqlResultSet.getInt("endIdx") + "," + sqlResultSet.getDouble("chainScore1") + "," + sqlResultSet.getInt("chainRank1") + "," + sqlResultSet.getDouble("chainScore2") + "," + sqlResultSet.getInt("chainRank2") + "\n");
                        }
                    }
                } else {
                    if (sqlResultSet.getString("clType").contentEquals("intra_protein")) {
                        if (sqlResultSet.getInt("hitType") == 0) {
                            intraTargetWriter.write(sqlResultSet.getInt("scanNum") + "," + entry.scanId + "," + sqlResultSet.getDouble("precursorMz") + "," + sqlResultSet.getDouble("precursorMass") + "," + sqlResultSet.getDouble("theoMass") + "," + sqlResultSet.getInt("rt") + "," + sqlResultSet.getInt("isotopeCorrectionNum") + "," + sqlResultSet.getInt("precursorCharge") + "," + sqlResultSet.getDouble("score") + "," + sqlResultSet.getDouble("deltaC") + "," + sqlResultSet.getDouble("ppm") + "," + sqlResultSet.getString("seq1") + "-" + sqlResultSet.getInt("linkSite1") + "-" + sqlResultSet.getString("seq2") + "-" + sqlResultSet.getInt("linkSite2") + "," + sqlResultSet.getString("proId1") + "-" + sqlResultSet.getString("proId2") + ",\"" + String.join(";", proAnnotationList1).replaceAll("\"", "") + "\",\"" + String.join(";", proAnnotationList2).replaceAll("\"", "") + "\"," + (cal_evalue ? String.format(Locale.US, "%E", sqlResultSet.getDouble("eValue")) : "-") + "," + entry.qValue + ",\"" + sqlResultSet.getString("mgfTitle").replaceAll("\"", "") + "\"\n");
                        } else {
                            intraDecoyWriter.write(sqlResultSet.getInt("scanNum") + "," + entry.scanId + "," + sqlResultSet.getDouble("precursorMz") + "," + sqlResultSet.getDouble("precursorMass") + "," + sqlResultSet.getDouble("theoMass") + "," + sqlResultSet.getInt("rt") + "," + sqlResultSet.getInt("isotopeCorrectionNum") + "," + sqlResultSet.getInt("precursorCharge") + "," + sqlResultSet.getDouble("score") + "," + sqlResultSet.getDouble("deltaC") + "," + sqlResultSet.getDouble("ppm") + "," + sqlResultSet.getString("seq1") + "-" + sqlResultSet.getInt("linkSite1") + "-" + sqlResultSet.getString("seq2") + "-" + sqlResultSet.getInt("linkSite2") + "," + sqlResultSet.getString("proId1") + "-" + sqlResultSet.getString("proId2") + ",\"" + String.join(";", proAnnotationList1).replaceAll("\"", "") + "\",\"" + String.join(";", proAnnotationList2).replaceAll("\"", "") + "\"," + (cal_evalue ? String.format(Locale.US, "%E", sqlResultSet.getDouble("eValue")) : "-") + ",\"" + sqlResultSet.getString("mgfTitle").replaceAll("\"", "") + "\",\n");
                        }
                    } else {
                        if (sqlResultSet.getInt("hitType") == 0) {
                            interTargetWriter.write(sqlResultSet.getInt("scanNum") + "," + entry.scanId + "," + sqlResultSet.getDouble("precursorMz") + "," + sqlResultSet.getDouble("precursorMass") + "," + sqlResultSet.getDouble("theoMass") + "," + sqlResultSet.getInt("rt") + "," + sqlResultSet.getInt("isotopeCorrectionNum") + "," + sqlResultSet.getInt("precursorCharge") + "," + sqlResultSet.getDouble("score") + "," + sqlResultSet.getDouble("deltaC") + "," + sqlResultSet.getDouble("ppm") + "," + sqlResultSet.getString("seq1") + "-" + sqlResultSet.getInt("linkSite1") + "-" + sqlResultSet.getString("seq2") + "-" + sqlResultSet.getInt("linkSite2") + "," + sqlResultSet.getString("proId1") + "-" + sqlResultSet.getString("proId2") + ",\"" + String.join(";", proAnnotationList1).replaceAll("\"", "") + "\",\"" + String.join(";", proAnnotationList2).replaceAll("\"", "") + "\"," + (cal_evalue ? String.format(Locale.US, "%E", sqlResultSet.getDouble("eValue")) : "-") + "," + entry.qValue + ",\"" + sqlResultSet.getString("mgfTitle").replaceAll("\"", "") + "\",\n");
                        } else {
                            interDecoyWriter.write(sqlResultSet.getInt("scanNum") + "," + entry.scanId + "," + sqlResultSet.getDouble("precursorMz") + "," + sqlResultSet.getDouble("precursorMass") + "," + sqlResultSet.getDouble("theoMass") + "," + sqlResultSet.getInt("rt") + "," + sqlResultSet.getInt("isotopeCorrectionNum") + "," + sqlResultSet.getInt("precursorCharge") + "," + sqlResultSet.getDouble("score") + "," + sqlResultSet.getDouble("deltaC") + "," + sqlResultSet.getDouble("ppm") + "," + sqlResultSet.getString("seq1") + "-" + sqlResultSet.getInt("linkSite1") + "-" + sqlResultSet.getString("seq2") + "-" + sqlResultSet.getInt("linkSite2") + "," + sqlResultSet.getString("proId1") + "-" + sqlResultSet.getString("proId2") + ",\"" + String.join(";", proAnnotationList1).replaceAll("\"", "") + "\",\"" + String.join(";", proAnnotationList2).replaceAll("\"", "") + "\"," + (cal_evalue ? String.format(Locale.US, "%E", sqlResultSet.getDouble("eValue")) : "-") + ",\"" + sqlResultSet.getString("mgfTitle").replaceAll("\"", "") + "\",\n");
                        }
                    }
                }
            } else {
                throw new NullPointerException(String.format(Locale.US, "Cannot find the record of scan %s.", entry.scanId));
            }
        }
        intraTargetWriter.close();
        intraDecoyWriter.close();
        interTargetWriter.close();
        interDecoyWriter.close();
        if (sqlResultSet != null) {
            sqlResultSet.close();
        }
        sqlStatement.close();
        sqlConnection.close();
    }

    private static void help() {
        String help_str = "ECL2 version " + version + "\r\n"
                + "A cross-linked peptides identification tool with exhaustive searching and linear computational complexity\r\n"
                + "Author: Fengchao Yu\r\n"
                + "Email: fyuab@connect.ust.hk\r\n"
                + "ECL2 usage: java -Xmx25g -jar /path/to/ECL2.jar <parameter_file> <data_file>\r\n"
                + "\t<parameter_file>: parameter file. Can be download along with ECL2.\r\n"
                + "\t<data_file>: spectra data file (mzXML or MGF)\r\n"
                + "\texample: java -Xmx32g -jar ECL2.jar parameter.def data.mzxml\r\n";
        System.out.print(help_str);
        System.exit(1);
    }
}
