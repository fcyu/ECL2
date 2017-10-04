package proteomics;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Index.BuildIndex;
import proteomics.Parameter.Parameter;
import proteomics.Search.SearchWrap;
import proteomics.Types.FinalResultEntry;
import proteomics.Search.Search;
import proteomics.Spectrum.PreSpectra;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.SpectrumEntry;
import proteomics.Validation.CalFDR;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;

import java.io.*;
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.util.*;
import java.util.concurrent.*;

public class ECL2 {

    public static final int score_point_t = 15000;

    private static final Logger logger = LoggerFactory.getLogger(ECL2.class);
    public static final String version = "2.1.5";

    public static boolean debug;
    public static boolean dev;

    public static void main(String[] args) {
        // Process inputs
        if (args.length != 2) {
            help();
        }

        // Set parameters
        String parameter_path = args[0].trim();
        String spectra_path = args[1].trim();

        logger.info("ECL2 Version {}", version);
        logger.info("Author: Fengchao Yu. Email: fyuab@connect.ust.hk");

        try {
            String hostName = InetAddress.getLocalHost().getHostName();
            logger.info("Computer: {}", hostName);
        } catch (UnknownHostException ex) {
            logger.warn("Cannot get the computer's name.");
        }

        logger.info("Parameter file: {}.", parameter_path);
        logger.info("Spectra file: {}.", spectra_path);

        // Get the parameter map
        Parameter parameter = new Parameter(parameter_path);
        Map<String, String> parameter_map = parameter.returnParameterMap();
        int max_common_ion_charge = Integer.valueOf(parameter_map.get("max_common_ion_charge"));

        boolean cal_evalue = true;
        if (parameter_map.containsKey("cal_evalue") && parameter_map.get("cal_evalue").trim().contentEquals("0")) {
            cal_evalue = false;
        }

        float delta_c_t = 0;
        if (parameter_map.containsKey("delta_c_t")) {
            delta_c_t = Float.valueOf(parameter_map.get("delta_c_t"));
        }

        boolean flankingPeaks = true;
        if (parameter_map.containsKey("flanking_peaks") && parameter_map.get("flanking_peaks").contentEquals("0")) {
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
        logger.info("Protein database: {}", parameter_map.get("db"));
        logger.info("Indexing protein database...");
        BuildIndex build_index_obj = new BuildIndex(parameter_map);
        MassTool mass_tool_obj = build_index_obj.returnMassTool();

        logger.info("Peptide chains number: {}", build_index_obj.getSeqEntryMap().size());

        logger.info("Reading spectra...");
        JMzReader spectra_parser = null;
        String ext = "";
        try {
            File spectra_file = new File(spectra_path);
            if ((!spectra_file.exists() || (spectra_file.isDirectory()))) {
                throw new FileNotFoundException("The spectra file not found.");
            }
            String[] temp = spectra_path.split("\\.");
            ext = temp[temp.length - 1];
            if (ext.contentEquals("mzXML")) {
                spectra_parser = new MzXMLFile(spectra_file);
            } else if (ext.toLowerCase().contentEquals("mgf")) {
                spectra_parser = new MgfFile(spectra_file);
            } else {
                logger.error("Unsupported data format {}. ECL2 only support mzXML and MGF.", ext);
                System.exit(1);
            }
        } catch (Exception ex) {
            logger.error(ex.getMessage());
            ex.printStackTrace();
            System.exit(1);
        }

        PreSpectra pre_spectra_obj = new PreSpectra(spectra_parser, build_index_obj, parameter_map, ext);
        Map<Integer, SpectrumEntry> num_spectrum_map = pre_spectra_obj.getNumSpectrumMap();
        Integer[] scanNumArray = num_spectrum_map.keySet().toArray(new Integer[num_spectrum_map.size()]);
        Arrays.sort(scanNumArray);

        logger.info("Useful MS/MS spectra number: {}", num_spectrum_map.size());

        logger.info("Start searching...");

        int thread_num = Integer.valueOf(parameter_map.get("thread_num"));
        if (thread_num == 0) {
            thread_num = 1 + Runtime.getRuntime().availableProcessors();
        }
        ExecutorService thread_pool = Executors.newFixedThreadPool(thread_num);
        Search search_obj = new Search(build_index_obj, parameter_map);
        List<Future<FinalResultEntry>> taskList = new LinkedList<>();
        for (int scanNum : scanNumArray) {
            taskList.add(thread_pool.submit(new SearchWrap(search_obj, num_spectrum_map.get(scanNum), build_index_obj, mass_tool_obj, max_common_ion_charge, build_index_obj.getSeqProMap(), cal_evalue, delta_c_t, flankingPeaks)));
        }

        // check progress every minute, record results,and delete finished tasks.
        int lastProgress = 0;
        Set<FinalResultEntry> final_search_results = new HashSet<>(scanNumArray.length + 1, 1);
        try {
            while (!taskList.isEmpty()) {
                // record search results and delete finished ones.
                List<Future<FinalResultEntry>> toBeDeleteTaskList = new LinkedList<>();
                for (Future<FinalResultEntry> task : taskList) {
                    if (task.isDone()) {
                        if (task.get() != null) {
                            final_search_results.add(task.get());
                        }
                        toBeDeleteTaskList.add(task);
                    } else if (task.isCancelled()) {
                        toBeDeleteTaskList.add(task);
                    }
                }
                taskList.removeAll(toBeDeleteTaskList);

                int progress = (scanNumArray.length - taskList.size()) * 20 / scanNumArray.length;
                if (progress != lastProgress) {
                    logger.info("Searching {}%...", progress * 5);
                    lastProgress = progress;
                }
                if (debug) {
                    Thread.sleep(1000);
                } else {
                    Thread.sleep(6000);
                }
            }
        } catch (InterruptedException | ExecutionException ex) {
            ex.printStackTrace();
            logger.error(ex.toString());
            System.exit(1);
        }

        // shutdown threads
        thread_pool.shutdown();
        try {
            if (!thread_pool.awaitTermination(60, TimeUnit.SECONDS)) {
                thread_pool.shutdownNow();
                if (!thread_pool.awaitTermination(60, TimeUnit.SECONDS))
                    System.err.println("Pool did not terminate");
            }
        } catch (InterruptedException ie) {
            thread_pool.shutdownNow();
            Thread.currentThread().interrupt();
            logger.error("Threads were not finished normally.");
            System.exit(1);
        }

        if (final_search_results.isEmpty()) {
            logger.warn("There is no useful PSM.");
        } else {
            // save result
            logger.info("Estimating q-value...");
            List<List<FinalResultEntry>> picked_result = pickResult(final_search_results);
            List<FinalResultEntry> intra_result = new LinkedList<>();
            List<FinalResultEntry> inter_result = new LinkedList<>();
            if (!picked_result.get(0).isEmpty()) {
                CalFDR cal_fdr_obj = new CalFDR(picked_result.get(0), cal_evalue);
                intra_result = cal_fdr_obj.includeStats(cal_evalue);
                intra_result.sort(Collections.reverseOrder());
            }
            if (!picked_result.get(1).isEmpty()) {
                CalFDR cal_fdr_obj = new CalFDR(picked_result.get(1), cal_evalue);
                inter_result = cal_fdr_obj.includeStats(cal_evalue);
                inter_result.sort(Collections.reverseOrder());
            }

            logger.info("Saving results...");
            saveTargetResult(intra_result, build_index_obj.getProAnnotateMap(), spectra_path, true, cal_evalue);
            saveTargetResult(inter_result, build_index_obj.getProAnnotateMap(), spectra_path, false, cal_evalue);
            saveDecoyResult(intra_result, build_index_obj.getProAnnotateMap(), spectra_path, true, cal_evalue);
            saveDecoyResult(inter_result, build_index_obj.getProAnnotateMap(), spectra_path, false, cal_evalue);
        }
        logger.info("Done.");
    }

    private static void saveTargetResult(List<FinalResultEntry> result, Map<String, String> pro_annotate_map, String id_file_name, boolean is_intra, boolean cal_evalue) {
        try {
            BufferedWriter writer;
            if (is_intra) {
                writer = new BufferedWriter(new FileWriter(id_file_name + ".intra.target.csv"));
            } else {
                writer = new BufferedWriter(new FileWriter(id_file_name + ".inter.target.csv"));
            }

            if (dev) {
                writer.write("scan_num,spectrum_id,spectrum_mz,spectrum_mass,peptide_mass,rt,C13_correction,charge,score,delta_C,ppm,peptide,protein,protein_annotation_1,protein_annotation_2,e_value,q_value,mgf_title,,candidate_num,point_num,r_square,slope,intercept,start_idx,end_idx,chain_score_1,chain_rank_1,chain_score_2,chain_rank_2\n");
            } else {
                writer.write("scan_num,spectrum_id,spectrum_mz,spectrum_mass,peptide_mass,rt,C13_correction,charge,score,delta_C,ppm,peptide,protein,protein_annotation_1,protein_annotation_2,e_value,q_value,mgf_title,\n");
            }
            for (FinalResultEntry re : result) {
                if (re.hit_type == 0) {
                    int link_site_1 = re.link_site_1;
                    int link_site_2 = re.link_site_2;
                    String pro_1 = re.pro_id_1;
                    if (re.pro_id_1.contains(";")) {
                        pro_1 = re.pro_id_1.split(";")[0];
                    }
                    String pro_2 = re.pro_id_2;
                    if (re.pro_id_2.contains(";")) {
                        pro_2 = re.pro_id_2.split(";")[0];
                    }

                    String annotate_1 = pro_annotate_map.get(pro_1).replace(',', ';');
                    String annotate_2 = pro_annotate_map.get(pro_2).replace(',', ';');

                    if (dev) {
                        writer.write(re.scan_num + "," + re.spectrum_id + "," + re.spectrum_mz + "," + re.spectrum_mass + "," + re.peptide_mass + "," + re.rt + "," + re.C13_correction + "," + re.charge + "," + String.format(Locale.US, "%.4f", re.score) + "," + re.delta_c + "," + String.format(Locale.US, "%.2f", re.ppm) + "," + re.seq_1 + "-" + link_site_1 + "-" + re.seq_2 + "-" + link_site_2 + "," + re.pro_id_1 + "-" + re.pro_id_2 + ",\"" + annotate_1 + "\",\"" + annotate_2 + "\"," + (cal_evalue ? String.format(Locale.US, "%E", re.e_value) : "-") + "," + String.format(Locale.US, "%.4f", re.qvalue) + ",\"" + re.mgfTitle + "\",," + re.candidate_num + "," + re.point_count + "," + String.format(Locale.US, "%.4f", re.r_square) + "," + String.format(Locale.US, "%.4f", re.slope) + "," + String.format(Locale.US, "%.4f", re.intercept) + "," + re.start_idx + "," + re.end_idx + "," + String.format(Locale.US, "%.4f", re.chain_score_1) + "," + re.chain_rank_1 + "," + String.format(Locale.US, "%.4f", re.chain_score_2) + "," + re.chain_rank_2 + "\n");
                    } else {
                        writer.write(re.scan_num + "," + re.spectrum_id + "," + re.spectrum_mz + "," + re.spectrum_mass + "," + re.peptide_mass + "," + re.rt + "," +  re.C13_correction + "," + re.charge + "," + String.format(Locale.US, "%.4f", re.score) + "," + re.delta_c + "," + String.format(Locale.US, "%.2f", re.ppm) + "," + re.seq_1 + "-" + link_site_1 + "-" + re.seq_2 + "-" + link_site_2 + "," + re.pro_id_1 + "-" + re.pro_id_2 + ",\"" + annotate_1 + "\",\"" + annotate_2 + "\"," + (cal_evalue ? String.format(Locale.US, "%E", re.e_value) : "-") + "," + String.format(Locale.US, "%.4f", re.qvalue) + ",\"" + re.mgfTitle + "\"\n");
                    }
                }
            }
            writer.close();
        } catch (IOException ex) {
            logger.error(ex.getMessage());
            System.exit(1);
        }
    }

    private static void saveDecoyResult(List<FinalResultEntry> result, Map<String, String> pro_annotate_map, String id_file_name, boolean is_intra, boolean cal_evalue) {
        try {
            BufferedWriter writer;
            if (is_intra) {
                writer = new BufferedWriter(new FileWriter(id_file_name + ".intra.decoy.csv"));
            } else {
                writer = new BufferedWriter(new FileWriter(id_file_name + ".inter.decoy.csv"));
            }

            if (dev) {
                writer.write("scan_num,spectrum_id,spectrum_mz,spectrum_mass,peptide_mass,rt,C13_correction,charge,score,delta_C,ppm,peptide,protein,protein_annotation_1,protein_annotation_2,e_value,mgf_title,,candidate_num,point_num,r_square,slope,intercept,start_idx,end_idx,chain_score_1,chain_rank_1,chain_score_2,chain_rank_2\n");
            } else {
                writer.write("scan_num,spectrum_id,spectrum_mz,spectrum_mass,peptide_mass,rt,C13_correction,charge,score,delta_C,ppm,peptide,protein,protein_annotation_1,protein_annotation_2,e_value,mgf_title\n");
            }
            for (FinalResultEntry re : result) {
                if ((re.hit_type == 1) || (re.hit_type == 2)) {
                    int link_site_1 = re.link_site_1;
                    int link_site_2 = re.link_site_2;

                    String annotate_1 = "DECOY";
                    String annotate_2 = "DECOY";

                    if (!re.pro_id_1.startsWith("DECOY")) {
                        String pro_1 = re.pro_id_1.split(";")[0];
                        annotate_1 = pro_annotate_map.get(pro_1).replace(',', ';');
                    }

                    if (!re.pro_id_2.startsWith("DECOY")) {
                        String pro_2 = re.pro_id_2.split(";")[0];
                        annotate_2 = pro_annotate_map.get(pro_2).replace(',', ';');
                    }

                    if (dev) {
                        writer.write(re.scan_num + "," + re.spectrum_id + "," + re.spectrum_mz + "," + re.spectrum_mass + "," + re.peptide_mass + "," + re.rt + "," + re.C13_correction + "," + re.charge + "," + String.format(Locale.US, "%.4f", re.score) + "," + re.delta_c + "," + String.format(Locale.US, "%.2f", re.ppm) + "," + re.seq_1 + "-" + link_site_1 + "-" + re.seq_2 + "-" + link_site_2 + "," + re.pro_id_1 + "-" + re.pro_id_2 + ",\"" + annotate_1 + "\",\"" + annotate_2 + "\"," + (cal_evalue ? String.format(Locale.US, "%E", re.e_value) : "-") + ",\"" + re.mgfTitle + "\",," + re.candidate_num + "," + re.point_count + "," + String.format(Locale.US, "%.4f", re.r_square) + "," + String.format(Locale.US, "%.4f", re.slope) + "," + String.format(Locale.US, "%.4f", re.intercept) + "," + re.start_idx + "," + re.end_idx + "," + String.format(Locale.US, "%.4f", re.chain_score_1) + "," + re.chain_rank_1 + "," + String.format(Locale.US, "%.4f", re.chain_score_2) + "," + re.chain_rank_2 + "\n");
                    } else {
                        writer.write(re.scan_num + "," + re.spectrum_id + "," + re.spectrum_mz + "," + re.spectrum_mass + "," + re.peptide_mass + "," + re.rt + "," + re.C13_correction + "," + re.charge + "," + String.format(Locale.US, "%.4f", re.score) + "," + re.delta_c + "," + String.format(Locale.US, "%.2f", re.ppm) + "," + re.seq_1 + "-" + link_site_1 + "-" + re.seq_2 + "-" + link_site_2 + "," + re.pro_id_1 + "-" + re.pro_id_2 + ",\"" + annotate_1 + "\",\"" + annotate_2 + "\"," +  (cal_evalue ? String.format(Locale.US, "%E", re.e_value) : "-") + ",\"" + re.mgfTitle + "\"\n");
                    }
                }
            }
            writer.close();
        } catch (IOException ex) {
            logger.error(ex.getMessage());
            System.exit(1);
        }
    }

    private static List<List<FinalResultEntry>> pickResult(Set<FinalResultEntry> search_result) {
        List<List<FinalResultEntry>> picked_result = new LinkedList<>();
        List<FinalResultEntry> inter_protein_result = new LinkedList<>();
        List<FinalResultEntry> intra_protein_result = new LinkedList<>();

        for (FinalResultEntry result_entry : search_result) {
            if (result_entry.cl_type.contentEquals("intra_protein")) {
                intra_protein_result.add(result_entry);
            } else {
                inter_protein_result.add(result_entry);
            }
        }

        picked_result.add(intra_protein_result);
        picked_result.add(inter_protein_result);

        return picked_result;
    }

    private static int finishedFutureNum(List<Future<FinalResultEntry>> futureList) {
        int count = 0;
        for (Future<FinalResultEntry> future : futureList) {
            if (future.isDone()) {
                ++count;
            }
        }
        return count;
    }

    private static boolean hasUnfinishedFuture(List<Future<FinalResultEntry>> futureList) {
        for (Future<FinalResultEntry> future : futureList) {
            if (!future.isDone()) {
                return true;
            }
        }
        return false;
    }

    private static void help() {
        String help_str = "ECL2 version " + version + "\r\n"
                + "A cross-linked peptides identification tool with exhaustive searching and linear computational complexity\r\n"
                + "Author: Fengchao Yu\r\n"
                + "Email: fyuab@connect.ust.hk\r\n"
                + "ECL2 usage: java -Xmx25g -jar /path/to/ECL2.jar <parameter_file> <data_file>\r\n"
                + "\t<parameter_file>: parameter file. Can be download along with ECL2.\r\n"
                + "\t<data_file>: spectra data file (mzXML)\r\n"
                + "\texample: java -Xmx32g -jar ECL2.jar parameter.def data.mzxml\r\n";
        System.out.print(help_str);
        System.exit(1);
    }
}
