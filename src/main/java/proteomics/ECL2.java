package proteomics;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Index.BuildIndex;
import proteomics.Parameter.Parameter;
import proteomics.Search.SearchWrap;
import proteomics.Spectrum.FilterSpectra;
import proteomics.Types.FinalResultEntry;
import proteomics.Search.Search;
import proteomics.Spectrum.PreSpectra;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.SpectrumEntry;
import proteomics.Validation.CalFDR;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLParsingException;

import java.io.*;
import java.util.*;
import java.util.concurrent.*;

public class ECL2 {

    public static final boolean cal_evalue = true;
    public static final double delta_c_t = 0;
    public static final int score_point_t = 6000;
    public static final int chain_score_point_t = 10;
    public static final int random_score_point_t = 2;

    private static final Logger logger = LoggerFactory.getLogger(ECL2.class);
    private static final String version = "2.0.0";

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
        logger.info("Parameter file: {}.", parameter_path);
        logger.info("Spectra file: {}.", spectra_path);

        // Get the parameter map
        Parameter parameter = new Parameter(parameter_path);
        Map<String, String> parameter_map = parameter.returnParameterMap();
        int batch_size = Integer.valueOf(parameter_map.get("batch_size"));
        int max_common_ion_charge = Integer.valueOf(parameter_map.get("max_common_ion_charge"));

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
        MzXMLFile spectra_parser = null;
        try {
            File spectra_file = new File(spectra_path);
            if ((!spectra_file.exists() || (spectra_file.isDirectory()))) {
                throw new FileNotFoundException("The spectra file not found.");
            }
            String[] temp = spectra_path.split("\\.");
            String ext = temp[temp.length - 1];
            if (ext.contentEquals("mzXML")) {
                spectra_parser = new MzXMLFile(spectra_file);
            } else {
                logger.error("Unsupported data format {}. ECL 2.0 only support mzXML.", ext);
                System.exit(1);
            }
        } catch (FileNotFoundException | MzXMLParsingException | NullPointerException ex) {
            logger.error(ex.getMessage());
            ex.printStackTrace();
            System.exit(1);
        }

        FilterSpectra filter_spectra_obj = new FilterSpectra(spectra_parser, parameter_map, mass_tool_obj.getMassTable());
        FilterSpectra.MassScan[] scan_num_array = filter_spectra_obj.getScanNumArray();

        logger.info("Useful MS/MS spectra number: {}", scan_num_array.length);

        int thread_num = Integer.valueOf(parameter_map.get("thread_num"));
        if (thread_num == 0) {
            thread_num = 1 + Runtime.getRuntime().availableProcessors();
        }
        ExecutorService thread_pool = Executors.newFixedThreadPool(thread_num);

        List<FinalResultEntry> final_search_results = new LinkedList<>();

        Search search_obj = new Search(build_index_obj, parameter_map);
        int start_idx;
        int end_idx = 0;
        while (end_idx < scan_num_array.length) {
            start_idx = end_idx;
            end_idx = Math.min(start_idx + batch_size, scan_num_array.length);

            logger.info("Searching batch {} - {} ({}%)...", start_idx, end_idx, String.format("%.1f", (float) end_idx * 100 / (float) scan_num_array.length));
            PreSpectra pre_spectra_obj = new PreSpectra(spectra_parser, build_index_obj, parameter_map, mass_tool_obj, scan_num_array, start_idx, end_idx);
            TreeMap<Integer, Set<SpectrumEntry>> bin_spectra_map = pre_spectra_obj.returnMassNumMap();
            Map<Integer, SpectrumEntry> num_spectrum_map = pre_spectra_obj.getNumSpectrumMap();

            if (!bin_spectra_map.isEmpty()) {
                // divide and run in parallel
                int batch_size_2 = (bin_spectra_map.size() / thread_num) + 1;
                Integer[] bin_array = bin_spectra_map.keySet().toArray(new Integer[bin_spectra_map.size()]);
                Collection<SearchWrap> task_list = new LinkedList<>();
                for (int i = 0; i < thread_num; ++i) {
                    int left_idx = i * batch_size_2;
                    int right_idx = Math.min((i + 1) * batch_size_2, bin_array.length - 1);
                    if (left_idx > bin_array.length - 1) {
                        break;
                    }

                    if (right_idx < bin_array.length - 1) {
                        task_list.add(new SearchWrap(search_obj, bin_spectra_map.subMap(bin_array[left_idx], true, bin_array[right_idx], false), num_spectrum_map, build_index_obj, mass_tool_obj, max_common_ion_charge));
                    } else {
                        task_list.add(new SearchWrap(search_obj, bin_spectra_map.subMap(bin_array[left_idx], true, bin_array[right_idx], true), num_spectrum_map, build_index_obj, mass_tool_obj, max_common_ion_charge));
                    }
                }

                // record search results and save them to disk for the sake of memory.
                try {
                    List<Future<List<FinalResultEntry>>> temp_result_list = thread_pool.invokeAll(task_list);
                    for (Future<List<FinalResultEntry>> temp_result : temp_result_list) {
                        if (temp_result.isDone() && !temp_result.isCancelled()) {
                            if (temp_result.get() != null) {
                                final_search_results.addAll(temp_result.get());
                            }
                        } else {
                            logger.error("Threads were not finished normally.");
                            System.exit(1);
                        }
                    }
                } catch (Exception ex) {
                    logger.error(ex.getMessage());
                    ex.printStackTrace();
                    System.exit(1);
                }
            }

            System.gc();
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
            CalFDR cal_fdr_obj = new CalFDR(picked_result.get(0));
            List<FinalResultEntry> intra_result = cal_fdr_obj.includeStats();
            Collections.sort(intra_result, Collections.<FinalResultEntry>reverseOrder());
            cal_fdr_obj = new CalFDR(picked_result.get(1));
            List<FinalResultEntry> inter_result = cal_fdr_obj.includeStats();
            Collections.sort(inter_result, Collections.<FinalResultEntry>reverseOrder());
            logger.info("Saving results...");
            saveTargetResult(intra_result, build_index_obj.getProAnnotateMap(), spectra_path, true);
            saveTargetResult(inter_result, build_index_obj.getProAnnotateMap(), spectra_path, false);
            saveDecoyResult(intra_result, build_index_obj.getProAnnotateMap(), spectra_path, true);
            saveDecoyResult(inter_result, build_index_obj.getProAnnotateMap(), spectra_path, false);
        }
        logger.info("Done.");
    }

    private static void saveTargetResult(List<FinalResultEntry> result, Map<String, String> pro_annotate_map, String id_file_name, boolean is_intra) {
        try {
            BufferedWriter writer;
            if (is_intra) {
                writer = new BufferedWriter(new FileWriter(id_file_name + ".intra.target.csv"));
            } else {
                writer = new BufferedWriter(new FileWriter(id_file_name + ".inter.target.csv"));
            }

            if (dev) {
                writer.write("scan_num,spectrum_id,spectrum_mz,spectrum_mass,peptide_mass,rt,C13_correction,charge,score,delta_C,ppm,peptide,protein,protein_annotation_1,protein_annotation_2,e_value,q_value,,candidate_num,point_num,r_square,slope,intercept,start_idx,end_idx,chain_score_1,chain_rank_1,chain_score_2,chain_rank_2\n");
            } else {
                writer.write("scan_num,spectrum_id,spectrum_mz,spectrum_mass,peptide_mass,rt,C13_correction,charge,score,delta_C,ppm,peptide,protein,protein_annotation_1,protein_annotation_2,e_value,q_value\n");
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

                    String annotate_1 = pro_annotate_map.get(pro_1).replace(",", ";");
                    String annotate_2 = pro_annotate_map.get(pro_2).replace(",", ";");

                    if (dev) {
                        writer.write(re.scan_num + "," + re.spectrum_id + "," + re.spectrum_mz + "," + re.spectrum_mass + "," + re.peptide_mass + "," + re.rt + "," + re.C13_correction + "," + re.charge + "," + String.format("%.4f", re.score) + "," + re.delta_c + "," + String.format("%.2f", re.ppm) + "," + re.seq_1 + "-" + link_site_1 + "-" + re.seq_2 + "-" + link_site_2 + "," + re.pro_id_1 + "-" + re.pro_id_2 + ',' + annotate_1 + "," + annotate_2 + "," + String.format("%E", re.e_value) + "," + String.format("%.4f", re.qvalue) + ",," + re.candidate_num + "," + re.point_count + "," + String.format("%.4f", re.r_square) + "," + String.format("%.4f", re.slope) + "," + String.format("%.4f", re.intercept) + "," + re.start_idx + "," + re.end_idx + "," + String.format("%.4f", re.chain_score_1) + "," + re.chain_rank_1 + "," + String.format("%.4f", re.chain_score_2) + "," + re.chain_rank_2 + "\n");
                    } else {
                        writer.write(re.scan_num + "," + re.spectrum_id + "," + re.spectrum_mz + "," + re.spectrum_mass + "," + re.peptide_mass + "," + re.rt + "," +  re.C13_correction + "," + re.charge + "," + String.format("%.4f", re.score) + "," + re.delta_c + "," + String.format("%.2f", re.ppm) + "," + re.seq_1 + "-" + link_site_1 + "-" + re.seq_2 + "-" + link_site_2 + "," + re.pro_id_1 + "-" + re.pro_id_2 + ',' + annotate_1 + "," + annotate_2 + "," + String.format("%E", re.e_value) + "," + String.format("%.4f", re.qvalue) + "\n");
                    }
                }
            }
            writer.close();
        } catch (IOException ex) {
            logger.error(ex.getMessage());
            System.exit(1);
        }
    }

    private static void saveDecoyResult(List<FinalResultEntry> result, Map<String, String> pro_annotate_map, String id_file_name, boolean is_intra) {
        try {
            BufferedWriter writer;
            if (is_intra) {
                writer = new BufferedWriter(new FileWriter(id_file_name + ".intra.decoy.csv"));
            } else {
                writer = new BufferedWriter(new FileWriter(id_file_name + ".inter.decoy.csv"));
            }

            if (dev) {
                writer.write("scan_num,spectrum_id,spectrum_mz,spectrum_mass,peptide_mass,rt,C13_correction,charge,score,delta_C,ppm,peptide,protein,protein_annotation_1,protein_annotation_2,e_value,,candidate_num,point_num,r_square,slope,intercept,start_idx,end_idx,chain_score_1,chain_rank_1,chain_score_2,chain_rank_2\n");
            } else {
                writer.write("scan_num,spectrum_id,spectrum_mz,spectrum_mass,peptide_mass,rt,C13_correction,charge,score,delta_C,ppm,peptide,protein,protein_annotation_1,protein_annotation_2,e_value\n");
            }
            for (FinalResultEntry re : result) {
                if ((re.hit_type == 1) || (re.hit_type == 2)) {
                    int link_site_1 = re.link_site_1;
                    int link_site_2 = re.link_site_2;

                    String annotate_1 = "DECOY";
                    String annotate_2 = "DECOY";

                    if (!re.pro_id_1.startsWith("DECOY")) {
                        String pro_1 = re.pro_id_1.split(";")[0];
                        annotate_1 = pro_annotate_map.get(pro_1).replace(",", ";");
                    }

                    if (!re.pro_id_2.startsWith("DECOY")) {
                        String pro_2 = re.pro_id_2.split(";")[0];
                        annotate_2 = pro_annotate_map.get(pro_2).replace(",", ";");
                    }

                    if (dev) {
                        writer.write(re.scan_num + "," + re.spectrum_id + "," + re.spectrum_mz + "," + re.spectrum_mass + "," + re.peptide_mass + "," + re.rt + "," + re.C13_correction + "," + re.charge + "," + String.format("%.4f", re.score) + "," + re.delta_c + "," + String.format("%.2f", re.ppm) + "," + re.seq_1 + "-" + link_site_1 + "-" + re.seq_2 + "-" + link_site_2 + "," + re.pro_id_1 + "-" + re.pro_id_2 + ',' + annotate_1 + "," + annotate_2 + "," + String.format("%E", re.e_value) + ",," + re.candidate_num + "," + re.point_count + "," + String.format("%.4f", re.r_square) + "," + String.format("%.4f", re.slope) + "," + String.format("%.4f", re.intercept) + "," + re.start_idx + "," + re.end_idx + "," + String.format("%.4f", re.chain_score_1) + "," + re.chain_rank_1 + "," + String.format("%.4f", re.chain_score_2) + "," + re.chain_rank_2 + "\n");
                    } else {
                        writer.write(re.scan_num + "," + re.spectrum_id + "," + re.spectrum_mz + "," + re.spectrum_mass + "," + re.peptide_mass + "," + re.rt + "," + re.C13_correction + "," + re.charge + "," + String.format("%.4f", re.score) + "," + re.delta_c + "," + String.format("%.2f", re.ppm) + "," + re.seq_1 + "-" + link_site_1 + "-" + re.seq_2 + "-" + link_site_2 + "," + re.pro_id_1 + "-" + re.pro_id_2 + ',' + annotate_1 + "," + annotate_2 + "," + String.format("%E", re.e_value) + "," + String.format("%.4f", re.qvalue) + "\n");
                    }
                }
            }
            writer.close();
        } catch (IOException ex) {
            logger.error(ex.getMessage());
            System.exit(1);
        }
    }

    private static List<List<FinalResultEntry>> pickResult(List<FinalResultEntry> search_result) {
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

    private static void help() {
        String help_str = "ECL 2.0 version " + version + "\r\n"
                + "A cross-linked peptides identification tool with exhaustive searching and linear computational complexity\r\n"
                + "Author: Fengchao Yu\r\n"
                + "Email: fyuab@connect.ust.hk\r\n"
                + "ECL2 usage: java -Xmx25g -jar /path/to/ECL2.jar <parameter_file> <data_file>\r\n"
                + "\t<parameter_file>: parameter file. Can be download along with ECL 2.0.\r\n"
                + "\t<data_file>: spectra data file (mzXML)\r\n"
                + "\texample: java -Xmx32g -jar ECL2.jar parameter.def data.mzxml\r\n";
        System.out.print(help_str);
        System.exit(1);
    }
}
