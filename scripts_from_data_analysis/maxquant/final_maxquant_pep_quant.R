
#

source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/maxquant/exp3/MQ_data_analysis_basic_imputation_edit_ggplot_coloring_without_FAIMS_final_version214_consider_ISOREF_fp.R")


final_maxquant_214_pep_quant_analysis_bio(file_path="D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_3/MQ_214_optimize_param_noFAIMS/",
                                              file_name="evidences.txt",
                                              #selected_spcies="MOUSE",
                                              exp_design=experimental_design,
                                              acquisiton_type="Exploris DDA no FAIMS",
                                              exp_id=3,
                                              software_name="MaxQuant v2.1.4",
                                              num_reps=3,
                                              test_type="limma")

final_maxquant_214_pep_quant_analysis_syn(file_path =file_path_exp3_tims,
                                     file_name = file_name_exp3_tims,
                                     sheet_name = sheet_name,
                                     theo_file_path = theo_file_path,
                                     theo_file_name = theo_file_name,
                                     sheet_theo_name = sheet_theo_name,
                                     selected_spcies = "_MOUSE",
                                     sample_size = 5)


file_path <- "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_3/MQ_214_optimize_param_noFAIMS/"
#"D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/MQ_214_optimize_param_noFAIMS/"
#file_path <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/MQ_data_analysis/with_FAIMS_MBR/min_mod_score_0/"
#file_path <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/MQ_data_analysis/exp2_wo_FAIMS/rerun_using_1634/"

acquisiton_type <- c("DDA Exploris with FAIMS","DDA Exploris no FAIMS","DDA TIMS-TOF")
software_name <- c("Proline", "MaxQuant", "PD")
experimental_design <- c("E3-A1-R1","E3-A1-R2",
                "E3-A1-R3",
                "E3-A2-R1",
                "E3-A2-R2",
                "E3-A2-R3",
                "E3-A3-R1",
                "E3-A3-R2",
                "E3-A3-R3",
                "E3-A4-R1",
                "E3-A4-R2",
                "E3-A4-R3",
                "E3-A5-R1",
                "E3-A5-R2",
                "E3-A5-R3")

output_path <- "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_3/MQ_214_optimize_param_noFAIMS/figures/"
# experiment_name <- c("E2-M1-R1",
#                      "E2-M1-R2",
#                      "E2-M1-R3",
#                      "E2-M2-R1",
#                      "E2-M2-R2",
#                      "E2-M2-R3",
#                      "E2-M3-R1",
#                      "E2-M3-R2",
#                      "E2-M3-R3",
#                      "E2-M4-R1",
#                      "E2-M4-R2",
#                      "E2-M4-R3",
#                      "E2-M5-R1",
#                      "E2-M5-R2",
#                      "E2-M5-R3")

