#
library(stringr)
library(dplyr)
library(data.table)
library(openxlsx)
library(tidyr)

#
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/proline/exp2/proline_exp2_with_faims_data_processing_basic_imputation_edit_coloring_ggplot_current_consider_ISO-REF_fp.R")

source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/proline/exp3/merging_quant_pep_and_t_cells_theme_edit_29032023.R")

file_paths <- c(paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_3/Proline_data_analysis/", c("quant_pep_06022023/","/with_FAIMS/")),
                "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_3/Proline_timsdata/",
                paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/", c("exp2_re_injection/","with_FAIMS/")),
                "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_timsdata/")


file_names <- c("PAL _T_cell_Exp3_( 5 conc 3reps)_NoFAIMS_DDA_with_cont_230206_2023-02-07_0947.xlsx",
                "PAL _T_cell_Exp3_( 5 conc 3reps)_withFAIMS_DDA_25052023_noDesign_2023-07-10_1134.xlsx",
                "PAL TimsTOF_data_conversion_all_exp_06062023 E3 corrected_2023-06-28_1522.xlsx",
                "PAL _Phosphopeptides exp2 ( 5 conc 3reps) DDA_230117 with Design_2023-06-07_1003.xlsx",
                "PAL _Exp2_( 5 conc 3reps)_withFAIMS_DDA_25052023_noDesign - correct_2023-06-13_1514.xlsx",
                "PAL TimsTOF_data_conversion_all_exp_06062023 E2 corrected_2023-06-28_1522.xlsx")




acquisiton_types <- c("DDA Exploris no FAIMS","DDA Exploris with FAIMS","DDA TIMS-TOF",
                      "DDA Exploris no FAIMS","DDA Exploris with FAIMS","DDA TIMS-TOF")
#test_type <- c("t.test","wilcoxon","limma")
#software_name <- c("Proline", "MaxQuant", "PD")
#exp_id <- 1:3
experimental_design <- c("E3_A1_R1","E3_A1_R2",
                         "E3_A1_R3",
                         "E3_A2_R1",
                         "E3_A2_R2",
                         "E3_A2_R3",
                         "E3_A3_R1",
                         "E3_A3_R2",
                         "E3_A3_R3",
                         "E3_A4_R1",
                         "E3_A4_R2",
                         "E3_A4_R3",
                         "E3_A5_R1",
                         "E3_A5_R2",
                         "E3_A5_R3")

#suppressPackageStartupMessages() = To Eliminate Package Startup Messages

for ( i in 1:3){
  final_proline_pep_quant_analysis_bio(file_path =file_paths[i],
                                       file_name = file_names[i],
                                       sheet_name = "Best PSM from protein sets",
                                       #sheet_theo_name <- "ISO-ref and OTHER with FC",
                                       selected_spcies = "_MOUSE",
                                       exp_design = experimental_design,
                                       acquisiton_type= acquisiton_types[i],
                                       num_reps = 3,
                                       exp_id=3,
                                       software_name = "Proline",
                                       test_type="limma",
                                       subtitle=file_names[i])
}



# Correct syn. peptide list for Experiment 2


# Experiment 2

experiment_name <- c("E2_A1_R1","E2_A1_R2",
                "E2_A1_R3",
                "E2_A2_R1",
                "E2_A2_R2",
                "E2_A2_R3",
                "E2_A3_R1",
                "E2_A3_R2",
                "E2_A3_R3",
                "E2_A4_R1",
                "E2_A4_R2",
                "E2_A4_R3",
                "E2_A5_R1",
                "E2_A5_R2",
                "E2_A5_R3")


for (i in 4:length(file_paths)){
    final_proline_pep_quant_analysis_syn(file_path =file_paths[i],
                                         file_name = file_names[i],
                                         sheet_name = "Best PSM from protein sets",
                                         theo_file_path= "D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/",
                                         theo_file_name="Synthetic peptides list_theo_conc_corrected_pool_id_iso_count_final.xlsx",
                                         sheet_theo_name = "ISO-ref and OTHER with FC",
                                         acquisiton_type= acquisiton_types[i],
                                         software_name = "Proline",
                                         test_type="limma",
                                         selected_spcies = "_HUMAN",
                                         background_species = "ECOLI",
                                         exp_id=2,
                                         num_reps=3,
                                         exp_design=experiment_name,
                                         fdr_threshold=0.05,
                                         actual_ratio=c(2,10,20,100),
                                         subtitle=file_names[i])
}








