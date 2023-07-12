#
library(stringr)
library(dplyr)
library(data.table)
library(openxlsx)
library(tidyr)

# Experiment 2
file_path_exp2 <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/"
file_name_exp2 <- "PAL - 2023 01 04 - Phosphopeptides exp2 ( 5 conc 3 reps) DDA_2023-01-06_1032.xlsx"
selected_spcies = "_HUMAN"

# Correct syn. peptide list for Experiment 2
theo_file_path <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/"
theo_file_name <- "Synthetic peptides list_theo_conc_added_080122.xlsx"
sheet_theo_name <- "ISO-ref and OTHER with FC"

# Experiment 3
file_path_exp3 <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_3/Proline_data_analysis/quant_pep_06022023/"
file_name_exp3 <- "PAL _Tcell_Exp3_( 5 conc 3reps)_NoFAIMS_DDA_with_cont_230206.xlsx"
selected_spcies = "_MOUSE"


# Common constant objects
sample_size <- 5
sheet_name <- "Quantified peptide ions"


#source("D:/dev/Desktop_copy/PHD/data_analysis/scripts/ggplot_functions.R")
#source("D:/dev/Desktop_copy/PHD/data_analysis/scripts/with_tidyr_merging_quant_pep_and_pep_list_w_theo_quant_24032023.R")

source("D:/dev/Desktop_copy/PHD/data_analysis/scripts/merging_quant_pep_and_t_cells_theme_edit_29032023.R")
final_proline_pep_quant_analysis_syn(file_path =file_path_exp3,
                                 file_name = file_name_exp3,
                                 sheet_name = sheet_name,
                                 theo_file_path = theo_file_path,
                                 theo_file_name = theo_file_name,
                                 sheet_theo_name = sheet_theo_name,
                                 selected_spcies = "_HUMAN",
                                 sample_size = 5)

final_proline_pep_quant_analysis_bio(file_path =file_path_exp3,
                                 file_name = file_name_exp3,
                                 sheet_name = sheet_name,
                                 theo_file_path = theo_file_path,
                                 theo_file_name = theo_file_name,
                                 sheet_theo_name = sheet_theo_name,
                                 selected_spcies = "_MOUSE",
                                 sample_size = 5)
