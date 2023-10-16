library(stringr)
library(dplyr)
library(data.table)
library(openxlsx)
library(tidyr)
library(ggplot2)
library(tidyverse)
###############################################
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/ggplot/ggplot_functions.R")
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/get_modification_func/getModificationPosition_func_for_all_mods.R")

file_path <- "D:/dev/Pinar/PHD/wet_lab_experiments/DIA_data_analysis/experiment_2/DIANN_software_with_monitor-mod_UniMod21/"
file_name <- "exp2_unimod_report.tsv"
mapping <- "D:/dev/Pinar/PHD/data_analysis/DIA_data_processing/mapping_files/exp2_mapping_btw_rawfile_exp_design.txt"



selected_spcies="HUMAN"
background_species= "ECOLI"
theo_file_path= "D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/"
theo_file_name="Synthetic peptides list_theo_conc_corrected_pool_id_iso_count_final.xlsx"
sheet_theo_name = "ISO-ref and OTHER with FC"
acquisiton_type="DIA no FAIMS Exploris"
subtitle = ""
fdr_threshold = 0.05
actual_ratio = c(2,10,20,100)
#exp_design=experiment_name
exp_id=2
software_name="DIANN"
num_reps=3
test_type="limma"

