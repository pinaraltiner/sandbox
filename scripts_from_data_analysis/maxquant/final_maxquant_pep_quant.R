
#

source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/maxquant/exp2/MQ_data_analysis_basic_imputation_edit_ggplot_coloring_without_FAIMS_final_version214_consider_ISOREF_fp.R")
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/maxquant/exp3/MQ_data_analysis_basic_imputation_edit_ggplot_coloring_without_FAIMS_final_version214_consider_ISOREF_fp.R")


file_paths <- c(paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_3/", c("MQ_214_optimize_noFAIMS/","MQ_214_optimize_withFAIMS/","MQ_timstof/")),
                paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/",c("MQ_214_optimize_noFAIMS/","MQ_214_optimize_withFAIMS/","MQ_timstof/")),
                "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/redesigned_exp_2/MaxQuant/")

acquisiton_types <- c("DDA Exploris no FAIMS",
                      "DDA Exploris with FAIMS",
                      "DDA TIMS-TOF",
                      "DDA Exploris no FAIMS",
                      "DDA Exploris with FAIMS",
                      "DDA TIMS-TOF",
                      "DDA Exploris no FAIMS")

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

for ( i in 1:3){ #1:(length(file_paths))[1:3]
  final_maxquant_214_pep_quant_analysis_bio(file_path=file_paths[i],
                                              file_name="evidences.txt",
                                              #selected_spcies="MOUSE",
                                              exp_design=experimental_design,
                                              acquisiton_type=acquisiton_types[i],
                                              exp_id=3,
                                              software_name="MaxQuant v2.1.4",
                                              num_reps=3,
                                              test_type="limma")}


# Experiment 2
experiment_name <- c("E2-A1-R1","E2-A1-R2",
                     "E2-A1-R3",
                     "E2-A2-R1",
                     "E2-A2-R2",
                     "E2-A2-R3",
                     "E2-A3-R1",
                     "E2-A3-R2",
                     "E2-A3-R3",
                     "E2-A4-R1",
                     "E2-A4-R2",
                     "E2-A4-R3",
                     "E2-A5-R1",
                     "E2-A5-R2",
                     "E2-A5-R3")

for ( i in 4:6){#(length(file_paths))){ 
  final_mq_pep_quant_analysis_syn(file_path=file_paths[i],
                                              file_name="evidence.txt",
                                              selected_spcies="HUMAN",
                                              background_species= "ECOLI",
                                              theo_file_path= "D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/",
                                              theo_file_name="Synthetic peptides list_theo_conc_corrected_isomericity_new_with_plates.xlsx",#"Synthetic peptides list_theo_conc_corrected_pool_id_iso_count_final.xlsx",
                                              sheet_theo_name = "ISO-refOTHER with FC_correct", #"ISO-ref and OTHER with FC",
                                              acquisiton_type=acquisiton_types[i],
                                              subtitle = "",
                                              loc_filter_opt=TRUE,
                                              loc_filter=0.75,
                                              fdr_threshold = 0.05,
                                              actual_ratio = c(2,10,20,100),
                                              exp_design=experiment_name,
                                              exp_id=2,
                                              size_variying_pep_size =110,
                                              software_name="MaxQuant v2.1.4",
                                              num_reps=3,
                                              test_type="limma",
                                              numerator=1)
}


##################### EXPERIMENT 2 REDESIGNED DATASET ###########################

#for ( i in 4:(length(file_paths))
i <- 7
experiment_name <-c("E2-A1-R1",
                    "E2-A1-R2",
                    "E2-A1-R3",
                    "E2-A2-R1",
                    "E2-A2-R2",
                    "E2-A2-R3",
                    "E2-A3-R1",
                    "E2-A3-R2",
                    "E2-A3-R3",
                    "E2-A4-R1",
                    "E2-A4-R2",
                    "E2-A4-R3",
                    "E2-A5-R1",
                    "E2-A5-R2",
                    "E2-A5-R3",
                    "E2-A6-R1",
                    "E2-A6-R2",
                    "E2-A6-R3")

final_mq_pep_quant_analysis_syn(file_path=file_paths[i],
                                file_name="evidence.txt",
                                selected_spcies="HUMAN",
                                background_species= "ECOLI",
                                theo_file_path= "D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/",
                                theo_file_name="Synthetic peptides list_theo_conc_corrected_isomericity_new_with_plates.xlsx", #"Synthetic peptides list_theo_conc_corrected_isomericity.xlsx",
                                sheet_theo_name="ISOREF_REF2_Others", #"ISO-ref and OTHER with FC",
                                acquisiton_type=acquisiton_types[i],
                                subtitle = "",
                                loc_filter_opt=TRUE,
                                loc_filter=0.75,
                                fdr_threshold = 0.05,
                                actual_ratio=c(2,10,20,60,301),#300001),##actual_ratio = c(2,10,20,100)
                                size_variying_pep_size =110,
                                exp_design=experiment_name,
                                exp_id=2,
                                software_name="MaxQuant v2.1.4",
                                num_reps=3,
                                test_type="limma",
                                numerator=1
                                
)



##################### EXPERIMENT 1 ###########################
  ### TODO: FIND LOGICAL WAY TO RUN THE FUNCTION ###


source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/maxquant/MQ_analysis_for_experiment_1.R")
file_path = "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/Exp_1/MQ2.1.4_def_params/"
all_dirs <- list.files(file_path)



final_MQ_pep_quant_analysis_exp1(file_path = "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/Exp_1/MQ2.1.4_def_params/",
                                 #file_name,
                                 curr_dir=all_dirs[1],
                                 theo_file_path="D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/",
                                 theo_file_name="Synthetic peptides list_theo_conc_corrected_isomericity_new_with_plates.xlsx",
                                 sheet_theo_name ="ISO-refOTHER with FC_correct",
                                 selected_species = "HUMAN",
                                 background_species ="ECOLI",
                                 exp_id=1,
                                 #exp_design,
                                 software_name="MaxQuant",
                                 acquisiton_type= "DDA_with_FAIMS")

final_MQ_pep_quant_analysis_exp1(file_path = "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/Exp_1/MQ2.1.4_def_params/",
                                 #file_name,
                                 curr_dir=all_dirs[2],
                                 theo_file_path="D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/",
                                 theo_file_name="Synthetic peptides list_theo_conc_corrected_isomericity_new_with_plates.xlsx",
                                 sheet_theo_name ="ISO-refOTHER with FC_correct",
                                 selected_species = "HUMAN",
                                 background_species ="",
                                 exp_id=1,
                                 #exp_design,
                                 software_name="MaxQuant",
                                 acquisiton_type= "DDA_TIMSTOF")


final_MQ_pep_quant_analysis_exp1(file_path = "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/Exp_1/MQ2.1.4_def_params/",
                                 #file_name,
                                 curr_dir=all_dirs[3],
                                 theo_file_path="D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/",
                                 theo_file_name="Synthetic peptides list_theo_conc_corrected_isomericity_new_with_plates.xlsx",
                                 sheet_theo_name ="ISO-refOTHER with FC_correct",
                                 selected_species = "HUMAN",
                                 background_species ="ECOLI",
                                 exp_id=1,
                                 #exp_design,
                                 software_name="MaxQuant",
                                 acquisiton_type= "DDA_with_FAIMS")



final_MQ_pep_quant_analysis_exp1(file_path = "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/Exp_1/MQ2.1.4_def_params/",
                                 #file_name,
                                 curr_dir=all_dirs[4],
                                 theo_file_path="D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/",
                                 theo_file_name="Synthetic peptides list_theo_conc_corrected_isomericity_new_with_plates.xlsx",
                                 sheet_theo_name ="ISO-refOTHER with FC_correct",
                                 selected_species = "HUMAN",
                                 background_species ="",
                                 exp_id=1,
                                 #exp_design,
                                 software_name="MaxQuant",
                                 acquisiton_type= "DDA_with_FAIMS")


final_MQ_pep_quant_analysis_exp1(file_path = "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/Exp_1/MQ2.1.4_def_params/",
                                 #file_name,
                                 curr_dir=all_dirs[5],
                                 theo_file_path="D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/",
                                 theo_file_name="Synthetic peptides list_theo_conc_corrected_isomericity_new_with_plates.xlsx",
                                 sheet_theo_name ="ISO-refOTHER with FC_correct",
                                 selected_species = "HUMAN",
                                 background_species ="ECOLI",
                                 exp_id=1,
                                 #exp_design,
                                 software_name="MaxQuant",
                                 acquisiton_type= c("DDA_no_FAIMS"))

final_MQ_pep_quant_analysis_exp1(file_path = "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/Exp_1/MQ2.1.4_def_params/",
                                 #file_name,
                                 curr_dir=all_dirs[6],
                                 theo_file_path="D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/",
                                 theo_file_name="Synthetic peptides list_theo_conc_corrected_isomericity_new_with_plates.xlsx",
                                 sheet_theo_name ="ISO-refOTHER with FC_correct",
                                 selected_species = "HUMAN",
                                 background_species ="",
                                 exp_id=1,
                                 #exp_design,
                                 software_name="MaxQuant",
                                 acquisiton_type= c("DDA_no_FAIMS"))




