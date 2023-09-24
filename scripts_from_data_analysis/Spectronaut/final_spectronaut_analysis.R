
#
#TODO: EXP3 for SPECTRONAUT
#source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/Spectronaut/exp3_spectronaut_DIA_data_analysis.R")
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/Spectronaut/exp2_spectronaut_DIA_data_analysis.R")

file_paths <- c(paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DIA_data_analysis/experiment_2/Spectronaut/",
                       c("correct_norm/",
                         "TIMS-TOF")),
                paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DIA_data_analysis/experiment_3/Spectronaut/",
                       c("",
                         "TIMS-TOF"))
                )

#"20230704_162303_OXPAL230118_Exp2woFAIMS_DIA_trypsinP_Report.tsv"
file_names <- c("20230915_Exp2_DIA_noFAIMS_ptm_workflow_correct_norm_output_ptm_peptide_quantification_Report.tsv",
                "20230705_125301_OXPAL230127_Exp3_woFAIMS_DIA_Report.tsv","","")


acquisiton_types <- c("Exploris_no_FAIMS","DIA TIMS-TOF","Exploris_no_FAIMS","DIA TIMS-TOF")

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

##TODO: Create Experimental design in a clever way

experiment_name <- c("E2-A1-R1",
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
                     "E2-A5-R3")
## TODO: DOES not EXIST yet
# 
# for (i in 2){ 
#   final_spectronaut_pep_quant_analysis_bio(file_path=file_paths[i],
#                                   file_name=file_names[i],
#                                   #selected_spcies="MOUSE",
#                                   exp_design=experimental_design,
#                                   acquisiton_type=acquisiton_types[i],
#                                   exp_id=3,
#                                   software_name= "Spectronaut",
#                                   num_reps=3,
#                                   test_type="t.test",
#                                   subtitle=file_names[i]) #limma
# }

fdr_threshold <- 0.05

for (i in 1:(length(file_paths)-3)){ 
  final_spectronaut_pep_quant_analysis_syn(file_path=file_paths[i],
                              file_name=file_names[i],
                              #sheet_name,
                              theo_file_path= "D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/",
                              theo_file_name="Synthetic peptides list_theo_conc_corrected_pool_id_iso_count_final.xlsx",
                              sheet_theo_name="ISO-ref and OTHER with FC",
                              background_species ="ECOLI",
                              selected_spcies ="HUMAN",
                              exp_design=experiment_name,
                              acquisiton_type=acquisiton_types[i],
                              exp_id=2,
                              loc_filter_opt=TRUE,
                              loc_filter=0.75,
                              software_name= "Spectronaut",
                              num_reps=3,
                              test_type="limma",
                              fdr_threshold=0.05,
                              actual_ratio=c(2,10,20,100),
                              subtitle=file_names[i]) #t.test
}
