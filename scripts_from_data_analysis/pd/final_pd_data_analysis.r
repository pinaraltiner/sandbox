
#

source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/pd/exp3/merging_quant_pep_and_t_cells__for_pd_09032023.R")
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/pd/exp2/pd_exp2_data_processing_new_consider_ISO-REF_fp.R")

file_paths <- c(paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_3/PD_data_analysis/",
                     c(#"percolator_no_FAIMS/",
                       #"percolator_with_FAIMS/",
                       "target_decoy_no_FAIMS/",
                       "target_decoy_with_FAIMS/")),
                paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/PD_data_analysis/",
                       c("target_decoy_no_FAIMS/","percolator_no_FAIMS/",
                         "target_decoy_with_FAIMS/","percolator_with_FAIMS/"))
                )

file_names <- c("Multiconsensus_Exp3_TCells_TargetDecoy_woFAIMS_230531_PeptideIsoforms.txt",
                "Multiconsensus_Exp3_TCells_TargetDecoy_withFAIMS_230531_PeptideIsoforms.txt",
                "Multiconsensus_Exp2_TargetDecoy_woFAIMS_230612_FerriesPhosphoMarkers_PeptideIsoforms.txt",
                "Multiconsensus_Exp2_Percolator_woFAIMS_230612_FerriesPhosphoMarkers_PeptideIsoforms.txt",
                "Multiconsensus_Exp2_TargetDecoy_withFAIMS_230531_PeptideIsoforms.txt",
                "Multiconsensus_Exp2_Percolator_withFAIMS_230601_PeptideIsoforms.txt")


acquisiton_types <- c("DDA Exploris no FAIMS","DDA Exploris with FAIMS")#"DDA TIMS-TOF")

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

experiment_name <- c("E2_A1_R1",
                     "E2_A1_R2",
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


for (i in 1:2){ 
    final_pd_pep_quant_analysis_bio(file_path=file_paths[i],
                                              file_name=file_names[i],
                                              #selected_spcies="MOUSE",
                                              exp_design=experimental_design,
                                              acquisiton_type=acquisiton_types[i],
                                              exp_id=3,
                                              software_name= "Proteome Discoverer",
                                              num_reps=3,
                                              test_type="t.test",
                                              subtitle=file_names[i]) #limma
}

acquisiton_types_exp2 <- c(NA,NA,"DDA Exploris no FAIMS",
                      "DDA Exploris no FAIMS",
                      "DDA Exploris with FAIMS",
                      "DDA Exploris with FAIMS")
fdr_threshold <- 0.05

for (i in 3:length(file_paths)){ 
    final_pd_pep_quant_analysis(file_path=file_paths[i],
                                    file_name=file_names[i],
                                    #sheet_name,
                                    theo_file_path= "D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/",
                                    theo_file_name="Synthetic peptides list_theo_conc_corrected_pool_id_iso_count_final.xlsx",
                                    sheet_theo_name="ISO-ref and OTHER with FC",
                                    background_species <-"ECOLI",
                                    selected_spcies ="PhosphoFerries",
                                    exp_design=experiment_name,
                                    acquisiton_type=acquisiton_types_exp2[i],
                                    exp_id=2,
                                    software_name= "Proteome Discoverer",
                                    num_reps=3,
                                    test_type="limma",
                                fdr_threshold=0.05,
                                actual_ratio=c(2,10,20,100),
                                subtitle=file_names[i]) #t.test
}


