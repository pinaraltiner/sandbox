




source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/DIANN/exp2_DIANN_data_analysis.R")


file_paths <- paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DIA_data_analysis/",c("Experiment_2_redesign/DIANN/",
                                                                                 "experiment_2/DIANN_software_with_monitor-mod_UniMod21/data_analysis_using_fasta_combine_fasta/combine_fasta/"))


mapping_files <- paste0("D:/dev/Pinar/PHD/data_analysis/DIA_mapping/",c("exp2_redesigned_mapping_btw_rawfile_exp_design.txt","exp2_mapping_btw_rawfile_exp_design.txt"))


actual_ratios <- list(new_design = c(2,10,20,60,301),exp2 =c(2,10,20,100))
i<-2

size_variying_pep_size <- c(110,141)

for( i in 1:length(file_paths)){
  final_diann_pep_quant_analysis_syn(file_path =file_paths[i],
                                     file_name = "report.tsv", #"exp2_unimod_report.tsv"
                                     mapping =mapping_files[i],
                                     #exp2_mapping_btw_rawfile_exp_design.txt
                                     #any_LC
                                     #mod_phospho_carb_only
                                     
                                     selected_spcies="HUMAN",
                                     background_species= "ECOLI",
                                     #theo_file_path= "D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/"
                                     #theo_file_name="Synthetic peptides list_theo_conc_corrected_pool_id_iso_count_final.xlsx"
                                     #sheet_theo_name = "ISO-ref and OTHER with FC"
                                     theo_file_path="D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/",
                                     theo_file_name="Synthetic peptides list_theo_conc_corrected_isomericity_new_with_plates.xlsx",
                                     sheet_theo_name ="ISOREF_REF2_Others",#"ISO-refOTHER with FC_correct"
                                     
                                     acquisiton_type="DIA no FAIMS Exploris",
                                     subtitle = "",
                                     fdr_threshold = 0.05,          ### THIS 301 was obtained by introducing pseudo count = (30 + 0.1) / (0 + 0.1)
                                     actual_ratio=as.numeric(unlist(actual_ratios[i])),     #c(2,10,20,60,301),#actual_ratio = c(2,10,20,100) #c(0.5,5,10,30,0)
                                     size_variying_pep =size_variying_pep_size[i],
                                     #exp_design=experiment_name
                                     exp_id=2,
                                     software_name="DIANN",
                                     num_reps=3,
                                     test_type="limma",
                                     numerator = 1)
}


