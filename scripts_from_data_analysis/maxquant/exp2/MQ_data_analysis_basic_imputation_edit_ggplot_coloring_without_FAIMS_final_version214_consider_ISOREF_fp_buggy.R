library(stringr)
library(dplyr)
#library(data.table) # REMOVED: unused
#library(openxlsx)
library(ggplot2)
library(reshape2)
library(readr)
library(tibble)
library(purrr)
library(readr)
library(tidyr)
library(openxlsx)
library(patchwork)
###############################################
## MAXQUANT DATA ANALYSIS
### Source code was taken from here: https://rdrr.io/github/singjc/mstools/src/R/getModificationPosition.R
## The code was modified based on what I want and based on software input tyoe (DIANN-and MaxQuant)
#source("D:/dev/Desktop_copy/PHD/data_analysis/scripts/getModificationPositionMQ_func_edit_v1_1634.R")
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/get_modification_func/getModificationPosition_general_change_condition_current_mod_sequence_MQ_Spectronaut.R")
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/ggplot/ggplot_functions.R") #ggplot_functions.R
#source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/roc_curve/roc_curve_generation_proline_edit.R")
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/roc_curve/new_roc_curve_generation_with_custom_threshold.R")

final_mq_pep_quant_analysis_syn <- function(file_path,
                                                 file_name,
                                                 #sheet_name,
                                                 theo_file_path,
                                                 theo_file_name,
                                                 sheet_theo_name,
                                                 background_species,
                                                 selected_species,
                                            norm_type,
                                                 exp_id,
                                                 size_variying_pep_size,
                                                 numerator,
                                                 loc_filter_opt,
                                                 loc_filter,
                                                 exp_design,
                                                 fdr_threshold,
                                                 acquisiton_type,
                                                 software_name,
                                                 test_type,
                                                 num_reps,
                                                 actual_ratio,
                                                 subtitle,
                                            dir_name){
  
  
  intended_dir <-paste0(file_path,dir_name,"/")#"outputs_imputed_5thpercent_one_sample_nonNA_fin")
  
  if(dir.exists(intended_dir)){
    new_path <- intended_dir
    
  }else{
    dir.create(intended_dir)
    new_path <- list.dirs(intended_dir)
    
  }
  
  
  sample_size <- length(exp_design) / num_reps
  sample_names <- paste0("A",1:sample_size)
  comparisons <- NULL
  
  for (i in 1:sample_size){
    tmp <- paste0(sample_names[numerator], "/",sample_names[i])
    comparisons[i] <- tmp
    rm(tmp)
  }
  #comparisons <- comparisons[-1]
  comparisons <- comparisons[-numerator]
  
  quant_peptides <- read_tsv(paste0(file_path,file_name))
  
  # pep_list_w_theo <- read.xlsx(paste0(theo_file_path, theo_file_name), sheet = sheet_theo_name)
  # pep_list_w_theo_quant <- pep_list_w_theo[,-1]
  # 
  # common_col_theo_quant <- as.data.frame(paste(pep_list_w_theo_quant$Phosphopeptide.sequence,
  #                                              pep_list_w_theo_quant$modified.position.in.peptide, sep = "_"))
  # colnames(common_col_theo_quant) <- "pep_with_pos"
  # pep_list_w_theo_quant_new <- cbind(common_col_theo_quant,pep_list_w_theo_quant)
  # 
  
  ## THEORETICAL PEPTIDE LIST
  pep_list_w_theo <- read.xlsx(paste0(theo_file_path, theo_file_name), sheet = sheet_theo_name)
  #pep_list_w_theo_quant <- pep_list_w_theo[,-1]
  
  if(sheet_theo_name == "ISO-refOTHER with FC_correct") {
    common_col_theo_quant <- as.data.frame(paste(pep_list_w_theo$Phosphopeptide.sequence,
                                                 pep_list_w_theo$modified.position.in.peptide, sep = "_"))
    colnames(common_col_theo_quant) <- "pep_with_pos"
    pep_list_w_theo <- cbind(common_col_theo_quant,pep_list_w_theo)
    
  } else if (sheet_theo_name == "ISOREF_REF2_Others"){
    
  }else{
    print("Please check the sheet name of the theo. peptide list.")
  }
  
  
  #################################################
  pep_list_w_theo_unique <- pep_list_w_theo %>% 
    select(Phosphopeptide.sequence,Pool) %>%
    distinct(Phosphopeptide.sequence,.keep_all = TRUE) %>%
    rename(Sequence = Phosphopeptide.sequence) %>%
    rename(Pool_for_seq_merge=Pool) %>%
    mutate(situation="Correct")
  
  
  pep_list_w_theo_seq_map <- pep_list_w_theo %>% select(Sequence,Neutral.mass.SH.Cys) %>%
    rename(Neutral_mass=Neutral.mass.SH.Cys) %>%
    mutate(Neutral_mass = round(Neutral_mass, 0)) %>%
    mutate(Neutral_mass_theo=Neutral_mass)
  
  pep_list_w_theo_pep_map <- pep_list_w_theo %>%
    select(pep_with_pos,Pool,isomericity)
  
  ### IMPUTATION
  # Calculate 1 percent quantile of each sample
  
   imputed_values <- quant_peptides %>% 
     group_by(Experiment) %>% 
     summarise(fifth_quantile=quantile(Intensity,probs=0.05,na.rm=TRUE))
   
   imputed_values_vec <- as.vector(imputed_values$fifth_quantile)
  
  # # # Impute missing values
  # 
  # abundances_for_impute <- quant_peptides %>% select(Experiment,Intensity) %>%
  #     group_by(Experiment) %>%
  #     mutate(imputed_intensity=case_when(grepl("A1-R1",Experiment) ~imputed_values_vec[1],
  #                                        grepl("A1-R2",Experiment) ~imputed_values_vec[2],
  #                                        grepl("A1-R3",Experiment) ~imputed_values_vec[3],
  #                                        grepl("A2-R1",Experiment) ~imputed_values_vec[4],
  #                                        grepl("A2-R2",Experiment) ~imputed_values_vec[5],
  #                                        grepl("A2-R3",Experiment) ~imputed_values_vec[6],
  #                                        grepl("A3-R1",Experiment) ~imputed_values_vec[7],
  #                                        grepl("A3-R2",Experiment) ~imputed_values_vec[8],
  #                                        grepl("A3-R3",Experiment) ~imputed_values_vec[9],
  #                                        grepl("A4-R1",Experiment) ~imputed_values_vec[10],
  #                                        grepl("A4-R2",Experiment) ~imputed_values_vec[11],
  #                                        grepl("A4-R3",Experiment) ~imputed_values_vec[12],
  #                                        grepl("A5-R1",Experiment) ~imputed_values_vec[13],
  #                                        grepl("A5-R2",Experiment) ~imputed_values_vec[14],
  #                                        grepl("A5-R3",Experiment) ~imputed_values_vec[15],
  #                                        TRUE ~ Intensity))
  # 
  ### LOCALIZATION THRESHOLD ###
  # quant_phospho <- quant_peptides %>%
  #   #select(-Experiment,Intensity) %>% 
  #   #mutate(abundances_for_impute) %>% 
  #   #filter(grepl(9606,`Taxonomy IDs`))%>%
  #   filter(grepl(selected_species,Proteins)) %>% 
  #   filter(grepl("Phospho",Modifications)) %>% drop_na(`Modified sequence`)
   
  # # IN THIS STEP PHOSPHO-PROBABILITIES WERE EXTRACTED FROM ITS ASSOCIATED COLUMN 
  # # AND SELECTED ONLY PHOSPHO-HUMAN PEPTIDES
  quant_phospho <- quant_peptides %>% 
    filter(grepl("Phospho",Modifications)) %>%
    filter(grepl(selected_species,Proteins) & !grepl("CON__",Proteins)) %>%
    mutate(species=selected_species) %>%
    mutate(
      extracted_values = sapply(str_extract_all(`Phospho (STY) Probabilities`, "\\(\\d+(\\.\\d+)?\\)"), function(x) {
        values <- as.numeric(str_extract_all(x, "\\d+(\\.\\d+)?"))
        if (length(values) == 0 || all(is.na(values))) NA else max(values)
      })
      
      #select(Sequence,Modifications,`Modified sequence`,Proteins,`Phospho (STY) Probabilities`) %>%
      #drop_na(`Phospho (STY) Probabilities`) %>%
      # mutate(
      #   extracted_values = sapply(str_extract_all(`Phospho (STY) Probabilities`, "\\(\\d+\\.\\d+\\)"), function(x) {
      #     values <- as.numeric(str_extract_all(x, "\\d+")) #.\\d+
      #     if (length(values) == 0) NA else max(values)
      #   })
    )%>% relocate(extracted_values,.after = Sequence) %>% #drop_na(extracted_values) %>%
    filter(grepl(selected_species,Proteins) & grepl("Phospho",Modifications)) %>%
    drop_na(`Modified sequence`)
  
  ## PHOSPHO-POSITION EXTRACTION ##
  df2 <- apply(X = as.data.frame(quant_phospho[,"Modified sequence"]),1,function(x){getModificationPosition_general(mod_seq = x,software_name = "MQ_214")})
  
  results1 <- map_dfr(df2, ~ enframe(.x)) %>%
    filter(grepl("modification_",name)| grepl("pep_seq", name)) %>%
    mutate(value = map_chr(value, str_c, collapse="&")) %>%
    mutate(mods=case_when(grepl("Phospho (STY)",fixed = T,name) ~ "phospho",
                          grepl("Oxidation (M)",fixed = T,name) ~ "Oxidation",
                          grepl("(Acetyl (Protein N-term))",fixed = T,name) ~ "N-term_Acetyl",
                          TRUE ~ ""))
  
  ## Adding indeces to use as pep-seq info
  results_with_index <- results1 %>%
    mutate(id = cumsum(name == "pep_seq")) 
  
  ## Creating a new object to combine everything;
  reshaped_results <- results1 %>% 
    ## ADDING INDEX
    mutate(id = cumsum(name == "pep_seq")) %>%
    ## REMOE rows contains "PEP_SEQ"
    filter(name != "pep_seq") %>%
    ## GROUPING
    group_by(id) %>%
    ## MERGING ALL MODS, POSITIONS, and their unimod id 
    ## ADDING "name" IS OPTIONAL 
    mutate(mods = paste(value, mods, collapse = "__")) %>% #name
    ## USING INITIAL INDECES, JOINING WILL BE DONE
    left_join(filter(results_with_index, name == "pep_seq"), by = "id") %>% 
    ungroup() %>%
    ## SELECTING USEFUL COLUMNS
    select(c(name.x,value.x,mods.x,value.y))
  
  
  result_with_common_col <- reshaped_results %>%
    filter(grepl("modification_(Phospho (STY))",fixed = T,name.x)) %>%
    mutate(pep_with_pos = paste(value.y,value.x, sep = "_")) %>%
    select(!value.y) %>%
    rename("mods_id" = "name.x",
           "phospho_positions" ="value.x",
           "all_mods_with_mod_type" = "mods.x")
  #getModificationPosition_MQ(mod_seq = id_syn_pep$`Modified sequence`,F)
  
  final_results_with_common_col <- mutate(result_with_common_col,quant_phospho)
  
  barplt_df <- final_results_with_common_col %>%
    select(Sequence,Experiment,Intensity,pep_with_pos,`Modified sequence`,Proteins,extracted_values) %>%
    separate(Experiment, into = c("Exp_id","Sample_id","Rep_id"),sep = "-",remove = F) %>%
    #mutate(sample_rep_id_seq = paste(pep_with_pos, Sample_id,Rep_id, sep = "_")) %>%
    group_by(pep_with_pos,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    slice(which.max(Intensity)) %>%
    ungroup() %>%
    separate(pep_with_pos,into = c("seq","Positions"),sep = "_",remove = F)
  
  ### Mapping by precursor mass + phospho  and sequence ####
  amino_acid_table <- read.delim("D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/amino_acid_table.txt")
  mapping_df <- as.data.frame(cbind(amino_acid_table$X1.letter.code,amino_acid_table$Monoisotopic.Mass))
  colnames(mapping_df) <- c("letters","mono_isotopic")
  
  mapping_df$mono_isotopic <- as.numeric(mapping_df$mono_isotopic)
  # Apply the function to each sequence in the sequences data frame
  barplt_df$Sum <- mapply(calculate_sum, barplt_df$Sequence, barplt_df$Positions, MoreArgs = list(mapping_df = mapping_df))
  
  barplt_df_mapped_seq <- barplt_df %>% 
    mutate(Neutral_mass = round(Sum, 0)) %>%
    mutate(Neutral_mass_res=Neutral_mass) %>%
    full_join(pep_list_w_theo_seq_map,by=c("Neutral_mass","Sequence")) %>% # #pep_with_pos
    
    mutate(map_seq = ifelse(is.na(Neutral_mass_res), "missing", NA)) %>%
    mutate(map_seq = ifelse(is.na(Neutral_mass_theo), "wrong seq", map_seq)) %>%
    mutate(map_seq = ifelse(!is.na(Neutral_mass_res) & !is.na(Neutral_mass_theo) & Neutral_mass_res == Neutral_mass_theo, 
                            "correct seq", 
                            map_seq)) %>%
    filter(grepl("correct seq",map_seq)) %>%
    select(!c(Positions,seq,Sum)) %>%
    group_by(pep_with_pos,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    slice(which.max(Intensity)) %>%
    ungroup() 
  
  barplt_df_mapped_pep <- barplt_df %>% 
    mutate(Neutral_mass = round(Sum, 0)) %>%
    mutate(Neutral_mass_res=Neutral_mass) %>%
    full_join(pep_list_w_theo_seq_map,by=c("Neutral_mass","Sequence")) %>% # #pep_with_pos
    
    mutate(map_seq = ifelse(is.na(Neutral_mass_res), "missing", NA)) %>%
    mutate(map_seq = ifelse(is.na(Neutral_mass_theo), "wrong seq", map_seq)) %>%
    mutate(map_seq = ifelse(!is.na(Neutral_mass_res) & !is.na(Neutral_mass_theo) & Neutral_mass_res == Neutral_mass_theo, 
                            "correct seq", 
                            map_seq)) %>%
    filter(grepl("correct seq",map_seq)) %>%
    left_join(pep_list_w_theo_pep_map,by="pep_with_pos") %>%
    mutate(Pool=ifelse(Pool=="Diluted",paste0(Pool,"_",isomericity),Pool)) %>%
    mutate(map_loc=ifelse(is.na(Pool),"Wrong Localization","Correct")) %>%
    mutate(Pool=ifelse(Pool=="Diluted_nonisomeric","Diluted_non-isomeric",Pool)) %>%
    select(!c(Positions,seq,Sum)) %>%
    group_by(pep_with_pos,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    slice(which.max(Intensity)) %>%
    ungroup() 
  #distinct(pep_with_pos,.keep_all = T) %>%
  #count(map_loc)
  
  # site_prob <- barplt_df_mapped_pep %>% 
  #   select(extracted_values, pep_with_pos,Experiment) %>% 
  #   pivot_wider(names_from = "Experiment",
  #               values_from = "extracted_values") %>% 
  #   relocate(exp_design) %>%
  #   rowwise() %>%
  #   mutate(max_value = max(c_across(all_of(exp_design)), na.rm = TRUE)) %>%
  #   select(!exp_design)
  #mutate(max_value=ifelse(!is.numeric(max_value),NA,max_value))
  #filter(!grepl(-Inf,max_value))

 #"/MQ_transposed_data_all_pep_list.txt")
  
  #######################################################################
  ecoli_seq <- quant_peptides %>% 
    select(Sequence,Modifications,Proteins) %>%
    filter(!grepl(selected_species, Proteins) & !grepl("CON__", Proteins)) %>%
    #pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
    #filter(grepl(background_species,Proteins)) %>%
    #distinct(Sequence,.keep_all = T) %>%
    mutate(species=background_species) 
  
  ecoli_seq_dist <- quant_peptides %>% 
    select(Sequence,Modifications,Modifications,Proteins) %>%
    filter(!grepl(selected_species, Proteins) & !grepl("CON__", Proteins)) %>%
    #pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
    #filter(grepl(background_species,Proteins)) %>%
    distinct(Sequence,.keep_all = T) %>%
    mutate(species=background_species) 
  
   #################################################
  all_seq <- barplt_df_mapped_seq %>% 
    select(Sequence,
           pep_with_pos,
           `Modified sequence`,map_seq,
           extracted_values,
           Experiment,Sample_id,Rep_id) %>%
    
    #full_join(pep_list_w_theo_unique,by="Sequence") %>%
    #mutate_at("Pool_for_seq_merge", ~replace_na(.,"Unexpected")) %>%
    #mutate(Pool_for_seq_merge= ifelse(is.na(species),"missing",Pool_for_seq_merge)) %>%
    #filter(!grepl("Unexpected",Pool_for_seq_merge)) %>%
    bind_rows(ecoli_seq) %>% 
    rename(Pool_for_seq_merge=map_seq) %>%
    mutate(Pool_for_seq_merge= ifelse(is.na(Pool_for_seq_merge),background_species,Pool_for_seq_merge)) %>%
    mutate(species=is.na(species),selected_species,species)
  
  
  write.table(all_seq, file=paste0(new_path,"/Experiment2",software_name,"_number_of_unique_sequence_for_each_species.txt"),sep = "\t",col.names = T,row.names = F)
  
  ########## ########## ########## ########## ########## ########## ########## ##########
  all_seq_syn <- all_seq %>%
    filter(!grepl(background_species,Pool_for_seq_merge)) %>%
    select(Sequence,`Modified sequence`,Proteins,Pool_for_seq_merge,extracted_values,Experiment,Sample_id,Rep_id) %>% 
    
    #bind_rows(ecoli_seq_dist) %>%
    #mutate(Pool_for_seq_merge= ifelse(is.na(Pool_for_seq_merge),background_species,Pool_for_seq_merge)) %>%
    mutate(acq_type=acquisiton_type) %>%
    mutate(soft_name=software_name) %>%
    group_by(Experiment) %>%
    distinct(Sequence, .keep_all = T)
  
  write.table(all_seq_syn, file=paste0(new_path,"/Number_of_human_phospho_sequences_",
                                       software_name,"_Experiment",exp_id,".txt"),sep = "\t",col.names = T,row.names = F)
  #paste0(new_path,"/",software_name,"_num_unique_seq_each_species.txt")
  
  all_pep_syn <- barplt_df_mapped_pep %>%
    #filter(!grepl(background_species,Pool_for_seq_merge)) %>%
    select(pep_with_pos, `Modified sequence`,map_loc,extracted_values,Experiment,Sample_id,Rep_id) %>%
    rename(Pool_for_pep_merge=map_loc) %>%
    #distinct(Sequence, .keep_all = T) %>%
    #bind_rows(ecoli_seq_dist) %>%
    #mutate(Pool_for_seq_merge= ifelse(is.na(Pool_for_seq_merge),background_species,Pool_for_seq_merge)) %>%
    mutate(acq_type=acquisiton_type) %>%
    mutate(soft_name=software_name) %>% 
    group_by(Experiment) %>%
    distinct(pep_with_pos,.keep_all = T)
  
  all_corr_pep_syn <- all_pep_syn %>% filter(grepl("Correct",Pool_for_pep_merge))
  
  # all_corr_seq <- quant_peptides %>% 
  #   filter(grepl("Phospho",Modifications)) %>%
  #   filter(grepl(selected_species,Proteins) & !grepl("CON__",Proteins)) %>%
  #   mutate(species=selected_species) %>% 
  #   group_by(Experiment) %>%
  #   distinct(Sequence,.keep_all = T) %>%
  #   full_join(pep_list_w_theo_unique,by="Sequence") %>%
  #   mutate_at("situation", ~replace_na(.,"Unexpected")) %>%
  #   mutate(situation= ifelse(is.na(species),"missing",situation)) %>%
  #   filter(!grepl("Unexpected",situation) & !grepl("missing",situation) ) %>%
  #   select(Sequence,Modifications,Experiment,situation,Pool_for_seq_merge) %>%
  # separate(Experiment, into = c("exp_id","samp_id","rep_id"),sep = "-")
  
  plot15 <- gg_barplt_id_pep_count(data_set = all_seq_syn,
                                   x_df = all_seq_syn$Sample_id,
                                   fill_df = all_seq_syn$Rep_id,
                                   ymax = nrow(all_seq_syn),
                                   size_num=10,
                                   header = paste("Total number of correctly identified ",selected_species,"phospho-sequence","across each sample",sep=" "),
                                   caption_lab = "Mapping was done without considering phospho-positions.",
                                   x_lab = "Sample id",
                                   fill_lab =  "Sample id",
                                   y_lab = "Number of identified peptides",
                                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) +
    theme(axis.text.x = element_text(angle = 90))
  
  write.table(all_seq_syn, file=paste0(new_path,"/Number_of_human_phospho_sequences_",
                                       software_name,"_Experiment",exp_id,".txt"),
              sep = "\t",col.names = T,row.names = F)
  ########## ########## ########## ########## ########## ########## ########## ##########
  
  plot12 <- gg_barplt_id_pep_count(data_set = all_corr_pep_syn,
                                   x_df = all_corr_pep_syn$Sample_id,
                                   fill_df = all_corr_pep_syn$Rep_id,
                                   ymax = length(unique(all_corr_pep_syn$pep_with_pos)),
                                   size_num=10,
                                   header = paste("Total number of identified & localized phosphorylated", selected_species,"across each sample",sep=" "),
                                   caption_lab = "NA values and multiple sequences are removed. \n In the case of muultiple PSMs, max. intensity was selected.",
                                   x_lab = "Sample id",
                                   fill_lab =  "Sample id",
                                   y_lab = "Number of identified peptides",
                                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) 
  #+ scale_fill_manual(values = c("#6baed6","#6baed6","#6baed6"))
 
  ########## ########## ########## ########## ########## ########## ########## ##########
  # plot13 <- gg_barplt_id_pep_count(data_set = all_seq_syn,
  #                                  x_df = all_seq_syn$Pool_for_seq_merge,
  #                                  fill_df = all_seq_syn$Pool_for_seq_merge,
  #                                  ymax = nrow(all_seq_syn),
  #                                  size_num=10,
  #                                  header = paste("Total number of identified phosphorylated", selected_species,"and", background_species,"sequences across each sample",sep=" "),
  #                                  caption_lab = "NA values are removed.",
  #                                  x_lab = "Sample id",
  #                                  fill_lab =  "Sample id",
  #                                  y_lab = "Number of identified sequences",
  #                                  subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  # 
  #################################################
  
  
# 
#   barplt_phospho_seq <- final_results_with_common_col %>%
#     select(Sequence,Experiment,Intensity,Proteins,extracted_values) %>%
#     separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F) %>%
#     mutate(sample_rep_id_seq = paste(Sequence, Sample_id,Rep_id, sep = "_")) %>%
#     group_by(sample_rep_id_seq) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
#     slice(which.max(Intensity)) %>%
#     ungroup() %>% 
#     mutate(Software_name=software_name) %>%
#     mutate(Acquisition_type=acquisiton_type)
#   
#   write.table(barplt_phospho_seq, file = paste0(new_path,"/Number_of_human_phospho_sequences_",
#                                                 software_name,"_Experiment",exp_id,".txt"),
#               sep = "\t",row.names = F)
#   
#   
  ######################################################################
  
  
  #############################################################################
  ### COMPARED TO pd AND proline, MQ DOES NOT HAVE ANY QUERY THAT CONTAINS BACKGROUND SPECIES.
  ## THUS, WE FILTERED ECOLI DATA USING NOT SELECTING HUMAN AND CONTAMINANTS
  ## EXAMPLE ID WAS CONFRIMED ON UNIPROT = P0AAX3 -> YBIJ_ECOLI
  
  barplt_df_ecoli <- quant_peptides %>% 
    filter(!grepl(selected_species, Proteins) & !grepl("CON__", Proteins)) %>%
    select(Sequence,Experiment,Intensity,Proteins) %>%
    separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F) %>%
    mutate(sample_rep_id_seq = paste(Sequence, Sample_id,Rep_id, sep = "_")) %>%
    group_by(sample_rep_id_seq) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    slice(which.max(Intensity)) %>%
    ungroup()
  
  
  plot1 <- gg_barplt_id_pep_count(data_set = barplt_df,
                               x_df = barplt_df$Sample_id,
                               fill_df = barplt_df$Rep_id,
                               ymax = 20000,
                               size_num=10,
                               header = "Total number of phospho-peptides across each sample",
                               caption_lab = "NA values are removed. \n Mapping with the theoretical list was not done yet.",
                               x_lab = "Sample id",
                               fill_lab =  "Sample id",
                               y_lab = "Number of identified peptides",
                               subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) #+
    #scale_fill_manual(values = c("#6baed6","#6baed6","#6baed6"))
  
  
  plot2 <- gg_barplt_id_pep_count(data_set = barplt_df_ecoli,
                               x_df = barplt_df_ecoli$Sample_id,
                               fill_df = barplt_df_ecoli$Rep_id,
                               ymax = 20500,
                               size_num=10,
                               header = "Total number of quantified Ecoli across each sample",
                               caption_lab = "NA values are removed.",
                               x_lab = "Sample id",
                               fill_lab =  "Sample id",
                               y_lab = "Number of identified peptides",
                               subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) 
  #+#scale_fill_manual(values =c("#969696","#969696","#969696"))
  
  barplt_prot_ecoli <- quant_peptides %>% 
    filter(!grepl(selected_species, Proteins) & !grepl("CON__", Proteins)) %>%
    select(Sequence,Experiment,Intensity,Proteins) %>%
    separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F) %>%
    #mutate(sample_rep_id_seq = paste(Sequence, Sample_id,Rep_id, sep = "_")) %>%
    group_by(Proteins,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    slice(which.max(Intensity)) %>%
    ungroup() %>%
    mutate(type=acquisiton_type)
  
  plot16 <- gg_barplt_id_pep_count(data_set = barplt_prot_ecoli,
                                   x_df = barplt_prot_ecoli$Sample_id,
                                   fill_df = barplt_prot_ecoli$Rep_id,
                                   ymax = 20500,
                                   size_num = 10,
                                   header = "Total number of identified Ecoli proteins across each sample",
                                   caption_lab = "NA values are removed.",
                                   x_lab = "Sample id",
                                   fill_lab =  "Sample id",
                                   y_lab = "Number of identified proteins",
                                   subtitle_txt = paste("Experiment - ", exp_id, 
                                                        acquisiton_type, " data processed by ", software_name)) 
    #scale_fill_manual(values =c("#969696","#969696","#969696"))
  
  
  write.table(barplt_prot_ecoli, file=paste0(new_path,"/Exp2",software_name,"_number_of_unique_",background_species,"proteins_",".txt"),sep = "\t",col.names = T,row.names = F)
  
 
  ## THIS RESHAPING IS ONLY FOR ELIMINATION OF MULTIPLE PHOSPHO-SITES and ECOLI PEPTIDES
  ## ELIMINATION STEP IS NOT NECESSARY FOR ECOLI, 1st STRATEGY can be used only (this will decrease lines of code)

  site_prob <- barplt_df_mapped_pep %>%
    select(extracted_values, pep_with_pos,Experiment) %>%
    pivot_wider(names_from = "Experiment",
                values_from = "extracted_values") %>%
    relocate(exp_design) %>%
    rowwise() %>%
    mutate(max_value = max(c_across(all_of(exp_design)), na.rm = TRUE)) %>%
    #mutate(mean_value = mean(c_across(all_of(exp_design)), na.rm = TRUE)) %>%
    select(!exp_design)
  
  
  barplt_df_wide_all_col <- barplt_df_mapped_pep %>%  ## If you select "charge" column, it will bring multiple rows for one seq
    select(Sequence,Experiment,Intensity,pep_with_pos,Proteins) %>%
    pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
    relocate(exp_design,.after = where(is.character)) %>%
    left_join(site_prob, by="pep_with_pos") %>%
    mutate(species=selected_species) %>%
    relocate(pep_with_pos, .after = Sequence) 
  
  
  # barplt_df_wide_all_col <- barplt_df_mapped_pep %>%  ## If you select "charge" column, it will bring multiple rows for one seq
  #   select(Sequence,Experiment,Intensity,pep_with_pos,extracted_values) %>% #Proteins
  #   pivot_wider(names_from = "Experiment",values_from = "Intensity") %>% #extracted_values
  #   relocate(exp_design,.after = where(is.character)) %>%
  #   #relocate(exp_design,.after = where(is.character)) %>%
  #   #left_join(site_prob, by="pep_with_pos") %>%
  #   mutate(species=selected_species) %>%
  #   relocate(pep_with_pos, .after = Sequence) 
  # 
  # 
  
  write.table(barplt_df_wide_all_col, file=paste0(new_path,"/Count_of_missing_unexpected_correct_phospho-sites_with_all_col_",
                                                  software_name,"_Experiment",exp_id,".txt"),
              sep = "\t",col.names = T,row.names = F)
  
  if(loc_filter_opt == TRUE){
    # 
    # barplt_df_wide <- barplt_df_mapped_pep %>%  ## If you select "charge" column, it will bring multiple rows for one seq
    #   select(Sequence,Experiment,Intensity,pep_with_pos,extracted_values,Proteins) %>%
    #   group_by(pep_with_pos,Experiment) %>%
    #   filter(extracted_values >= loc_filter) %>%
    #   ungroup() %>%
    #   #relocate(exp_design,.after = where(is.character)) %>%
    #   #left_join(site_prob, by="pep_with_pos") %>%
    #   pivot_wider(names_from = "Experiment",values_from =c( "Intensity")) %>% #extracted_values
    #   mutate(species=selected_species) #%>%
    #   #select(!contains("extracted_values")) %>%  rename_with(~ exp_design, matches("^Intensity_")) #%>%
 
    barplt_df_wide <- barplt_df_mapped_pep %>%  ## If you select "charge" column, it will bring multiple rows for one seq
      select(Sequence,Experiment,Intensity,pep_with_pos,Proteins) %>%
      pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
      relocate(exp_design,.after = where(is.character)) %>%
      left_join(site_prob, by="pep_with_pos") %>%
      mutate(species=selected_species) %>%
      filter(max_value >= loc_filter) %>%
      relocate(pep_with_pos, .after = Sequence) 

  }else{
    
    barplt_df_wide <- barplt_df_mapped_pep %>%  ## If you select "charge" column, it will bring multiple rows for one seq
      select(Sequence,Experiment,Intensity,pep_with_pos,Proteins) %>%
      pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
      left_join(site_prob, by="pep_with_pos") %>%
      mutate(species=selected_species) %>%
      #filter(max_value >= loc_filter) %>%
      relocate(pep_with_pos, .after = Sequence) 
    
    # barplt_df_wide <- barplt_df_mapped_pep %>%  ## If you select "charge" column, it will bring multiple rows for one seq
    #   select(Sequence,Experiment,Intensity,pep_with_pos,extracted_values,Proteins) %>%
    #   pivot_wider(names_from = "Experiment",values_from = "Intensity") %>% #extracted_values
    #   #select(!contains("extracted_values")) %>%  rename_with(~ exp_design, matches("^Intensity_"))
    #   #left_join(site_prob, by="pep_with_pos") %>%
    #   mutate(species=selected_species) #%>%
    #   #filter(max_value >= loc_filter) %>%
    #   #relocate(pep_with_pos, .after = Sequence) 
  }
  
  # # Function to merge rows based on the described conditions
  # merged_df <- barplt_df_wide %>%
  #   #select(pep_with_pos,ptm_score,contains(exp_design)) %>%
  #   group_by(pep_with_pos) %>%  # Group by identical 'pep' values
  #   mutate(
  #     # Calculate the number of non-NA values in columns containing "_r" for each row
  #     non_na_count = rowSums(!is.na(across(all_of(exp_design))))
  #   ) %>%
  #   mutate(ptm_score_num=as.numeric(extracted_values)) %>%
  #   arrange(desc(ptm_score_num), desc(non_na_count)) %>%  # Arrange to have the row with higher ptm and more non-NA values on top
  #   summarise(
  #     merge_ptm_score =  mean(ptm_score_num, na.rm = TRUE),#first(ptm_score_num),  # Keep the higher ptm value for each group
  #     # For each column with "-r", if the first row has NA, fill it from other rows
  #     across(all_of(exp_design), ~ ifelse(is.na(first(.)), mean(., na.rm = TRUE), first(.)), .names = "{.col}")
  #   ) %>%
  #   ungroup() %>%
  #   mutate(across(where(is.numeric), ~ na_if(., NaN)))
  # 
  
  barplt_df_ecoli_wide <- barplt_df_ecoli %>% 
    select(Sequence,Experiment, Proteins,Intensity) %>%
    pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
    mutate(species=background_species) %>% relocate(exp_design,.after = "Proteins")
  
  write.table(barplt_df_wide,file = paste0(new_path,"/barplt_df_wide",
                                           software_name,"_Experiment",exp_id,".txt"),
              sep = "\t",row.names = F)
  ####### ADDITIONAL PLOT TO DISPLAY MISSING and UNEXPECTED PEPTIDES ########
  # df_merge_syn <- barplt_df_wide %>%
  #   full_join(pep_list_w_theo,by="pep_with_pos") %>% 
  #   mutate_at("Pool", ~replace_na(.,"Unexpected")) %>%
  #   mutate(Pool= ifelse(is.na(Proteins),"missing",Pool)) %>%
  #   select(pep_with_pos,starts_with(exp_design),Pool) %>%
  #   mutate(soft_name=software_name, Acquisition_type=acquisiton_type)
  
  #filtered_abundances_ecoli <- barplt_df_ecoli_wide
  #filtered_abundances <- merged_df
  
  # Nothing is changed
  #filtered_abundances<-barplt_df_wide[rowSums(!is.na(select(barplt_df_wide,starts_with(exp_design))))>0,]
  #filtered_abundances_ecoli <-barplt_df_ecoli_wide[rowSums(!is.na(select(barplt_df_ecoli_wide,starts_with(exp_design))))>0,]
  ### DATA FILTERING
  
  filtered_abundances <-  filter_NA(df = barplt_df_wide,samp_names = sample_names,num_allowed_NA = 1,num_expected_nonNA = 3)
  filtered_abundances_ecoli <-  filter_NA(df = barplt_df_ecoli_wide,samp_names = sample_names,num_allowed_NA = 5,num_expected_nonNA = 1)
  
  write.table(filtered_abundances,file = paste0(new_path,"/filtered_abundances",
                                                software_name,"_Experiment",exp_id,".txt"),
              sep = "\t",row.names = F)
  ## REMOVE SEQUENCE COLUMN AFTER NA TABLE
  filtered_abundances <- filtered_abundances %>% mutate(species=selected_species) #%>% select(!Sequence)
  
  all_data<- filtered_abundances_ecoli%>% rename(pep_with_pos=Sequence) %>%
    bind_rows(filtered_abundances) 
  
  all_na_count <- NULL
  for (i in 1:length(exp_design)){
    na_val <- as.data.frame(sum(is.na(select(all_data, contains(exp_design[i])))))
    all_na_count <- bind_rows(all_na_count,na_val)
    rm(na_val)
  }
  
  ref_col_indx <- which(all_na_count == min(all_na_count$`sum(is.na(select(all_data, contains(exp_design[i]))))`))
  ref_col <- all_data %>% select(contains(exp_design[ref_col_indx]))
  reference_int <- all_data %>% select(contains(exp_design)) %>% select(ref_col_indx) 
  colnames(reference_int) <- "reference_int"
  
  divide_int_cols <- all_data %>% 
    ungroup() %>%
    select(contains(exp_design)) %>% 
    rename(reference_int=ref_col_indx) %>% 
    mutate(across(everything(), ~ . / reference_int))
  
  exp_design_mods <- exp_design
  exp_design_mods[ref_col_indx] <- "reference_int"
  
#library(DescTools) # REMOVED: unused
  norm_factors <- NULL
  for (i in 1:length(exp_design)){
    
    ratios= as.numeric(na.omit(unlist(divide_int_cols[,exp_design_mods[i]])))
    
    if(norm_type=="mod-based"){
      
      norm_factors[i] <- 2^ModEstM::ModEstM(log2(ratios))[[1]][1]
      
      all_data_exp <- all_data %>% select(contains(exp_design)) 
      after_norm <- all_data_exp %>%
        mutate(across(everything(), ~ . / norm_factors[which(names(all_data_exp) == cur_column())]))
      
      final_data_aft_norm <- all_data %>% 
        select(!contains(exp_design)) %>%
        bind_cols(after_norm) %>%
        relocate(exp_design,.after = Proteins)
      
      filtered_abundances <- final_data_aft_norm %>% filter(grepl(selected_species,species))
      filtered_abundances_ecoli <- final_data_aft_norm %>% filter(!grepl(selected_species,species))
      
    }else if(norm_type=="median-based"){
      
      norm_factors[i] <- median(x=as.numeric(na.omit(unlist(divide_int_cols[,exp_design_mods[i]]))))
      
      all_data_exp <- all_data %>% select(contains(exp_design)) 
      after_norm <- all_data_exp %>%
        mutate(across(everything(), ~ . / norm_factors[which(names(all_data_exp) == cur_column())]))
      
      final_data_aft_norm <- all_data %>% 
        select(!contains(exp_design)) %>%
        bind_cols(after_norm) %>%
        relocate(exp_design,.after = Proteins)
      
      filtered_abundances <- final_data_aft_norm %>% filter(grepl(selected_species,species))
      filtered_abundances_ecoli <- final_data_aft_norm %>% filter(!grepl(selected_species,species))
      
    }else{
      print("No normalization method was applied.")
    }
  }
  
  filtered_abundances_ecoli_bfr_impt <- filtered_abundances_ecoli
  filtered_abundances_bfr_impt <- filtered_abundances
  
  
  ######################################
  df_merge_syn <- filtered_abundances %>%
    #select(pep_with_pos,Experiment,Intensity,Proteins) %>% 
    #pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
    full_join(pep_list_w_theo,by="pep_with_pos") %>% 
    mutate_at("Pool", ~replace_na(.,"Unexpected")) %>%
    mutate(Pool= ifelse(is.na(species),"missing",Pool)) %>%
    mutate(Pool=ifelse(Pool=="Diluted",paste0(Pool,'_',isomericity),Pool)) %>%
    select(pep_with_pos,starts_with(exp_design),Pool) %>%
    mutate(soft_name=software_name,ion_mobility=acquisiton_type)
  
  plot11 <- gg_barplt_id_pep_count(data_set = df_merge_syn,
                                   x_df = df_merge_syn$Pool,
                                   fill_df = df_merge_syn$Pool,
                                   ymax = 20000,
                                   size_num=10,
                                   header = "Total number of quantified phospho-site across each sample",
                                   caption_lab = "NA values are removed. \n No filtering based on localization",
                                   x_lab = "Sample id",
                                   fill_lab =  "Sample id",
                                   y_lab = "Number of identified peptides",
                                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
  write.table(df_merge_syn,file = paste0(new_path,"/Count_of_miss_unexpctd_corr_phosphosites_",
                                         software_name,".txt"),
              sep = "\t",row.names = F)
  

  
  #extract_df <- filtered_abundances %>% full_join(pep_list_w_theo,by="pep_with_pos")
  #write.table(extract_df,file=paste0(new_path,"/MQ_transposed_data_corr_pep_list.txt"),sep = "\t",col.names = T,row.names = F)
  
  ###############################################
  
  filtered_abundances[,"na_val"] <- apply(X = !is.na(select(filtered_abundances,contains(exp_design))), MARGIN = 1, FUN = sum)
  
  phospho_completeness <- filtered_abundances %>% 
    dplyr::count(na_val) %>%
    mutate(data_complete=((n/dim(filtered_abundances)[1])*100)) %>% mutate(species=selected_species)
  
  filtered_abundances_ecoli[,"na_val"] <- apply(X = !is.na(select(filtered_abundances_ecoli,contains(exp_design))), MARGIN = 1, FUN = sum)
  
  completeness <- filtered_abundances_ecoli %>% 
    dplyr::count(na_val) %>%
    mutate(data_complete=((n/dim(filtered_abundances_ecoli)[1])*100)) %>%
    mutate(species=background_species) %>%
    bind_rows(phospho_completeness)
  
  write.table(completeness,file = paste0(new_path,"/data_completeness",acquisiton_type,software_name,".txt"))
  
  plot23 <- ggplot(completeness, aes(x=na_val,y=data_complete,color=species)) + geom_point(size=2.5) +
    geom_line(size=2)+
    scale_x_reverse(limits=c(18,0),breaks=seq(0, sample_size*num_reps, by = 2)) +
    ## If you look for is.na() in apply function, you should use the one below:
    #### scale_x_continuous(limits=c(0,18),breaks=seq(0, 18, by = 2)) +
    scale_y_continuous(limits = c(0,100), breaks = seq(from =0, to=100,by=10)) +
    theme_minimal() +
    theme(legend.text = element_text(size=30), 
          axis.title.x = element_text(size=30),
          axis.title.y = element_text(size=30),
          plot.title = element_text(size=35),
          plot.subtitle = element_text(size = 20),
          legend.title=element_text(size=30),
          axis.text=element_text(size=30),
          axis.title=element_text(size=30),
          strip.text.x = element_text(
            size = 15
          )
    ) + scale_color_manual(values = c("#636363","#3182bd"))+
    ggtitle(label = paste("Data completeness of",background_species,"and",selected_species)) +
    labs(x="n Sample", y="% of peptides",subtitle = paste("Experiment",exp_id,software_name,acquisiton_type,"\n",file_name))
  
  
  # na_phospho_selected <- apply(X = is.na(filtered_abundances %>% select(Sequence, Proteins,starts_with("E2"))), MARGIN = 2, FUN = sum)
  # na_ecoli <- apply(X = is.na(filtered_abundances_ecoli %>% select(Sequence,Proteins,starts_with("E2"))), MARGIN = 2, FUN = sum)
  # 
  # na_table <- bind_rows(na_phospho_selected,na_ecoli)
  # na_table$species <- c(selected_species,background_species)
  # na_table$total <- c(dim(filtered_abundances)[1],dim(filtered_abundances_ecoli)[1])
  # 
  # na_table %>% select(-Sequence) %>%
  #   kbl(caption = paste("Number of NA values across all samples \n Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) %>%
  #   kable_material(c("striped", "hover")) %>%
  #   kable_styling(bootstrap_options = "striped", full_width = F, position = "left", font_size = 12) %>%
  #   kable_minimal(full_width = F) %>%
  #   footnote(general =  paste("This table was created after elimination of multiple charges by selecting either phospho-sites of", selected_species, "and", background_species, "sequences \n that has the highest abundace."),
  #            # number = c("Footnote 1; ", "Footnote 2; "),
  #            # alphabet = c("Footnote A; ", "Footnote B; "),
  #            # symbol = c("Footnote Symbol 1; ", "Footnote Symbol 2")
  #            footnote_as_chunk = T, title_format = c("italic", "underline")) %>%
  #   #as_image(width = 8) %>%
  #   save_kable(paste0(file_path,"outputs_with_new_script/table1.png"))
  # 
 
  
  abundances_rowMeans <- NULL
  abundances_ecoli_rowMeans <- NULL
  for (k in 1:sample_size){
    # If separate version of row means is not needed, it can be commented later.
    # Separate row Means can be collected in temp object to merge in "log_10_filtered_abundances_rowMeans"
    assign(paste0("abundances_A",k),as.data.frame(rowMeans(filtered_abundances  %>%
                                                             select(contains(paste0("A",k)))))) #%>%
    #select(starts_with("abundance_")))))
    abundances_rowMeans<- bind_cols(abundances_rowMeans,get(paste0("abundances_A",k)))
    
    
    assign(paste0("ecoli_abundances_A",k),as.data.frame(rowMeans(filtered_abundances_ecoli  %>%
                                                                   select(contains(paste0("A",k))))))# %>%
    #select(starts_with("abundance_")))))  
    
    abundances_ecoli_rowMeans<- bind_cols(abundances_ecoli_rowMeans,get(paste0("ecoli_abundances_A",k)))
    
  }
  
  quant_peptides_ECOLI_density_plt <- filtered_abundances_ecoli %>%
    select(!starts_with("E")) %>%
    bind_cols(abundances_ecoli_rowMeans) %>%
    tibble() %>%
    rename_with(~ paste0("mean abundance",1:sample_size), matches("^row")) %>%
    pivot_longer(cols = starts_with("mean"), 
                 values_to = "Intensity",
                 names_to = "sample_ids",
                 values_drop_na = T) %>%
    mutate(sample_id_seq = paste(pep_with_pos, sample_ids, sep = "_"))
  
  ecoli_density_plt<- quant_peptides_ECOLI_density_plt %>%
    select(contains(c("sample_ids","Intensity","species"))) 
  
  
  colnames(abundances_rowMeans) <- paste0("mean abundance",1:sample_size)
  
  #### MEAN ABUNDANCE RATIO WITH  DENSITY PLOT ####
  ### BEFORE IMPUTATION ###
  quant_phospho_density_plt <- filtered_abundances %>%
    select(!starts_with("E")) %>%
    bind_cols(abundances_rowMeans) %>% 
    rename(ptm_value=max_value) %>%
    #rename_with(~ paste0("mean_abun",1:5), matches("^row")) %>%
    tibble() %>% #mutate(pep_with_pos = sequence) %>% ###  At this stage, no need for phospho-position#   
    pivot_longer(cols = starts_with("mean"),
                 names_to = "sample_ids",
                 values_to = "Intensity",
                 values_drop_na = T) %>%
    mutate(sample_id_seq = paste(pep_with_pos, sample_ids, sep = "_"))
  
  density_df <-quant_phospho_density_plt %>%
    select(c(sample_ids,Intensity,species)) %>%
    bind_rows(ecoli_density_plt) #%>%
  #separate(accession, into = c("prot_id","species","position"),sep = "_",remove = F)
  
  
  plot3 <- gg_density(data_set = density_df, 
                   x_df = density_df$Intensity,
                   fill_df = density_df$species,
                   color_df = NULL,
                   header="Distribution of mean abundance of every sample before imputation",
                   facet_df = "sample_ids",
                   x_lab = "log10(intensities)",
                   color_lab= "",
                   fill_lab = "species",
                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))+scale_fill_manual(values = c("#636363","#3182bd"))
  
  
  
  
  ## ADDITIONAL IMPUTATION METHOD with MICE()
  #library(tidyverse)
  #library(tidyr)
  #library(mice)

  is.imputed_df_syn <- filtered_abundances %>%  pivot_longer(cols = starts_with(exp_design),
                                                        names_to = "Experiment",
                                                        values_to = "Intensity",
                                                        values_drop_na = F) %>%
    mutate(is.imputed=FALSE) %>%
    mutate(is.imputed=ifelse(is.na(Intensity), TRUE,is.imputed)) 
  
  
  #filtered_abundances <- as.data.frame(filtered_abundances)
  #rownames(filtered_abundances) <- paste0(filtered_abundances$pep_with_pos,"@",filtered_abundances$species,"@",1:nrow(filtered_abundances))
  
  #intensities <- filtered_abundances %>%
    #select(starts_with(exp_design))
  
  #intensities_short <- intensities[1:10,]
  #filtered_abundances_ecoli <- as.data.frame(filtered_abundances_ecoli)
  #rownames(filtered_abundances_ecoli) <- paste0(filtered_abundances_ecoli$Sequence,"@",filtered_abundances_ecoli$species,"@",1:nrow(filtered_abundances_ecoli))
  
  is.imputed_df_ecoli <- filtered_abundances_ecoli %>%  pivot_longer(cols = starts_with(exp_design),
                                                                names_to = "Experiment",
                                                                values_to = "Intensity",
                                                                values_drop_na = F) %>%
    mutate(is.imputed=FALSE) %>%
    mutate(is.imputed=ifelse(is.na(Intensity), TRUE,is.imputed)) 
  
  intensitiesECOLI <- filtered_abundances_ecoli %>%
    select(starts_with(exp_design))
  
  # Impute missing values
   for (j in 1:length(imputed_values_vec)){
  #   # Number NA
  #   #num_NA <- length(abundances_for_impute_all[,j+2][is.na(abundances_for_impute_all[,j+2])])
  #   
     filtered_abundances[,j+2][is.na(filtered_abundances[,j+2])] <- imputed_values_vec[j]
     filtered_abundances_ecoli[,j+2][is.na(filtered_abundances_ecoli[,j+2])] <- imputed_values_vec[j]
  #   
  #   #abundances_for_impute_all[,j+2][is.na(abundances_for_impute_all)[,j+2]] <- impute_values[j]
  #   # After imputation number of imputed values
  #   #num_imp <-length(abundances_for_impute_all[,j+2][(abundances_for_impute_all[,j+2]==impute_values[j])])
  #   
  #   # This is verification of imputation is done successfully
  #   # Because we expect to see that number of imputed values should be the same amount as number of NA
  #   #print(setequal(num_NA,num_imp))
  #   #print(num_NA)
  #   #print(num_imp)
   }
  abundances_all_aft_imputation <- filtered_abundances_ecoli %>%
    bind_rows(filtered_abundances)
  # new_experiment_name <-c("E2_A1_R1",
  #                         "E2_A1_R2",
  #                         "E2_A1_R3",
  #                         "E2_A2_R1",
  #                         "E2_A2_R2",
  #                         "E2_A2_R3",
  #                         "E2_A3_R1",
  #                         "E2_A3_R2",
  #                         "E2_A3_R3",
  #                         "E2_A4_R1",
  #                         "E2_A4_R2",
  #                         "E2_A4_R3",
  #                         "E2_A5_R1",
  #                         "E2_A5_R2",
  #                         "E2_A5_R3",
  #                         "E2_A6_R1",
  #                         "E2_A6_R2",
  #                         "E2_A6_R3")
  # 
  # colnames(intensities) <- new_experiment_name
  # colnames(intensitiesECOLI) <- new_experiment_name
  # 
  # imp_intensities <-  mice(intensities,m=5,maxit=10,meth='cart',seed=500)
  # #imp_intensities_norm <-  mice(intensities,m=5,maxit=10,method ='norm.nob',seed=500)
  # 
  # imp_intensitiesECOLI <- mice(intensitiesECOLI,m=5,maxit=5,meth='cart',seed=500) #maxit=10
  # 
  # completeData <- complete(imp_intensities,2)
  # colnames(completeData) <- exp_design
  # 
  # completeDataECOLI <- complete(imp_intensitiesECOLI,2)
  # colnames(completeDataECOLI) <- exp_design
  # #pattern <- md.pattern(select(quant_phospho_peptides,starts_with(exp_design)))
  # #library(VIM)
  # #aggr_plot <- aggr(intensities, col=c('navyblue','red'),
  # #numbers=TRUE, sortVars=TRUE,
  # #labels=names(intensities), cex.axis=.7,
  # #gap=3, ylab=c("Histogram of missing data","Pattern"))
  # 
  # completeDataECOLIs <-completeDataECOLI %>%
  #   mutate(tmp=rownames(completeDataECOLI)) %>%
  #   separate(tmp,into=c("pep_with_pos","species","indx"),sep="@") %>%
  #   select(!indx)
  # 
  # 
  # # UNNECESSARY TO KEEP DATA BEFORE IMPUTATION
  # abundances_all_aft_imputation <- completeData %>% 
  #   mutate(tmp=rownames(completeData)) %>%
  #   separate(tmp, into = c("pep_with_pos","species","indx"),sep = "@") %>%
  #   select(!indx) %>% 
  #   #separate(pep_with_pos, into = c("pep","pos","species"),sep = "_") %>%
  #   #mutate(pep_with_pos=paste0(pep,"_",pos)) %>%
  #   #select(!c(pep,pos)) %>%
  #   bind_rows(completeDataECOLIs)
  
  ### DISTRIBUTION OF IMPUTED VALUES ACROSS non-NA values
  is.imputed_df <- is.imputed_df_ecoli %>% 
    #rename(pep_with_pos=Sequence) %>% 
    bind_rows(is.imputed_df_syn)
  
  imputed_dataset <- abundances_all_aft_imputation %>% 
    pivot_longer(cols = starts_with(exp_design),
                 names_to = "Experiment",
                 values_to = "Intensity",
                 values_drop_na = F) %>%
    left_join(is.imputed_df,by=c("pep_with_pos","Experiment")) %>%
    mutate(Intensity.y=ifelse(is.na(Intensity.y),0,Intensity.y))
  
  p21  <- imputed_dataset %>% filter(grepl(selected_species,species.x)) %>% 
    ggplot(aes(x=log2(Intensity.x),fill=is.imputed)) + geom_histogram(bins = 30) +
    theme_minimal()+ scale_fill_manual(values = c("#7fbf7b","#af8dc3")) +
    theme(legend.text = element_text(size = 45), #aspect.ratio=6.5/11, 
          axis.title.x = element_text(size = 45),
          axis.title.y = element_text(size = 45),
          plot.title = element_text(size = 55),
          legend.title = element_text(size = 45),
          axis.text.x = element_text(size = 45),
          axis.title = element_text(size = 45),
          axis.text.y = element_text(size = 45),
          plot.subtitle = element_text(size = 45)) +
    labs( y= "Count of Intensity", x="Intenisity",
          title = paste("Distribution of imputed values \n",selected_species), 
          subtitle = paste('Experiment 2 ', acquisiton_type, " data processed by ", software_name))
  
  p22 <- imputed_dataset %>% filter(grepl(background_species,species.x)) %>% 
    ggplot(aes(x=log2(Intensity.x),fill=is.imputed)) + geom_histogram(bins = 30) +
    theme_minimal() + scale_fill_manual(values = c("#7fbf7b","#af8dc3")) +
    theme(legend.text = element_text(size = 45), #aspect.ratio=6.5/11, 
          axis.title.x = element_text(size = 45),
          axis.title.y = element_text(size = 45),
          plot.title = element_text(size = 55),
          legend.title = element_text(size = 45),
          axis.text.x = element_text(size = 45),
          axis.title = element_text(size = 45),
          axis.text.y = element_text(size = 45),
          plot.subtitle = element_text(size = 45)) +
    labs( y= "Count of Intensity", x="Intenisity",
          title = paste(background_species)) 
  #subtitle = paste('Experiment 2 ', acquisiton_type, " data processed by ", software_name))
  
  plot22 <- p21/p22
  
  # UNNECESSARY TO KEEP DATA BEFORE IMPUTATION
  # abundances_all_aft_imputation <- filtered_abundances_ecoli %>%
  #   rename_with(~ paste0("pep_with_pos"), matches("^Seq")) %>%
  #   bind_rows(filtered_abundances) 

  filtered_abundances_rowMeans <- NULL
  filtered_abundances_log10 <- NULL
  filtered_abundances_log10_rowMeans <- NULL
  
  for (k in 1:sample_size){
    # If separate version of row means is not needed, it can be commented later.
    # Separate row Means can be collected in temp object to merge in "log_10_filtered_abundances_rowMeans"
    
    #### COMMENTED CODES ARE CORRESPOND TO IMPUTATION FOR ONLY MOUSE ####
    
    # assign(paste0("aft_imp_abundances_A",k),as.data.frame(rowMeans(abundances_for_impute_mouse  %>% select(contains(paste0("A",k))))))
    # filtered_abundances_rowMeans<- bind_cols(filtered_abundances_rowMeans,get(paste0("aft_imp_abundances_A",k)))
    # 
    # assign(paste0("log10_abundances_A",k),as.data.frame(log10(abundances_for_impute_mouse  %>% select(contains(paste0("A",k)))))) 
    # filtered_abundances_log10 <- bind_cols(filtered_abundances_log10, get(paste0("log10_abundances_A",k)))
    # 
    # assign(paste0("rowMean_log10_abundances_A",k),as.data.frame(rowMeans(get(paste0("log10_abundances_A",k)) %>% select(contains(paste0("A",k))))))
    # filtered_abundances_log10_rowMeans<- bind_cols(filtered_abundances_log10_rowMeans,get(paste0("rowMean_log10_abundances_A",k)))
    
    assign(paste0("aft_imp_abundances_A",k),as.data.frame(rowMeans(abundances_all_aft_imputation  %>% select(contains(paste0("A",k))))))
    filtered_abundances_rowMeans<- bind_cols(filtered_abundances_rowMeans,get(paste0("aft_imp_abundances_A",k)))
    
    assign(paste0("log10_abundances_A",k),as.data.frame(log10(abundances_all_aft_imputation  %>% select(contains(paste0("A",k))))))
    filtered_abundances_log10 <- bind_cols(filtered_abundances_log10, get(paste0("log10_abundances_A",k)))
    
    assign(paste0("rowMean_log10_abundances_A",k),as.data.frame(rowMeans(get(paste0("log10_abundances_A",k)) %>% select(contains(paste0("A",k))))))
    filtered_abundances_log10_rowMeans<- bind_cols(filtered_abundances_log10_rowMeans,get(paste0("rowMean_log10_abundances_A",k)))
    
  }
  
  colnames(filtered_abundances_rowMeans) <- paste0("mean_abundances_aft_imp_A",1:sample_size)
  colnames(filtered_abundances_log10_rowMeans) <- paste0("mean_log10_abundances_A",1:sample_size)
  filtered_abundances_log10 <- filtered_abundances_log10 %>%
    rename_with(~ paste0("log10_", .x), everything())
  
  # Calculate Fold Change by keeping A1 constant (mean(S1)/mean(S2), etc.)
  cols <- ncol(filtered_abundances_rowMeans)
  for(An in 1:cols){
    filtered_abundances_rowMeans[,paste0("exp_FC_A",numerator,"/A",An)] <- filtered_abundances_rowMeans[,numerator]/filtered_abundances_rowMeans[,An]
    
  }
  rmv_col <- paste0("exp_FC_A",numerator,"/A",numerator)
  filtered_abundances_rowMeans <- filtered_abundances_rowMeans %>% select(!rmv_col)
  
  # To calculate all binary combination in the data frame
  #mat <- do.call(cbind, lapply(cols, function(xj) 
  #  sapply(cols, function(xi) (filtered_abundances_rowMeans[, xj]/(filtered_abundances_rowMeans[, xj])))))
  #colnames(mat) <-  outer(names(filtered_abundances_rowMeans), names(filtered_abundances_rowMeans), paste0)
  
  final_imputed_data <- cbind(abundances_all_aft_imputation,
                              filtered_abundances_rowMeans,
                              filtered_abundances_log10,
                              filtered_abundances_log10_rowMeans) #filtered_abundances
  
  final_imputed_data_syn <- final_imputed_data %>% filter(grepl(selected_species,species))
  
  final_imputed_data_ecoli <- final_imputed_data  %>% filter(!grepl(selected_species, species))
  
  write.table(final_imputed_data, 
              file =paste0(new_path,"/final_imputed_data_MQ_",
                           exp_id, 
                           acquisiton_type,".txt"),sep = "\t",row.names = F)
  
  df_merge_syn_after_impt <- final_imputed_data_syn %>%
    left_join(pep_list_w_theo,by="pep_with_pos") %>% 
    mutate_at("Pool", ~replace_na(.,"Unexpected")) %>%
    mutate(Pool=ifelse(Pool=="Diluted",paste0(Pool,"_",isomericity),Pool))
  
  df_merge <- df_merge_syn_after_impt %>%
    bind_rows(final_imputed_data_ecoli) %>%
    mutate_at("Pool", ~replace_na(.,background_species)) 
 
   # df_merge <- final_imputed_data_syn %>%
   #  left_join(pep_list_w_theo,by="pep_with_pos") %>% 
   #  mutate_at("Pool", ~replace_na(.,"Unexpected")) %>%
   #  bind_rows(final_imputed_data_ecoli) %>%
   #  mutate_at("Pool", ~replace_na(.,background_species)) 
   # 
  #write.table(final_imputed_data, file = "final_imputed_normalized_data_PAL _T_cell_Exp3_( 5 conc 3reps)_NoFAIMS_DDA_with_cont_230206_2023-02-07_0947.txt",sep = "\t",row.names = F)
  
  df_mean_ab_after_impt <- df_merge %>% 
    select(contains("aft_imp") | contains("species"),Pool) %>%
    tibble() %>% 
    #separate(accession, into = c("uniprot_id", "species", "position"), remove = F) %>%
    pivot_longer(cols = contains("aft_imp"),
                 names_to = "Mean_abundance",
                 values_to = "values")# %>%
  #mutate_at("Pool", ~replace_na(.,background_species))
  
  # df_FC_ratio_after_impt <- df_merge %>% 
  #   select(starts_with("exp_") | contains("species"),Pool) %>%
  #   #separate(accession, into = c("uniprot_id", "species", "position"), remove = F) %>%
  #   tibble() %>% 
  #   pivot_longer(cols = starts_with("exp_"),
  #                names_to = "exp_FC",
  #                values_to = "values")# %>%
  # #mutate_at("Pool", ~replace_na(.,background_species))
  # 
  ##########################################################################
  df_FC_ratio_after_impt <- df_merge %>% 
    select(starts_with("exp_") | contains("species"),Pool) %>% #accession
    #separate(accession, into = c("uniprot_id", "species", "position"), remove = F) %>%
    tibble() %>% 
    pivot_longer(cols = starts_with("exp_"),
                 names_to = "exp_FC",
                 values_to = "values") %>%
    mutate_at("Pool", ~replace_na(.,background_species)) %>%
    filter(!grepl("A6",exp_FC)) %>%
    mutate(actual_ratio_val= case_when(grepl(comparisons[1],exp_FC) ~ actual_ratio[1],
                                       grepl(comparisons[2],exp_FC) ~actual_ratio[2],
                                       grepl(comparisons[3],exp_FC) ~actual_ratio[3],
                                       grepl(comparisons[4],exp_FC) ~actual_ratio[4])) %>%
    mutate(log2_act_val=log2(actual_ratio_val)) %>%
    mutate(log2_exp_val=log2(values)) %>%
    filter(grepl(selected_species,species))
  
  median_val <- df_FC_ratio_after_impt %>% group_by(exp_FC) %>%
    summarise(exp_median=median(log2_exp_val))
  
  df_FC_ratio_after_impt_final <- df_FC_ratio_after_impt %>% mutate(acquisition=acquisiton_type) %>% mutate(software_name=software_name)
  
  write.table(df_FC_ratio_after_impt_final,file = paste0(new_path,"/df_FC_ratio_after_impt",software_name,".txt"),sep = "\t",row.names =F )
  #test <- df_FC_ratio_after_impt %>% group_by(exp_FC) %>% summarise(min_val=min(log2_exp_val),max_val=max(log2_exp_val))
  
  plot21 <-gg_quant_ratio_acc(data_set = df_FC_ratio_after_impt,
                              x_df = df_FC_ratio_after_impt$log2_act_val,
                              y_df = df_FC_ratio_after_impt$log2_exp_val,
                              color_df = df_FC_ratio_after_impt$exp_FC,
                              median_col = "red",
                              header="Quantitative Ratio Assessment",
                              x_lab="log2(Actual Ratio)",
                              y_lab="log2(Experimental Ratio",
                              color_lab="Comparisons",
                              subtitle_txt=paste("Experiment", exp_id,"data acquired from", acquisiton_type,"processed by ",software_name)
  )+ scale_color_manual(values = c("#08519c","#3182bd","#6baed6","#a6bddb"))
  #+ geom_errorbar(data = test,aes(ymin = test$min_val, ymax=  test$max_val,color=test$exp_FC), width=0.5) 
  
  ###########################################################
  plot4 <- gg_density(data_set = df_mean_ab_after_impt, 
                   x_df = df_mean_ab_after_impt$values,
                   fill_df = df_mean_ab_after_impt$species, #Mean_abundance
                   color_df = NULL,#df_mean_ab_after_impt$Pool,
                   header="Distribution of mean abundance of every sample after imputation",
                   facet_df = "Mean_abundance",
                   x_lab = "log10(values)",
                   color_lab= "",
                   fill_lab = "Sample Names",
                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))# +
    #scale_fill_manual(values = c("#636363","#3182bd"))
  
  ############################## ############################
  pep_list_w_theo_sel <- pep_list_w_theo %>%
    select(pep_with_pos,isomericity,Pool,pool_id) 
  
  #################" ###################
  ############################## ############################
  ############      ############    CV ASSESSMENT BEFORE & AFTER IMPUTATION      ############      ############  
  ### BEFEORE  ### 
  filtered_abundances_ecoli_bef_impt <- filtered_abundances_ecoli_bfr_impt %>% 
    #separate(accession, into = c("uniprot_id", "species", "position"), remove = T) %>%
    rename(Pool=species) 
    #rename(pep_with_pos=Sequence)# %>%
    #select(!c(uniprot_id,position)) %>% 
    #drop_na()
  
  pep_list_w_theo_iso <- pep_list_w_theo %>% select(pep_with_pos,isomericity)
  
  df_merge_before_imp <- df_merge_syn %>% 
    left_join(pep_list_w_theo_iso,by="pep_with_pos") %>%
    mutate(Pool=ifelse(Pool=="Diluted",paste0(Pool,"_",isomericity),Pool)) %>% 
    #pivot_longer(cols = contains("E2"),
    #names_to = "Experiment",
    #values_to = "Abundances") %>% 
    #drop_na() %>%
    select(!isomericity) %>%
    #pivot_wider(names_from = "Experiment",values_from = "Abundances")%>% drop_na() 
    bind_rows(filtered_abundances_ecoli_bef_impt)
  
  abundances_all_before_impt <- NULL
  
  for (k in 1:(sample_size-1)){
    assign(paste0("df_merge_before_imp_A",k),as.data.frame(rowMeans(df_merge_before_imp  %>% select(contains(paste0("A",k))))))
    
    abundances_all_before_impt<- bind_cols(abundances_all_before_impt,get(paste0("df_merge_before_imp_A",k)))
  }
  colnames(abundances_all_before_impt) <- paste0("mean_abundance_A",1:(sample_size-1))
  
  df_merge_before_imp_all <- df_merge_before_imp %>% bind_cols(abundances_all_before_impt) %>% select(!contains("A6"))
  
  df_before_impt_CV <- CV_calculator(data = df_merge_before_imp_all,
                                     abundance_col = "mean_abundance_A",
                                     iter = 5,replicate = num_reps,
                                     exp_id = 2,
                                     acquisiton_type = acquisiton_type,
                                     software_name = software_name)
  
  write.table(df_before_impt_CV,file = paste0(new_path,"/CV_calculation_before_imputation",acquisiton_type,software_name,".txt"))
  
  ############      ############                     ############      ############  
  ### AFTER  ### 
  
  df_merge_aft_impt <- df_merge_syn_after_impt %>% #mutate(Pool=ifelse(Pool=="Diluted",paste0(Pool,"_",isomericity),Pool)) %>%
    select(!c(species)) %>% select(!contains("A6") & !contains("log10"))
  
  df_aft_impt_CV <- CV_calculator(data = df_merge_aft_impt,abundance_col = "mean_abundances_aft_imp_A",
                                  iter = 5,replicate = num_reps,
                                  exp_id = 2,
                                  acquisiton_type = acquisiton_type,
                                  software_name = software_name)
  
  write.table(df_aft_impt_CV,file = paste0(new_path,"/CV_calculation_after_imputation",acquisiton_type,software_name,".txt"))
  
  library(ggridges)
#library(paletteer) # REMOVED: unused 
  
  
  p18<- df_before_impt_CV %>% filter(grepl("Diluted",Pool) | grepl("Fixed",Pool)) %>%
    ggplot(aes(x=CV_values,y=Pool,
               color=CV_samples
    )) + 
    geom_density_ridges(scale = 0.9,
                        jittered_points = TRUE,
                        position = position_points_jitter(width = 0.05, height = 0),
                        point_shape = '*',
                        point_size = 3,
                        point_alpha = 1,
                        alpha = 0,linewidth=1.5)+ 
    #scale_fill_brewer(palette ="PuOr")+
    scale_color_manual(values = c("#08519c","#3182bd","#6baed6","#a6bddb","#d0d1e6"))+
    #scale_color_manual(values = c("#7fc97f","#beaed4","#fb8072","#fdc086","#386cb0"))+
    #scale_color_manual(values = c("#6A51A3","#6A51A3","#FD8D3C","#FD8D3C"))+ #"grey68"
    scale_linetype_manual(values = c("solid","dashed","solid","dashed"))+
    theme_minimal() +
    scale_y_discrete(expand = expand_scale(mult = c(0, 0)))+
    #facet_wrap(~Pool) +
    theme(legend.text = element_text(size = 45), #aspect.ratio=6.5/11, 
          axis.title.x = element_text(size = 45),
          axis.title.y = element_text(size = 45),
          plot.title = element_text(size = 55),
          legend.title = element_text(size = 45),
          axis.text.x = element_text(size = 45),
          axis.title = element_text(size = 45),
          axis.text.y = element_text(size = 45),
          plot.subtitle = element_text(size = 45)) +
    labs( y= "Intensity", x="Percentage of CVs",
          #title = paste("Experiment 2 Percentage of CVs", acquisiton_type, " data processed by ", software_name), 
          subtitle = paste("Before imputation")) +
    xlim(0,round(max(df_before_impt_CV$CV_values))) 
  # scale_y_continuous(
  # limits = c(0,100),
  # breaks = seq(0, 100,10))
  
  p19 <- df_aft_impt_CV %>% filter(grepl("Diluted",Pool) | grepl("Fixed",Pool)) %>%
    ggplot(aes(x=CV_values,y=Pool,
               color=CV_samples
    )) + 
    geom_density_ridges(scale = 0.9,
                        jittered_points = TRUE,
                        position = position_points_jitter(width = 0.05, height = 0),
                        point_shape = '*',
                        point_size = 3,
                        point_alpha = 1,
                        alpha = 0,linewidth=1.5)+ 
    #scale_fill_brewer(palette ="PuOr")+
    scale_color_manual(values = c("#08519c","#3182bd","#6baed6","#a6bddb","#d0d1e6"))+
    #scale_color_manual(values = c("#6A51A3","#6A51A3","#FD8D3C","#FD8D3C"))+ #"grey68"
    scale_linetype_manual(values = c("solid","dashed","solid","dashed"))+
    theme_minimal() +
    scale_y_discrete(expand = expand_scale(mult = c(0, 0)))+
    #facet_wrap(~Pool) +
    theme(legend.text = element_text(size = 45), #aspect.ratio=6.5/11, 
          axis.title.x = element_text(size = 45),
          axis.title.y = element_text(size = 45),
          plot.title = element_text(size = 55),
          legend.title = element_text(size = 45),
          axis.text.x = element_text(size = 45),
          axis.title = element_text(size = 45),
          axis.text.y = element_text(size = 45),
          plot.subtitle = element_text(size = 45)) +
    labs( y= "Intensity", x="Percentage of CVs",
          #title = paste("Experiment 2 Percentage of CVs", acquisiton_type, " data processed by ", software_name), 
          subtitle = paste("After imputation")) +
    xlim(0,round(max(df_before_impt_CV$CV_values)))
  #scale_y_continuous(
  #limits = c(0,max(df_aft_impt_CV$CV_values)), 
  #breaks = seq(0, max(df_aft_impt_CV$CV_values),10))
  
  plot18 <- p18/p19
  
  # plot19 <- df_aft_impt_CV %>% filter(grepl("Diluted",Pool) | grepl("Fixed",Pool)) %>%
  #   ggplot(aes(x=CV_values,y=Pool,fill=CV_samples)) + 
  #   geom_density_ridges(jittered_points = TRUE,position = position_points_jitter(width = 0.05, height = 0),point_shape = '*', point_size = 3, point_alpha = 1, alpha = 0.75)+ 
  #   scale_color_manual(values =c("#803233FF","#8A4D00FF","#A16928FF","#B0986CFF","#DCBE9BFF"))+ #"#EFDDCFFF",,"#EDEAC2FF"
  #   #scale_fill_grey()+
  #   theme_bw() +
  #   #facet_wrap(~Pool) +
  #   theme(aspect.ratio=6.5/11, 
  #         legend.text = element_text(size = 45),
  #         axis.title.x = element_text(size = 45),
  #         axis.title.y = element_text(size = 45),
  #         plot.title = element_text(size = 55),
  #         legend.title = element_text(size = 45),
  #         axis.text.x = element_text(size = 45),
  #         axis.title = element_text(size = 45),
  #         axis.text.y = element_text(size = 30),
  #         plot.subtitle = element_text(size = 45)) +
  #   labs( y= "Intensity", x="Percentage of CVs",
  #         title = paste("Experiment 2 Percentage of CVs", acquisiton_type, " data processed by ", software_name), 
  #         subtitle = paste("After imputation")) 
  # 
  library(ggpattern)
  p20 <-  df_before_impt_CV %>% filter(grepl("Diluted",Pool) | grepl("Fixed",Pool)) %>%
    separate(Pool,into = c("Pool_id","isomer"),sep = "_",remove = F) %>%
    mutate(isomer=ifelse(is.na(isomer),"Wrong Localization",isomer)) %>%
    mutate(Pool_id=ifelse(Pool_id=="Unexpected","Wrong Localization",Pool_id)) %>%
    
    ggplot(aes(x=CV_samples,y=CV_values ,fill=Pool)) + 
    #geom_boxplot(aes(x=CV_samples,y=CV_values ,fill=Pool)) + theme_bw() +
    #geom_half_violin(side = "r")+
    #facet_wrap(~Pool) +
    geom_boxplot_pattern(aes(fill=Pool_id,pattern=isomer),
                         pattern_fill = "grey90" ,#pattern = isomer,
                         pattern_color = "grey90",
                         pattern_density = 0.1,
                         pattern_spacing = 0.025,
                         pattern_key_scale_factor = 0.6 #outlier.shape = NA
    ) +
    scale_pattern_manual(values= c("isomeric" = "stripe", 
                                   "nonisomeric" = "none",
                                   "Wrong Localization"="none")) +
    theme_minimal() +
    labs( x= "Sample ID", y="Percentage of CVs",
          title = paste("Experiment 2 Percentage of CVs", acquisiton_type, " data processed by ", software_name), 
          subtitle = paste("Before imputation"),caption = "Outliers above 100% were removed from the figure.") +
    scale_fill_manual(values = c("#6A51A3","#E6550D","cyan3")) + #"grey68" cyan3,"#6A51A3","#E6550D",
    scale_y_continuous(
      limits = c(0,100), #max(df_aft_impt_CV$CV_values)
      breaks = seq(0, 100,10)) + #max(df_aft_impt_CV$CV_values
    theme(#legend.text = element_text(size=30), 
      strip.text = element_text(size=50,colour = "black"),
      strip.background = element_rect(color="black",fill="white", size=0.65),
      axis.title = element_text(size=50),
      #legend.position='bottom',
      #axis.title.y = element_text(size=30),
      plot.title = element_text(size=50),
      legend.title=element_text(size=50),
      legend.text = element_text(size=50),
      axis.text=element_text(size=50),
      plot.subtitle = element_text(size = 50),
      plot.caption = element_text(size = 20)
    ) 
  
  p21 <-  df_aft_impt_CV %>% filter(grepl("Diluted",Pool) | grepl("Fixed",Pool)) %>% 
    separate(Pool,into = c("Pool_id","isomer"),sep = "_",remove = F) %>%
    mutate(isomer=ifelse(is.na(isomer),"Wrong Localization",isomer)) %>%
    mutate(Pool_id=ifelse(Pool_id=="Unexpected","Wrong Localization",Pool_id)) %>%
    
    ggplot(aes(x=CV_samples,y=CV_values ,fill=Pool)) +
    #geom_boxplot(aes(x=CV_samples,y=CV_values ,fill=Pool)) + theme_bw() +
    
    geom_boxplot_pattern(aes(fill=Pool_id,pattern=isomer),
                         pattern_fill = "grey90" ,#pattern = isomer,
                         pattern_color = "grey90",
                         pattern_density = 0.1,
                         pattern_spacing = 0.025,
                         pattern_key_scale_factor = 0.6 #outlier.shape = NA
    ) +
    scale_pattern_manual(values= c("isomeric" = "stripe", 
                                   "nonisomeric" = "none",
                                   "Wrong Localization"="none")) +
    theme_minimal() +
    labs( x= "Sample ID", y="Percentage of CVs",
          title = paste("Experiment 2 Percentage of CVs", acquisiton_type, " data processed by ", software_name), 
          subtitle = paste("After imputation"),caption = "Outliers above 100% were removed from the figure.") +
    scale_fill_manual(values = c("#6A51A3","#E6550D","cyan3")) + #"grey68" cyan3,"#6A51A3","#E6550D",
    scale_y_continuous(
      limits = c(0,100), #max(df_aft_impt_CV$CV_values)
      breaks = seq(0, 100,10)) + #max(df_aft_impt_CV$CV_values
    theme(#legend.text = element_text(size=30), 
      strip.text = element_text(size=50,colour = "black"),
      strip.background = element_rect(color="black",fill="white", size=0.65),
      axis.title = element_text(size=50),
      #legend.position='bottom',
      #axis.title.y = element_text(size=30),
      plot.title = element_text(size=50),
      legend.title=element_text(size=50),
      legend.text = element_text(size=50),
      axis.text=element_text(size=50),
      plot.subtitle = element_text(size = 50),
      plot.caption = element_text(size = 20)
    ) 

  plot20 <- p20/p21
  
  ############################## ############################
  ratio_supp_FC <-df_merge %>% 
    select(starts_with("exp_") | contains("accession"),pep_with_pos) %>%
    left_join(pep_list_w_theo_sel,by="pep_with_pos") %>%
    pivot_longer(cols = starts_with("exp_"),
                 names_to = "exp_FC",
                 values_to = "values")
  
  ratio_supp_FC_Fixed <- ratio_supp_FC %>%
    filter(grepl("Fixed",Pool)) %>%
    drop_na(Pool)
  
  ratio_supp_FC_Spiked <- ratio_supp_FC %>%
    filter(!grepl("Fixed",Pool)) %>%
    drop_na(Pool) %>%
    mutate(Pool=paste(Pool,isomericity,sep = "_"))
  
  p15 <- gg_raincloud(data_set = ratio_supp_FC_Spiked,
                      x_df = ratio_supp_FC_Spiked$exp_FC,
                      y_df = ratio_supp_FC_Spiked$values,
                      fill_df = ratio_supp_FC_Spiked$Pool,
                      
                      header = "Distribution of mean abundance of every sample after imputation",
                      x_lab = "Sample Names",
                      y_lab = " Density of log10(Mean Abundance)",
                      fill_lab = "Sample Names",
                      caption_lab = "",
                      subtitle_txt = "Spiked Pool") + 
    scale_fill_manual(values = c("#6A51A3","#544082"))
  
  p16 <-  gg_raincloud(data_set = ratio_supp_FC_Fixed,
                       x_df = ratio_supp_FC_Fixed$exp_FC,
                       y_df = ratio_supp_FC_Fixed$values,
                       fill_df = ratio_supp_FC_Fixed$Pool,
                       
                       header = "Distribution of mean abundance of every sample after imputation",
                       x_lab = "Sample Names",
                       y_lab = " Density of log10(Mean Abundance)",
                       fill_lab = "Sample Names",
                       caption_lab = "",
                       subtitle_txt = "Fixed Pool") + 
    scale_fill_manual(values = c("#E6550D","#b8440a"))
  
  
  #library(patchwork)
  plot17 <- p15/p16
  ##########################################################
  ### BOX-PLOT: Experimental Quantity Ratio of Phospho Peptides  
  
  # df_FC_ratio_absErr_after_impt <- df_merge %>% 
  #   select(pep_with_pos,starts_with("exp_") | contains("species"),Pool,isomericity) %>%
  #   #separate(accession, into = c("uniprot_id", "species", "position"), remove = F) %>%
  #   tibble() %>% 
  #   pivot_longer(cols = starts_with("exp_"),
  #                names_to = "exp_FC",
  #                values_to = "values") %>%
  #   mutate_at("Pool", ~replace_na(.,background_species)) %>%
  #   mutate(actual_ratio_val= case_when(grepl(comparisons[1],exp_FC) ~ as.numeric(actual_ratio[1]),
  #                                      grepl(comparisons[2],exp_FC) ~as.numeric(actual_ratio[2]),
  #                                      grepl(comparisons[3],exp_FC) ~as.numeric(actual_ratio[3]),
  #                                      grepl(comparisons[4],exp_FC) ~as.numeric(actual_ratio[4]))) %>%
  #   drop_na(actual_ratio_val) %>%
  #   mutate(actual_ratio_val=ifelse(Pool=="ECOLI",1,actual_ratio_val)) %>%
  #   mutate(actual_ratio_val=ifelse((Pool=="Fixed_isomeric" | Pool=="Fixed_non-isomeric"),1,actual_ratio_val)) %>%
  #   mutate(log2_fc_values =log2(values)) %>%
  #   filter(!grepl("Unexpected",Pool)) %>%
  #   #filter(!grepl("ECOLI",species)) %>%
  #   mutate(Pool=ifelse(Pool=="Diluted",paste0(Pool,"_",isomericity),Pool)) %>%
  #   #mutate(lower_bound = quantile(log2_fc_values, 0.25) - 1.5 * IQR(log2_fc_values),
  #   # upper_bound = quantile(log2_fc_values, 0.75) + 1.5 * IQR(log2_fc_values)) %>%
  #   mutate(acq_type=acquisiton_type) %>%
  #   mutate(soft_name=software_name)
  # 
  # ranges <- df_FC_ratio_absErr_after_impt %>% 
  #   group_by(exp_FC,Pool) %>%
  #   summarise(firstQ=(quantile(log2_fc_values,probs = 0.25)),
  #             iqr_val=IQR(log2_fc_values),
  #             thirdQ=(quantile(log2_fc_values,probs = 0.75))) %>%
  #   mutate(lower_bound=firstQ - 1.5*iqr_val) %>%
  #   mutate(upper_bound=thirdQ + 1.5*iqr_val) %>%
  #   ungroup() %>%
  #   mutate(Pool_FC=paste0(Pool,"_",exp_FC)) %>%
  #   select(Pool_FC,lower_bound,upper_bound)
  # 
  # df_FC_ratio_absErr_after_impt_filt <- df_FC_ratio_absErr_after_impt %>% #filter(log2_fc_values >= lower_bound & log2_fc_values <= upper_bound) %>%
  #   mutate(Pool_FC = paste0(Pool,"_",exp_FC)) %>%
  #   left_join(ranges,by = "Pool_FC") %>%
  #   group_by(Pool_FC) %>%
  #   filter(log2_fc_values >= lower_bound & log2_fc_values <= upper_bound) %>%
  #   mutate(abs_err = abs(log2(actual_ratio_val)-log2_fc_values)) %>%
  #   mutate(rel_err=(abs_err/abs(actual_ratio_val))*100) 
  # 
  # write.table(df_FC_ratio_absErr_after_impt_filt,file = paste0(file_path,"AbsError_FC",software_name,"_",acquisiton_type,".txt"),sep = 
  #               "\t",col.names = T,row.names = F)
  # 
  # 
  # library(gghalves)
  # plot18 <- gg_half_boxplt_exp_ratio_nolog(data_set = df_FC_ratio_absErr_after_impt_filt , 
  #                                          x_df = df_FC_ratio_absErr_after_impt_filt$exp_FC,
  #                                          y_df = df_FC_ratio_absErr_after_impt_filt$rel_err,
  #                                          fill_df = df_FC_ratio_absErr_after_impt_filt$Pool,
  #                                          header="Relative Absolute Error (%) using Experimental Quantity Ratio of Phospho Peptides",
  #                                          x_lab="Sample Names",
  #                                          y_lab="Relative Absolute Error (%)",
  #                                          fill_lab = "Pool Names",
  #                                          
  #                                          subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) + 
  #   #scale_y_continuous(limits = c(0,200)) + #breaks = seq(from =0, to=100,by=10)) +
  #   
  #   scale_fill_manual(values = c("#6A51A3","#6A51A3","grey68","#E6550D","#E6550D"))+ #"grey68"
  #   scale_color_manual(values =  c("#6A51A3","#6A51A3","grey68","#E6550D","#E6550D")) #"grey68"
  # #geom_half_point(alpha = 1, show.legend = TRUE, aes(color=df_FC_ratio_absErr_after_impt_filt$Pool,shape=df_FC_ratio_absErr_after_impt_filt$Pool))#+ scale_shape_manual(values = c(15,24,8,15,24)) 
  # 
  # plot19 <- gg_violin_exp_ratio_nolog(data_set = df_FC_ratio_absErr_after_impt_filt, 
  #                                     x_df = df_FC_ratio_absErr_after_impt_filt$exp_FC,
  #                                     y_df = df_FC_ratio_absErr_after_impt_filt$log2_fc_values,
  #                                     fill_df = df_FC_ratio_absErr_after_impt_filt$Pool,
  #                                     trim=TRUE,
  #                                     header="Fold change Ratio of every sample based on pool names after imputation",
  #                                     x_lab="Sample Names",
  #                                     y_lab="Fold change values",
  #                                     fill_lab = "Pool Names",
  #                                     subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name,"\n after filtering outliers")
  # ) + scale_fill_manual(values = c("#6A51A3","#6A51A3","grey68","#E6550D","#E6550D"))+ #"grey68"
  #   scale_color_manual(values =  c("#6A51A3","#6A51A3","grey68","#E6550D","#E6550D")) 
  # #scale_y_continuous(limits = c(0,200)) + #breaks = seq(from =0, to=100,by=10)) 
  # 
 
  ############################## ###########################
  
  # px <- gg_raincloud(data_set = df_mean_ab_after_impt,
  #                    x_df = df_mean_ab_after_impt$Mean_abundance,
  #                    y_df = df_mean_ab_after_impt$values,
  #                    fill_df = df_mean_ab_after_impt$Mean_abundance,
  #                    header = "Distribution of mean abundance of every sample after imputation",
  #                    x_lab = "Sample Names",
  #                    y_lab = " Density of log10(Mean Abundance)",
  #                    fill_lab = "Sample Names",
  #                    caption_lab = "",
  #                    subtitle_txt = "")
  
  
  plot5 <- gg_density(data_set = df_FC_ratio_after_impt,
                   x_df = df_FC_ratio_after_impt$values,
                   fill_df = df_FC_ratio_after_impt$exp_FC,
                   color_df = df_FC_ratio_after_impt$Pool,
                   header="Distribution of Fold change Ratio of every sample after imputation",
                   facet_df = "exp_FC",
                   x_lab = "log10(values)",
                   color_lab= "",
                   fill_lab = "Sample Names",
                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))#+scale_fill_manual(values = c("#636363","#3182bd"))
  # 
  
  ### BOX-PLOT: Experimental Quantity Ratio of Phospho Peptides 
  
  df_FC_ratio_plt6 <- df_FC_ratio_after_impt %>%
    separate(Pool,into = c("Pool_id","isomer"),sep = "_",remove = F) %>%
    mutate(isomer=ifelse(is.na(isomer),"Wrong Localization",isomer)) %>%
    mutate(Pool_id=ifelse(Pool_id=="Unexpected","Wrong Localization",Pool_id))
  
  hline_df <- data.frame(exp_FC =unique(df_FC_ratio_plt6$exp_FC),log2_act_val =log2(actual_ratio[-length(comparisons)]))
  #hline_df <- hline_df %>% separate(exp_FC,into = c("tmp","tmp2","FC"),sep = "_",remove = T) %>%
    #mutate(exp_FC=paste("FC isomeric",FC))
  
  library(ggpattern)
  plot6 <-  ggplot(df_FC_ratio_plt6,aes(x=exp_FC,y=log2_exp_val)) +
    geom_boxplot_pattern(aes(fill=Pool,pattern=isomer),
                         pattern_fill = "black" ,#pattern = isomer,
                         pattern_density = 0.1,
                         pattern_spacing = 0.025,
                         pattern_key_scale_factor = 0.6 
                         ) + #outlier.shape = NA
    #geom_line(data = hline_df,aes(x=as.factor(hline_df$exp_FC),y=hline_df$log2_act_val,group=1)) +
    geom_hline(data = hline_df,aes(yintercept=log2_act_val),color="grey",linetype="dashed",size=1) +
    geom_text(data = hline_df,color="grey40",aes(0,log2_act_val,label = paste(hline_df$exp_FC,"=",round(hline_df$log2_act_val)),hjust = -0.05, vjust = -1),size=8)+
    scale_y_continuous(breaks = seq(0,10,by=2),limits = c(0,10))+
    scale_pattern_manual(values= c("isomeric" = "stripe", 
                                   "nonisomeric" = "none",
                                   "Wrong Localization"="none")) + # manually assign pattern
    #scale_pattern_fill_manual(values=c("isomeric" = "black", "nonisomeric" = "grey90","Wrong Localization"= "cyan3")) + # manually assign colors
    scale_fill_manual(values = c("#6A51A3","#6A51A3","#E6550D","#E6550D","cyan3")) +
    theme_minimal() +
    theme(#legend.text = element_text(size=30), 
      strip.text = element_text(size=30,colour = "black"),
      strip.background = element_rect(color="black",fill="white", size=0.65),
      axis.title = element_text(size=30),
      #axis.title.y = element_text(size=30),
      plot.title = element_text(size=35),
      legend.title=element_text(size=30),
      legend.text = element_text(size=25),
      axis.text=element_text(size=30),
      
      plot.subtitle = element_text(size = 25)
      
    ) +
    labs(title = "Experimental Quantity Ratio of Phospho Peptides",
         x="Sample Names",
         y="Abundance Ratios",subtitle = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
  
     #gg_boxplt_exp_ratio(data_set = df_FC_ratio_after_impt, 
  #                           x_df = df_FC_ratio_after_impt$exp_FC,
  #                           y_df = df_FC_ratio_after_impt$values,
  #                           fill_df = df_FC_ratio_after_impt$Pool,
  #                           header="Experimental Quantity Ratio of Phospho Peptides",
  #                           x_lab="Sample Names",
  #                           y_lab="Abundance Ratios",
  #                           fill_lab = "Sample Names",
  #                           subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  # 
  
  ### HALF-BOX-PLOT & HALF-SCATTER-PLOT: Experimental Quantity Ratio of Synthetic Peptides  
  library(gghalves)
  
  plot7 <- gg_half_boxplt_exp_ratio(data_set = df_FC_ratio_after_impt, 
                                 x_df = df_FC_ratio_after_impt$exp_FC,
                                 y_df = df_FC_ratio_after_impt$values,
                                 fill_df = df_FC_ratio_after_impt$Pool,
                                 header="Experimental Quantity Ratio of Phospho Peptides with Background",
                                 x_lab="Sample Names",
                                 y_lab="Abundance Ratios",
                                 fill_lab = "Sample Names",
                                 subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
  ### VIOLIN-PLOT: Experimental Quantity Ratio of Synthetic Peptides   
  
  ### TODO: fix y scaling without trimming 
  plot8 <- gg_violin_exp_ratio(data_set = df_FC_ratio_after_impt, 
                            x_df = df_FC_ratio_after_impt$exp_FC,
                            y_df = df_FC_ratio_after_impt$values,
                            fill_df = df_FC_ratio_after_impt$Pool,
                            header="Experimental Quantity Ratio of Phospho Peptides with Background",
                            x_lab="Sample Names",
                            y_lab="Abundance Ratios",
                            fill_lab = "Sample Names",
                            trim=TRUE,
                            subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
  
  ttest_func <- function(x, y) {
    # if (sum(!is.na(x)) < 2 | sum(!is.na(y)) < 2) {
    #   return(NA)
    # }else{
    #   
    # }
    t.test(x, y,alternative = c("two.sided"))$p.value
  }
  
  wilcox_func <- function(x, y) {
    # if (sum(!is.na(x)) < 2 | sum(!is.na(y)) < 2) {
    #   return(NA)
    # }else{
    #   
    # }
    wilcox.test(as.numeric(x), as.numeric(y),alternative = c("two.sided"))$p.value
    
    
  }
  ## DATA FOR STATISTICAL ANALYSIS
  stat_analysis <- df_merge %>%
    select(pep_with_pos,Pool,species,isomericity,starts_with("log10_") | starts_with("mean_log10_") | starts_with("exp_FC")) %>%
    filter(!grepl(background_species, Pool))
  #filter(grepl("HUMAN",accession)) #spectrum_title
  
  rownames(stat_analysis) <- paste0(stat_analysis$pep_with_pos,"@",stat_analysis$Pool,"@",(1:nrow(stat_analysis))) #stat_analysis$spectrum_title,"@"
#library(multtest) # REMOVED: replaced by base R p.adjust()
  if(test_type== "t.test" | test_type== "wilcoxon"){
    # Pre-allocate memory for results
    num_iterations <- sample_size - 1
    all_pvalues <- matrix(NA, nrow(stat_analysis), num_iterations)
    all_adjust_pval <- matrix(NA, nrow(stat_analysis), num_iterations)
    
    # Loop through iterations using lapply
    for (i in 2:sample_size) {
      p_values_tmp <- lapply(1:dim(stat_analysis)[1], function(j) {
        if (test_type == "t.test") {
          ttest_func(select(stat_analysis, contains(paste0("A",numerator,"-")) & contains("log10_"))[j,],
                     select(stat_analysis, contains(paste0("A", i, "-")) & contains("log10_"))[j,])
        } else if (test_type == "wilcoxon") {
          wilcox.test(select(stat_analysis, contains(paste0("A",numerator,"-")) & contains("log10_"))[j,],
                      select(stat_analysis, contains(paste0("A", i, "-")) & contains("log10_"))[j,])
        }
      })
      
      # Extract p-values from the list
      p_values_tmp <- sapply(p_values_tmp, function(x) x)
      
      # Store p-values in the pre-allocated matrix
      all_pvalues[, i - 1] <- p_values_tmp
      
      # Perform adjustment
      qval <- p.adjust(p_values_tmp, method = "BH")
      
      # Store adjusted p-values in the pre-allocated matrix
      all_adjust_pval[, i - 1] <- qval
    }
    
    # Create data frames from matrices
    all_pvalues <- as.data.frame(all_pvalues)
    colnames(all_pvalues) <- paste0("pvalues_A",numerator,"/", "A", 2:sample_size)
    rownames(all_pvalues) <- row.names(stat_analysis)
    
    all_adjust_pval <- as.data.frame(all_adjust_pval)
    colnames(all_adjust_pval) <- paste0("pvalues_A",numerator,"/", "A", 2:sample_size)
    rownames(all_adjust_pval) <- row.names(stat_analysis)
    
    
    all_pvalues_common_col <- stat_analysis %>% 
      select(pep_with_pos, isomericity, Pool) %>% #spectrum_title
      bind_cols(all_pvalues) %>%  #all_adjust_pval
      pivot_longer(cols = starts_with("pvalues"), values_to = "pval", names_to ="pratios") %>% #adj_pval and adj_pratios
      separate(pratios, into = c("tmp","fold_change_comp","tmp1"),sep = "_") %>%
      select(!c(tmp,tmp1)) %>%
      mutate(common_col = paste(pep_with_pos,Pool,fold_change_comp,1:((sample_size-1)*nrow(stat_analysis)),sep="@")) #spectrum_title
    
    ## THE BEST WAY TO DO is this:
    merge_stat_df <- stat_analysis %>%
      select(pep_with_pos, Pool, isomericity,starts_with("exp_FC")) %>%
      pivot_longer(cols = starts_with("exp_FC"), values_to = "fold_change_values", names_to ="fold_change_ratios") %>%
      separate(fold_change_ratios, into = c("tmp","tmp1","ratio"),sep = "_") %>%
      select(!c(tmp,tmp1)) %>%
      #mutate(common_col = paste(pep_with_pos,Pool,1:((sample_size-1)*nrow(stat_analysis)),sep="@")) %>%
      bind_cols(all_pvalues_common_col$fold_change_comp,all_pvalues_common_col$pval) %>%
      rename_with(.col =6 , ~"fold_change_comp") %>%
      rename_with(.col=7, ~ "P.Value") %>%
      #filter(!grepl("ECOLI",Pool))
      #separate(accession, into = c("prot_id","species"),sep = "_")
      mutate(isomericity = ifelse(is.na(isomericity), "False Positive", isomericity)) %>%
      unite(Pool_new, Pool, isomericity,sep = "_",remove = FALSE) %>%
      unite('new_col_coloring',Pool_new,ratio,sep = "_",remove = FALSE) %>%
      mutate(new_col_coloring = if_else(grepl("Fixed_mono", new_col_coloring), "Fixed_mono", new_col_coloring)) %>%
      mutate(new_col_coloring = if_else(grepl("Fixed_multi", new_col_coloring), "Fixed_multi", new_col_coloring)) %>%
      mutate(new_col_coloring = if_else(grepl("Unexpected", new_col_coloring), "Unexpected", new_col_coloring))
    
    
    ###############################################################################
    merge_stat_df_final <- merge_stat_df #%>% filter(!grepl("A1/A6",fold_change_ratios))
    
    
  }else if(test_type=="limma"){
    # library(limma)
    # design_matrix <- model.matrix(~factor(c(rep(2,num_reps),rep(1,num_reps))))
    # merge_stat_df <-NULL
    # for ( i in 2:sample_size){
    #   # Change only the colname iteratively makes fit to every comparison
    #   colnames(design_matrix) <- c("Intercept", paste0("A1-A",i))
    #   #print(colnames(design_matrix))
    #   # Col selection for each comparison
    #   assign(paste0("df_A1vsA",i),stat_analysis %>% select(1:2 | contains("A1-") & contains("log10_") | contains(paste0("A",i,"-")) & contains("log10_")))
    #   # First, linear model was built
    #   assign(paste0("fit",i) ,lmFit(get(paste0("df_A1vsA",i))[,3:8], design_matrix))
    #   assign(paste0("fit",i), eBayes(get(paste0("fit",i))))
    #   # Readable dataframe format was generated 
    #   assign(paste0("alllimma",i), topTable(get(paste0("fit",i)), coef=2,adjust.method="BH",p.value=1,"P"))
    #   # Colnames were labeled in each comparison to make easier data merging
    #   #assign(paste0("alllimma",i),get(paste0("alllimma",i)) %>% rename_with(~ paste0(colnames(design_matrix)[2],"_", .x), everything()))
    #   
    #   assign(paste0("alllimma",i),get(paste0("alllimma",i)) %>% bind_cols(colnames(design_matrix)[2])) 
    #   # Collect everything into one object
    #   merge_stat_df <- bind_rows(merge_stat_df,get(paste0("alllimma",i)))
    # }
    # 
    
    library(limma)
    design_matrix <- model.matrix(~factor(c(rep(2,num_reps),rep(1,num_reps))))
    merge_stat_df <-NULL
    for (i in 2:sample_size){
      # Change only the colname iteratively makes fit to every comparison
      colnames(design_matrix) <- c("Intercept", paste0("A",numerator,"-A",i))
      #print(colnames(design_matrix))
      # Col selection for each comparison
      if(numerator != i){
        assign(paste0("df_A",numerator,"vsA",i),stat_analysis %>% 
                 select(1:2 | contains(paste0("A",numerator,"-")) & contains("log10_") | contains(paste0("A",i,"-")) & contains("log10_")))
        # First, linear model was built
        assign(paste0("fit",i) ,lmFit(get(paste0("df_A",numerator,"vsA",i))[,3:8], design_matrix))
        assign(paste0("fit",i), eBayes(get(paste0("fit",i))))
        # Readable dataframe format was generated 
        assign(paste0("alllimma",i), topTable(get(paste0("fit",i)), coef=2,adjust.method="BH",p.value=1,"P"))
        # Colnames were labeled in each comparison to make easier data merging
        #assign(paste0("alllimma",i),get(paste0("alllimma",i)) %>% rename_with(~ paste0(colnames(design_matrix)[2],"_", .x), everything()))
        
        assign(paste0("alllimma",i),get(paste0("alllimma",i)) %>% bind_cols(colnames(design_matrix)[2])) 
        # Collect everything into one object
        merge_stat_df <- bind_rows(merge_stat_df,get(paste0("alllimma",i)))
      }else{}
    }
    
    merge_stat_df <- bind_cols(rownames(merge_stat_df),merge_stat_df) 
    colnames(merge_stat_df)[1] <- "common_col"
    
    merge_stat_df1 <- merge_stat_df %>% 
      separate(common_col, into = c("pep_with_pos","Pool","tmp"),sep = "@") %>%
      select(!tmp) %>%
      rename_with(.col=9, ~ "fold_change_comp") %>%
      separate(fold_change_comp, into = c("first","second"),sep = "-") %>%
      mutate(fold_change_comp = paste(first,second,sep = "/")) %>%
      mutate(common_col = paste(pep_with_pos,Pool,fold_change_comp,sep = "@")) %>%
      select(!c(first,second))
    
    merge_stat_df_final <- stat_analysis %>%
      select(pep_with_pos, Pool, isomericity,starts_with("exp_FC")) %>%
      pivot_longer(cols = starts_with("exp_FC"), values_to = "fold_change_values", names_to ="fold_change_ratios") %>%
      separate(fold_change_ratios, into = c("tmp","tmp1","ratio"),sep = "_") %>%
      select(!c(tmp,tmp1)) %>%
      #mutate(common_col = paste(common_col_for_merging,Pool,1:((sample_size-1)*nrow(stat_analysis)),sep="@")) %>%
      #bind_cols(merge_stat_df$logFC,merge_stat_df$P.Value,merge_stat_df$adj.P.Val) %>%
      mutate(common_col = paste(pep_with_pos,Pool,ratio,sep = "@")) %>%
      #rename_with(.col =8 , ~"common_col") %>%
      left_join(merge_stat_df1,by="common_col") %>%
      # rename_with(.col=7, ~ "pvalues") %>%
      #rename_with(.col=8, ~ "adj_pvalues") %>%
      
      #filter(!grepl("ECOLI",Pool))
      #separate(accession, into = c("prot_id","species"),sep = "_")
      mutate(isomericity = ifelse(is.na(isomericity), "False Positive", isomericity)) %>%
      
      #unite('Pool_new',Pool.x, isomericity,sep = "_",remove = FALSE) %>%
      #mutate(Pool_new=ifelse(Pool_new=="Fixed_isomeric_isomeric","Fixed_isomeric",Pool_new)) %>%
      #mutate(Pool_new=ifelse(Pool_new=="Fixed_non-isomeric_nonisomeric","Fixed_non-isomeric",Pool_new)) %>%
      separate(Pool.x,into = c("Pool_name","isomer_type"),sep = "_",remove = F) %>%
      mutate(coloring_comp=paste0(Pool_name,"_",ratio)) %>%
      unite('new_col_coloring',Pool.x,ratio,sep = "_",remove = FALSE) %>%
      
      #mutate(new_col_coloring = if_else(grepl("Fixed_isomeric", new_col_coloring), "Fixed_isomeric", new_col_coloring)) %>%
      #mutate(new_col_coloring = if_else(grepl("Fixed_non-isomeric_nonisomeric", new_col_coloring), "Fixed_nonisomeric", new_col_coloring)) %>%
      mutate(new_col_coloring = if_else(grepl("Unexpected", new_col_coloring), "Unexpected", new_col_coloring)) %>%
      rename(pep_with_pos=pep_with_pos.x) %>% rename(Pool=Pool.x) %>%
      #unite(Pool_new, Pool.x, isomericity,sep = "_",remove = FALSE) %>%
      #unite('new_col_coloring',Pool_new,ratio,sep = "_",remove = FALSE) %>%
      #mutate(new_col_coloring = if_else(grepl("Fixed_multi", new_col_coloring), "Fixed_multi", new_col_coloring)) %>%
      #mutate(new_col_coloring = if_else(grepl("Fixed_mono", new_col_coloring), "Fixed_mono", new_col_coloring)) %>%
      #mutate(new_col_coloring = if_else(grepl("Unexpected", new_col_coloring), "Unexpected", new_col_coloring)) %>%
      #rename(pep_with_pos=pep_with_pos.x) %>% rename(Pool=Pool.x) %>%
      filter(!grepl("A1/A6",fold_change_comp))
    
    # merge_stat_df_final <- stat_analysis %>%
    #   select(pep_with_pos, Pool, isomericity,starts_with("exp_FC")) %>%
    #   pivot_longer(cols = starts_with("exp_FC"), values_to = "fold_change_values", names_to ="fold_change_ratios") %>%
    #   separate(fold_change_ratios, into = c("tmp","tmp1","ratio"),sep = "_") %>%
    #   select(!c(tmp,tmp1)) %>%
    #   #mutate(common_col = paste(common_col_for_merging,Pool,1:((sample_size-1)*nrow(stat_analysis)),sep="@")) %>%
    #   #bind_cols(merge_stat_df$logFC,merge_stat_df$P.Value,merge_stat_df$adj.P.Val) %>%
    #   mutate(common_col = paste(pep_with_pos,Pool,ratio,sep = "@")) %>%
    #   #rename_with(.col =8 , ~"common_col") %>%
    #   left_join(merge_stat_df1,by="common_col") %>%
    #   # rename_with(.col=7, ~ "pvalues") %>%
    #   #rename_with(.col=8, ~ "adj_pvalues") %>%
    #   
    #   #filter(!grepl("ECOLI",Pool))
    #   #separate(accession, into = c("prot_id","species"),sep = "_")
    #   mutate(isomericity = ifelse(is.na(isomericity), "False Positive", isomericity)) %>%
    #   
    #   unite('Pool_new',Pool.x, isomericity,sep = "_",remove = FALSE) %>%
    #   mutate(Pool_new=ifelse(Pool_new=="Fixed_isomeric_isomeric","Fixed_isomeric",Pool_new)) %>%
    #   mutate(Pool_new=ifelse(Pool_new=="Fixed_non-isomeric_nonisomeric","Fixed_non-isomeric",Pool_new)) %>%
    #   unite('new_col_coloring',Pool_new,ratio,sep = "_",remove = FALSE) %>%
    #   #mutate(new_col_coloring = if_else(grepl("Fixed_isomeric", new_col_coloring), "Fixed_isomeric", new_col_coloring)) %>%
    #   #mutate(new_col_coloring = if_else(grepl("Fixed_non-isomeric_nonisomeric", new_col_coloring), "Fixed_nonisomeric", new_col_coloring)) %>%
    #   mutate(new_col_coloring = if_else(grepl("Unexpected", new_col_coloring), "Unexpected", new_col_coloring)) %>%
    #   rename(pep_with_pos=pep_with_pos.x) %>% rename(Pool=Pool.x) %>%
    #   #unite(Pool_new, Pool.x, isomericity,sep = "_",remove = FALSE) %>%
    #   #unite('new_col_coloring',Pool_new,ratio,sep = "_",remove = FALSE) %>%
    #   #mutate(new_col_coloring = if_else(grepl("Fixed_multi", new_col_coloring), "Fixed_multi", new_col_coloring)) %>%
    #   #mutate(new_col_coloring = if_else(grepl("Fixed_mono", new_col_coloring), "Fixed_mono", new_col_coloring)) %>%
    #   #mutate(new_col_coloring = if_else(grepl("Unexpected", new_col_coloring), "Unexpected", new_col_coloring)) %>%
    #   #rename(pep_with_pos=pep_with_pos.x) %>% rename(Pool=Pool.x) #%>%
    #   filter(!grepl("A1/A6",fold_change_comp))
    
  }else{
    print("Statistical test could not be assessed. Check the input files!")
  }

  
  ##TODO: Find more logical way to map these values!!!
  comparisons <- comparisons[-5]
  
  if(length(comparisons) == 4){
    actual_ratio_col <- merge_stat_df_final %>%
      select(fold_change_comp) %>% distinct() %>%
      mutate(actual_ratio_val= case_when(grepl(comparisons[1],fold_change_comp) ~ actual_ratio[1],
                                         grepl(comparisons[2],fold_change_comp) ~actual_ratio[2],
                                         grepl(comparisons[3],fold_change_comp) ~actual_ratio[3],
                                         grepl(comparisons[4],fold_change_comp) ~actual_ratio[4]))
    
  }else if(length(comparisons) == 5){
    actual_ratio_col <- merge_stat_df_final %>%
      select(fold_change_comp) %>% distinct() %>%
      mutate(actual_ratio_val= case_when(grepl(comparisons[1],fold_change_comp) ~ actual_ratio[1],
                                         grepl(comparisons[2],fold_change_comp) ~actual_ratio[2],
                                         grepl(comparisons[3],fold_change_comp) ~actual_ratio[3],
                                         grepl(comparisons[4],fold_change_comp) ~actual_ratio[4],
                                         grepl(comparisons[5],fold_change_comp) ~actual_ratio[5]))
  }else{
    print("Mapping between theoretical ratio and comparison cannot be done. Please make sure that you have either 4 or 5 comparisons overall.")
  }
  
  ## Generation of df -> expected abundance ratio for volcano plot
  # actual_ratio_col <- merge_stat_df_final %>%
  #   select(A1vs_Ai) %>% distinct() %>%
  #   mutate(actual_ratio_val = case_when(grepl(comparisons[1],A1vs_Ai) ~actual_ratio[1],
  #                                       grepl(comparisons[2],A1vs_Ai) ~actual_ratio[2],
  #                                       grepl(comparisons[3],A1vs_Ai) ~actual_ratio[3],
  #                                       grepl(comparisons[4],A1vs_Ai) ~actual_ratio[4]))
  # 
  point_count_y_axis <- merge_stat_df_final %>%
    group_by(fold_change_comp, new_col_coloring,coloring_comp) %>%
    filter(P.Value < 0.05) %>% 
    dplyr::count(new_col_coloring) %>% left_join(actual_ratio_col)
  
  ymax <-11 + 0.5 # max(-log10(merge_stat_df_final$P.Value)) 
  y_decrement <- 0.5
  
  
  calculate_y_pos <- function(group) {
    group_length <- length(group)
    y_pos <- ymax - seq(0, by = y_decrement, length.out = group_length)
    return(y_pos)
  }
  
  #comparisons <- comparisons[-5]
  
  # Apply the function to calculate y_pos within each group
  point_count_y_axis$y_pos <- unlist(by(point_count_y_axis$fold_change_comp, 
                                        point_count_y_axis$fold_change_comp, calculate_y_pos))
  
  ######### COLORING STRATEGY BASED ON ONLY RATIO and SAMPLE NAME (coloring_comp) #########
  col_vline <- rep("#000000",4)
  #"#FD8D3C" ligther orange ,"#9E9AC8"=lighter purple instead of #807DBA
  
  
  if(length(unique(merge_stat_df_final$coloring_comp)) == 8){
    cols <- c(rep("#E6550D" ,4),"#3F007D","#6A51A3","#807DBA","#BCBDDC") #"#FD8D3C" ligther orange ,"#9E9AC8"=lighter purple instead of #807DBA
    
    corr_level <-  c("Fixed_A1/A2","Fixed_A1/A3","Fixed_A1/A4","Fixed_A1/A5",
                     "Diluted_A1/A2","Diluted_A1/A3","Diluted_A1/A4","Diluted_A1/A5"
    )
  }else if(length(unique(merge_stat_df_final$coloring_comp)) == 12){
    cols <- c(rep("#E6550D" ,4),"#3F007D","#6A51A3","#807DBA","#BCBDDC",rep("cyan3",4))
    corr_level <-  c("Fixed_A1/A2","Fixed_A1/A3","Fixed_A1/A4","Fixed_A1/A5",
                     "Diluted_A1/A2","Diluted_A1/A3","Diluted_A1/A4","Diluted_A1/A5",
                     "Unexpected_A1/A2", "Unexpected_A1/A3", "Unexpected_A1/A4", "Unexpected_A1/A5"
    )
  }else{
    print("There are some Unexpected peptides found! Then, they don't appear in all samples!")
  }
  
  labels <- factor(unique(merge_stat_df_final$coloring_comp),levels=corr_level)
  
  sorted_labels <- sort(labels)
  
  actual_ratio_col <- as.data.frame(actual_ratio_col %>% bind_cols(col_vline)) %>% rename(col=3)
  #actual_ratio_col_filt <- as.data.frame(actual_ratio_col_filt %>% bind_cols(col_vline)) %>% rename(col=3)
  
  mapped_coloring_dat <- as.data.frame(cols)
  
  point_count_y_axis <- mapped_coloring_dat %>% 
    bind_cols(sorted_labels) %>%
    rename(colors=1,coloring_comp=2) %>%
    right_join(point_count_y_axis,by="coloring_comp")
  
  plot9 <- ggplot(merge_stat_df_final, aes(x =log2(fold_change_values), y = -log10(P.Value))) +
    geom_point(size = 9, aes(color = coloring_comp,shape=isomericity)) + # Pool_new might be use to 
    #geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 5)
    scale_color_manual(values =setNames(cols, sorted_labels)) +
    #scale_shape_identity() +                                        # differentiate peptides are found as Unexpected and reference 
    # in their associated concentration 
    # Pool can be used to show only difference btw pools same as coloring 'less complex visualization)
    #geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "red") +
    #scale_fill_manual(values=setNames(mapped_coloring, labels)) + 
    #geom_line(aes(color =setNames(mapped_coloring, labels)), size = 1) +  # Add color aesthetic to geom_line()
    ##scale_color_manual(values =setNames(mapped_coloring, labels)) +
    #scale_y_continuous(limits = c(0, max(-log10(merge_stat_df_final$adj.P.Val))), breaks = seq(0, max(-log10(merge_stat_df_final$adj.P.Val)), by = 0.8)) +
    #scale_x_continuous(limits = c(min(log2(merge_stat_df_final$fold_change_values)),max(log2(merge_stat_df_final$fold_change_values)))) +#facet_wrap(~ratio) +
  #scale_x_continuous(breaks = seq(from =round(min(log2(merge_stat_df_final$fold_change_values))), to=(round(max(log2(merge_stat_df_final$fold_change_values)))+2),by=1)) +
  ##scale_y_continuous(breaks = seq(from =round(min(-log10(merge_stat_df_final$P.Value))), to=(round(max(-log10(merge_stat_df_final$P.Value)))+2),by=1)) +
  #scale_y_continuous(breaks = seq(0, max(-log10(volcano_final1$pvalues_value)), length.out = 21)) +
  theme_bw() +
    scale_y_continuous(
      limits = c(0,12), 
      breaks = seq(0, 12,2)
    ) +
    scale_x_continuous(
      limits = c(-9,11), 
      breaks = seq(-9, 11,2)
    ) +
    theme(legend.text = element_text(size = 45),
          axis.title.x = element_text(size = 45),
          axis.title.y = element_text(size = 45),
          plot.title = element_text(size = 55),
          legend.title = element_text(size = 45),
          axis.text.x = element_text(size = 45),
          axis.title = element_text(size = 45),
          axis.text.y = element_text(size = 45),
          plot.subtitle = element_text(size = 45)) +
    labs( y= "-log10(p values)", x="log2(fold change)",title = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name), subtitle = paste("Limma was used \n")) + #plt_data_type[j], , subtitle,"\n"
    geom_vline(data = actual_ratio_col, aes(xintercept = log2(actual_ratio_val), show.legend = FALSE,color=col),size=2) +
    geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "#006D2C",size=2) + 
    geom_label(data =point_count_y_axis, aes(x = log2(actual_ratio_val), y = y_pos, fill=coloring_comp ,label = n),color="white",size=18,show.legend = FALSE) + #coloring_comp
    scale_fill_manual(values =setNames(cols, sorted_labels))

  merge_stat_df_final_text <- merge_stat_df_final %>% mutate(soft_name=paste0(software_name)) %>% mutate(acq_type=paste0(acquisiton_type))
  
  write.table(merge_stat_df_final_text,file = paste0(new_path,"/volcano_plot_",software_name,"_",acquisiton_type,".txt"),sep = 
                "\t",col.names = T,row.names = F)
  
  df_roc <- merge_stat_df_final %>%
    select(pep_with_pos, Pool,P.Value) %>%
    #select(pep_with_pos, Pool,P.Value) %>%
    mutate(Pool = if_else(grepl("Diluted_isomeric", Pool), "Diluted", Pool)) %>%
    mutate(Pool = if_else(grepl("Diluted_nonisomeric", Pool), "Diluted", Pool))
  #filter(!grepl("Unexpected",Pool))
  
  ### ROC analysis custom func
  df_roc_order <- df_roc[order(df_roc$P.Value),]
  
  df_roc_func <- compute_roc_curve(df=df_roc_order, flag = "Diluted",expected = (length(comparisons)*size_variying_pep_size))
  
  write.table(df_roc_func, file = paste0(new_path,"/new_custom_Roc_analysis_",exp_id,"_",software_name,"_",".txt"),sep = "\t",row.names = F)
  
  plot14 <- ggplot(df_roc_func, aes(y=tpr, x = fdr)) +
    geom_path(size=1.5) +
    #geom_vline(aes(xintercept=fdr)) +
    #geom_text(data=as.data.frame(result),aes(label=fdr)) +
    #scale_x_reverse() + 
    theme_bw() +
    theme(legend.text = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          plot.title = element_text(size = 25),
          legend.title = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title = element_text(size = 20),
          axis.text.y = element_text(size = 20)) +
    scale_color_brewer(palette = "Dark2") +
    labs(y="True Positive Rate \n (Sensitivity)", x="False Positive Rate \n (Specificity)",
         title =  paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name), 
         subtitle = paste(subtitle,"including unexpected"), color="Pool Type")
  
  write.table(df_roc_func, file = paste0(new_path,"/new_custom_Roc_analysis_",exp_id,"_",software_name,"_",".txt"),sep = "\t",row.names = F)
  
  #### ROC Analysis using pROC 
  
  df_roc$variant <- ifelse(df_roc$Pool == "Diluted", TRUE, FALSE)
  #df_roc$non_var <- ifelse(df_roc$Pool == "ISO-REF", TRUE, FALSE)
  
  library(pROC)
  # Calculate ROC curve for raw p-values
  roc_raw_variant <- roc(df_roc$variant, df_roc$P.Value)
  #tpr_and_fpr_variant  <- cbind(roc_raw_variant$sensitivities,
  #                             roc_raw_variant$specificities,
  #                              "Variant Pool")
  
  fpr <- as.data.frame(1 - roc_raw_variant$specificities)
  tpr_and_fpr_variant  <- cbind(roc_raw_variant$sensitivities,
                                fpr,#roc_raw_variant$specificities,
                                "Diluted Pool")
  
  
  
  #roc_raw_non_var <- roc(df_roc$non_var, df_roc$P.Value)
  #tpr_and_fpr_non_var  <- cbind(roc_raw_non_var$sensitivities,
  #roc_raw_non_var$specificities,
  #"Non-variant Pool")
  
  roc_plt_df <- as.data.frame(tpr_and_fpr_variant) %>% 
    #bind_rows(as.data.frame(tpr_and_fpr_non_var)) 
    bind_cols(software_name)
  
  colnames(roc_plt_df) <- c("sensitivity", "fpr","Pool_type","Software_name")
  
  plot10 <- roc_plt_df %>% #group_by(Pool_type) %>% 
    ggplot( aes(y=as.numeric(sensitivity), x = as.numeric(fpr), color=Pool_type)) +
    geom_path(size=1.5)  + theme_bw() + #+  scale_x_reverse()
    theme(legend.text = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          plot.title = element_text(size = 25),
          legend.title = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title = element_text(size = 20),
          axis.text.y = element_text(size = 20)) +
    scale_color_brewer(palette = "Dark2") +
    labs(y="True Positive Rate \n (Sensitivity)", x="False Positive Rate \n (Specificity)",
         title =  paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name), 
         subtitle = paste(subtitle), color="Pool Type")
  
  write.table(roc_plt_df, file = paste0(new_path,"/pRoc_analysis_",exp_id,"_",software_name,"_",".txt"),sep = "\t",row.names = F)
  
  # sapply(1:14,function(x) ggsave(filename = paste0("p",x,".tiff"),
  #                                width = 60, height = 45, 
  #                                path = paste0(file_path,"/outputs_with_new_script/"),
  #                                units = "cm",
  #                                get(paste0("p",x)),
  #                                device = "tiff", #".svg"
  # ))
  # 
  plt_obj <- ls(pattern="plot")
  plt_obj <- plt_obj[!is.na(plt_obj)]
  sapply(1:length(plt_obj),function(x) ggsave(filename = paste0("p",x,".png"),
                                               width = 90, height = 60, 
                                               path = paste0(file_path,"/",dir_name,"/"), #outputs_imputed_5thpercent_one_sample_nonNA_fin/
                                               units = "cm",
                                               get(plt_obj[x]),
                                               device = "png", #".svg"
  ))
  
  
    
   
    
    
}


  
  
  
  
  
  