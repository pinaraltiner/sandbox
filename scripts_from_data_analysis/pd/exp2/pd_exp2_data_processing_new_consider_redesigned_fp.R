library(stringr)
library(dplyr)
library(data.table)
#library(openxlsx)
library(ggplot2)
library(reshape2)
library(readr)
library(tibble)
library(purrr)
library(readr)
library(tidyr)
library(openxlsx)
###############################################
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/ggplot/ggplot_functions.R")
#source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/roc_curve/roc_curve_generation_proline_edit.R")
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/roc_curve/new_roc_curve_generation_with_custom_threshold.R")

# Experiment 2

# 
# file_path <- "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/PD_data_analysis/target_decoy_no_FAIMS/new/"
# file_name <- "Multiconsensus_Exp2_TargetDecoy_woFAIMS_230612_FerriesPhosphoMarkers_PeptideIsoforms.txt"
# 
# # Correct syn. peptide list for Experiment 2
# theo_file_path <- "D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/"
# theo_file_name <- "Synthetic peptides list_theo_conc_corrected_pool_id_iso_count_final.xlsx" ### DONE!
# #theo_file_name <- "Synthetic peptides list_theo_conc_added_pool_id_iso_count.xlsx"  ### TODO CHANGE me
# sheet_theo_name <- "ISO-ref and OTHER with FC"
# 
# 

final_pd_pep_quant_analysis_redesigned <- function(file_path,
                                             file_name,
                                             #sheet_name,
                                             theo_file_path,
                                             theo_file_name,
                                             sheet_theo_name,
                                             selected_spcies,
                                             background_species,
                                             acquisiton_type,
                                             exp_id,
                                             loc_filter_opt,
                                             size_variying_pep_size,
                                             numerator,
                                             loc_filter,
                                             software_name,
                                             test_type,
                                             exp_design,
                                             num_reps,
                                        subtitle,
                                        fdr_threshold,
                                        actual_ratio
                                        ){
 

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
  
  #### PHOSPHO-POSITION EXTRACTION FUNCTION ####
  ################
  extract_phospho_numbers <- function(input_string) {
      phospho_part <- regmatches(input_string, gregexpr("Phospho \\[[^]]+\\]", input_string))
      
      if (length(phospho_part) == 0) {
          return(list(NULL, NULL))
      }
      
      phospho_text <- phospho_part[[1]]
      letter_numbers <- gregexpr("[A-Z](\\d+)", phospho_text)
      extracted_letters <- regmatches(phospho_text, letter_numbers)[[1]]
      extracted_values <- gregexpr("\\((\\d+(?:\\.\\d+)?)\\)", phospho_text)
      extracted_values <- regmatches(phospho_text, extracted_values)[[1]]
      
      extracted_values <- gsub("\\(|\\)", "", extracted_values)  # Remove parentheses
      extracted_letters <- gsub("[A-Z]", "", extracted_letters)  # Remove letters
                                                                                ## Selecting the max is necessary for applying 
                                                                                  ##           filtering at localization score
      combined_results <- list(paste(extracted_letters, collapse = "&"), paste(max(extracted_values), collapse = "&"))
      return(combined_results)
  }
  ######################
  
  setwd(file_path)
  quant_peptides <- read.table(file_name, sep = "\t", header = T)
  
  
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
  #################################################
  pep_list_w_theo_unique <- pep_list_w_theo %>% 
    distinct(Phosphopeptide.sequence,.keep_all = TRUE) %>%
    #rename(Sequence = Phosphopeptide.sequence) %>% 
    select(Sequence,Pool) %>%
    rename(Pool_for_seq_merge=Pool) %>%
    mutate(situation="Correct")
  
  ecoli_seq <- quant_peptides %>% 
    select(Sequence, Modifications, Master.Protein.Descriptions) %>%
    filter(grepl(background_species,Master.Protein.Descriptions) & !grepl("CON__",Master.Protein.Descriptions)) %>%
    #distinct(Sequence,.keep_all = T) %>%
    mutate(species=background_species) 
  
  ecoli_seq_dist <- quant_peptides %>% 
    select(Sequence, Modifications, Master.Protein.Descriptions) %>%
    filter(grepl(background_species,Master.Protein.Descriptions) & !grepl("CON__",Master.Protein.Descriptions)) %>%
    distinct(Sequence,.keep_all = T) %>%
    mutate(species=background_species) 
  
  all_seq <- quant_peptides %>% 
    filter(grepl("Homo sapiens",Master.Protein.Descriptions) & !grepl("CON__",Master.Protein.Descriptions)) %>%
    filter(grepl("Phospho",Modifications)) %>%
    mutate(species="Homo sapiens") %>%
    full_join(pep_list_w_theo_unique,by="Sequence") %>%
    mutate_at("Pool_for_seq_merge", ~replace_na(.,"Unexpected")) %>%
    mutate(Pool_for_seq_merge= ifelse(is.na(species),"missing",Pool_for_seq_merge)) %>%
    filter(!grepl("Unexpected",Pool_for_seq_merge)) %>%
    bind_rows(ecoli_seq) %>% 
    mutate(Pool_for_seq_merge= ifelse(is.na(Pool_for_seq_merge),background_species,Pool_for_seq_merge)) 
  
  ########## ########## ########## ########## ########## ########## ########## ##########
  
  all_corr_seq <-  quant_peptides %>%
    filter(grepl("Phospho",Modifications)) %>%
    select(contains("Abundances.Normalized."),Sequence, Modifications,Master.Protein.Descriptions) %>%
    filter(grepl("Homo sapiens",Master.Protein.Descriptions) & !grepl("CON__",Master.Protein.Descriptions)) %>%
    rename_with(~ exp_design, starts_with("Abundances.Normalized")) %>%
    pivot_longer(cols = starts_with("E2-"), 
                 values_to = "Intensity",
                 names_to = "Experiment",
                 values_drop_na = T) %>%
    group_by(Experiment) %>%
    distinct(Sequence,.keep_all = T) %>%
    mutate(species=selected_spcies) %>% 
    full_join(pep_list_w_theo_unique,by="Sequence") %>%
    mutate_at("situation", ~replace_na(.,"Unexpected")) %>%
    mutate(situation= ifelse(is.na(species),"missing",situation)) %>%
    filter(!grepl("Unexpected",situation) & !grepl("missing",situation) ) %>%
    select(Sequence, Modifications,Experiment,situation,Pool_for_seq_merge) %>%
    separate(Experiment, into = c("exp_id","samp_id","rep_id"),sep = "-")
  
  
  plot15 <- gg_barplt_id_pep_count(data_set = all_corr_seq,
                                   x_df = all_corr_seq$samp_id,
                                   fill_df = all_corr_seq$rep_id,
                                   ymax = nrow(all_corr_seq),
                                   size_num=10,
                                   header = paste("Total number of correctly identified ",selected_spcies,"phospho-sequence","across each sample",sep=" "),
                                   caption_lab = "Mapping was done without considering phospho-positions. \n In the case of muultiple PSMs, max. intensity was selected.",
                                   x_lab = "Sample id",
                                   fill_lab =  "Sample id",
                                   y_lab = "Number of identified peptides",
                                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) +
    theme(axis.text.x = element_text(angle = 90))
  
  
  ########## ########## ########## ########## ########## ########## ########## ##########  
  all_seq_syn <- all_seq %>% 
    select(Sequence, Modifications, Pool_for_seq_merge) %>% #Master.Protein.Descriptions
    distinct(Sequence,.keep_all = T) %>% 
    bind_rows(ecoli_seq_dist) %>%
    mutate(Pool_for_seq_merge= ifelse(is.na(Pool_for_seq_merge),background_species,Pool_for_seq_merge)) %>%
    mutate(acq_type=acquisiton_type) %>%
    mutate(soft_name=software_name)
    
  
  plot13 <- gg_barplt_id_pep_count(data_set = all_seq_syn,
                                x_df = all_seq_syn$Pool_for_seq_merge,
                                fill_df = all_seq_syn$Pool_for_seq_merge,
                                ymax = 20000,
                                size_num = 10,
                                header = paste("Total number of identified phosphorylated", selected_spcies,"and", background_species,"sequences across each sample",sep=" "),
                                caption_lab = "NA values are removed.",
                                x_lab = "Sample id",
                                fill_lab =  "Sample id",
                                y_lab = "Number of identified sequences",
                                subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
  write.table(all_seq_syn, file=paste0(file_path,"Experiment2",software_name,"_number_of_unique_sequence_for_each_species.txt"),sep = "\t",col.names = T,row.names = F)
  
  #################################################
  #mapping_df <- read.delim(paste0(file_path,mapping_file))
  
  ### IMPUTATION
  
  ### HERE MAPPING FILE WAS CREATED TO OBTAIN AND SORT EXPERIMENT NAME
  ## THIS IS THE SAFEST WAY TO GET EXPERIMENT NAME; 
  ## HOWEVER, IT IS NOT NECESSARY TO USE FOR THIS DATASET, SINCE THIS FUNCTION IS SPECIFIC TO THIS DATASET.
  ## ORDER OF ABUNDANCES IS CORRECT
  ## 
  
  # abundances_for_impute <- quant_peptides %>% 
  #     select(starts_with("Abundances.Normalized")) %>% #### BEFORE RENAME IT BE SURED THAT COLUMNS ARE THE SAME ORDER AS EXP_DESIGN
  #     #rename_with(~ exp_design, starts_with("Abundances.Normalized"))
  #     pivot_longer(cols = starts_with("Abundances.Normalized"), 
  #                values_to = "Intensity",
  #                names_to = "Files",
  #                values_drop_na = T) %>%
  #   mutate(Study.File.ID = str_extract(pattern = "F\\d+",Files)) %>%
  #   right_join(mapping_df,by="Study.File.ID") 
  # 
  # ordered_abun_cols <- abundances_for_impute_names[order(abundances_for_impute_names$Experiment),]
  # 
  
  abundances_for_impute <- quant_peptides %>% 
         select(starts_with("Abundances.Normalized")) %>% #### BEFORE RENAME IT BE SURED THAT COLUMNS ARE THE SAME ORDER AS EXP_DESIGN
         rename_with(~ exp_design, starts_with("Abundances.Normalized"))
  
  #################################################  
  
  
  
  if(loc_filter_opt == TRUE){
    quant_phospho_peptides <-all_seq %>%
      filter(grepl("Homo sapiens",Master.Protein.Descriptions) & 
               grepl("Phospho",Modifications) & 
               !grepl("positions not distinguishable", Modification.Pattern)) %>%
      select(Sequence,
             Modifications,
             Modification.Pattern,
             Number.of.PSMs,
             Master.Protein.Descriptions,
             Protein.Accessions, #Marked.as
             starts_with("Abundances.Normalized"),
             starts_with("Abundance.Ratio.P.Value"),
             starts_with("Abundance.Ratio.log2")) %>%
      rowwise() %>%
      mutate(species=selected_spcies) %>%
      mutate(results = list(extract_phospho_numbers(Modifications)),
             phospho_pos = results[[1]],
             phospho_score = results[[2]]) %>%
      select(!results) %>% 
      filter(as.numeric(phospho_score) >= loc_filter)  %>%  
      #### BEFORE RENAME IT BE SURED THAT COLUMNS ARE THE SAME ORDER AS EXP_DESIGN
      rename_with(~ exp_design, starts_with("Abundances.Normalized"))
  }else{
    quant_phospho_peptides <-all_seq %>%
      filter(grepl("Homo sapiens",Master.Protein.Descriptions) & 
               grepl("Phospho",Modifications) & 
               !grepl("positions not distinguishable", Modification.Pattern)) %>%
      select(Sequence,
             Modifications,
             Modification.Pattern,
             Number.of.PSMs,
             Master.Protein.Descriptions,
             Protein.Accessions, #Marked.as
             starts_with("Abundances.Normalized"),
             starts_with("Abundance.Ratio.P.Value"),
             starts_with("Abundance.Ratio.log2")) %>%
      mutate(species=selected_spcies) %>%
      rowwise() %>%
      mutate(results = list(extract_phospho_numbers(Modifications)),
             phospho_pos = results[[1]],
             phospho_score = results[[2]]) %>%
      select(!results) %>% 
      #filter(as.numeric(phospho_score) >= 75)  %>%  
      #### BEFORE RENAME IT BE SURED THAT COLUMNS ARE THE SAME ORDER AS EXP_DESIGN
      rename_with(~ exp_design, starts_with("Abundances.Normalized"))
  }
 
  
  #quant_phospho_peptides$Marked.as <-"HUMAN"
  
  quant_peptides_ECOLI <- quant_peptides %>%
      filter(grepl(background_species,Master.Protein.Descriptions) & !grepl("CON__",Master.Protein.Descriptions)) %>%
      #!grepl("positions not distinguishable", Modification.Pattern)) 
      select(Sequence,
             Modifications,Number.of.PSMs,
             Master.Protein.Descriptions,
             Protein.Accessions,
             starts_with("Abundances.Normalized"), #Marked.as, 
             starts_with("Abundance.Ratio.P.Value"),
             starts_with("Abundance.Ratio.log2")) %>%
      rename_with(~ exp_design, starts_with("Abundances.Normalized")) %>%
    mutate(species=background_species)
  
  filtered_abundances<-quant_phospho_peptides[rowSums(!is.na(select(quant_phospho_peptides,starts_with(exp_design))))>0,]
  filtered_abundances_ecoli <-quant_peptides_ECOLI[rowSums(!is.na(select(quant_peptides_ECOLI,starts_with(exp_design))))>0,]
  
  # 
  # # Calculate 1 percent quantile of each sample
  # triplicate_indx <- c(1,3,4,6,7,9,10,12,13,15)
  # impute_values_2NA <- NULL
  # impute_values_1NA <- NULL
  # tmp_2NA <- NULL
  # tmp_1NA <- NULL
  # #imputed_abundances <- NULL
  # for (i in 1:5){
  #   tmp_1NA <- quantile (abundances_for_impute[,experiment_name[triplicate_indx[1]:triplicate_indx[2]]],
  #                              probs = 0.01 , na.rm = TRUE )
  #   tmp_2NA <- quantile (abundances_for_impute[,experiment_name[triplicate_indx[1]:triplicate_indx[2]]],
  #                        probs = 0.001 , na.rm = TRUE )
  #   impute_values_2NA[i] <- as.numeric(tmp_2NA)
  #   impute_values_1NA[i] <- as.numeric(tmp_1NA)
  #   
  #   triplicate_indx <- triplicate_indx[-c(1:2)]
  #   
  #   # tmp1 <- abundances_for_impute %>% select(contains(paste0("A",i))) %>%
  #   #   mutate(across(contains(paste0("A",i)), ~ifelse(is.na(.), impute_values[i], .)))
  #   # imputed_abundances <- bind_cols(imputed_abundances, tmp1)
  # }
  # 
  
  df_id_pep <- filtered_abundances %>% 
      select(Sequence,Modifications, phospho_pos,phospho_score, starts_with(exp_design),species) %>%
      tibble() %>%
      mutate(pep_with_pos = paste0(Sequence,"_",phospho_pos)) %>%
      # rename(A1_R1= 4, # Using column index to rename the colnames
      #        A1_R2= 5,
      #        A1_R3= 6,
      #        A2_R1= 7,
      #        A2_R2= 8,
      #        A2_R3= 9,
      #        A3_R1= 10,
      #        A3_R2= 11,
      #        A3_R3= 12,
      #        A4_R1= 13,
      #        A4_R2= 14,
  #        A4_R3= 15,
  #        A5_R1= 16,
  #        A5_R2= 17,
  #        A5_R3= 18) %>%
  pivot_longer(cols = starts_with("E2"), 
               values_to = "intensity",
               names_to = "sample_ids",
               values_drop_na = T) %>%
      separate(sample_ids, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F)
  
  barplt_df <- df_id_pep %>%
      #mutate(sample_rep_id_seq = paste(pep_with_pos, Sample_id,Rep_id, sep = "_"))%>%
      group_by(pep_with_pos,sample_ids) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
      slice(which.max(intensity)) %>% ## ELIMINATE MULTIPLE CHARGES
      ungroup()
  
  ## SAME STRATEGIES ABOVE (3rd) WAS APPLIED TO BACKGROUND AS WELL
  barplt_df_ecoli <- filtered_abundances_ecoli %>% 
      select(Sequence,Modifications, starts_with(exp_design),species) %>% #Marked.as
      pivot_longer(cols = starts_with("E2"), 
                   values_to = "intensity",
                   names_to = "sample_ids",
                   values_drop_na = T) %>%
      separate(sample_ids, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F) %>%
      mutate(sample_rep_id_seq = paste(Sequence, Sample_id,Rep_id, sep = "_")) %>%
      group_by(sample_rep_id_seq,sample_ids) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
      slice(which.max(intensity)) %>% ## ELIMINATE MULTIPLE CHARGES
      ungroup()
  
  barplt_prot_ecoli <- filtered_abundances_ecoli %>% 
    filter(grepl(background_species,Master.Protein.Descriptions) & !grepl("CON__",Master.Protein.Descriptions)) %>%
    select(Sequence,Modifications,Protein.Accessions, starts_with(exp_design),species) %>% #Marked.as
    pivot_longer(cols = starts_with("E2"), 
                 values_to = "Intensity",
                 names_to = "Experiment",
                 values_drop_na = T) %>%
    separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F) %>%
    #mutate(sample_rep_id_seq = paste(Sequence, Sample_id,Rep_id, sep = "_")) %>%
    group_by(Protein.Accessions,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    slice(which.max(Intensity)) %>%
    ungroup() %>%
    mutate(type=subtitle)
  
  plot16 <- gg_barplt_id_pep_count(data_set = barplt_prot_ecoli,
                                   x_df = barplt_prot_ecoli$Sample_id,
                                   fill_df = barplt_prot_ecoli$Rep_id,
                                   ymax = 20500,
                                   size_num = 10,
                                   header = "Total number of quantified Ecoli across each sample",
                                   caption_lab = "NA values are removed.",
                                   x_lab = "Sample id",
                                   fill_lab =  "Sample id",
                                   y_lab = "Number of identified peptides",
                                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
  write.table(barplt_prot_ecoli, file=paste0(file_path,"Exp2",software_name,"_number_of_unique_",background_species,"proteins_",".txt"),sep = "\t",col.names = T,row.names = F)
  
  
  
  barplt_phospho_seq <- filtered_abundances %>%  select(Sequence,Modifications, starts_with(exp_design),species) %>% #Marked.as
    pivot_longer(cols = starts_with("E2"), 
                 values_to = "intensity",
                 names_to = "sample_ids",
                 values_drop_na = T) %>%
    separate(sample_ids, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F) %>%
    mutate(sample_rep_id_seq = paste(Sequence, Sample_id,Rep_id, sep = "_")) %>%
    group_by(sample_rep_id_seq,sample_ids) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    slice(which.max(intensity)) %>% ## ELIMINATE MULTIPLE CHARGES
    ungroup() %>% mutate(Software_name=software_name) %>%
    mutate(Acquisition_type=acquisiton_type)
    
  write.table(barplt_phospho_seq, file = paste0(file_path,"Number_of_human_phospho_sequences_",
                                                software_name,"_Experiment",exp_id,".txt"),
              sep = "\t",row.names = F)
  
  
  ####### ADDITIONAL PLOT TO DISPLAY MISSING and UNEXPECTED PEPTIDES ########
  df_merge_syn <- barplt_df %>%
    select(pep_with_pos,sample_ids,intensity, species) %>% #Marked.as
    pivot_wider(names_from = "sample_ids",values_from = "intensity") %>%
    full_join(pep_list_w_theo,by="pep_with_pos") %>% 
    mutate_at("Pool", ~replace_na(.,"Unexpected")) %>%
    mutate(Pool= ifelse(is.na(species),"missing",Pool)) %>%
    select(pep_with_pos,starts_with(exp_design),Pool) %>%
    mutate(soft_name=software_name,ion_mobility=acquisiton_type)
  
  ##############################################################################
  ####### GATHERING ALL COLUMNS OF MAIN OUTPUT FROM PROLINE WITH THE CORRECT RESULTS ########
  ### This is necessary only for Proline and PD additionally to compare 
  ## the missing peptides with their scan number.
  
  merge_phospho_peptides <- quant_peptides %>% 
    filter(grepl("Homo sapiens",Master.Protein.Descriptions) & 
             grepl("Phospho",Modifications) & !grepl("CON__",Master.Protein.Descriptions)) %>%
    rowwise() %>%
    mutate(results = list(extract_phospho_numbers(Modifications)),
           phospho_pos = results[[1]],
           phospho_score = results[[2]]) %>%
    select(!results) %>%
    mutate(pep_with_pos = paste0(Sequence,"_",phospho_pos))
 
  df_merge_all_col <- merge_phospho_peptides %>% tibble() %>%
    full_join(pep_list_w_theo,by="pep_with_pos") %>%
    mutate_at("Pool", ~replace_na(.,"Unexpected")) %>%
    mutate(Pool= ifelse(is.na(Protein.Accessions),"missing",Pool)) %>%
    #select(pep_with_pos,starts_with(exp_design),Pool) %>%
    mutate(soft_name=software_name,ion_mobility=acquisiton_type)
  
  
  plot11 <- gg_barplt_id_pep_count(data_set = df_merge_syn,
                                x_df = df_merge_syn$Pool,
                                fill_df = df_merge_syn$Pool,
                                ymax = 20000,
                                size_num = 10,
                                header = "Total number of quantified phospho-site across each sample",
                                caption_lab = "NA values are removed.",
                                x_lab = "Sample id",
                                fill_lab =  "Sample id",
                                y_lab = "Number of identified peptides",
                                subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
  write.table(df_merge_all_col,file = paste0(file_path,"Count_of_missing_unexpected_correct_phospho-sites_with_all_col_",
                                         software_name,"_Experiment",exp_id,".txt"),
              sep = "\t",row.names = F)
   
  write.table(df_merge_syn,file = paste0(file_path,"Count_of_missing_unexpected_correct_phospho-sites_",
                                             software_name,"_Experiment",exp_id,".txt"),
              sep = "\t",row.names = F)
  #############################################################################
  
  plot1 <- gg_barplt_id_pep_count(data_set = barplt_df,
                               x_df = barplt_df$Sample_id,
                               fill_df = barplt_df$Rep_id,
                               ymax = 20000,
                               size_num = 10,
                               header = "Total number of quantified phospho-site across each sample",
                               caption_lab = "NA values are removed.",
                               x_lab = "Sample id",
                               fill_lab =  "Sample id",
                               y_lab = "Number of identified peptides",
                               subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name, subtitle))
  
  
  plot2 <- gg_barplt_id_pep_count(data_set = barplt_df_ecoli,
                               x_df = barplt_df_ecoli$Sample_id,
                               fill_df = barplt_df_ecoli$Rep_id,
                               ymax = 20000,
                               size_num = 10,
                               header = "Total number of quantified Ecoli across each sample",
                               caption_lab = "NA values are removed.",
                               x_lab = "Sample id",
                               fill_lab =  "Sample id",
                               y_lab = "Number of identified peptides",
                               subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name, subtitle))
  
  
  
  plot12 <- gg_barplt_id_pep_count(data_set = barplt_phospho_seq,
                               x_df = barplt_phospho_seq$Sample_id,
                               fill_df = barplt_phospho_seq$Rep_id,
                               ymax = 20000,
                               size_num = 10,
                               header = "Total number of quantified phospho-sequence across each sample",
                               caption_lab = "NA values and multiple sequences are removed.",
                               x_lab = "Sample id",
                               fill_lab =  "Sample id",
                               y_lab = "Number of identified peptides",
                               subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name, subtitle))
  
  
  
  ## Since Multiple Charges were eliminated, number of rows are not the same as before applying pivot_longer()
  barplt_df_wide <- barplt_df %>% 
      select(Sequence,pep_with_pos, species, sample_ids,intensity) %>% #Marked.as
      pivot_wider(names_from = "sample_ids",values_from = "intensity")
  
  barplt_df_ecoli_wide <- barplt_df_ecoli %>% 
      select(Sequence, species, sample_ids,intensity) %>% #Marked.as
      pivot_wider(names_from = "sample_ids",values_from = "intensity")
  
  ## DENSITY PLOT OF BEFORE IMPUTATION 

  abundances_rowMeans <- NULL
  abundances_ecoli_rowMeans <- NULL
  for (k in 1:sample_size){
      # If separate version of row means is not needed, it can be commented later.
      # Separate row Means can be collected in temp object to merge in "log_10_filtered_abundances_rowMeans"
      assign(paste0("abundances_A",k),as.data.frame(rowMeans(barplt_df_wide  %>%
                                                                 select(contains(paste0("A",k)))))) #%>%
      #select(starts_with("abundance_")))))
      abundances_rowMeans<- bind_cols(abundances_rowMeans,get(paste0("abundances_A",k)))
      
      
      assign(paste0("ecoli_abundances_A",k),as.data.frame(rowMeans(barplt_df_ecoli_wide  %>%
                                                                       select(contains(paste0("A",k))))))# %>%
      #select(starts_with("abundance_")))))  
      
      abundances_ecoli_rowMeans<- bind_cols(abundances_ecoli_rowMeans,get(paste0("ecoli_abundances_A",k)))
      
  }
  
  quant_peptides_ECOLI_density_plt <- barplt_df_ecoli_wide %>%
      select(!starts_with("E")) %>%
      bind_cols(abundances_ecoli_rowMeans) %>%
      tibble() %>%
      rename_with(~ paste0("mean_abun",1:sample_size), matches("^row")) %>%
      pivot_longer(cols = starts_with("mean"), 
                   values_to = "intensity",
                   names_to = "sample_ids",
                   values_drop_na = T) %>%
      mutate(sample_id_seq = paste(Sequence, sample_ids, sep = "_"))
  
  ecoli_density_plt<- quant_peptides_ECOLI_density_plt %>%
      select(contains(c("sample_ids","intensity","species"))) #Marked.as
  
  
  colnames(abundances_rowMeans) <- paste0("mean_abun",1:sample_size)
  
  #### MEAN ABUNDANCE RATIO WITH  DENSITY PLOT ####
  ### BEFORE IMPUTATION ###
  quant_phospho_density_plt <- barplt_df_wide %>%
      select(!starts_with("E")) %>%
      bind_cols(abundances_rowMeans) %>% 
      #rename_with(~ paste0("mean_abun",1:sample_size), matches("^row")) %>%
      #rename_with(~ paste0("mean_abun",1:5), matches("^row")) %>%
      tibble() %>% #mutate(pep_with_pos = sequence) %>% ###  At this stage, no need for phospho-position#   
      pivot_longer(cols = starts_with("mean"),
                   names_to = "sample_ids",
                   values_to = "intensity",
                   values_drop_na = T) %>%
      mutate(sample_id_seq = paste(pep_with_pos, sample_ids, sep = "_"))
  
  density_df <-quant_phospho_density_plt %>%
      select(c(sample_ids,intensity,species)) %>% #Marked.as
      bind_rows(ecoli_density_plt) # %>%
  #separate(accession, into = c("prot_id","species"),remove = F) ### Go back to beginning of inital data creation and continue with "Mark.as" column to use species column here 
  
  
  plot3 <- gg_density(data_set = density_df, 
                   x_df = density_df$intensity,
                   fill_df = density_df$species,
                   color_df = NULL,
                   header="Distribution of mean abundance of every sample before imputation",
                   facet_df = "sample_ids",
                   x_lab = "log10(intensities)",
                   color_lab= "",
                   fill_lab = "species",
                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name, subtitle))
  
  library(kableExtra)
  na_phospho_selected <- apply(X = is.na(barplt_df_wide %>% select(Sequence, species,starts_with("E2"))), MARGIN = 2, FUN = sum)
  na_ecoli <- apply(X = is.na(barplt_df_ecoli_wide %>% select(Sequence, species,starts_with("E2"))), MARGIN = 2, FUN = sum)
  
  na_table <- bind_rows(na_phospho_selected,na_ecoli)
  na_table$species <- c(selected_spcies,background_species)
  na_table$total <- c(dim(barplt_df_wide)[1],dim(barplt_df_ecoli_wide)[1])
  
  na_table %>% select(exp_design,total,species) %>%
      kbl(caption = paste("Number of NA values across all samples \n Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) %>%
      kable_material(c("striped", "hover")) %>%
      kable_styling(bootstrap_options = "striped", full_width = F, position = "left", font_size = 12) %>%
      kable_minimal(full_width = F) %>%
      footnote(general = paste("This table was created after elimination of multiple charges by selecting either phospho-sites of Homo sapiens and", background_species, "sequences \n that has the highest abundace."),
               # number = c("Footnote 1; ", "Footnote 2; "),
               # alphabet = c("Footnote A; ", "Footnote B; "),
               # symbol = c("Footnote Symbol 1; ", "Footnote Symbol 2")
               footnote_as_chunk = T, title_format = c("italic", "underline")) %>%
      #as_image(width = 8) %>%
      save_kable(paste0(file_path,"outputs_with_new_script/table1.png"))
  
  ## REMOVE SEQUENCE COLUMN AFTER NA TABLE
  barplt_df_wide <- barplt_df_wide %>% select(!Sequence)
  
  
  abundances_for_before_impt <- barplt_df_wide %>%
      bind_rows(barplt_df_ecoli_wide) %>% select(exp_design)
  
  # Calculate 1 percent quantile of each sample
  impute_values <- apply(abundances_for_impute, 2 , quantile , probs = 0.01 , na.rm = TRUE,type=6 ) # Type =6 is to obtain the same result as Excel.
  
  # Impute missing values
  for (j in 1:length(impute_values)){
      # Number NA
      #num_NA <- length(abundances_for_impute_all[,j+2][is.na(abundances_for_impute_all[,j+2])])
      
      barplt_df_wide[,j+2][is.na(barplt_df_wide[,j+2])] <- impute_values[j]
      barplt_df_ecoli_wide[,j+2][is.na(barplt_df_ecoli_wide[,j+2])] <- impute_values[j]
      
      #abundances_for_impute_all[,j+2][is.na(abundances_for_impute_all)[,j+2]] <- impute_values[j]
      # After imputation number of imputed values
      #num_imp <-length(abundances_for_impute_all[,j+2][(abundances_for_impute_all[,j+2]==impute_values[j])])
      
      # This is verification of imputation is done successfully
      # Because we expect to see that number of imputed values should be the same amount as number of NA
      #print(setequal(num_NA,num_imp))
      #print(num_NA)
      #print(num_imp)
  }
  
  # UNNECESSARY TO KEEP DATA BEFORE IMPUTATION
  abundances_all_aft_imputation <- barplt_df_ecoli_wide %>%
      rename_with(~ paste0("pep_with_pos"), matches("^Seq")) %>%
      bind_rows(barplt_df_wide) 
  
  
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
  
  final_imputed_data <- cbind(abundances_all_aft_imputation, filtered_abundances_rowMeans,filtered_abundances_log10,filtered_abundances_log10_rowMeans) #filtered_abundances
  
  final_imputed_data_syn <- final_imputed_data %>% filter(grepl(selected_spcies, species))
  
  final_imputed_data_ecoli <- final_imputed_data  %>% filter(!grepl(selected_spcies, species))
  
  df_merge <- final_imputed_data_syn %>%
      left_join(pep_list_w_theo,by="pep_with_pos") %>% 
      mutate_at("Pool", ~replace_na(.,"Unexpected")) %>%
      bind_rows(final_imputed_data_ecoli) %>%
      mutate_at("Pool", ~replace_na(.,background_species)) 
  #write.table(final_imputed_data, file = "final_imputed_normalized_data_PAL _T_cell_Exp3_( 5 conc 3reps)_NoFAIMS_DDA_with_cont_230206_2023-02-07_0947.txt",sep = "\t",row.names = F)
  
  df_mean_ab_after_impt <- df_merge %>% 
      select(contains("aft_imp") | contains("species"),Pool) %>%
      tibble() %>% 
      #separate(accession, into = c("uniprot_id", "species"), remove = F) %>%
      pivot_longer(cols = contains("aft_imp"),
                   names_to = "Mean_abundance",
                   values_to = "values")
  
  
  df_FC_ratio_after_impt <- df_merge %>% 
      select(starts_with("exp_")| contains("species"),Pool) %>%
      #separate(accession, into = c("uniprot_id", "species"), remove = F) %>%
      tibble() %>% 
      pivot_longer(cols = starts_with("exp_"),
                   names_to = "exp_FC",
                   values_to = "values")
  
  plot4 <- gg_density(data_set = df_mean_ab_after_impt, 
                   x_df = df_mean_ab_after_impt$values,
                   fill_df = df_mean_ab_after_impt$Mean_abundance,
                   color_df = df_mean_ab_after_impt$Pool,
                   header="Distribution of mean abundance of every sample after imputation",
                   facet_df = "Mean_abundance",
                   x_lab = "log10(values)",
                   color_lab= "",
                   fill_lab = "Sample Names",
                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name, subtitle))
  
  write.table(final_imputed_data, 
              file =paste0(file_path,"/outputs_with_new_script/final_imputed_data_PD_",
                           exp_id, 
                           acquisiton_type,".txt"),sep = "\t",row.names = F)
  
  
  ############################## ############################
  
  ### RATIO SUPPRESION ASSESSMENT ### 
  ### For this calculation, wrong localizations were removed. 
  
  pep_list_w_theo_sel <- pep_list_w_theo %>%
    select(pep_with_pos,isomericity,Pool,pool_id) 
  
  ratio_supp_FC <-df_merge %>% 
    select(starts_with("exp_") | contains("species"),pep_with_pos) %>%
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
                      subtitle_txt = "Spiked Pool")
  
  p16 <-  gg_raincloud(data_set = ratio_supp_FC_Fixed,
                       x_df = ratio_supp_FC_Fixed$exp_FC,
                       y_df = ratio_supp_FC_Fixed$values,
                       fill_df = ratio_supp_FC_Fixed$Pool,
                       header = "Distribution of mean abundance of every sample after imputation",
                       x_lab = "Sample Names",
                       y_lab = " Density of log10(Mean Abundance)",
                       fill_lab = "Sample Names",
                       caption_lab = "",
                       subtitle_txt = "Fixed Pool")
  
  
  library(patchwork)
  plot17 <- p15/p16
  
  ############################## ############################
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
                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name, subtitle))
  ### BOX-PLOT: Experimental Quantity Ratio of Phospho Peptides  
  
  plot6 <- gg_boxplt_exp_ratio(data_set = df_FC_ratio_after_impt, 
                            x_df = df_FC_ratio_after_impt$exp_FC,
                            y_df = df_FC_ratio_after_impt$values,
                            fill_df = df_FC_ratio_after_impt$Pool,
                            header="Experimental Quantity Ratio of Phospho Peptides",
                            x_lab="Sample Names",
                            y_lab="Abundance Ratios",
                            fill_lab = "Sample Names",
                            subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name, subtitle))
  
  
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
                                 subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name, subtitle))
  
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
                            subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name, subtitle))
  
 
  # filter_at(vars(Pool), all_vars(is.na(.)))
  
  #output_path <- "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/comparision_DDA_noFAIMS_mq_proline_pd/"
  #write.xlsx(unexpectedly_id_peps,file = paste0(output_path,"Corrected_imputation_proline_unexpectedly_identified_phosphopeptides_exp2_noFAIMS_DDA.xlsx"))
  
  #test_data_one_roc <- final_imputed_data %>% select(common_col_for_merging,pool_id)
  # 
  # stat_analysis <- final_imputed_data %>%
  #   select(common_col_for_merging | starts_with("log10_E2"))
  
  #### STATISTICAL PART: ####
  ## T-TEST
  
  ttest_func <- function(x, y) {
      # if (sum(!is.na(x)) < 2 | sum(!is.na(y)) < 2) {
      #   return(NA)
      # }else{
      #   
      # }
      t.test(x, y,alternative = c("two.sided"),var.equal=TRUE)$p.value
  }
  
  wilcox_func <- function(x, y) {
      # if (sum(!is.na(x)) < 2 | sum(!is.na(y)) < 2) {
      #   return(NA)
      # }else{
      #   
      # }
      wilcox.test(as.numeric(x), as.numeric(y),alternative = c("two.sided"))$p.value
  }
  
  ## TODO: ADD LIMMA
  stat_analysis <- df_merge %>%
      select(pep_with_pos,Pool,isomericity,starts_with("log10_") | starts_with("mean_log10_") | starts_with("exp_FC")) %>%
      filter(!grepl(background_species, Pool))
  #filter(grepl("HUMAN",accession)) #spectrum_title
  
  rownames(stat_analysis) <- paste0(stat_analysis$pep_with_pos,"@",stat_analysis$Pool,"@",(1:nrow(stat_analysis))) #stat_analysis$spectrum_title,"@"
  
  library(multtest)
  
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
                      select(stat_analysis, contains(paste0("A", i, "_")) & contains("log10_"))[j,])
        }
      })
      
      # Extract p-values from the list
      p_values_tmp <- sapply(p_values_tmp, function(x) x)
      
      # Store p-values in the pre-allocated matrix
      all_pvalues[, i - 1] <- p_values_tmp
      
      # Perform adjustment
      adjust_pval_tmp <- mt.rawp2adjp(p_values_tmp, proc = "BH", alpha = 0.05)
      qval <- data.frame(adjust_pval_tmp$adjp, adjust_pval_tmp$index)[order(adjust_pval_tmp$index), 2]
      
      # Store adjusted p-values in the pre-allocated matrix
      all_adjust_pval[, i - 1] <- qval
    }
    
    # Create data frames from matrices
    all_pvalues <- as.data.frame(all_pvalues)
    colnames(all_pvalues) <- paste0("pvalues_A",numerator,"/", "A", 2:sample_size)
    rownames(all_pvalues) <- row.names(stat_analysis)
    
    all_adjust_pval <- as.data.frame(all_adjust_pval)
    colnames(all_adjust_pval) <- paste0("adjust_pval_A",numerator,"/", "A", 2:sample_size)
    rownames(all_adjust_pval) <- row.names(stat_analysis)
    
    
    all_pvalues_common_col <- stat_analysis %>% 
      select(pep_with_pos, isomericity, Pool) %>% #spectrum_title
      bind_cols(all_pvalues) %>%  #all_adjust_pval
      pivot_longer(cols = starts_with("pvalues"), values_to = "pval", names_to ="pratios") %>% #adj_pval and adj_pratios
      separate(pratios, into = c("tmp","ratio","tmp1"),sep = "_") %>%
      select(!c(tmp,tmp1)) %>%
      mutate(common_col = paste(pep_with_pos,Pool,ratio,1:((sample_size-1)*nrow(stat_analysis)),sep="@")) #spectrum_title
    
    ## THE BEST WAY TO DO is this:
    merge_stat_df <- stat_analysis %>%
      select(pep_with_pos, Pool, isomericity,starts_with("exp_FC")) %>%
      pivot_longer(cols = starts_with("exp_FC"), values_to = "fold_change_values", names_to ="fold_change_ratios") %>%
      separate(fold_change_ratios, into = c("tmp","tmp1","ratio"),sep = "_") %>%
      select(!c(tmp,tmp1)) %>%
      #mutate(common_col = paste(pep_with_pos,Pool,1:((sample_size-1)*nrow(stat_analysis)),sep="@")) %>%
      bind_cols(all_pvalues_common_col$ratio,all_pvalues_common_col$pval) %>%
      rename_with(.col =6 , ~"fold_change_comp") %>%
      rename_with(.col=7, ~ "P.Value") %>%
      #filter(!grepl("ECOLI",Pool))
      #separate(accession, into = c("prot_id","species"),sep = "_")
      mutate(isomericity = ifelse(is.na(isomericity), "False Positive", isomericity)) %>%
      unite(Pool_new, Pool, isomericity,sep = "_",remove = FALSE) %>%
      unite('new_col_coloring',Pool_new,ratio,sep = "_",remove = FALSE) %>%
      mutate(new_col_coloring = if_else(grepl("Fixed_mono", new_col_coloring), "Fixed_mono", new_col_coloring)) %>%
      mutate(new_col_coloring = if_else(grepl("Fixed-multi", new_col_coloring), "Fixed-multi", new_col_coloring)) %>%
      mutate(new_col_coloring = if_else(grepl("Unexpected", new_col_coloring), "Unexpected", new_col_coloring))
    
    
    
    ###############################################################################
    ## Generation of df -> expected abundance ratio for volcano plot
    merge_stat_df_final <- merge_stat_df
    
  }else if(test_type=="limma"){
    
    library(limma)
    design_matrix <- model.matrix(~factor(c(rep(2,num_reps),rep(1,num_reps))))
    merge_stat_df <-NULL
    for ( i in 2:sample_size){
      # Change only the colname iteratively makes fit to every comparison
      colnames(design_matrix) <- c("Intercept", paste0("A",numerator,"-A",i))
      #print(colnames(design_matrix))
      # Col selection for each comparison
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
      unite(Pool_new, Pool.x, isomericity,sep = "_",remove = FALSE) %>%
      unite('new_col_coloring',Pool_new,ratio,sep = "_",remove = FALSE) %>%
      mutate(new_col_coloring = if_else(grepl("Fixed_multi", new_col_coloring), "Fixed_multi", new_col_coloring)) %>%
      mutate(new_col_coloring = if_else(grepl("Fixed_mono", new_col_coloring), "Fixed_mono", new_col_coloring)) %>%
      mutate(new_col_coloring = if_else(grepl("Unexpected", new_col_coloring), "Unexpected", new_col_coloring)) %>%
      rename(pep_with_pos=pep_with_pos.x) %>% rename(Pool=Pool.x)
    
  }else{
    print("Statistical test could not be assessed. Check the input files!")
  }
  
  
  ## Generation of df -> expected abundance ratio for volcano plot
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
  point_count_y_axis <- merge_stat_df_final %>%
    group_by(fold_change_comp, new_col_coloring) %>%
    filter(P.Value < 0.05) %>% 
    count(new_col_coloring) %>% left_join(actual_ratio_col)
  
  ymax <- max(-log10(merge_stat_df_final$P.Value)) + 0.5
  y_decrement <- 0.5
  
  
  calculate_y_pos <- function(group) {
    group_length <- length(group)
    y_pos <- ymax - seq(0, by = y_decrement, length.out = group_length)
    return(y_pos)
  }
  
  
  # Apply the function to calculate y_pos within each group
  #point_count_y_axis$y_pos <- unlist(by(point_count_y_axis$A1vs_Ai, point_count_y_axis$A1vs_Ai, calculate_y_pos))
  
  point_count_y_axis$y_pos <- unlist(by(point_count_y_axis$fold_change_comp, 
                                        point_count_y_axis$fold_change_comp, calculate_y_pos))
  
  
  col <- RColorBrewer::brewer.pal(n=length(comparisons),name = "Dark2")
  
  labels <- unique(merge_stat_df_final$new_col_coloring)
  
  if (sheet_theo_name == "ISOREF_REF2_Others"){
    
    colors <- c(RColorBrewer::brewer.pal(n=length(comparisons),name = "Dark2"),"#2171b5","#999999")
    new_comparisons <- c(comparisons,"Unexpected","Fixed_mono")
    
  }else if (sheet_theo_name == "ISO-refOTHER with FC_correct"){
    
    new_comparisons <- c(comparisons,"Unexpected")
    colors <- c(RColorBrewer::brewer.pal(n=length(comparisons),name = "Dark2"),"#999999")
    
  }else{
    print("Sheet_theo_name could not be found, please make sure that you selected the correct sheet_name.")
  }
  
  mapped_coloring <- rep("#000000",length(labels))
  
  for (i in 1:length(labels)) {
    # Check if the color_element contains any of the comparisons
    if (any(new_comparisons %in% str_extract_all(labels[i], paste(new_comparisons, collapse = "|"))[[1]])) {
      # Find the index of the matching comparison in the comparisons list
      comp_index <- match(TRUE, sapply(new_comparisons, function(comp) comp %in% str_extract_all(labels[i], comp)))
      
      # Assign the corresponding color to the data frame
      mapped_coloring[i]<- paste0(colors[comp_index])
      #mapped_coloring[i]<- paste(paste0(labels[i],'"'),paste0('"',colors[comp_index]),sep = "=")
    }
  }
  
  merge_stat_df_final <- as.data.frame(merge_stat_df_final)
  actual_ratio_col <- as.data.frame(actual_ratio_col %>% bind_cols(col))
  
  mapped_coloring_dat <- as.data.frame(mapped_coloring)
  
  point_count_y_axis <- mapped_coloring_dat %>% 
    bind_cols(labels) %>%
    rename(colors=1,new_col_coloring=2) %>%
    right_join(point_count_y_axis,by="new_col_coloring")
  
  plot9 <- ggplot(merge_stat_df_final,aes(x =log2(merge_stat_df_final$fold_change_values), y = -log10(merge_stat_df_final$P.Value))) +
    geom_point(size = 3, aes(color = new_col_coloring,shape=Pool)) + # Pool_new might be use to 
    scale_color_manual(values =setNames(mapped_coloring, labels)) +
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
    
    #scale_x_continuous(breaks = seq(from =-7, to=11,by=1)) +
    #scale_y_continuous(breaks = seq(from =0,to=16,by=2)) +
    scale_x_continuous(breaks = seq(from =round(min(log2(merge_stat_df_final$fold_change_values))), to=(round(max(log2(merge_stat_df_final$fold_change_values)))+2),by=2)) +
    scale_y_continuous(breaks = seq(from =round(min(-log10(merge_stat_df_final$P.Value))), to=(round(max(-log10(merge_stat_df_final$P.Value)))+2),by=2)) +
    
    theme(legend.text = element_text(size = 30),
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30),
          plot.title = element_text(size = 35),
          legend.title = element_text(size = 30),
          axis.text.x = element_text(size = 30),
          axis.title = element_text(size = 30),
          axis.text.y = element_text(size = 30),
          plot.subtitle = element_text(size = 30)) +
    labs( y= "-log10(p values)", x="log2(fold change)",title = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name), subtitle = paste("Limma was used \n", subtitle)) +
    geom_vline(data = actual_ratio_col, aes(xintercept = log2(actual_ratio_val), show.legend = FALSE),color=col,size=1.5) +
    geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "red",size=1.5) + 
    geom_label(data = point_count_y_axis, aes(x = log2(actual_ratio_val), y = y_pos, fill=new_col_coloring,label = n),color="white",size=12,show.legend = FALSE) +
    scale_fill_manual(values =setNames(mapped_coloring, labels)) 
  
  
  
  # p9 <- ggplot(merge_stat_df_final,aes(x =log2(merge_stat_df_final$fold_change_values), y = -log10(merge_stat_df_final$P.Value))) +
  #     geom_point(aes(color = new_col_coloring,shape=Pool_new), size = 4) +
  #     #geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "red") +
  #     scale_fill_manual(values = c("ISO-REF" = "#000000",
  #                                  "Unexpected_False Positive_A1/A2"="#999999",
  #                                  "Unexpected_False Positive_A1/A3" ="#999999",
  #                                  "Unexpected_False Positive_A1/A4"="#999999",
  #                                  "Unexpected_False Positive_A1/A5"="#999999",
  #                                  "Others_multi_A1/A2" = "#CC79A7",
  #                                  "Others_mono_A1/A2" = "#CC79A7",
  #                                  "Others_multi_A1/A3" = "#E69F00",
  #                                  "Others_mono_A1/A3" = "#E69F00",
  #                                  "Others_multi_A1/A4" = "#56B4E9",
  #                                  "Others_mono_A1/A4" = "#56B4E9",
  #                                  "Others_multi_A1/A5" = "#009E73",
  #                                  "Others_mono_A1/A5" = "#009E73")) + 
  #     #geom_line(aes(color = new_col_coloring), size = 1) +  # Add color aesthetic to geom_line()
  #     scale_color_manual(values = c("ISO-REF" = "#000000",
  #                                   "Unexpected_False Positive_A1/A2"="#999999",
  #                                   "Unexpected_False Positive_A1/A3" ="#999999",
  #                                   "Unexpected_False Positive_A1/A4"="#999999",
  #                                   "Unexpected_False Positive_A1/A5"="#999999",
  #                                   "Others_multi_A1/A2" = "#CC79A7",
  #                                   "Others_mono_A1/A2" = "#CC79A7",
  #                                   "Others_multi_A1/A3" = "#E69F00",
  #                                   "Others_mono_A1/A3" = "#E69F00",
  #                                   "Others_multi_A1/A4" = "#56B4E9",
  #                                   "Others_mono_A1/A4" = "#56B4E9",
  #                                   "Others_multi_A1/A5" = "#009E73",
  #                                   "Others_mono_A1/A5" = "#009E73"),
  #                        
  #                        labels = c('Non-variant', 'Variant non-isomeric A1 vs A2',
  #                                   'Variant non-isomeric A1 vs A3',
  #                                   'Variant non-isomeric A1 vs A4',
  #                                   'Variant non-isomeric A1 vs A5',
  #                                   'Variant isomeric A1 vs A2',
  #                                   'Variant isomeric A1 vs A3',
  #                                   'Variant isomeric A1 vs A4',
  #                                   'Variant isomeric A1 vs A5',
  #                                   "Unexpected_False Positive A1/A2",
  #                                   "Unexpected_False Positive A1/A3",
  #                                   "Unexpected_False Positive A1/A4",
  #                                   "Unexpected_False Positive A1/A5")) +
  #     scale_shape_manual(values = c(16, 15, 12, 17),
  #                        labels = c('Non-variant', 'Variant non-isomeric', 'Variant isomeric', 'Unexpected')) +
  #     #scale_y_continuous(limits = c(0, max(-log10(merge_stat_df_final$adj.P.Val))), breaks = seq(0, max(-log10(merge_stat_df_final$adj.P.Val)), by = 0.8)) +
  #     #scale_x_continuous(limits = c(min(log2(merge_stat_df_final$fold_change_values)),max(log2(merge_stat_df_final$fold_change_values)))) +#facet_wrap(~ratio) +
  #     #scale_x_continuous(breaks = seq(from =round(min(log2(merge_stat_df_final$fold_change_values))), to=(round(max(log2(merge_stat_df_final$fold_change_values)))+2),by=1)) +
  #     #scale_y_continuous(breaks = seq(from =round(min(-log10(merge_stat_df_final$P.Value))), to=(round(max(-log10(merge_stat_df_final$P.Value)))+2),by=1)) +
  #     #scale_y_continuous(breaks = seq(0, max(-log10(volcano_final1$pvalues_value)), length.out = 21)) +
  #   scale_x_continuous(limits = c(-8, 10),breaks = seq(from = -8, to = 10, by = 2)) +  # Set the ticks for the x-axis
  #   scale_y_continuous(limits = c(0, 10),breaks = seq(from = 0, to = 10, by = 2))+  # Set the ticks for the y-axis
  #   
  #   theme_bw() +
  #     theme(legend.text = element_text(size = 30),
  #           axis.title.x = element_text(size = 30),
  #           axis.title.y = element_text(size = 30),
  #           plot.title = element_text(size = 35),
  #           legend.title = element_text(size = 30),
  #           axis.text.x = element_text(size = 30),
  #           axis.title = element_text(size = 30),
  #           axis.text.y = element_text(size = 30),
  #           plot.subtitle = element_text(size = 30)) +
  #     labs( y= "-log10(p values)", x="log2(fold change)",title = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name), subtitle = paste("Limma was used \n",subtitle)) +
  #     geom_vline(data = actual_ratio_col, aes(xintercept = log2(actual_ratio_val), show.legend = FALSE),color=c("#CC79A7","#E69F00","#56B4E9","#009E73"),size=2) +
  #     geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "red",size=2) + 
  #     geom_label(data = point_count_y_axis, aes(x = log2(actual_ratio_val), y = y_pos,fill=new_col_coloring, label = n),size=14, colour="white",show.legend = FALSE) 
  # 
  # 
  merge_stat_df_final_text <- merge_stat_df_final %>%
    mutate(soft_name=paste0(software_name)) %>%
    mutate(acq_type=paste0(acquisiton_type))
  
  write.table(merge_stat_df_final_text,file = paste0(file_path,"volcano_plot_",software_name,"_",acquisiton_type,".txt"),sep = 
                "\t",col.names = T,row.names = F)
  
  
  df_roc <- merge_stat_df_final %>%
    select(pep_with_pos, Pool,P.Value)
  
  ### ROC analysis custom func
  df_roc_order <- df_roc[order(df_roc$P.Value),]
  
  df_roc_func <-  compute_roc_curve(df=df_roc_order, flag = "Spiked",expected =length(comparisons)*size_variying_pep_size)
  
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
  
  write.table(df_roc_func, file = paste0(file_path,"outputs_with_new_script/new_custom_Roc_analysis_",exp_id,"_",software_name,"_",".txt"),sep = "\t",row.names = F)

  #### ROC Analysis using pROC 
  
  df_roc$variant <- ifelse(df_roc$Pool == "Spiked", TRUE, FALSE)
  #df_roc$non_var <- ifelse(df_roc$Pool == "ISO-REF", TRUE, FALSE)
  
  library(pROC)
  # Calculate ROC curve for raw p-values
  roc_raw_variant <- roc(df_roc$variant, df_roc$P.Value)
  #tpr_and_fpr_variant  <- cbind(roc_raw_variant$sensitivities,
  #                              roc_raw_variant$specificities,
  #                              "Variant Pool")

  fpr <- as.data.frame(1 - roc_raw_variant$specificities)
  tpr_and_fpr_variant  <- cbind(roc_raw_variant$sensitivities,
                                fpr,#roc_raw_variant$specificities,
                                "Variant Pool")
  
  
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
      geom_line(size=1.5) + theme_bw() + # +  scale_x_reverse()
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
      
  write.table(roc_plt_df, file = paste0(file_path,"outputs_with_new_script/pRoc_analysis_",exp_id,"_",software_name,"_",".txt"),sep = "\t",row.names = F)
  
  plt_obj <- ls(pattern="plot")
  plot_obj <- plt_obj[!is.na(plt_obj)]
  sapply(1:length(plot_obj),function(x) ggsave(filename = paste0("p",x,".png"),
                                               width = 80, height = 60, 
                                               path = paste0(file_path,"/outputs_with_new_script/"),
                                               units = "cm",
                                               get(plot_obj[x]),
                                               device = "png", #".svg"
  ))
  
  
  
  # 
  # 
  # 
  # roc_curve_data_exp2_noFAIMS_MQ <- read.delim("D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/MQ_data_analysis/exp2_wo_FAIMS/exp2_wo_FAIMSwith_MBR/roc_curve_data_exp2_noFAIMS_MQ.txt")
  # roc_proline_mq <- rbind(roc_curve_data_exp2_noFAIMS_MQ,roc_final)
  # 
  # 
  # roc_data_unexpected121 <- cbind(roc_data_unexpected12, "A1 vs A2","unexcepted")
  # roc_data_unexpected131 <- cbind(roc_data_unexpected13, "A1 vs A3","unexcepted")
  # roc_data_unexpected141 <- cbind(roc_data_unexpected14, "A1 vs A4","unexcepted")
  # roc_data_unexpected151 <- cbind(roc_data_unexpected15, "A1 vs A5","unexcepted")
  # 
  # colnames(roc_data_unexpected121)[c(4,5)] <- c("Comparison","Pool")
  # colnames(roc_data_unexpected131)[c(4,5)] <- c("Comparison","Pool")
  # colnames(roc_data_unexpected141)[c(4,5)] <- c("Comparison","Pool")
  # colnames(roc_data_unexpected151)[c(4,5)] <- c("Comparison","Pool")
  # 
  # 
  # 
  # roc_final2 <- rbind(roc_data_unexpected121,roc_data_unexpected131,
  #                     roc_data_unexpected141,roc_data_unexpected151)
  # 
  # 
  # 
  # roc_final <- rbind(roc_data121,roc_data131,roc_data141,roc_data151)
  # write.table(roc_final,file = paste0(file_path,"/","ROC_curve_data_PD_target_decoy_noFAIMS_consider_isoref_false.txt"),sep = "\t",row.names = FALSE)
  # line_types <- c("12345678","dotted","solid","dashed" )
  # 
  # ##TODO: make it more professional
  # ggplot(roc_final,aes(x=fdp,y=tpr,color=Comparison,linetype=Software)) + 
  #   geom_line(size=1,values=line_types) +
  #   scale_color_manual(values = c("#CC79A7", "#E69F00", "#56B4E9", "#009E73"),
  #                      labels=c('A1 vs A2','A1 vs A3','A1 vs A4','A1 vs A5')) +
  #   scale_linetype_manual(values = line_types)+
  #   theme_bw() +
  #   theme(legend.text = element_text(size=15), 
  #         axis.title.x = element_text(size = 15),
  #         axis.title.y = element_text(size = 15),
  #         plot.title = element_text(size=30),
  #         legend.title=element_text(size=15),
  #         axis.text.x=element_text(size=15),
  #         axis.title=element_text(size=15),
  #         axis.text.y = element_text(size = 15)) + 
  #   expand_limits(x = 0, y = 0) +
  #   expand_limits(x = 0, y = 0) +
  #   scale_y_continuous(limits = c(0,100)) +
  #   #facet_wrap(~Ratio_col,scales = "free_x") +
  #   labs(title = "Experiment 2 - DDA no FAIMS - \n Comparision of Softwares")
  # 
  # ##TODO: make it more professional
  # ggplot(roc_data14,aes(x=fdp,y=tpr)) + geom_line()
  # 
  # 
  # roc_curve_combine_proline_mq_score_0_40 <- read.delim("D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/comparision_DDA_noFAIMS_mq_proline_pd/comparison_wrong_localized_pep_MQ_Proline_with_out_FAIMS/final_figures/roc_curve_combine_proline_mq_score_0_40.txt")
  # 
  # 
  # 
  # 
  # 
  # 
  # merging_all_pvalues <- final_imputed_data %>% 
  #   select(col_sel_stat_analysis) %>%
  #   left_join(p_values_12,by="common_col_for_merging") %>%
  #   left_join(p_values_13,by="common_col_for_merging") %>%
  #   left_join(p_values_14,by="common_col_for_merging") %>%
  #   left_join(p_values_15,by="common_col_for_merging")
  # 
  # merging_all_pvalues <- merging_all_pvalues %>%
  #   mutate(Pool=replace_na(Pool,"unexpected"))
  # 
  # 
  # 
  # 
  # 
  #   
  # # 
  # # 
  # # 
  # # 
  # # 
  # # #adjusted_p_values <- as.data.frame(p.adjust(p_values_A1_vs_A5, method = "BH", n = length(p_values_A1_vs_A5)))
  # # adjusted_p_values <- lapply(stat_analysis, function(x) p.adjust (x, method = "BH", n = length(x)))
  # # p_values_12[,"adj_pvalue12"] <- p.adjust(p_values_12$p_values_12, method = "BH")
  # # 
  # # all_adjusted_p_values <- cbind(unlist(adjusted_p_values[["p_values_12"]]),
  # #                                unlist(adjusted_p_values[["p_values_13"]]),
  # #                                unlist(adjusted_p_values[["p_values_14"]]),
  # #                                unlist(adjusted_p_values[["p_values_15"]]))
  # # colnames(all_adjusted_p_values) <- c("adj_pvalues_1_2","adj_pvalues_1_3","adj_pvalues_1_4","adj_pvalues_1_5")
  # # 
  # # 
  # 
  # ##### p_values_for_all_ratio changed it later to adjusted!!!!
  # df_volcano <- final_imputed_data %>% 
  #   select(!contains("mean_abundance")) %>%
  #   mutate(as.data.frame(p_values_for_all_ratio)) %>%
  #   #mutate(across(everything(), ~replace(., is.infinite(.), NA)))
  #   rowwise() %>%
  #   filter(across(where(is.numeric), ~!is.infinite(.)))
  # 
  # 
  # pvalues_df <- merging_all_pvalues %>%
  #   pivot_longer(cols = starts_with("p_values"), names_to = "pvalues_col", values_to = "pvalues_value") %>%
  #   select(common_col_for_merging,Pool,isomericity,pvalues_col,pvalues_value)
  # 
  # # Pivot the columns containing "Ratio"
  # ratio_df <- merging_all_pvalues %>%
  #   pivot_longer(cols = ends_with("_Ratio"), names_to = "Ratio_col", values_to = "Ratio_value") %>%
  #   select(Ratio_col,Ratio_value)
  # 
  # # Pivot the columns containing "exp_FC"
  # exp_FC_df <- merging_all_pvalues %>%
  #   pivot_longer(cols = starts_with("exp_FC"), names_to = "exp_FC_col", values_to = "exp_FC_value") %>%
  #   select(exp_FC_col,exp_FC_value)
  # 
  # volcano_final <- bind_cols(pvalues_df,ratio_df,exp_FC_df)
  # 
  # actual_ratio <- c(2,10,20,100)
  # actual_ratio <- data.frame(Ratio_col=unique(volcano_final$Ratio_col),log2(c(2,10,20,100)))
  # 
  # log10_p_thresholds <- data.frame(ratio_col=unique(volcano_final$Ratio_col), -log10(p_thresholds))
  # 
  # 
  # #### DETERMINE THRESHOLD OF EACH P VALUE COMPARISON BASED on 
  # ####    WHERE WE CAN SEE THE CHANGE FROM TRUE to FALSE
  # 
  # ## THIS WILL BE ADDED AS AN EXTRA LINE FOR EACH COMPARISION 
  # ##     WHEN I MERGED ALL VOLCANO PLOTS INTO ONE
  # 
  # volcano_final1 <- volcano_final %>% select(-Ratio_value) %>%
  #   mutate(isomericity = ifelse(is.na(isomericity), "False Positive", isomericity)) %>%
  #   unite(Pool_new, Pool, isomericity,sep = "_",remove = FALSE) %>%
  #   unite('new_col_coloring',Pool_new,Ratio_col,sep = "_",remove = FALSE) %>%
  #   mutate(new_col_coloring = if_else(grepl("ISO-REF", new_col_coloring), "ISO-REF", new_col_coloring)) %>%
  #   mutate(new_col_coloring = if_else(grepl("unexpected", new_col_coloring), "unexpected", new_col_coloring))
  # 
  # write.table(volcano_final1,file = "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/PD_data_analysis/target_decoy_no_FAIMS/new/volcano_plot_output_PD.tsv",
  #             sep = "\t",row.names = F,col.names = T)
  # 
  # ggplot(volcano_final1, aes(x = log2(exp_FC_value), y = -log10(pvalues_value))) +
  #   geom_point(aes(shape = Pool_new, color = new_col_coloring), size = 2.5) +
  #   #geom_line(aes(color = new_col_coloring), size = 1) +  # Add color aesthetic to geom_line()
  #   scale_color_manual(values = c("ISO-REF" = "#000000", "unexpected" = "#999999",
  #                                 "Others_multi_A1-A2_Ratio" = "#CC79A7",
  #                                 "Others_mono_A1-A2_Ratio" = "#CC79A7",
  #                                 "Others_multi_A1-A3_Ratio" = "#E69F00",
  #                                 "Others_mono_A1-A3_Ratio" = "#E69F00",
  #                                 "Others_multi_A1-A4_Ratio" = "#56B4E9",
  #                                 "Others_mono_A1-A4_Ratio" = "#56B4E9",
  #                                 "Others_multi_A1-A5_Ratio" = "#009E73",
  #                                 "Others_mono_A1-A5_Ratio" = "#009E73"),
  #                      labels = c('Non-variant', 'Variant non-isomeric A1 vs A2',
  #                                 'Variant non-isomeric A1 vs A3',
  #                                 'Variant non-isomeric A1 vs A4',
  #                                 'Variant non-isomeric A1 vs A5',
  #                                 'Variant isomeric A1 vs A2',
  #                                 'Variant isomeric A1 vs A3',
  #                                 'Variant isomeric A1 vs A4',
  #                                 'Variant isomeric A1 vs A5',
  #                                 'Unexpected')) +
  #   scale_shape_manual(values = c(16, 15, 12, 17),
  #                      labels = c('Non-variant', 'Variant isomeric', 'Variant non-isomeric', 'Unexpected')) +
  #   scale_y_continuous(limits = c(0, 7.2), breaks = seq(0, 7.2, by = 0.8)) +
  #   scale_x_continuous(limits = c(-7,7)) +
  #   #scale_y_continuous(breaks = seq(0, max(-log10(volcano_final1$pvalues_value)), length.out = 21)) +
  #   theme_bw() +
  #   theme(legend.text = element_text(size = 15),
  #         axis.title.x = element_text(size = 15),
  #         axis.title.y = element_text(size = 15),
  #         plot.title = element_text(size = 30),
  #         legend.title = element_text(size = 15),
  #         axis.text.x = element_text(size = 15),
  #         axis.title = element_text(size = 15),
  #         axis.text.y = element_text(size = 15)) +
  #   expand_limits(x = 0, y = 0) +
  #   geom_vline(data = actual_ratio, aes(xintercept = actual_ratio$log2.c.2..10..20..100..),color=c("#CC79A7","#E69F00","#56B4E9","#009E73"), size = 1, show.legend = FALSE) +
  #   geom_hline(data = log10_p_thresholds, aes(yintercept = log10_p_thresholds$X.log10.p_thresholds.),color=c("#CC79A7","#E69F00","#56B4E9","#009E73"), size = 1, linetype = 2, show.legend = FALSE)+ 
  #   labs(title = "Experiment 2 - DDA no FAIMS processed by Proteome Discoverer", color = "Classes", shape="Type")
  # 
  #   #coord_fixed()
  # 
  # final_imputed_data %>% 
  #   pivot_longer(cols = starts_with("mean_"),
  #                names_to = "Abundance_col",
  #                values_to = "Abundance_value") %>%
  # ggplot(aes(x = Abundance_col, y = common_col_for_merging, fill = Abundance_value)) +
  #   geom_tile(color = "white") +
  #   scale_fill_viridis() +
  #   labs(x = "Condition", y = "Peptide", fill = "Abundance") +
  #   ggtitle("Heatmap")
  # 
  # 
  # 
  # ggplot(df_volcano, aes(x=df_volcano$`exp_FC_A1/A4`, y=-log10(pvalues_1_4),color=Pool)) + 
  #   geom_point(aes(Pool)) + geom_hline(yintercept = pvalue_threshold, size=1) + 
  #   geom_vline(xintercept = actual_ratio[3], size=1)
  # #### Distribution of experimental quantitative ratios with two different pools
  # 
  # 
  # 
  # library(ggplot2)
  # 
  # ## This can be used as an object name
  # df_for_figure <- final_imputed_data %>% 
  #   select(contains(col_sel_stat_analysis)) %>%
  #   pivot_longer(cols = col_sel_stat_analysis[3:11], names_to = "Theo_Exp", values_to = "Ratios",values_drop_na = T) %>%
  #   #select(!contains(c("common_col_for_merging"))) %>%
  #   filter(grepl("exp_FC",Theo_Exp)) %>%
  #   filter_at(vars(Pool), all_vars(!is.na(.))) %>%
  #   filter_at(vars(Ratios), all_vars(!is.infinite(.)))
  # 
  # p1 <- gg_density(data_set = df_for_figure, 
  #                  x_df = df_for_figure$Ratios,
  #                  fill_df = df_for_figure$Theo_Exp,
  #                  color_df = df_for_figure$Pool,
  #                  header="Distribution of experimental quantitative ratios with two different pools",
  #                  facet_df = "Theo_Exp",
  #                  x_lab = "log10(Ratios)",
  #                  color_lab= "Pool",
  #                  fill_lab = "Ratios")
  # 
  # 
  # #### Distribution of mean abundance of every sample with two different pools
  # 
  # figure_with_mean_abundance <- final_imputed_data %>%
  #   select(starts_with("mean_") | contains(c("Pool", "common_col_for_merging"))) %>%
  #   pivot_longer(cols = starts_with("mean"),
  #                names_to = "Theo_Exp",
  #                values_to = "Mean_abundance",
  #                values_drop_na = T) %>%
  #   filter_at(vars(Pool), all_vars(!is.na(.)))
  # 
  # 
  # p2 <- gg_density(data_set = figure_with_mean_abundance, 
  #                  x_df = figure_with_mean_abundance$Mean_abundance,
  #                  fill_df = figure_with_mean_abundance$Theo_Exp,
  #                  color_df = figure_with_mean_abundance$Pool,
  #                  header="Distribution of mean abundance of every sample with two different pools",
  #                  facet_df = "Theo_Exp",
  #                  x_lab = "log10(Ratios)",
  #                  color_lab= "Pool",
  #                  fill_lab = "Ratios")
  # 
  # # final_imputed_normalized_data %>%
  # #   select(starts_with("mean_") | contains(c("Pool", "common_col_for_merging"))) %>%  
  # #   filter_at(vars(Pool), all_vars(!is.na(.))) %>%
  # 
  # #### MANUAL PLOTTING   
  # #
  # #   ggplot(aes(x =log10(mean_abundances_A1) , y = log10(mean_abundances_A5), color=Pool)) +
  # #   geom_point()+
  # #   #facet_wrap(vars(df_for_figure$common_col_for_merging))  + 
  # #   geom_smooth(formula = y ~ x,method = "loess", colour = "green", fill = "green") +
  # #   theme_minimal() +
  # #   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
  # #         axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
  # #         plot.title = element_text(size=20),
  # #         legend.title=element_text(size=15),
  # #         axis.text=element_text(size=15),
  # #         axis.title=element_text(size=15)
  # #   )
  # #### PLOTTING WITH FUNCTION
  # #
  # # gg_density_mean_abun <- function(data_set = final_imputed_normalized_data,
  # #                                  x_df = final_imputed_normalized_data$mean_abundances_A1,
  # #                                  y_df = final_imputed_normalized_data$mean_abundances_A5,
  # #                                  color_df =final_imputed_normalized_data$Pool,
  # #                                  #header,
  # #                                  facet_df= "common_col_for_merging")
  # 
  # ### BOX-PLOT: Experimental Quantity Ratio of Synthetic Peptides  
  # 
  # p3 <- gg_boxplt_exp_ratio(data_set = df_for_figure, 
  #                           x_df = df_for_figure$Theo_Exp,
  #                           y_df = df_for_figure$Ratios,
  #                           fill_df = df_for_figure$Pool,
  #                           header="Experimental Quantity Ratio of Synthetic Peptides",
  #                           x_lab="Sample Names",
  #                           y_lab="Abundance Ratios",
  #                           fill_lab = "Pool")
  # 
  # ### HALF-BOX-PLOT & HALF-SCATTER-PLOT: Experimental Quantity Ratio of Synthetic Peptides  
  # library(gghalves)
  # 
  # p4 <- gg_half_boxplt_exp_ratio(data_set = df_for_figure, 
  #                                x_df = df_for_figure$Theo_Exp,
  #                                y_df = df_for_figure$Ratios,
  #                                fill_df = df_for_figure$Pool,
  #                                header="Experimental Quantity Ratio of Synthetic Peptides",
  #                                x_lab="Sample Names",
  #                                y_lab="Abundance Ratios",
  #                                fill_lab = "Pool")
  # 
  # 
  # 
  # ### VIOLIN-PLOT: Experimental Quantity Ratio of Synthetic Peptides   
  # 
  # ### TODO: fix y scaling without trimming 
  # p5 <- gg_violin_exp_ratio(data_set = df_for_figure, 
  #                           x_df = df_for_figure$Theo_Exp,
  #                           y_df = df_for_figure$Ratios,
  #                           fill_df = df_for_figure$Pool,
  #                           header="Experimental Quantity Ratio of Synthetic Peptides",
  #                           x_lab="Sample Names",
  #                           y_lab="Abundance Ratios",
  #                           fill_lab = "Pool",
  #                           trim=TRUE)
  # 
  # ### SAVE ALL PLOT AUTOMATICALLY (without giving a custom file name)
  # sapply(1:5,function(x) ggsave(filename = paste0("p",x,".tiff"),
  #                               width = 50, height = 40, 
  #                               path = file_path,
  #                               units = "cm",
  #                               get(paste0("p",x)),
  #                               device = "tiff", #".svg"
  # ))
  # 
  # 
  # 
  
}

