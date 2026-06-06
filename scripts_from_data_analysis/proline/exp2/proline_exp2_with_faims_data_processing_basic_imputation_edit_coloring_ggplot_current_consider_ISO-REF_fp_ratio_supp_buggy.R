#library(PhosR)
library(stringr)
library(dplyr)
#library(data.table) # REMOVED: unused
library(openxlsx)
library(tidyr)
library(ggplot2)
library(patchwork)
###############################################
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/ggplot/ggplot_functions.R")
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/roc_curve/new_roc_curve_generation_with_custom_threshold.R")
       #roc_curve_generation_proline_edit.R
# Experiment 2
 #file_path <- "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/exp2_re_injection/"
 #file_name <- "PAL _Phosphopeptides exp2 ( 5 conc 3reps) DDA_230117 with Design_2023-06-07_1003.xlsx"


#file_path <- "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/with_FAIMS/"
#file_name <- "PAL _Exp2_( 5 conc 3reps)_withFAIMS_DDA_25052023_noDesign - correct_2023-06-13_1514.xlsx"
#selected_species = "_HUMAN"


# # Common constant objects
# sample_size <- 5
# sheet_name <- "Best PSM from protein sets"
# #sheet_name <- "Quantified peptide ions"
# exp_id <- 3
# background_species <- "ECOLI"

final_proline_pep_quant_analysis_syn_buggy <- function(file_path,
                                             file_name,
                                             sheet_name,
                                             theo_file_path,
                                             theo_file_name,
                                             sheet_theo_name,
                                             background_species,
                                             selected_species,
                                             loc_filter_opt,
                                             loc_filter,
                                             exp_id,
                                             exp_design,
                                             fdr_threshold,
                                             acquisiton_type,
                                             size_variying_pep_size,
                                             software_name,
                                             test_type,
                                             num_reps,
                                             numerator,
                                             actual_ratio,
                                             subtitle,
                                             dir_name){
  
  
  
  intended_dir <-paste0(file_path,dir_name)
  
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
  

    proline_phospho_pos_extraction <- function(df){
      # Extraction of phospho positions from quant peptides object
      phospho_ptm_pos <- lapply(df, function(each_ptm_protein_positions) {
        
        ptm_list <- as.list(strsplit(each_ptm_protein_positions,"; ", fixed=TRUE)[[1]]) # Split ptm_protein_position depending on ";"
        ptm_list <- ptm_list[grepl("Phospho", ptm_list, fixed = TRUE)] # Extract only which contains "Phospho"
        phospho_positions <- lapply(ptm_list, function(ptm) { 
          #sub('Phospho \\(([A-Z]\\d+)\\)', "\\d+", ptm) #then, remove "Phospho" and remain only positions
          as.character(str_extract(ptm, "\\d+"))
          
        })
        
        phospho_positions_as_str <- paste(phospho_positions, collapse="&") #combine each position with "|"
        
      })
      # Data conversion 
      phospho_ptm_pos_df <- t(as.data.frame(phospho_ptm_pos))
      rownames(phospho_ptm_pos_df) <- 1:length(phospho_ptm_pos_df)
      
      return(phospho_ptm_pos_df)
      
    }
    
    quant_peptides <- read.xlsx(paste0(file_path,file_name), sheet = sheet_name,sep.names = " ")
    # 
    # if (acquisiton_type == "DDA Exploris with FAIMS" | acquisiton_type== "DDA TIMS-TOF"){
    #   
    #   search_info <-  read.xlsx(paste0(file_path,file_name), sheet = "Search settings and infos",sep.names = ".")
    # 
    #   tr_search_info_df <- search_info %>% t 
    #   colnames(tr_search_info_df) <- tr_search_info_df[1,]
    #     
    #   tr_search_info_df <- as.data.frame(tr_search_info_df[-1,])
    #   
    #   mapping_abund_cols <- tr_search_info_df %>%
    #     select(quant_channel_name, result_set_name) %>%
    #     mutate(tmp = str_split(result_set_name, "_") %>% sapply(tail, n = 1)) %>%
    #     mutate(raw_file_names=paste(quant_channel_name,tmp,sep = "_")) %>%
    #     select(!c(result_set_name,tmp)) 
    #   
    #   quant_peptides_longer <- quant_peptides %>% select(sequence,
    #                                                      modifications,
    #                                                      accession,
    #                                                      ptm_score,
    #                                                      ptm_sites_confidence,
    #                                                      ptm_protein_positions,
    #                                                      charge,
    #                                                      spectrum_title,
    #                                                      starts_with("abundance_")) %>%
    #     pivot_longer(cols = starts_with("abundance_"),
    #                  values_to = "intensity",
    #                  names_to = "quant_channel_name") %>%
    #     mutate(tmp=str_remove("abudance",quant_channel_name))
    #     left_join(mapping_abund_cols,by="quant_channel_name")
    #                                                             
    #   
    # }else{}
    
    
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
      distinct(Phosphopeptide.sequence,.keep_all = TRUE) %>%
      select(Sequence,Pool) %>%
      rename(Pool_for_seq_merge=Pool) %>%
      mutate(situation="Correct")
    
    pep_list_w_theo_seq_map <- pep_list_w_theo %>% select(Sequence,Neutral.mass.SH.Cys) %>%
      rename(Neutral_mass=Neutral.mass.SH.Cys) %>%
      mutate(Neutral_mass = round(Neutral_mass, 0)) %>%
      mutate(Neutral_mass_theo=Neutral_mass)
    
    pep_list_w_theo_pep_map <- pep_list_w_theo %>%
      select(pep_with_pos,Pool,isomericity)
    
    pep_list_w_theo_sel <- pep_list_w_theo %>% select(pep_with_pos,Pool,isomericity)
    
    #rechecking_df <- exp2_noFAIMS_DDA_proline %>% filter(grepl("HUMAN",accession)) %>% filter(grepl("Phospho",modifications))
    #write.table(rechecking_df,file="exp2_no_faims_Proline_output_filtered_human_phospho.tsv",sep="\t",row.names = FALSE)
    
    ##### REMOVE REDUNDANCY OF ECOLI SEQUENCE
    # ecoli_pep_max_absum <- exp2_noFAIMS_DDA_proline %>% filter(grepl("ECOLI", accession)) %>%
    #   mutate(across(all_of(sort_col), ~ ifelse(is.na(.), 0, .))) %>%
    #   rowwise() %>%
    #   mutate(row_sum = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
    #   group_by(sequence) %>%
    #   #summarise(max_row_sum = max(row_sum, na.rm = TRUE)) %>%
    #   slice(which.max(row_sum))
    
    ### IMPUTATION
    abundances_for_impute_names <- quant_peptides %>% 
      select(starts_with("abundance"))  ## spectrum_title remove it because it was not make it as rownames (has duplicates)
    
    ordered_abun_cols <- names(abundances_for_impute_names)[order(names(abundances_for_impute_names), decreasing = FALSE)]
    
    abundances_for_impute <- quant_peptides %>%
      select(ordered_abun_cols) %>%
      rename_with(~exp_design,matches("^abundance"))
      
    quant_peptides_cor_abun <- quant_peptides %>% 
      relocate(ordered_abun_cols) %>%
      rename_with(~exp_design,matches("^abundance")) %>%
      rename(Sequence=sequence)
    
    #write.table(quant_peptides_cor_abun,file = paste0(new_path,"/test.tsv"),sep = "\t",col.names = T)
#################################################  
   
    quant_phospho_peptides_tmp <-quant_peptides_cor_abun %>%
      filter(grepl(selected_species,
                   accession) & !grepl("CON__",accession)) %>%
      filter(grepl("Phospho",modifications)) %>% # & !grepl("positions not distinguishable", Modification.Pattern)
      select(Sequence,
             modifications,
             accession,
             ptm_score,
             ptm_sites_confidence,
             ptm_protein_positions,
             charge,
             spectrum_title,
             starts_with(exp_design)) %>%
      rowwise() %>%
      mutate(species=selected_species) %>%
      rowwise() %>%
      mutate(phospho_pos = proline_phospho_pos_extraction(modifications)) %>%
      mutate(Positions=phospho_pos) %>%
      mutate(pep_with_pos=paste0(Sequence,"_",phospho_pos))
      
    
    ######################################   #######################################
    ######################################   #######################################
    df_id_pep <- quant_phospho_peptides_tmp %>% 
      select(Sequence,starts_with("E2"),pep_with_pos,
             accession,modifications,species) %>%
      #rename_with(~ exp_design, starts_with("Abundances.Normalized")) %>%
      pivot_longer(cols = starts_with("E2"), 
                   values_to = "Intensity",
                   names_to = "Experiment",
                   values_drop_na = T) %>%
      separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F)
    
    barplt_df <- df_id_pep %>%
      #mutate(sample_rep_id_seq = paste(pep_with_pos, Sample_id,Rep_id, sep = "_"))%>%
      group_by(pep_with_pos,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
      slice(which.max(Intensity)) %>% ## ELIMINATE MULTIPLE CHARGES
      ungroup()
    
    plot1 <- gg_barplt_id_pep_count(data_set = barplt_df,
                                    x_df = barplt_df$Sample_id,
                                    fill_df = barplt_df$Rep_id,
                                    ymax = 20000,
                                    size_num = 10,
                                    header = "Total number of phospho-peptides across each sample",
                                    caption_lab = "NA values are removed. \n Mapping with the theoretical list was not done yet. \n wrong sequence may include !",
                                    x_lab = "Sample id",
                                    fill_lab =  "Sample id",
                                    y_lab = "Number of identified peptides",
                                    subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name, subtitle))
    
    #######################################  #######################################
    #######################################   #######################################

    ecoli_seq <- quant_peptides_cor_abun %>% 
      #select(sequence, modifications, accession) %>%
      filter(grepl(background_species,accession) & !grepl("CON__",accession)) %>%
      mutate(species=background_species) 
    
    
    ecoli_seq_dist <- quant_peptides_cor_abun %>% 
      select(Sequence, modifications, accession) %>%
      filter(grepl(background_species,accession) & !grepl("CON__",accession)) %>%
      distinct(Sequence,.keep_all = T) %>%
      mutate(species=background_species) 
    
    
    quant_peptides_ECOLI <- quant_peptides_cor_abun %>% 
      filter(grepl(background_species,accession)) %>% 
      select(Sequence,
             modifications,
             accession,
             spectrum_title,
             starts_with(exp_design))
    
    ## SAME STRATEGIES ABOVE (3rd) WAS APPLIED TO BACKGROUND AS WELL
    barplt_df_ecoli <- quant_peptides_ECOLI %>% 
      select(Sequence,modifications, starts_with(exp_design),accession) %>%
      pivot_longer(cols = starts_with("E2"), 
                   values_to = "Intensity",
                   names_to = "Experiment",
                   values_drop_na = T) %>%
      mutate(species=background_species) %>%
      separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F) %>%
      mutate(sample_rep_id_seq = paste(Sequence, Sample_id,Rep_id, sep = "@")) %>%
      group_by(sample_rep_id_seq,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
      slice(which.max(Intensity)) %>%
      ungroup()
    
    plot2 <- gg_barplt_id_pep_count(data_set = barplt_df_ecoli,
                                    x_df = barplt_df_ecoli$Sample_id,
                                    fill_df = barplt_df_ecoli$Rep_id,
                                    ymax = 20500,
                                    size_num = 8,
                                    header = "Total number of quantified Ecoli across each sample",
                                    caption_lab = "NA values are removed.",
                                    x_lab = "Sample id",
                                    fill_lab =  "Sample id",
                                    y_lab = "Number of identified peptides",
                                    subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
    
    
    barplt_prot_ecoli <- quant_peptides_ECOLI %>% 
      filter(grepl(background_species, accession) & !grepl("CON__",accession)) %>%
      select(Sequence,starts_with(exp_design),accession) %>%
      pivot_longer(cols = starts_with("E2"), 
                   values_to = "Intensity",
                   names_to = "Experiment",
                   values_drop_na = T) %>%
      separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F) %>%
      #mutate(sample_rep_id_seq = paste(Sequence, Sample_id,Rep_id, sep = "_")) %>%
      group_by(accession,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
      slice(which.max(Intensity)) %>%
      ungroup() %>%
      mutate(type=subtitle)
    
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
                                     subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
    
    write.table(barplt_prot_ecoli, file=paste0(new_path,"/Exp2",software_name,"_number_of_unique_",background_species,"proteins_",".txt"),sep = "\t",col.names = T,row.names = F)
    
    ####################################################################
    
    amino_acid_table <- read.delim("D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/amino_acid_table.txt")
    mapping_df <- as.data.frame(cbind(amino_acid_table$X1.letter.code,amino_acid_table$Monoisotopic.Mass))
    colnames(mapping_df) <- c("letters","mono_isotopic")
    
    mapping_df$mono_isotopic <- as.numeric(mapping_df$mono_isotopic)
    # Apply the function to each sequence in the sequences data frame
    
    
    quant_phospho_peptides_tmp$Sum <- mapply(calculate_sum, quant_phospho_peptides_tmp$Sequence, 
                                             quant_phospho_peptides_tmp$Positions,
                                             MoreArgs = list(mapping_df = mapping_df))
    
    ######################################################################
    quant_phospho_map_seq <- quant_phospho_peptides_tmp %>%
      mutate(Neutral_mass = round(Sum, 0)) %>%
      mutate(Neutral_mass_res=Neutral_mass) %>%
      pivot_longer(cols = starts_with("E2"), 
                   values_to = "Intensity",
                   names_to = "Experiment",
                   values_drop_na = T) %>%
      full_join(pep_list_w_theo_seq_map,by=c("Neutral_mass","Sequence")) %>% # #pep_with_pos
      
      mutate(map_seq = ifelse(is.na(Neutral_mass_res), "missing", NA)) %>%
      mutate(map_seq = ifelse(is.na(Neutral_mass_theo), "wrong seq", map_seq)) %>%
      mutate(map_seq = ifelse(!is.na(Neutral_mass_res) & !is.na(Neutral_mass_theo) & Neutral_mass_res == Neutral_mass_theo, 
                              "correct seq", 
                              map_seq)) %>%
      filter(grepl("correct seq",map_seq)) %>%
      select(!c(Positions,Sum)) %>%
      rename(Pool_for_seq_merge=map_seq) %>%
      #distinct(Sequence,.keep_all = T) %>%
      group_by(Sequence,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
      slice(which.max(Intensity)) %>%
      ungroup() %>%
      separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F)
    
    
    all_seq <- quant_phospho_map_seq %>% 
      # filter(grepl("Phospho",modifications)) %>%
      # filter(grepl(selected_species,accession) & !grepl("CON__",accession)) %>%
      # #distinct(sequence, .keep_all = T) %>%
      # full_join(pep_list_w_theo_unique,by="sequence") %>%
      # mutate_at("Pool_for_seq_merge", ~replace_na(.,"Unexpected")) %>%
      # mutate(Pool_for_seq_merge= ifelse(is.na(accession),"missing",Pool_for_seq_merge)) %>%
      # filter(!grepl("Unexpected",Pool_for_seq_merge)) %>%
      bind_rows(barplt_df_ecoli) %>% 
      mutate(Pool_for_seq_merge= ifelse(is.na(Pool_for_seq_merge),background_species,Pool_for_seq_merge)) %>%
      separate(accession, into = c("protein","species"),remove = F,sep="_") 
    #mutate(acq_type=acquisiton_type) %>%
    #mutate(soft_name=software_name)
    
    write.table(all_seq, file=paste0(new_path,"/Experiment2",software_name,"_number_of_unique_sequence_for_each_species.txt"),sep = "\t",col.names = T,row.names = F)
    
    
    all_seq_syn <- all_seq %>%
      filter(!grepl(background_species,Pool_for_seq_merge)) %>%
      #select(Sequence, Modified.Sequence, Pool_for_seq_merge,PTM.Site.Confidence,Experiment,Sample_id,Rep_id) %>%
      #distinct(Sequence, .keep_all = T) %>%
      #bind_rows(ecoli_seq_dist) %>%
      #mutate(Pool_for_seq_merge= ifelse(is.na(Pool_for_seq_merge),background_species,Pool_for_seq_merge)) %>%
      mutate(acq_type=acquisiton_type) %>%
      mutate(soft_name=software_name) %>% group_by(Experiment) %>%
      distinct(Sequence,.keep_all = T)
    
    
    write.table(all_seq_syn, file = paste0(new_path,"/Number_of_human_phospho_sequences_",
                                           software_name,"_Experiment",exp_id,".txt"),
                sep = "\t",row.names = F)
    
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
                                     subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) 
    
    # all_corr_seq <- quant_peptides_cor_abun %>% 
    #   filter(grepl("Phospho",modifications)) %>%
    #   select(sequence, modifications, starts_with("E2"),accession) %>%
    #   filter(grepl(selected_species,accession) & !grepl("CON__",accession)) %>%
    #   pivot_longer(cols = starts_with("E2"), 
    #                values_to = "intensity",
    #                names_to = "Experiment",
    #                values_drop_na = T) %>%
    #   group_by(Experiment) %>%
    #   distinct(sequence,.keep_all = T) %>%
    #   mutate(species=selected_species) %>% 
    #   full_join(pep_list_w_theo_unique,by="sequence") %>%
    #   mutate_at("situation", ~replace_na(.,"Unexpected")) %>%
    #   mutate(situation= ifelse(is.na(species),"missing",situation)) %>%
    #   filter(!grepl("Unexpected",situation) & !grepl("missing",situation)) %>%
    #   select(sequence,modifications,Experiment,situation,Pool_for_seq_merge) %>%
    #   separate(Experiment, into = c("exp_id","samp_id","rep_id"),sep = "-")
    #   
    #   
    # plot15 <- gg_barplt_id_pep_count(data_set = all_corr_seq,
    #                                  x_df = all_corr_seq$samp_id,
    #                                  fill_df = all_corr_seq$rep_id,
    #                                  ymax = nrow(all_corr_seq),
    #                                  size_num=10,
    #                                  header = paste("Total number of correctly identified ",selected_species,"phospho-sequence","across each sample",sep=" "),
    #                                  caption_lab = "Mapping was done without considering phospho-positions. \n In the case of multiple PSMs, max. intensity was selected.",
    #                                  x_lab = "Sample id",
    #                                  fill_lab =  "Sample id",
    #                                  y_lab = "Number of identified peptides",
    #                                  subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) #+
    #   #theme(axis.text.x = element_text(angle = 90))
    #   
    # all_seq_dist <- all_seq %>%
    #   select(sequence, modifications, Pool_for_seq_merge) %>%
    #   distinct(sequence, .keep_all = T) %>%
    #   bind_rows(ecoli_seq_dist) %>%
    #   mutate(Pool_for_seq_merge= ifelse(is.na(Pool_for_seq_merge),background_species,Pool_for_seq_merge)) %>%
    #   mutate(acq_type=acquisiton_type) %>%
    #   mutate(soft_name=software_name)
    # 
    # plot13 <- gg_barplt_id_pep_count(data_set = all_seq_dist,
    #                        x_df = all_seq_dist$Pool_for_seq_merge,
    #                        fill_df = all_seq_dist$Pool_for_seq_merge,
    #                        ymax = nrow(all_seq_dist),
    #                        size_num=10,
    #                        header = paste("Total number of identified phosphorylated", selected_species,"and", background_species,"sequences across each sample",sep=" "),
    #                        caption_lab = "NA values are removed.",
    #                        x_lab = "Sample id",
    #                        fill_lab =  "Sample id",
    #                        y_lab = "Number of identified sequences",
    #                        subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
    
    #################################################  
    ######################################################################
    
    quant_phospho_map_pep <- quant_phospho_peptides_tmp %>% 
      mutate(Neutral_mass = round(Sum, 0)) %>%
      mutate(Neutral_mass_res=Neutral_mass) %>%
     
      pivot_longer(cols = starts_with("E2"), 
                   values_to = "Intensity",
                   names_to = "Experiment",
                   values_drop_na = T) %>%
      full_join(pep_list_w_theo_seq_map,by=c("Neutral_mass","Sequence")) %>%
      mutate(map_seq = ifelse(is.na(Neutral_mass_res), "missing", NA)) %>%
      mutate(map_seq = ifelse(is.na(Neutral_mass_theo), "wrong seq", map_seq)) %>%
      mutate(map_seq = ifelse(!is.na(Neutral_mass_res) & !is.na(Neutral_mass_theo) & Neutral_mass_res == Neutral_mass_theo, 
                              "correct seq", 
                              map_seq)) %>%
      filter(grepl("correct seq",map_seq)) %>%
      
      full_join(pep_list_w_theo_pep_map,by="pep_with_pos") %>%
      mutate(map_loc=NA) %>%
      mutate(map_loc = case_when(
        !is.na(Pool) & !is.na(species) ~ "Correct",                
        is.na(Pool) & !is.na(species) ~ "Wrong Localization",      
        !is.na(Pool) & is.na(species) ~ "Missing",               
        TRUE ~ NA_character_  # Default case if none of the above match
      )) %>%
      mutate(Pool=ifelse(Pool=="Diluted",paste0(Pool,"_",isomericity),Pool)) %>%
      mutate(Pool=ifelse(Pool=="Diluted_nonisomeric","Diluted_non-isomeric",Pool)) %>%
      # full_join(pep_list_w_theo_pep_map,by="pep_with_pos") %>%
      # mutate(map_loc=NA) %>%
      # mutate(map_loc = case_when(
      #   !is.na(Pool) & !is.na(species) ~ "Correct",                
      #   is.na(Pool) & !is.na(species) ~ "Wrong Localization",      
      #   !is.na(Pool) & is.na(species) ~ "Missing",               
      #   TRUE ~ NA_character_  # Default case if none of the above match
      # )) %>%
      select(!c(Positions,Sum)) 
    
   
      
      # #pep_with_pos
     
    quant_phospho_aft_inner_map_pep <- quant_phospho_map_pep %>% 
      group_by(pep_with_pos,ptm_score) %>%
      summarize(
        sum_int = sum(Intensity), 
        .groups = 'drop'
      ) %>%
      # Join the summarized Intensity with the original dataframe to get the row with the highest ptm_score
      inner_join(quant_phospho_map_pep, by =c("pep_with_pos","ptm_score")) %>%
      # Select the row with the highest ptm_score
      group_by(pep_with_pos) %>%
      #slice_max(sum_int, n = 1) %>%
      slice_max(ptm_score, n = 1) %>%
      ungroup() %>%
      rename(Pool_for_pep_merge=map_loc) %>%
      group_by(pep_with_pos,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
      slice(which.max(Intensity)) %>%
      ungroup() 
    
     
    
    # df_merge_all_col_wide <-df_merge_all_col %>%
    #   mutate(Pool=ifelse(is.na(Pool),"Wrong localization",Pool)) %>%
    #   pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
    #   rowwise() %>%
    #   mutate(row_sum = sum(c_across(all_of(exp_design)), na.rm = TRUE)) %>%
    #   group_by(pep_with_pos) %>%
    #   filter(row_sum == max(row_sum)) %>%
    #   ungroup() %>%
    #   mutate(soft_name=software_name,ion_mobility=acquisiton_type)
   
    # quant_phospho_aft_inner_map_pep <- quant_phospho_map_pep %>% 
    #   group_by(pep_with_pos,ptm_score) %>%
    #   summarize(
    #     sum_int = sum(Intensity), 
    #     .groups = 'drop'
    #   ) %>%
    #   # Join the summarized Intensity with the original dataframe to get the row with the highest ptm_score
    #   inner_join(quant_phospho_map_pep, by =c("pep_with_pos","ptm_score")) %>%
    #   # Select the row with the highest ptm_score
    #   group_by(pep_with_pos) %>%
    #   #slice_max(sum_int, n = 1) %>%
    #   slice_max(ptm_score, n = 1) %>%
    #   ungroup() %>%
    #   rename(Pool_for_pep_merge=map_loc) %>%
    #   group_by(pep_with_pos,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    #   slice(which.max(Intensity)) %>%
    #   ungroup() 
    
    ########## ########## ########## ########## ########## ########## ########## ##########
    ########## ########## ########## ########## ########## ########## ########## ##########  
 
    all_pep_syn <- quant_phospho_aft_inner_map_pep %>%
    
      #filter(!grepl(background_species,Pool_for_seq_merge)) %>%
      #select(pep_with_pos, Modified.Sequence, map_loc,PTM.Site.Confidence,Experiment,Sample_id,Rep_id) %>%
      
      #distinct(Sequence, .keep_all = T) %>%
      #bind_rows(ecoli_seq_dist) %>%
      #mutate(Pool_for_seq_merge= ifelse(is.na(Pool_for_seq_merge),background_species,Pool_for_seq_merge)) %>%
      mutate(acq_type=acquisiton_type) %>%
      mutate(soft_name=software_name) %>% 
      group_by(Experiment) %>%
      distinct(pep_with_pos,.keep_all = T)
    
    all_corr_pep_syn <- all_pep_syn %>% 
      filter(grepl("Correct",Pool_for_pep_merge)) %>%
      separate(Experiment,into = c("Exp_id","Sample_id","Rep_id"),sep = "-",remove = F)
    
    plot12 <- gg_barplt_id_pep_count(data_set = all_corr_pep_syn,
                                     x_df = all_corr_pep_syn$Sample_id,
                                     fill_df = all_corr_pep_syn$Rep_id,
                                     ymax = 20000,
                                     size_num = 10,
                                     header = paste("Total number of correctly identified & localized phosphorylated", selected_species,"across each sample",sep=" "),
                                     #,"and", background_species,"sequence 
                                     caption_lab = "NA values are removed. (p13)",
                                     x_lab = "Sample id",
                                     fill_lab =  "Sample id",
                                     y_lab = "Number of identified sequence",
                                     subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
    
    ##############################################################################
    
    df_merge_all_col <- quant_phospho_aft_inner_map_pep %>%
      select(pep_with_pos,Experiment,Intensity,ptm_score,species) %>% #Marked.as,species
      #rename(Pool_leftjoin=Pool) %>%
      full_join(pep_list_w_theo_sel,by="pep_with_pos") %>%
      # mutate(map_loc = case_when(
      #   !is.na(Pool) & !is.na(species) ~ "Correct",                
      #   is.na(Pool) & !is.na(species) ~ "Wrong Localization",      
      #   !is.na(Pool) & is.na(species) ~ "Missing",               
      #   TRUE ~ NA_character_  # Default case if none of the above match
      # )) %>%
      mutate_at("Pool", ~replace_na(.,"Wrong Localization")) %>%
      mutate(Pool= ifelse(is.na(species),"Missing",Pool)) %>%
      mutate(Pool=ifelse(Pool=="Diluted",paste0(Pool,"_",isomericity),Pool)) %>%
      mutate(Pool=ifelse(Pool=="Diluted_nonisomeric","Diluted_non-isomeric",Pool))
    
    df_merge_all_col_wide <-df_merge_all_col %>%
      filter(!grepl("Missing",Pool)) %>%
      pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
      rowwise() %>%
      mutate(row_sum = sum(c_across(all_of(exp_design)), na.rm = TRUE)) %>%
      group_by(pep_with_pos) %>%
      filter(row_sum == max(row_sum)) %>%
      ungroup() %>%
      mutate(soft_name=software_name,ion_mobility=acquisiton_type) 
    
    write.table(df_merge_all_col_wide,file = paste0(new_path,"/Count_of_missing_unexpected_correct_phospho-sites_with_all_col_",
                                                    software_name,"_Experiment",exp_id,".txt"),
                sep = "\t",row.names = F)
    #############################################################################
    
    if(loc_filter_opt == TRUE){
      quant_phospho_peptides <- df_merge_all_col_wide %>%
        # filter(grepl(selected_species,accession) & 
        #          grepl("Phospho",modifications)) %>%
        # select(sequence,
        #        modifications,
        #        accession,
        #        ptm_score,
        #        ptm_sites_confidence,
        #        ptm_protein_positions,
        #        charge,
        #        spectrum_title,
        #        starts_with(exp_design)) %>%
      # rowwise() %>%
      # mutate(species=selected_species) %>%
      filter(as.numeric(ptm_score) >= loc_filter)  #%>%  
       
      #mutate(phospho_pos = proline_phospho_pos_extraction(modifications)) %>%
      #mutate(pep_with_pos=paste0(sequence,"_",phospho_pos))
      
    }else{
      quant_phospho_peptides <- df_merge_all_col_wide
      # filter(grepl(selected_species,accession) & 
      #          grepl("Phospho",modifications)) %>%
      # select(sequence,
      #        modifications,
      #        accession,
      #        ptm_score,
      #        ptm_sites_confidence,
      #        ptm_protein_positions,
      #        charge,
      #        spectrum_title,
      #        starts_with(exp_design)) %>%
      # rowwise() %>%
      # mutate(species=selected_species) %>%
      # #filter(as.numeric(ptm_score) >= loc_filter)  %>%  
      # mutate(phospho_pos = proline_phospho_pos_extraction(modifications)) %>%
      # mutate(pep_with_pos=paste0(sequence,"_",phospho_pos))
      
    }
    
    #################################################
    #############################################################################
    # quant_phospho_peptides <- all_seq %>% 
    #     filter(grepl(selected_species,accession)) %>% 
    #     filter(grepl("Phospho",modifications)) %>%
    #     select(sequence,
    #            modifications,
    #            accession,
    #            ptm_protein_positions,
    #            charge,
    #            spectrum_title,
    #            starts_with(exp_design)
    
    # tmp <- all_seq %>%
    #   filter(!grepl("CON__",accession)) %>%
    #   select(sequence,
    #          modifications,
    #          accession,
    #          ptm_score,
    #          ptm_sites_confidence,
    #          ptm_protein_positions,
    #          spectrum_title,
    #          species,
    #          starts_with(exp_design)) %>%
    #   pivot_longer(cols = starts_with("E2"), 
    #                values_to = "intensity",
    #                names_to = "sample_ids",
    #                values_drop_na = T) 
    # 
    # 
    # pool_wise <- tmp %>% filter(grepl("Phospho",modifications)) %>%
    #   rowwise() %>%
    #   mutate(species=selected_species) %>%
    #   #filter(as.numeric(ptm_score) >= loc_filter)  %>%  
    #   mutate(phospho_pos = proline_phospho_pos_extraction(modifications)) %>%
    #   mutate(pep_with_pos=paste0(sequence,"_",phospho_pos)) 

    
    ######### THIS PART WAS MOVED INSIDE THE CODE WAS DONE PTM FILTERING (ABOVE) ######### 
    
    # Extraction of phospho positions from quant peptides object
    #phospho_ptm_pos <- proline_phospho_pos_extraction(quant_phospho_peptides$modifications)
    # Data conversion 
    #phospho_ptm_pos_df <- t(as.data.frame(phospho_ptm_pos))
    #rownames(phospho_ptm_pos_df) <- 1:length(phospho_ptm_pos_df)
    
    # Creation of common column merging peptide sequence and phospho positions -> experimental data
    #common_col_exp_quant <- as.data.frame(paste(quant_phospho_peptides$sequence, phospho_ptm_pos, sep = "_"))
    #colnames(common_col_exp_quant) <- "pep_with_pos"
    
    # Bind it to the quant data
    #syn_phospho_pep_proline_new <- cbind(common_col_exp_quant,quant_phospho_peptides)
    
    # Creation of common column merging peptide sequence and phospho positions -> theoretical data
    #common_col_theo_quant <- as.data.frame(paste(pep_list_w_theo_quant$Phosphopeptide.sequence,
                                                 #pep_list_w_theo_quant$modified.position.in.peptide, sep = "_"))
    #colnames(common_col_theo_quant) <- "pep_with_pos"
    #pep_list_w_theo_quant_new <- cbind(common_col_theo_quant,pep_list_w_theo_quant)
    
    ###########################################################################
    # Nothing is changed
    ### INPUT PARAM: num_allowed_na and num_expected_na are linked to each other!!
    
      # num_allowed_NA: total number of NA values overall dataset 
         ## It enables to control number of the number of allowed NAs !!
         ## (if it is set as number of sample, it will allow only what was defined in num_expected_nonNA.
         ## However, it can be modified by decreasing the number till zero to allow more NA values 
          ## ATTENTION: this time it will not be considered by sample
         ## ATTENTION: it cannot be higher than total number of sample
            
    
      # num_expected_nonNA: count of NAs per sample 
         ## (if it is set as 1 means that it will search at least 1 non-NA values)
         ## (if it is set as 2 means that it will search at least 2 non-NA values)
    filtered_abundances <-  filter_NA(df = quant_phospho_peptides,samp_names = sample_names,num_allowed_NA = 1,num_expected_nonNA = 3)
    filtered_abundances_ecoli <-  filter_NA(df = quant_peptides_ECOLI,samp_names = sample_names,num_allowed_NA = 5,num_expected_nonNA = 1)
    
    filtered_abundances_ecoli_bfr_impt <- filtered_abundances_ecoli
    #filtered_abundances<-quant_phospho_peptides[rowSums(!is.na(select(quant_phospho_peptides,starts_with(exp_design))))>3,] #0
    #filtered_abundances_ecoli <-quant_peptides_ECOLI[rowSums(!is.na(select(quant_peptides_ECOLI,starts_with(exp_design))))>3,] #0
    # 
    # df_id_pep <- filtered_abundances %>% 
    #     select(sequence,modifications,charge, pep_with_pos, starts_with(exp_design),accession) %>%
    #     tibble() %>%
    # pivot_longer(cols = starts_with("E2"), 
    #              values_to = "intensity",
    #              names_to = "sample_ids",
    #              values_drop_na = T) %>%
    #     separate(sample_ids, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F)
    # 
    ### ELIMINATION of MULTIPLE CHARGED PEPTIDES ### 
    
    #### 1st STRATEGY: without pivot_longer() func. 
    
    ### IT's an alternative way to remove duplicates from mouse dataset
    ### DOES NOT WORK BACKGROUND DF BECAUSE POSITION IS IMPORTANT FOR THE DISCRIMINATION
    
    # syn_phospho_pep_proline_unique_max_absum <- filtered_abundances %>% 
    #  #mutate(across(all_of(experiment_name), ~ ifelse(is.na(.), 0, .))) %>%
    #  rowwise() %>%
    #  mutate(row_sum = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
    #  group_by(pep_with_pos) %>%
    #  #summarise(max_row_sum = max(row_sum, na.rm = TRUE)) %>%
    #  slice(which.max(row_sum)) %>%
    #  ungroup()
    
    ### THE RESULT OF barplt_df is the same as syn_phospho_pep_proline_unique_max_absum
    
    ## 2nd STRATEGY: removing duplicates automatically
    
    # barplt_df <- df_id_pep %>%
    #   mutate(sample_rep_id_seq = paste(pep_with_pos, Sample_id,Rep_id, sep = "_"))%>%
    #   filter(duplicated(sample_rep_id_seq)==FALSE)
    
    ## 3rd STRATEGY: removing duplicates the one has lower abundances
    # barplt_df <- df_id_pep %>%
    #     #mutate(sample_rep_id_seq = paste(common_col_for_merging, Sample_id,Rep_id, sep = "@"))%>%
    #     group_by(pep_with_pos,sample_ids) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    #     slice(which.max(intensity)) %>%
    #     ungroup() 
    #   
    
    ####### ADDITIONAL PLOT TO DISPLAY MISSING and UNEXPECTED PEPTIDES ########
    df_merge_syn <- filtered_abundances %>%
      #select(pep_with_pos,sample_ids,Inten, species) %>% #Marked.as
      #pivot_wider(names_from = "sample_ids",values_from = "intensity") %>%
      #full_join(pep_list_w_theo,by="pep_with_pos") %>% 
      #mutate_at("Pool", ~replace_na(.,"Unexpected")) %>%
      #mutate(Pool= ifelse(is.na(species),"missing",Pool)) %>%
      #mutate(Pool=ifelse(Pool=="Diluted",paste0(Pool,'_',isomericity),Pool)) %>%
      select(pep_with_pos,starts_with(exp_design),Pool) %>%
      mutate(soft_name=software_name,ion_mobility=acquisiton_type)
    
    ################################################################################ 
    
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
   
    write.table(df_merge_syn,file = paste0(new_path,"/Count_of_missing_unexpected_correct_phospho-sites_",
                                           software_name,"_Experiment",exp_id,".txt"),
                sep = "\t",row.names = F)
    #############################################################################
    
    ##############################################################################
    ## 
    tmp_wide <- barplt_df %>% 
      select(pep_with_pos, species, Experiment,Intensity) %>% #Marked.as
      pivot_wider(names_from = "Experiment",values_from = "Intensity")
    
    barplt_df_wide <- filtered_abundances %>%
      select(pep_with_pos, species, starts_with("E2"))#Marked.as
      #select(colnames(tmp_wide))
    
    barplt_df_ecoli_wide <- barplt_df_ecoli %>% 
    
      select(Sequence, species, Experiment,Intensity) %>% #Marked.as
      pivot_wider(names_from = "Experiment",values_from = "Intensity")
    
    
    #merge_phospho_pos <- proline_phospho_pos_extraction(merge_phospho_peptides[,"ptm_protein_positions"])
    
    #pep_with_pos_merge <- cbind(merge_phospho_peptides, paste(merge_phospho_peptides$sequence,merge_phospho_pos,sep = "_"))
    #colnames(pep_with_pos_merge)[dim(pep_with_pos_merge)[2]] <- "pep_with_pos"
    # 
    # df_merge_all_col <- quant_peptides_cor_abun %>% 
    #   filter(grepl(selected_species,accession)) %>% 
    #   filter(grepl("Phospho",modifications)) %>%
    #   tibble() %>%
    #   mutate(phospho_pos =proline_phospho_pos_extraction(modifications)) %>%
    #   mutate(pep_with_pos=paste0(sequence,"_",phospho_pos)) %>%
    #   distinct(pep_with_pos,.keep_all = T) %>%
    #   full_join(pep_list_w_theo,by="pep_with_pos") %>%
    #   mutate_at("Pool", ~replace_na(.,"Unexpected")) %>%
    #   mutate(Pool= ifelse(is.na(accession),"missing",Pool)) %>%
    #   mutate(Pool=ifelse(Pool=="Diluted",paste0(Pool,'_',isomericity),Pool)) %>%
    #   select(pep_with_pos,starts_with(exp_design),Pool,ptm_score) %>%
    #   mutate(soft_name=software_name,ion_mobility=acquisiton_type) 
    # # 
    # barplt_phospho_seq <- filtered_abundances %>%
    #   select(sequence,modifications, sequence, starts_with(exp_design),accession) %>%
    #   pivot_longer(cols = starts_with("E2"),
    #                values_to = "intensity",
    #                names_to = "sample_ids",
    #                values_drop_na = T) %>%
    #   separate(sample_ids, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F) %>%
    #   mutate(sample_rep_id_seq = paste(sequence, Sample_id,Rep_id, sep = "@")) %>%
    #   group_by(sample_rep_id_seq,sample_ids) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    #   slice(which.max(intensity)) %>%
    #   ungroup() %>% mutate(Software_name=software_name) %>%
    #   mutate(Acquisition_type=acquisiton_type)
    # 
    # write.table(barplt_phospho_seq, file = paste0(new_path,"/Number_of_human_phospho_sequences_",
    #                                               software_name,"_Experiment",exp_id,".txt"),
    #             sep = "\t",row.names = F)

    #colnames(abundances_for_impute) <- experiment_name
    
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
    
    
    
    # plot1 <- gg_barplt_id_pep_count(data_set = barplt_df,
    #                              x_df = barplt_df$Sample_id,
    #                              fill_df = barplt_df$Rep_id,
    #                              ymax = 20000,
    #                              size_num = 10,
    #                              header = "Total number of quantified phospho-site across each sample",
    #                              caption_lab = "NA values are removed.",
    #                              x_lab = "Sample id",
    #                              fill_lab =  "Sample id",
    #                              y_lab = "Number of identified peptides",
    #                              subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
    # 
    
   
    # plot12 <- gg_barplt_id_pep_count(data_set = barplt_phospho_seq,
    #                               x_df = barplt_phospho_seq$Sample_id,
    #                               fill_df = barplt_phospho_seq$Rep_id,
    #                               ymax = 20000,
    #                               size_num = 10,
    #                               header = "Total number of quantified phospho-sequence across each sample",
    #                               caption_lab = "NA values and multiple sequences are removed.",
    #                               x_lab = "Sample id",
    #                               fill_lab =  "Sample id",
    #                               y_lab = "Number of identified peptides",
    #                               subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
    # 
    
    ## THIS RESHAPING IS ONLY FOR ELIMINATION OF MULTIPLE PHOSPHO-SITES and ECOLI PEPTIDES
    ## ELIMINATION STEP IS NOT NECESSARY FOR ECOLI, 1st STRATEGY can be used only (this will decrease lines of code)
    # 
    # barplt_df_wide <- barplt_df %>%  ## If you select "charge" column, it will bring multiple rows for one seq
    #   select(sequence,pep_with_pos,accession, sample_ids,intensity) %>%
    #   pivot_wider(names_from = "sample_ids",values_from = "intensity") %>%
    #   separate(accession,into = c("proteins","species","positions"),sep = "_") %>%
    #   select(!c(sequence,positions,proteins))
    #   
    # barplt_df_ecoli_wide <- barplt_df_ecoli %>% 
    #   select(sequence,sample_ids, accession,intensity) %>%
    #   pivot_wider(names_from = "sample_ids",values_from = "intensity") %>%
    #   separate(accession,into = c("proteins","species"),sep = "_") %>%
    #   select(!c(proteins))
    # 

    barplt_df_wide[,"na_val"] <- apply(X = !is.na(select(barplt_df_wide,contains(exp_design))), MARGIN = 1, FUN = sum)
    
    phospho_completeness <- barplt_df_wide %>% 
      count(na_val) %>%
      mutate(data_complete=((n/dim(barplt_df_wide)[1])*100)) %>% mutate(species=selected_species)
    
    barplt_df_ecoli_wide[,"na_val"] <- apply(X = !is.na(select(barplt_df_ecoli_wide,contains(exp_design))), MARGIN = 1, FUN = sum)
    
    completeness <- barplt_df_ecoli_wide %>% 
      count(na_val) %>%
      mutate(data_complete=((n/dim(barplt_df_ecoli_wide)[1])*100)) %>%
      mutate(species=background_species) %>%
      bind_rows(phospho_completeness)
    
    write.table(completeness,file = paste0(new_path,"/data_completeness",acquisiton_type,software_name,".txt"))
      
    plot23 <- ggplot(completeness, aes(x=na_val,y=data_complete,color=species)) + geom_point(size=2.5) +
       geom_line(size=2)+
       scale_x_reverse(limits=c(18,0),breaks=seq(0, 18, by = 2)) +
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
       ) + scale_color_manual(values = c("#3182bd","#a6bddb"))+
       ggtitle(label = paste("Data completeness of",background_species,"and",selected_species)) +
       labs(x="n Sample", y="% of peptides",subtitle = paste("Experiment",exp_id,software_name,acquisiton_type,"\n",file_name))
  
    
    #na_phospho_selected <- apply(X = is.na(barplt_df_wide %>% select(sequence, accession,starts_with("E2"))), MARGIN = 2, FUN = sum)
    #na_ecoli <- apply(X = is.na(barplt_df_ecoli_wide %>% select(sequence,accession,starts_with("E2"))), MARGIN = 2, FUN = sum)
    
    #na_table <- bind_rows(na_phospho_selected,na_ecoli)
    #na_table$species <- c(selected_species,background_species)
    #na_table$total <- c(dim(barplt_df_wide)[1],dim(barplt_df_ecoli_wide)[1])
    
    #na_table %>% select(-sequence) %>%
        #kbl(caption = paste("Number of NA values across all samples \n Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) %>%
        #kable_material(c("striped", "hover")) %>%
        #kable_styling(bootstrap_options = "striped", full_width = F, position = "left", font_size = 12) %>%
        #kable_minimal(full_width = F) %>%
        #footnote(general =  paste("This table was created after elimination of multiple charges by selecting either phospho-sites of", selected_species, "and", background_species, "sequences \n that has the highest abundace."),
                 # number = c("Footnote 1; ", "Footnote 2; "),
                 # alphabet = c("Footnote A; ", "Footnote B; "),
                 # symbol = c("Footnote Symbol 1; ", "Footnote Symbol 2")
                 #footnote_as_chunk = T, title_format = c("italic", "underline")) %>%
        #as_image(width = 8) %>%
        #save_kable(paste0(file_path,"/outputs/table1.png"))
    
    ## REMOVE SEQUENCE COLUMN AFTER NA TABLE
    #barplt_df_wide <- barplt_df_wide %>% select(!sequence)
    
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
                     values_to = "Intensity",
                     names_to = "Experiment",
                     values_drop_na = T) %>%
        mutate(sample_id_seq = paste(Sequence, Experiment, sep = "_"))
    
    ecoli_density_plt<- quant_peptides_ECOLI_density_plt %>%
        select(contains(c("sample_ids","intensity","species"))) 
    
    
    colnames(abundances_rowMeans) <- paste0("mean_abun",1:sample_size)
    
    #### MEAN ABUNDANCE RATIO WITH  DENSITY PLOT ####
    ### BEFORE IMPUTATION ###
    quant_phospho_density_plt <- barplt_df_wide %>%
        select(!starts_with("E")) %>%
        bind_cols(abundances_rowMeans) %>% 
        #rename_with(~ paste0("mean_abun",1:5), matches("^row")) %>%
        tibble() %>% #mutate(pep_with_pos = sequence) %>% ###  At this stage, no need for phospho-position#   
        pivot_longer(cols = starts_with("mean"),
                     names_to = "sample_ids",
                     values_to = "intensity",
                     values_drop_na = T) %>%
        mutate(sample_id_seq = paste(pep_with_pos, sample_ids, sep = "_"))
    
    density_df <-quant_phospho_density_plt %>%
        select(c(sample_ids,intensity,species)) %>%
        bind_rows(ecoli_density_plt) #%>%
        #separate(accession, into = c("prot_id","species","position"),sep = "_",remove = F)
    
    
    plot3 <- gg_density(data_set = density_df, 
                     x_df = density_df$intensity,
                     fill_df = density_df$species,
                     color_df = NULL,
                     header="Distribution of mean abundance of every sample before imputation",
                     facet_df = "sample_ids",
                     x_lab = "log10(intensities)",
                     color_lab= "",
                     fill_lab = "species",
                     subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
    
    ### RATIO SUPPRESION ASSESSMENT ### 
    ### For this calculation, wrong localizations were removed. 
    
    pep_list_w_theo_sel <- pep_list_w_theo %>%
      select(pep_with_pos,isomericity,Pool,pool_id) 
    
    ratio_supp_fixed <- barplt_df_wide %>%
      select(!starts_with("E")) %>%
      bind_cols(abundances_rowMeans) %>% 
      left_join(pep_list_w_theo_sel,by="pep_with_pos") %>% 
      filter(grepl("Fixed",Pool)) %>%
      rename_with(~ paste0("mean_abun_",sample_names), matches("^mean")) %>%
      #rename_with(~ paste0("mean_abun",1:5), matches("^row")) %>%
      tibble() %>% #mutate(pep_with_pos = sequence) %>% ###  At this stage, no need for phospho-position#   
      pivot_longer(cols = starts_with("mean"),
                   names_to = "sample_ids",
                   values_to = "Intensity",
                   values_drop_na = T) %>%
      drop_na(Pool)
    
    gg_raincloud(data_set = ratio_supp_fixed,
                 x_df = ratio_supp_fixed$sample_ids,
                 y_df = ratio_supp_fixed$Intensity,
                 fill_df = ratio_supp_fixed$Pool,
                 header = "Distribution of mean abundance of every sample after imputation",
                 x_lab = "Sample Names",
                 y_lab = " Density of log10(Mean Abundance)",
                 fill_lab = "Sample Names",
                 caption_lab = "",
                 subtitle_txt = "Fixed Pool")
    
    ratio_supp_spiked <- barplt_df_wide %>%
      select(!starts_with("E")) %>%
      bind_cols(abundances_rowMeans) %>% 
      left_join(pep_list_w_theo_sel,by="pep_with_pos") %>% 
      filter(!grepl("Fixed",Pool)) %>%
      drop_na(Pool) %>%
      rename_with(~ paste0("mean_abun_",sample_names), matches("^mean")) %>%
      mutate(Pool=paste(Pool,isomericity,sep = "_")) %>%
      #rename_with(~ paste0("mean_abun",1:5), matches("^row")) %>%
      tibble() %>% #mutate(pep_with_pos = sequence) %>% ###  At this stage, no need for phospho-position#   
      pivot_longer(cols = starts_with("mean"),
                   names_to = "sample_ids",
                   values_to = "Intensity",
                   values_drop_na = T) 
    
    gg_raincloud(data_set = ratio_supp_spiked,
                 x_df = ratio_supp_spiked$sample_ids,
                 y_df = ratio_supp_spiked$Intensity,
                 fill_df = ratio_supp_spiked$Pool,
                 header = "Distribution of mean abundance of every sample after imputation",
                 x_lab = "Sample Names",
                 y_lab = " Density of log10(Mean Abundance)",
                 fill_lab = "Sample Names",
                 caption_lab = "",
                 subtitle_txt = "Spiked Pool")
    

    # abundances_for_before_impt <- barplt_df_wide %>%
    #     bind_rows(barplt_df_ecoli_wide) %>% select(exp_design)
    # 
    # Calculate 5 percent quantile of each sample
    impute_values <- apply(abundances_for_impute, 2 , quantile , probs = 0.05 , na.rm = TRUE )
    
    ###########################################################################
    ## ADDITIONAL IMPUTATION METHOD with MICE()
    #library(tidyverse)
    #library(tidyr)
    #library(mice)
    
    
    is.imputed_df_syn <- barplt_df_wide %>%  pivot_longer(cols = starts_with(exp_design),
                                                          names_to = "Experiment",
                                                          values_to = "Intensity",
                                                          values_drop_na = F) %>%
      mutate(is.imputed=FALSE) %>%
      mutate(is.imputed=ifelse(is.na(Intensity), TRUE,is.imputed)) 
    
    
    #barplt_df_wide <- as.data.frame(barplt_df_wide)
    #rownames(barplt_df_wide) <- paste0(barplt_df_wide$pep_with_pos,"@",barplt_df_wide$species,"@",1:nrow(barplt_df_wide))
    
    #intensities <- barplt_df_wide %>%
    # select(starts_with(exp_design))
    
    #intensities_short <- intensities[1:10,]
    # barplt_df_ecoli_wide <- as.data.frame(barplt_df_ecoli_wide)
    #rownames(barplt_df_ecoli_wide) <- paste0(barplt_df_ecoli_wide$sequence,"@",barplt_df_ecoli_wide$species,"@",1:nrow(barplt_df_ecoli_wide))
    
    is.imputed_df_ecoli <- barplt_df_ecoli_wide %>%  pivot_longer(cols = starts_with(exp_design),
                                                                  names_to = "Experiment",
                                                                  values_to = "Intensity",
                                                                  values_drop_na = F) %>%
      mutate(is.imputed=FALSE) %>%
      mutate(is.imputed=ifelse(is.na(Intensity), TRUE,is.imputed)) 
    
    #intensitiesECOLI <- barplt_df_ecoli_wide %>%
    # select(starts_with(exp_design))
    
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
    # colnames(intensities) <- new_experiment_name
    # colnames(intensitiesECOLI) <- new_experiment_name
    
    #imp_intensities <-  mice(intensities,m=5,maxit=10,meth='cart',seed=500)
    #imp_intensities_norm <-  mice(intensities,m=5,maxit=10,method ='norm.nob',seed=500)
    
    #imp_intensitiesECOLI <- mice(intensitiesECOLI,m=5,maxit=5,meth='cart',seed=500) #maxit=10
    
    #completeData <- complete(imp_intensities,2)
    #colnames(completeData) <- exp_design
    
    #completeDataECOLI <- complete(imp_intensitiesECOLI,2)
    #colnames(completeDataECOLI) <- exp_design
    #pattern <- md.pattern(select(quant_phospho_peptides,starts_with(exp_design)))
    
    #library(VIM)
    #aggr_plot <- aggr(intensities, col=c('navyblue','red'),
    #numbers=TRUE, sortVars=TRUE,
    #labels=names(intensities), cex.axis=.7,
    #gap=3, ylab=c("Histogram of missing data","Pattern"))
    
    #completeDataECOLIs <-completeDataECOLI %>%
    #mutate(tmp=rownames(completeDataECOLI)) %>%
    #separate(tmp,into=c("pep_with_pos","species","indx"),sep="@") %>%
    #select(!indx)
    ##########################################################################
    
    # # Impute missing values
    for (j in 1:length(impute_values)){
      # # Number NA
      # #num_NA <- length(abundances_for_impute_all[,j+2][is.na(abundances_for_impute_all[,j+2])])
      # 
      barplt_df_wide[,j+2][is.na(barplt_df_wide[,j+2])] <- impute_values[j]
      barplt_df_ecoli_wide[,j+2][is.na(barplt_df_ecoli_wide[,j+2])] <- impute_values[j]
      # 
      #abundances_for_impute_all[,j+2][is.na(abundances_for_impute_all)[,j+2]] <- impute_values[j]
      # After imputation number of imputed values
      #num_imp <-length(abundances_for_impute_all[,j+2][(abundances_for_impute_all[,j+2]==impute_values[j])])
      
      # This is verification of imputation is done successfully
      # Because we expect to see that number of imputed values should be the same amount as number of NA
      #print(setequal(num_NA,num_imp))
      #print(num_NA)
      #print(num_imp)
    }
    
    #barplt_df_wide <- as.data.frame(barplt_df_wide)
    #rownames(barplt_df_wide) <- paste0(barplt_df_wide$pep_with_pos,"@",barplt_df_wide$species,"@",1:nrow(barplt_df_wide))
    
    #intensities <- barplt_df_wide %>%
    # select(starts_with(exp_design))
    
    #intensities_short <- intensities[1:10,]
    #barplt_df_ecoli_wide <- as.data.frame(barplt_df_ecoli_wide)
    #rownames(barplt_df_ecoli_wide) <- paste0(barplt_df_ecoli_wide$sequence,"@",barplt_df_ecoli_wide$species,"@",1:nrow(barplt_df_ecoli_wide))
    
    abundances_all_aft_imputation <- barplt_df_ecoli_wide %>%
      rename_with(~ paste0("pep_with_pos"), matches("^seq")) %>%
      bind_rows(barplt_df_wide) 
    
    # conditional_imputation <- function(df,impute_val,samp_names,num_NA_imputed){
    #   all_complete_df <- NULL
    #   for (i in 1:length(samp_names)){
    #     assign(paste0("df_",samp_names[i]), df %>% 
    #              mutate(impute=ifelse(rowSums(is.na(select(df,contains(samp_names[i]))))>=num_NA_imputed,TRUE,FALSE)) %>% 
    #              select(contains(samp_names[i]) | contains("impute")))
    #     
    #     sel_for_impt <- get(paste0("df_",samp_names[i])) %>% 
    #       filter(grepl(TRUE,impute)) %>% 
    #       select(!impute)
    #     
    #     nonimp_df <- get(paste0("df_",samp_names[i])) %>% 
    #       filter(!grepl(TRUE,impute)) %>% 
    #       select(!impute)
    #     
    #     #for(j in 0:5){
    #     sel_impute_values <- impute_val[(3*(i-1)+1):(3*(i-1)+3)]
    #     
    #     for(k in 1:length(sel_impute_values)){
    #       sel_for_impt[,k]<-sel_impute_values[k]
    #     }
    #     
    #     assign(paste0("df_complete",samp_names[i]),bind_rows(nonimp_df, sel_for_impt))
    #     all_complete_df <- bind_cols(all_complete_df, get(paste0("df_complete",samp_names[i])))
    #     #}
    #     
    #     rm(sel_for_impt,nonimp_df)
    #   }
    #   return(all_complete_df)
    # }
    # 
    # barplt_df_wide <- conditional_imputation(df=barplt_df_wide,samp_names = sample_names,num_NA_imputed = 2,impute_val = impute_values)
    # barplt_df_ecoli_wide <- conditional_imputation(df=barplt_df_ecoli_wide,samp_names = sample_names,num_NA_imputed = 2,impute_val = impute_values)
    # 
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
    # colnames(barplt_df_wide) <- new_experiment_name
    # colnames(barplt_df_ecoli_wide) <- new_experiment_name
    # 
    # library(mice)
    # imp_intensities <-  mice(barplt_df_wide,m=5,maxit=10,meth='cart',seed=500)
    # imp_intensitiesECOLI <- mice(barplt_df_ecoli_wide,m=5,maxit=1,meth='cart',seed=500) #maxit=10
    # 
    # completeData <- complete(imp_intensities,2)
    # colnames(completeData) <- exp_design
    # 
    # completeDataECOLI <- complete(imp_intensitiesECOLI,2)
    # colnames(completeDataECOLI) <- exp_design
    # 
    # 
    # barplt_df_wide_rowname <- barplt_df_wide %>%
    #   mutate(tmp=rownames(barplt_df_wide))
    #   
    # abundances_all_aft_imputation <- barplt_df_ecoli_wide %>%
    #   mutate(tmp=rownames(barplt_df_ecoli_wide)) %>%
    #   bind_rows(barplt_df_wide_rowname) %>%
    #   separate(tmp,into = c("pep_with_pos","species","indx"),sep = "@") %>%
    #   select(!indx)
    # #   
    # completeDataECOLIs <-completeDataECOLI %>%
    # mutate(tmp=rownames(completeDataECOLI)) %>%
    # separate(tmp,into=c("pep_with_pos","species","indx"),sep="@") %>%
    # select(!indx)
    #  
    # abundances_all_aft_imputation <- completeData %>%
    #   mutate(tmp=rownames(completeData)) %>%
    #   separate(tmp, into = c("pep_with_pos","species","indx"),sep = "@") %>%
    #   select(!indx) %>%
    #   #separate(pep_with_pos, into = c("pep","pos","species"),sep = "_") %>%
    #   #mutate(pep_with_pos=paste0(pep,"_",pos)) %>%
    #   #select(!c(pep,pos)) %>%
    #   bind_rows(completeDataECOLIs)
    ############################################################################
    
    ### DISTRIBUTION OF IMPUTED VALUES ACROSS non-NA values
    is.imputed_df <- is.imputed_df_ecoli %>% 
      rename(pep_with_pos=Sequence) %>% 
      bind_rows(is.imputed_df_syn)
    
    imputed_dataset <- abundances_all_aft_imputation %>% 
      pivot_longer(cols = starts_with(exp_design),
                   names_to = "Experiment",
                   values_to = "Intensity",
                   values_drop_na = F) %>%
      left_join(is.imputed_df,by=c("pep_with_pos","Experiment")) %>%
      mutate(Intensity.y=ifelse(is.na(Intensity.y),0,Intensity.y))
    
    write.table(imputed_dataset,file=paste0(new_path,"/imputed_dataset",software_name,".txt"),sep = "\t",row.names = F)
    
    p21  <- imputed_dataset %>% filter(grepl(selected_species,species.x)) %>% 
      ggplot(aes(x=log2(Intensity.x),fill=is.imputed)) + geom_histogram(bins = 30) +
      theme_minimal() + scale_fill_brewer(palette = "Set1",direction = -1) +
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
      theme_minimal() + scale_fill_brewer(palette = "Set1",direction = -1) +
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
    
    final_imputed_data_syn <- final_imputed_data %>% filter(grepl(selected_species,species)) #accession
    
    final_imputed_data_ecoli <- final_imputed_data  %>% filter(!grepl(selected_species, species))
    
    df_merge_syn_after_impt <- final_imputed_data_syn %>%
      left_join(pep_list_w_theo,by="pep_with_pos") %>% 
      mutate_at("Pool", ~replace_na(.,"Unexpected")) %>%
      mutate(Pool=ifelse(Pool=="Diluted",paste0(Pool,"_",isomericity),Pool))
    
    df_merge <- df_merge_syn_after_impt %>%
        bind_rows(final_imputed_data_ecoli) %>%
        mutate_at("Pool", ~replace_na(.,background_species)) 
    
    #write.table(final_imputed_data, file = "final_imputed_normalized_data_PAL _T_cell_Exp3_( 5 conc 3reps)_NoFAIMS_DDA_with_cont_230206_2023-02-07_0947.txt",sep = "\t",row.names = F)
    
    df_mean_ab_after_impt <- df_merge %>% 
        select(pep_with_pos,contains("aft_imp") | contains("species"),Pool) %>% #accession
        tibble() %>% 
        #separate(accession, into = c("uniprot_id", "species", "position"), remove = F) %>%
        pivot_longer(cols = contains("aft_imp"),
                     names_to = "Mean_abundance",
                     values_to = "values") %>%
        mutate_at("Pool", ~replace_na(.,background_species)) 
    
    
    
    # for (i in 1:(sample_size-1)){
    #   
    #   dftmp_mean <-  df_merge %>% select(pep_with_pos,contains(paste0("log10_E2-A",i))) %>%
    #     rowwise() %>%
    #     mutate(row_mean = mean(c_across(where(is.numeric)), na.rm = TRUE)) %>%
    #     select(row_mean)
    #   
    #   if (i==1){
    #     
    #     assign(paste0('log10_E2_A_CV',i), df_merge %>% 
    #              select(pep_with_pos,contains(paste0("log10_E2-A",i)),Pool,accession) %>%
    #              rowwise() %>%
    #              mutate(row_sd = sd(c_across(where(is.numeric)), na.rm = TRUE)) %>%
    #              bind_cols(dftmp_mean) %>%
    #              mutate(CV=row_sd/row_mean) %>%
    #              select(!contains("log10_")))
    #     }else{
    #       
    #       assign(paste0('log10_E2_A_CV',i), df_merge %>% 
    #              select(contains(paste0("log10_E2-A",i))) %>%
    #              rowwise() %>%
    #              mutate(row_sd = sd(c_across(where(is.numeric)), na.rm = TRUE)) %>%
    #              bind_cols(dftmp_mean) %>%
    #              mutate(CV=row_sd/row_mean) %>%
    #              select(CV))
    #       }
    #   print(i)
    # }
    # 
    # colnames(log10_E2_A_CV1)[6] <-"CV1" 
    # colnames(log10_E2_A_CV2) <- "CV2"
    # colnames(log10_E2_A_CV3) <- "CV3"
    # colnames(log10_E2_A_CV4) <- "CV4"
    # colnames(log10_E2_A_CV5) <- "CV5"
    # 
    # df_CV_aft_impt <- log10_E2_A_CV1 %>% 
    #   
    #   bind_cols(log10_E2_A_CV2,
    #             log10_E2_A_CV3,
    #             log10_E2_A_CV4,
    #             log10_E2_A_CV5) %>% 
    #   pivot_longer(cols = contains("CV"),
    #                names_to = "CV_samples",
    #                values_to = "CV_values") %>%
    #   filter(!grepl("Unexpected",Pool)) %>%
    #   mutate(CV_values=(CV_values*100)) %>%
    #   mutate(acq_type=acquisiton_type) %>%
    #   mutate(soft_name=software_name)
    #   
    
    
    # df_FC_ratio_after_impt <- df_merge %>% 
    #     select(starts_with("exp_") | contains("species"),Pool) %>% #accession
    #     #separate(accession, into = c("uniprot_id", "species", "position"), remove = F) %>%
    #     tibble() %>% 
    #     pivot_longer(cols = starts_with("exp_"),
    #                  names_to = "exp_FC",
    #                  values_to = "values") %>%
    #     mutate_at("Pool", ~replace_na(.,background_species))
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
    ) + scale_color_manual(values = c("#08519c","#3182bd","#6baed6","#a6bddb"))
      #+ geom_errorbar(data = test,aes(ymin = test$min_val, ymax=  test$max_val,color=test$exp_FC), width=0.5) 
    df_FC_ratio_after_impt_final <- df_FC_ratio_after_impt %>% mutate(acquisition=acquisiton_type) %>% mutate(software_name=software_name)
    
    write.table(df_FC_ratio_after_impt_final,file = paste0(new_path,"/df_FC_ratio_after_impt",software_name,".txt"),sep = "\t",row.names =F )
 #############################################################################
    plot4 <- gg_density(data_set = df_mean_ab_after_impt, 
                     x_df = df_mean_ab_after_impt$values,
                     fill_df = df_mean_ab_after_impt$species,
                     color_df = NULL, #df_mean_ab_after_impt$Pool,
                     header="Distribution of mean abundance of every sample after imputation",
                     facet_df = "Mean_abundance",
                     x_lab = "log10(values)",
                     color_lab= "",
                     fill_lab = "Sample Names",
                     subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
    
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
    
    ############################## ############################
    
    ##########################################################
    ### BOX-PLOT: Experimental Quantity Ratio of Phospho Peptides  
    
    # df_FC_ratio_absErr_after_impt <- df_merge %>% 
    #   select(pep_with_pos,starts_with("exp_") | contains("accession"),Pool,isomericity) %>%
    #   separate(accession, into = c("uniprot_id", "species", "position"), remove = F) %>%
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
    # 
    # 
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
  
    ############################## ############################
    ############      ############    CV ASSESSMENT BEFORE & AFTER IMPUTATION      ############      ############  
    ### BEFEORE  ### 
    filtered_abundances_ecoli_bef_impt <- barplt_df_ecoli_wide %>% 
      #separate(accession, into = c("uniprot_id", "species", "position"), remove = T) %>%
      rename(Pool=species) %>%
      rename(pep_with_pos=Sequence) %>%
      #select(!c(uniprot_id,position)) %>% 
      drop_na()
    
    pep_list_w_theo_iso <- pep_list_w_theo %>% select(pep_with_pos,isomericity)
    
    df_merge_before_imp <- df_merge_syn %>% 
      left_join(pep_list_w_theo_iso,by="pep_with_pos") %>%
      mutate(Pool=ifelse(Pool=="Diluted",paste0(Pool,"_",isomericity),Pool)) %>% 
      #pivot_longer(cols = contains("E2"),
      #names_to = "Experiment",
      #values_to = "Abundances") %>% 
      drop_na() %>% select(!isomericity) %>%
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
                                       exp_id = 2,replicate = num_reps,
                                       iter = 5,acquisiton_type = acquisiton_type,
                                       software_name = software_name)
    
    write.table(df_before_impt_CV,file = paste0(new_path,"/CV_calculation_before_imputation",acquisiton_type,software_name,".txt"))
    
    ############      ############                     ############      ############  
    ### AFTER  ### 
    
    df_merge_aft_impt <- df_merge_syn_after_impt %>% #mutate(Pool=ifelse(Pool=="Diluted",paste0(Pool,"_",isomericity),Pool)) %>%
      select(!c(species)) %>% select(!contains("A6") & !contains("log10"))
    
    df_aft_impt_CV <- CV_calculator(data = df_merge_aft_impt,abundance_col = "mean_abundances_aft_imp_A",
                                    iter = 5, exp_id = 2,replicate = num_reps,
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
    library(gghalves)
    p20 <-  df_before_impt_CV %>% filter(grepl("Diluted",Pool) | grepl("Fixed",Pool)) %>% 
      ggplot(aes(x=CV_samples,y=CV_values ,fill=Pool)) + 
      geom_boxplot(aes(x=CV_samples,y=CV_values ,fill=Pool)) + theme_bw() +
      #geom_half_violin(side = "r")+
      #facet_wrap(~Pool) +
      theme(legend.text = element_text(size = 45),
            axis.title.x = element_text(size = 45),
            axis.title.y = element_text(size = 45),
            plot.title = element_text(size = 55),
            legend.title = element_text(size = 45),
            axis.text.x = element_text(size = 45),
            axis.title = element_text(size = 45),
            axis.text.y = element_text(size = 45),
            plot.subtitle = element_text(size = 45)) +
      labs( x= "Sample ID", y="Percentage of CVs",
            title = paste("Experiment 2 Percentage of CVs", acquisiton_type, " data processed by ", software_name), 
            subtitle = paste("Before imputation")) +
      scale_fill_manual(values = c("#6A51A3","#6A51A3","#E6550D","#E6550D","cyan3")) + #"grey68"
      scale_y_continuous(
        limits = c(0,max(df_aft_impt_CV$CV_values)), 
        breaks = seq(0, max(df_aft_impt_CV$CV_values),10))
    
    p21 <-  df_aft_impt_CV %>% filter(grepl("Diluted",Pool) | grepl("Fixed",Pool)) %>% 
      ggplot(aes(x=CV_samples,y=CV_values ,fill=Pool)) +
      geom_boxplot(aes(x=CV_samples,y=CV_values ,fill=Pool)) + theme_bw() +
      #geom_half_violin(side = "r")+
      #facet_wrap(~Pool) +
      theme(legend.text = element_text(size = 45),
            axis.title.x = element_text(size = 45),
            axis.title.y = element_text(size = 45),
            plot.title = element_text(size = 55),
            legend.title = element_text(size = 45),
            axis.text.x = element_text(size = 45),
            axis.title = element_text(size = 45),
            axis.text.y = element_text(size = 45),
            plot.subtitle = element_text(size = 45)) +
      labs( x= "Sample ID", y="Percentage of CVs",
            #title = paste("Experiment 2 Percentage of CVs", acquisiton_type, " data processed by ", software_name), 
            subtitle = paste("After imputation")) +
      scale_fill_manual(values = c("#6A51A3","#6A51A3","#E6550D","#E6550D","cyan3")) + #"grey68"
      scale_y_continuous(
        limits = c(0,max(df_aft_impt_CV$CV_values)), 
        breaks = seq(0, max(df_aft_impt_CV$CV_values),10))
    
    
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
    
    ############################## ############################
    plot5 <- gg_density(data_set = df_FC_ratio_after_impt,
                     x_df = df_FC_ratio_after_impt$values,
                     fill_df = df_FC_ratio_after_impt$exp_FC,
                     color_df = df_FC_ratio_after_impt$Pool,
                     header="Distribution of Fold change Ratio of every sample after imputation",
                     facet_df = "exp_FC",
                     x_lab = "log10(values)",
                     color_lab= "",
                     fill_lab = "Sample Names",
                     subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
    # 
    
    ### BOX-PLOT: Experimental Quantity Ratio of Phospho Peptides  
    
    # plot6 <- gg_boxplt_exp_ratio(data_set = df_FC_ratio_after_impt, 
    #                           x_df = df_FC_ratio_after_impt$exp_FC,
    #                           y_df = df_FC_ratio_after_impt$values,
    #                           fill_df = df_FC_ratio_after_impt$Pool,
    #                           header="Experimental Quantity Ratio of Phospho Peptides",
    #                           x_lab="Sample Names",
    #                           y_lab="Abundance Ratios",
    #                           fill_lab = "Sample Names",
    #                           subtitle_txt = paste("Experiment - ",
    #                                                exp_id, acquisiton_type,
    #                                                " data processed by ",
    #                                                software_name)) +
    #   scale_fill_manual(values = c("#6A51A3","#6A51A3","#E6550D","#E6550D","cyan3",'grey60')) 
    #
    #actual_ratio_col
    library(ggpattern)
    df_FC_ratio_plt6 <- df_FC_ratio_after_impt %>%
      separate(Pool,into = c("Pool_id","isomer"),sep = "_",remove = F) %>%
      mutate(isomer=ifelse(is.na(isomer),"Wrong Localization",isomer)) %>%
      mutate(Pool_id=ifelse(Pool_id=="Unexpected","Wrong Localization",Pool_id))
    
    hline_df <- data.frame(exp_FC =unique(df_FC_ratio_plt6$exp_FC),log2_act_val =log2(actual_ratio[-length(comparisons)]))
    #hline_df <- hline_df %>% separate(exp_FC,into = c("tmp","FC"),sep = "_") %>%
      #mutate(exp_FC=paste("FC isomeric",FC))

     plot6 <-  ggplot(df_FC_ratio_plt6,aes(x=exp_FC,y=log2_exp_val)) +
      geom_boxplot_pattern(aes(fill=Pool,pattern=isomer),
                           pattern_fill = "black" ,#pattern = isomer,
                           pattern_density = 0.1,
                           pattern_spacing = 0.025,
                           pattern_key_scale_factor = 0.6,
                           outlier.shape = NA) +
       #geom_line(data = hline_df,aes(x=as.factor(hline_df$exp_FC),y=hline_df$log2_act_val,group=1)) +
       geom_hline(data = hline_df,aes(yintercept=log2_act_val),color="grey",linetype="dashed") +
       geom_text(data = hline_df,color="grey40",aes(0,log2_act_val,label = paste(hline_df$exp_FC,"=",round(hline_df$log2_act_val)),hjust = -0.05, vjust = -1),size=8)+
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
    
   
    
    ### BOX-PLOT: Experimental Quantity Ratio of Phospho Peptides  
    # df_FC_ratio_absErr_after_impt <- df_merge %>% 
    #   select(pep_with_pos,starts_with("exp_") | contains("accession"),Pool) %>%
    #   separate(accession, into = c("uniprot_id", "species", "position"), remove = F) %>%
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
    #   mutate(abs_err = abs(actual_ratio_val-values)) %>%
    #   filter(!grepl("Unexpected",Pool)) %>%
    #   mutate(acq_type=acquisiton_type) %>%
    #   mutate(soft_name=software_name)
    # 
    # write.table(df_FC_ratio_absErr_after_impt,file = paste0(file_path,"AbsError_FC",software_name,"_",acquisiton_type,".txt"),sep = 
    #               "\t",col.names = T,row.names = F)
    
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
    #scale_y_continuous(limits = c(0,200)) + #breaks = seq(from =0, to=100,by=10)) 
    
    # 
    # plot19 <-  gg_half_boxplt_exp_ratio(data_set = df_CV_aft_impt, 
    #                          x_df = df_CV_aft_impt$CV_samples,
    #                          y_df = df_CV_aft_impt$CV_values,
    #                          fill_df = df_CV_aft_impt$Pool,
    #                          header="Experimental Quantity Ratio of Phospho Peptides",
    #                          x_lab="Sample Names",
    #                          y_lab="Abundance Ratios",
    #                          fill_lab = "Sample Names",
    #                          
    #                          subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) + 
    #   scale_fill_manual(values = c("#6A51A3","#6A51A3","grey68","#E6550D","#E6550D"))+
    #   scale_color_manual(values = c("white","black","grey68","white","black")) + 
    #   theme_dark() 
    # 
    # write.table(df_CV_aft_impt,file = paste0(file_path,"CV_percentage",software_name,"_",acquisiton_type,".txt"),sep = 
    #               "\t",col.names = T,row.names = F)
    # 
    
    
    #cols <- c(rep("#E6550D" ,4),"#3F007D","#6A51A3","#807DBA","#BCBDDC",rep("cyan3",4))
    ### HALF-BOX-PLOT & HALF-SCATTER-PLOT: Experimental Quantity Ratio of Synthetic Peptides  
  
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
    
    #   selected_cols <- colnames(exp2_FAIMS_DDA_proline[c(1:4,6,28)])
    # #selected_cols <- colnames(exp2_noFAIMS_DDA_proline[c(1:4,15)])
    # 
    # syn_phospho_pep_proline <- exp2_FAIMS_DDA_proline %>% 
    #   select(-starts_with("abundance_")) %>% 
    #   mutate(abundances_for_impute) %>%
    #   filter(grepl("Phospho", ptm_protein_positions)) %>% 
    #   filter(grepl("_HUMAN", accession)) %>%
    #   select(contains(selected_cols) | starts_with("abundance"))
    # 
    # 
    # # Sorting columns A1 to A5 and R1 to R3
    # sorted_colnames <- c("psm_count","abundance") # Abundance contains both w/wo raw_abundance
    # for (i in 1:length(sorted_colnames)){
    #   assign(paste0("tmp",i), sort(colnames(select(syn_phospho_pep_proline_new,
    #                                                contains(sorted_colnames[i])))))
    # }
    # 
    # sort_col <- c(tmp1,tmp2)
    # 
    # sort_syn_phospho_pep_proline_new <- syn_phospho_pep_proline_new %>%
    #   relocate(all_of(sort_col), .after = 7) %>% rename_with(~ experiment_name, all_of(sort_col)) 
    # 
    # 
    ##### REMOVE REDUNDANCY OF PHOSPHO SEQUENCES DUE TO DIFFERENT CHARGES
    # syn_phospho_pep_proline_unique_max_absum <- sort_syn_phospho_pep_proline_new %>% 
    #   #mutate(across(all_of(experiment_name), ~ ifelse(is.na(.), 0, .))) %>%
    #   rowwise() %>%
    #   mutate(row_sum = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
    #   group_by(common_col_for_merging) %>%
    #   #summarise(max_row_sum = max(row_sum, na.rm = TRUE)) %>%
    #   slice(which.max(row_sum)) %>%
    #   ungroup()
    
    ####################################################
    # #### IMPUTATION BASED ON COMPLEX CONDITIONS ####
    # 
    # df_long <- syn_phospho_pep_proline_unique_max_absum %>% 
    #   select(common_col_for_merging,all_of(experiment_name)) %>%
    #   pivot_longer(cols = -common_col_for_merging, names_to = "column", values_to = "value") %>%
    #   separate(column, into = c("col_group", "col_index", "row_index"), sep = "_")
    # #mutate(col_index = str_remove(col_index, "A"))
    # 
    # # Count the number of NA values in each row for each column group
    # na_counts <- df_long %>%
    #   group_by(common_col_for_merging, col_index) %>%
    #   summarise(na_count = sum(is.na(value))) 
    #   #pivot_wider(names_from = col_index, values_from = na_count)
    #   
    #   # Replace the NA values based on the count
    #   df_result <- df_long %>%
    #   left_join(na_counts, by = c("common_col_for_merging", "col_index")) %>%
    #   mutate(value = case_when(
    #     na_count == 3 ~ 0,
    #     
    #     na_count == 2 & col_index == "A1" ~ impute_values_2NA[1],
    #     na_count == 1 & col_index == "A1" ~ impute_values_1NA[1],
    #     
    #     na_count == 2 & col_index == "A2" ~ impute_values_2NA[2],
    #     na_count == 1 & col_index == "A2" ~ impute_values_1NA[2],
    #     
    #     na_count == 2 & col_index == "A3" ~ impute_values_2NA[3],
    #     na_count == 1 & col_index == "A3" ~ impute_values_1NA[3],
    #     
    #     na_count == 2 & col_index == "A4" ~ impute_values_2NA[4],
    #     na_count == 1 & col_index == "A4" ~ impute_values_1NA[4],
    #     
    #     na_count == 2 & col_index == "A5" ~ impute_values_2NA[5],
    #     na_count == 1 & col_index == "A5" ~ impute_values_1NA[5],
    #     TRUE ~ value
    #   )) %>%
    #   select(-na_count) %>% unite(new_col,col_group,col_index, row_index,sep = "_") %>%
    #   pivot_wider(names_from = new_col, values_from = value) %>%
    #   #select(-col_group) %>% 
    #   filter(rowSums(across(where(is.numeric)))!=0)
    # 
    # 
    # ####################
    
    # Combine data frame (we will continue with this for further step)
    #df_merge <- df_result %>% left_join(syn_phospho_pep_proline_unique_max_absum, 
    # by="common_col_for_merging") 
    #df_merge <- syn_phospho_pep_proline_unique_max_absum %>% left_join(pep_list_w_theo_quant_new,by="common_col_for_merging")
    
    ## ALWAYS KEEP IT COMMENTED TO PREVENT OVERWRITE
    #write.xlsx(df_merge,file = "D:/dev/Pinar/PHD/wet_lab_experiments/experiment2_quantification_peptide_level/Correct_imputation_exp2proline_and_peplist_theo_before_stast_analysis.xlsx")
    
    
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
                                select(stat_analysis, contains(paste0("A", i, "_")) & contains("log10_"))[j,])
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
        merge_stat_df_final <- merge_stat_df %>% filter(!grepl("A1/A6",fold_change_comp))
        
        
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
        
    }else{
        print("Statistical test could not be assessed. Check the input files!")
    }
    ### VISUALIZATION OF P-VALUE DISTRIBUTION ###
    #ggplot(merge_stat_df_final,aes(x=-log10(P.Value),color=Pool_new)) + geom_density(linewidth=1) + facet_wrap(~ratio)
  
   
    comparisons <- comparisons[-5]
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
                                           grepl(comparisons[4],fold_change_comp) ~actual_ratio[4]))#,
                                           #grepl(comparisons[5],fold_change_comp) ~actual_ratio[5]))
                                          
    }else{
      print("Mapping between theoretical ratio and comparison cannot be done. Please make sure that you have either 4 or 5 comparisons overall.")
    }
    
    point_count_y_axis <- merge_stat_df_final %>%
        group_by(fold_change_comp, new_col_coloring, coloring_comp) %>% #new_col_coloring
        filter(P.Value < 0.05) %>% 
        count(new_col_coloring) %>% left_join(actual_ratio_col)
    
   
    ymax <- 11 + 0.5
    y_decrement <- 0.5
    
    #ymax <- 9 + 0.5 #max(-log10(merge_stat_df_final$P.Value))
    #y_decrement <- 0.50
    
    calculate_y_pos <- function(group) {
      group_length <- length(group)
      y_pos <- ymax - seq(0, by = y_decrement, length.out = group_length)
      return(y_pos)
    }
    
    # Apply the function to calculate y_pos within each group
    #point_count_y_axis$y_pos <- unlist(by(point_count_y_axis$A1vs_Ai, point_count_y_axis$A1vs_Ai, calculate_y_pos))
    
    # Apply the function to calculate y_pos within each group
    point_count_y_axis$y_pos <- unlist(by(point_count_y_axis$fold_change_comp, 
                                          point_count_y_axis$fold_change_comp, calculate_y_pos))
    
    # Apply the function to calculate y_pos within each group
    ######################################################################################################
    ######### COLORING STRATEGY BASED ON ISOMERICITY AND RATIO OF EACH SAMPLE (new_col_coloring) #########
    # col <- RColorBrewer::brewer.pal(n=length(comparisons),name = "Dark2")
    # 
    # labels <- unique(merge_stat_df_final$new_col_coloring)
    # 
    # if (sheet_theo_name == "ISOREF_REF2_Others"){
    #   
    #   colors <- c(RColorBrewer::brewer.pal(n=length(comparisons),name = "Dark2"),"#2171b5","#999999")
    #   new_comparisons <- c(comparisons,"Unexpected","Fixed_mono")
    #   
    # }else if (sheet_theo_name == "ISO-refOTHER with FC_correct"){
    #   
    #   new_comparisons <- c(comparisons,"Unexpected")
    #   colors <- c(RColorBrewer::brewer.pal(n=length(comparisons),name = "Dark2"),"#999999")
    #   
    # }else{
    #   print("Sheet_theo_name could not be found, please make sure that you selected the correct sheet_name.")
    # }
    # 
    # mapped_coloring <- rep("#000000",length(labels))
    # 
    # for (i in 1:length(labels)) {
    #   # Check if the color_element contains any of the comparisons
    #   if (any(new_comparisons %in% str_extract_all(labels[i], paste(new_comparisons, collapse = "|"))[[1]])) {
    #     # Find the index of the matching comparison in the comparisons list
    #     comp_index <- match(TRUE, sapply(new_comparisons, function(comp) comp %in% str_extract_all(labels[i], comp)))
    #     
    #     # Assign the corresponding color to the data frame
    #     mapped_coloring[i]<- paste0(colors[comp_index])
    #     #mapped_coloring[i]<- paste(paste0(labels[i],'"'),paste0('"',colors[comp_index]),sep = "=")
    #   }
    # }
    # 
    # actual_ratio_col <- actual_ratio_col %>% bind_cols(col)
    # 
    # mapped_coloring_dat <- as.data.frame(mapped_coloring)
    # 
    # point_count_y_axis <- mapped_coloring_dat %>% 
    #   bind_cols(labels) %>%
    #   rename(colors=1,new_col_coloring=2) %>%
    #   right_join(point_count_y_axis,by="new_col_coloring")
    
###############################################################################
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
    
    
    
  ########################################################################################  
    
    
    
  
    
    # plot9 <- ggplot(merge_stat_df_final,aes(x =log2(merge_stat_df_final$fold_change_values), y = -log10(merge_stat_df_final$P.Value))) +
    #   geom_point(size = 4, aes(color = new_col_coloring,shape=Pool)) + # Pool_new might be use to 
    #   #scale_shape_identity() +                                        # differentiate peptides are found as Unexpected and reference 
    #                                                                    # in their associated concentration 
    #                                                                    # Pool can be used to show only difference btw pools same as coloring 'less complex visualization)
    #   #geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "red") +
    #   #scale_fill_manual(values=setNames(mapped_coloring, labels)) + 
    #   #geom_line(aes(color =setNames(mapped_coloring, labels)), size = 1) +  # Add color aesthetic to geom_line()
    #   scale_color_manual(values =setNames(mapped_coloring, labels)) +
    #   #scale_x_continuous(breaks = seq(from =round(min(log2(merge_stat_df_final$fold_change_values))), to=(round(max(log2(merge_stat_df_final$fold_change_values)))+2),by=2)) +
    #   #scale_y_continuous(breaks = seq(from =round(min(-log10(merge_stat_df_final$P.Value))), to=(round(max(-log10(merge_stat_df_final$P.Value)))+2),by=2)) +
    #   scale_y_continuous(
    #     limits = c(0,12), 
    #     breaks = seq(0, 12,2)
    #   ) +
    #   scale_x_continuous(
    #     limits = c(-9,11), 
    #     breaks = seq(-9, 11,2)
    #   ) +
    #   
    #   #expand_limits(x=c(-7,11), y=c(0, 16))+
    #   #scale_x_continuous(limits = c(-7,11)) +
    #   #scale_y_continuous(limits = c(0,16)) +
    #   #scale_y_continuous(breaks = seq(0, max(-log10(volcano_final1$pvalues_value)), length.out = 21)) +
    #   theme_bw() +
    #   theme(legend.text = element_text(size = 45),
    #         axis.title.x = element_text(size = 45),
    #         axis.title.y = element_text(size = 45),
    #         plot.title = element_text(size = 55),
    #         legend.title = element_text(size = 45),
    #         axis.text.x = element_text(size = 45),
    #         axis.title = element_text(size = 45),
    #         axis.text.y = element_text(size = 45),
    #         plot.subtitle = element_text(size = 45)) +
    #   
    #   labs( y= "-log10(p values)", x="log2(fold change)",title = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name), subtitle = paste("Limma was used \n", subtitle)) +
    #   geom_vline(data = actual_ratio_col, aes(xintercept = log2(actual_ratio_val), show.legend = FALSE),color=col,size=1.5) +
    #   geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "red",size=1.5) + 
    #   geom_label(data = point_count_y_axis, aes(x = log2(actual_ratio_val), y = y_pos, fill=new_col_coloring,label = n),color="white",size=14,show.legend = FALSE) +
    #   scale_fill_manual(values =setNames(mapped_coloring, labels))
    # 

    merge_stat_df_final_text <- merge_stat_df_final %>% 
      mutate(soft_name=paste0(software_name)) %>%
      mutate(acq_type=paste0(acquisiton_type))
    
    write.table(merge_stat_df_final_text,file = paste0(new_path,"/volcano_plot_",software_name,"_",acquisiton_type,".txt"),sep = 
            "\t",col.names = T,row.names = F)
    
    df_roc <- merge_stat_df_final %>%
      select(pep_with_pos, Pool,P.Value) %>%
      mutate(Pool = if_else(grepl("Diluted_isomeric", Pool), "Diluted", Pool)) %>%
      mutate(Pool = if_else(grepl("Diluted_nonisomeric", Pool), "Diluted", Pool))
    
   
    ### ROC analysis custom func
    df_roc_order <- df_roc[order(df_roc$P.Value),]
    
    df_roc_func <- compute_roc_curve(df=df_roc_order, flag = "Diluted",expected =length(comparisons)*size_variying_pep_size)
   
    
    write.table(df_roc_func, file = paste0(new_path,"/new_custom_Roc_analysis_",exp_id,"_",software_name,"_",".txt"),sep = "\t",row.names = F)
    #write.table(df_roc_func_filt, file = paste0(new_path,"/new_custom_Roc_analysis_",exp_id,"_",software_name,"_","filtered.txt"),sep = "\t",row.names = F)
    
   
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
    
    
    #### ROC Analysis using pROC 
    
    df_roc$variant <- ifelse(df_roc$Pool == "Diluted", TRUE, FALSE)
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
                                  "Diluted Pool")
    
    
    roc_plt_df <- as.data.frame(tpr_and_fpr_variant) %>% 
        #bind_rows(as.data.frame(tpr_and_fpr_non_var)) %>% 
      bind_cols(software_name)
    
    colnames(roc_plt_df) <-    c("sensitivity", "fpr","Pool_type","Software_name")
    
    
    plot10 <- roc_plt_df %>% group_by(Pool_type) %>% 
        ggplot( aes(y=as.numeric(sensitivity), x = as.numeric(fpr), color=Pool_type)) +
        geom_path(size=1.5) + theme_bw() + #scale_x_reverse() 
        theme(legend.text = element_text(size = 20),
              axis.title.x = element_text(size = 20),
              axis.title.y = element_text(size = 20),
              plot.title = element_text(size = 25),
              legend.title = element_text(size = 20),
              axis.text.x = element_text(size = 20),
              axis.title = element_text(size = 20),
              axis.text.y = element_text(size = 20)) +
        scale_color_brewer(palette = "Dark2") +
        labs(y="True Positive Rate \n (Sensitivity)", x="False Discovery Rate",
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

    plt_obj <- ls(pattern="plot")
    plt_obj <- plt_obj[!is.na(plt_obj)]
    sapply(1:length(plt_obj),function(x) ggsave(filename = paste0("p",x,".png"),
                                                 width = 90, height = 60, 
                                                 path = paste0(file_path,"/",dir_name,"/"),
                                                 units = "cm",
                                                 get(plt_obj[x]),
                                                 device = "png", #".svg"
    ))
    
    
    
    #  
    #  
    #  
    #  
    #  
    #  
    #  df_roc <- merge_stat_df_final %>%
    #      select(pep_with_pos.x, Pool.x,P.Value) %>%
    #      
    #    
    #    roc_data_all <- compute_roc_curve(complete_pvalues_all,flag = "Others",expected = 145)
    #    
    #    
    #    ## ROC CURVE GENERATION ### 
    #    roc_data121 <- cbind(roc_data12, "A1 vs A2","Others", "Proline")
    #    roc_data131 <- cbind(roc_data13, "A1 vs A3","Others", "Proline")
    #    roc_data141 <- cbind(roc_data14, "A1 vs A4","Others", "Proline")
    #    roc_data151 <- cbind(roc_data15, "A1 vs A5","Others", "Proline")
    #    
    #    colnames(roc_data121)[c(4:6)] <- c("Comparison","Pool","Software")
    #    colnames(roc_data131)[c(4:6)] <- c("Comparison","Pool","Software")
    #    colnames(roc_data141)[c(4:6)] <- c("Comparison","Pool","Software")
    #    colnames(roc_data151)[c(4:6)] <- c("Comparison","Pool","Software")
    #    
    #    roc_final <- rbind(roc_data121,roc_data131,roc_data141,roc_data151)
    #    
    #    #roc_curve_data_exp2_noFAIMS_MQ <- read.delim("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/MQ_data_analysis/exp2_wo_FAIMS/exp2_wo_FAIMSwith_MBR/roc_curve_data_exp2_noFAIMS_MQ.txt")
    #    write.table(roc_final,file = "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/exp2_re_injection/exp2_no_FAIMS_proline_roc_curve_data_consider_isoref_fp.tsv",sep = "\t",row.names = F)
    #    #roc_curve_data_pro_faims <- read.delim("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/exp2_re_injection/exp2wo_FAIMS_DDA_Proline_volcano_plot_data.tsv",sep = "\t")
    #    
    #    
    #    roc_proline_mq <- rbind(roc_curve_data_exp2_noFAIMS_MQ,roc_final)
    #    
    #    
    #    roc_data_unexpected121 <- cbind(roc_data_unexpected12, "A1 vs A2","unexcepted")
    #    roc_data_unexpected131 <- cbind(roc_data_unexpected13, "A1 vs A3","unexcepted")
    #    roc_data_unexpected141 <- cbind(roc_data_unexpected14, "A1 vs A4","unexcepted")
    #    roc_data_unexpected151 <- cbind(roc_data_unexpected15, "A1 vs A5","unexcepted")
    #    
    #    colnames(roc_data_unexpected121)[c(4,5)] <- c("Comparison","Pool")
    #    colnames(roc_data_unexpected131)[c(4,5)] <- c("Comparison","Pool")
    #    colnames(roc_data_unexpected141)[c(4,5)] <- c("Comparison","Pool")
    #    colnames(roc_data_unexpected151)[c(4,5)] <- c("Comparison","Pool")
    #    
    #    
    #    
    #    roc_final2 <- rbind(roc_data_unexpected121,roc_data_unexpected131,
    #                        roc_data_unexpected141,roc_data_unexpected151)
    #    
    #    roc_curve_combine_proline_mq_score_40 <- read.delim("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/comparision_DDA_noFAIMS_mq_proline_pd/comparison_wrong_localized_pep_MQ_Proline_with_out_FAIMS/final_figures/roc_curve_combine_proline_mq_score_40.txt")
    #    
    #    line_types <- c("solid","dashed" )
    #    
    #    ##TODO: make it more professional
    #    ggplot(roc_final,aes(x=fdp,y=tpr)) + 
    #      geom_line(size=1.2,aes(color=Comparison,linetype=Software)) +
    #      scale_color_manual(values = c("#CC79A7", "#E69F00", "#56B4E9", "#009E73"),
    #                         labels=c('A1 vs A2','A1 vs A3','A1 vs A4','A1 vs A5')) +
    #      scale_linetype_manual(values = line_types)+
    #      theme_bw() +
    #      theme(legend.text = element_text(size=15), 
    #            axis.title.x = element_text(size = 15),
    #            axis.title.y = element_text(size = 15),
    #            plot.title = element_text(size=30),
    #            legend.title=element_text(size=15),
    #            axis.text.x=element_text(size=15),
    #            axis.title=element_text(size=15),
    #            axis.text.y = element_text(size = 15)) + 
    #      expand_limits(x = 0, y = 0) +
    #      scale_y_continuous(limits = c(0,100)) +
    #      #facet_wrap(~Ratio_col,scales = "free_x") +
    #      labs(title = "Experiment 2 - DDA without FAIMS - \n Comparision of Softwares")
    #    
    # #    
    # #  # 
    # #  # # iterate through each column of the data set and calculate t-tests
    # #  # triplicate_indx <- c(1,3,4,6,7,9,10,12,13,15) +1
    # #  # 
    # # 
    # #  p_values_12 <- NULL
    # #  p_values_13 <- NULL
    # #  p_values_14 <- NULL
    # #  p_values_15 <- NULL
    # # 
    # #  ##TODO put this into another for loop to reduce redundancy of the code
    # #    for(j in 1:dim(stat_analysis)[1]){
    # # 
    # #      #if ((sum(stat_analysis[j,2:4])/3) !=  stat_analysis[j,2:4][1] | (sum(stat_analysis[j,5:7])/3) !=  stat_analysis[j,5:7][1]){
    # # 
    # #      p_values_12[j] <- ttest_func(stat_analysis[j,3:5], stat_analysis[j,6:8])
    # #      #p_values_13[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,8:10])
    # #      #p_values_14[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,11:13])
    # #      #p_values_15[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,14:16])
    # #      #}else{
    # # 
    # #      #}
    # #      #if ((sum(stat_analysis[j,2:4])/3) !=  stat_analysis[j,2:4][1] | (sum(stat_analysis[j,8:10])/3) !=  stat_analysis[j,8:10][1]){
    # # 
    # #        p_values_13[j] <- ttest_func(stat_analysis[j,3:5], stat_analysis[j,9:11])
    # # 
    # #      #}else{
    # # 
    # #      #}
    # #      #if ((sum(stat_analysis[j,2:4])/3) !=  stat_analysis[j,2:4][1] | (sum(stat_analysis[j,11:13])/3) !=  stat_analysis[j,11:13][1]){
    # # 
    # # 
    # #        p_values_14[j] <- ttest_func(stat_analysis[j,3:5], stat_analysis[j,12:14])
    # # 
    # #      #}else{
    # # 
    # #     # }
    # #      #if ((sum(stat_analysis[j,2:4])/3) !=  stat_analysis[j,2:4][1] | (sum(stat_analysis[j,14:16])/3) !=  stat_analysis[j,14:16][1]){
    # # 
    # #        p_values_15[j] <- ttest_func(stat_analysis[j,3:5], stat_analysis[j,15:17])
    # #     # }else{
    # # 
    # #     # }
    # #      }
    # # 
    # # 
    # #  p_values_12 <- as.data.frame(p_values_12)
    # #  p_values_13 <- as.data.frame(p_values_13)
    # #  p_values_14 <- as.data.frame(p_values_14)
    # #  p_values_15 <- as.data.frame(p_values_15)
    # # 
    # # 
    # #    for (i in 2:sample_size){
    # # 
    # #      assign(paste0("p_values_1",i),
    # #             cbind(stat_analysis$common_col_for_merging,
    # #                   as.data.frame(get(paste0("p_values_1",i)))))
    # # 
    # #    }
    # # 
    # #  
    # #  test_num <- nrow(stat_analysis)
    # #  
    # #  for (j in 1:nrow(stat_analysis)){
    # #    
    # #    p_values_12 <- p_values_12[order(p_values_12$p_values_12),]
    # #    colnames(p_values_12)[1] <- "common_col_for_merging"
    # #    p_values_12["rank12"] <- 1:test_num
    # #    p_values_12[j,"test12"] <- (p_values_12[j,"rank12"]/test_num)*0.05
    # #    p_values_12[j,"test_bool"] <- p_values_12[j,"test12"] > p_values_12[j,"p_values_12"]
    # #  
    # #    
    # #    p_values_13 <- p_values_13[order(p_values_13$p_values_13),]
    # #    colnames(p_values_13)[1] <- "common_col_for_merging"
    # #    p_values_13["rank13"] <- 1:test_num
    # #    p_values_13[j,"test13"] <- (p_values_13[j,"rank13"]/test_num)*0.05
    # #    p_values_13[j,"test_bool"] <- p_values_13[j,"test13"] > p_values_13[j,"p_values_13"]
    # #    
    # #    
    # #    p_values_14 <- p_values_14[order(p_values_14$p_values_14),]
    # #    colnames(p_values_14)[1] <- "common_col_for_merging"
    # #    p_values_14["rank14"] <- 1:test_num
    # #    p_values_14[j,"test14"] <- (p_values_14[j,"rank14"]/test_num)*0.05
    # #    p_values_14[j,"test_bool"] <- p_values_14[j,"test14"] > p_values_14[j,"p_values_14"]
    # #    
    # #    
    # #    p_values_15 <- p_values_15[order(p_values_15$p_values_15),]
    # #    colnames(p_values_15)[1] <- "common_col_for_merging"
    # #    p_values_15["rank15"] <- 1:test_num
    # #    p_values_15[j,"test15"] <- (p_values_15[j,"rank15"]/test_num)*0.05
    # #    p_values_15[j,"test_bool"] <- p_values_15[j,"test15"] > p_values_15[j,"p_values_15"]
    # #    
    # #  }
    # # ### TODO: Extract where the first FALSE was generated to use for threshold of each comparison.
    # #  #p_thresholds <- c(2.100824e-02,2.381169e-02,3.119312e-02,2.567998e-02)
    # #  p_thresholds<- c(2.045695e-02,2.359024e-02,2.751729e-02,2.703404e-02)
    # #  col_sel_stat_analysis <- c("common_col_for_merging",
    # #                             #"modifications",
    # #                             #"accession",
    # #                             "isomericity",
    # #                             "isomeric_count",
    # #                             #"pool_id",
    # #                             "Pool",
    # #                             "A1-A2_Ratio",
    # #                             "A1-A3_Ratio",
    # #                             "A1-A3_Ratio",
    # #                             "A1-A4_Ratio",
    # #                             "A1-A5_Ratio",
    # #                             "exp_FC_A1/A2",
    # #                             "exp_FC_A1/A3",
    # #                             "exp_FC_A1/A4",
    # #                             "exp_FC_A1/A5")
    # #    
    # #  
    # 
    #  # for (i in 12:15){
    #  #   assign(paste0("complete_p_values_",i), 
    #  #           final_imputed_data %>% 
    #  #             select(col_sel_stat_analysis) %>%
    #  #             left_join(get(paste0("p_values_",i)),
    #  #                       by="common_col_for_merging") %>%
    #  #             mutate(ls_correctness = ifelse(test_bool == FALSE, 0, 1)))}
    #  #   
    #  #   assign(paste0("roc_data",i),compute_roc_curve(get(paste0("complete_p_values_",i)),flag = "Others",expected = 145 ))
    #  #   #assign(paste0("roc_data_unexpected",i),compute_roc_curve(get(paste0("complete_p_values_",i)),flag = "unexpected",expected = 100))
    #  #   #assign(paste0("roc_data_isoref",i),compute_roc_curve(get(paste0("complete_p_values_",i)),flag = "ISO-REF",expected = 37 ))
    #  #   
    #  # common_colnames <-c("common_col_for_merging", "p_values","rank","test","test_bool")
    #  # 
    #  # 
    #  # 
    #  # colnames(p_values_12) <- common_colnames
    #  # colnames(p_values_13) <- common_colnames
    #  # colnames(p_values_14) <- common_colnames
    #  # colnames(p_values_15) <- common_colnames
    #  # complete_pvalues_all <- rbind(p_values_12,p_values_13,p_values_14,p_values_15)
    # 
    #  
    #  ##TODO: make it more professional
    #  ggplot(roc_data12,aes(x=fdp,y=tpr)) + geom_line()
    #  
    #  
    #  num_wrong_local_pep_Proline_wo_FAIMS <- c("ARSRTPPSAPSQSR_3&11",
    #  "ARSRTPPSAPSQSR_3&13",
    #  "ATSLPSLDTPGELR_2",
    #  "FSDQAGPAIPTSNSYSK_14",
    #  "HTDDEMTGYVATR_12",
    #  "HTDDEMTGYVATR_7",
    #  "LMTGDTYTAHAGAK_8",
    #  "LPLTRSHNNFVAILDLPEGEHQYK_4",
    #  "NGSLKPGSSHR_9",
    #  "RLSSTSLASGHSVR_3&5",
    #  "RLSSTSLASGHSVR_4&5",
    #  "RLSSTSLASGHSVR_4&6",
    #  "RLSSTSLASGHSVR_5",
    #  "RLSSTSLASGHSVR_5&6",
    #  "RLSSTSLASGHSVR_6",
    #  "SFGSPNRAYTHQVVTR_1&10",
    #  "SNSTSSMSSGLPEQDR_1",
    #  "STGDPQGVIR_1",
    #  "STVASMMHR_2",
    #  "TAGTSFMMTPYVVTR_14",
    #  "TGMGSGSAGKEGGPFK_1",
    #  "TVSTSSQPEENVDR_4",
    #  "VIEDNEYTAR_8",
    #  "VSPSPTTYR_6",
    #  "VSPSPTTYR_7",
    #  "YATPQVIQAPGPR_3")
    #  
    #  num_wrong_local_pept_MQ_score_40 <- c("ARSRTPPSAPSQSR_3&13",
    #    "DIYSTDYYR_8",
    #    "ESKSSPRPTAEK_2",
    #    "FSDQAGPAIPTSNSYSK_14",
    #    "IQPAGNTSPR_7",
    #    "LPLTRSHNNFVAILDLPEGEHQYK_4",
    #    "LSYYEYDFER_3",
    #    "RLSSFVTK_7",
    #    "RLSSTSLASGHSVR_3&5",
    #    "RLSSTSLASGHSVR_3&6",
    #    "RLSSTSLASGHSVR_5",
    #    "RSMSPFRGPK_2",
    #    "SFGSPNRAYTHQVVTR_1&10",
    #    "STVASMMHR_2",
    #    "TAGTSFMMTPYVVTR_14",
    #    "TVSTSSQPEENVDR_4&5")
    #  
    #  png_width <- 1800
    #  png_height <- 1500
    #  
    #  venn.diagram(
    #    x = list(num_wrong_local_pep_Proline_wo_FAIMS, num_wrong_local_pept_MQ_score_40),
    #    #category.names = c("Proline_wo_FAIMS", "MQ_wo_FAIMS"),
    #    #fill = c("#00AFBB", "#D16103"),  # Specify custom colors for the sets
    #    alpha = 0.5,  # Set transparency level
    #    fontfamily = "sans",  # Specify font family
    #    fontface = "bold",  # Specify font face
    #    fontcolor = "black",  # Specify font color
    #    fontsize = 12,  # Specify font size
    #    cat.fontfamily = "sans",  # Specify category label font family
    #    cat.fontface = "bold",  # Specify category label font face
    #    cat.fontsize = 12,  # Specify category label font size
    #    cat.cex = 1.2,  # Specify category label expansion factor
    #    cex = 1,  # Specify overall expansion factor for the diagram
    #    filename = "venn.png",  # Optional: specify a filename to save the plot as an image
    #    width = png_width,
    #    height = png_height
    #  )
    #  
    #  
    #  
    #  
    #  merging_all_pvalues <- final_imputed_data %>% 
    #    select(col_sel_stat_analysis) %>%
    #    left_join(p_values_12,by="common_col_for_merging") %>%
    #    left_join(p_values_13,by="common_col_for_merging") %>%
    #    left_join(p_values_14,by="common_col_for_merging") %>%
    #    left_join(p_values_15,by="common_col_for_merging")
    #  
    #  merging_all_pvalues <- merging_all_pvalues %>%
    #    mutate(Pool=replace_na(Pool,"unexpected"))
    #  
    #  
    #  
    #  
    #  
    #    
    #  # 
    #  # 
    #  # 
    #  # 
    #  # 
    #  # #adjusted_p_values <- as.data.frame(p.adjust(p_values_A1_vs_A5, method = "BH", n = length(p_values_A1_vs_A5)))
    #  # adjusted_p_values <- lapply(stat_analysis, function(x) p.adjust (x, method = "BH", n = length(x)))
    #  # p_values_12[,"adj_pvalue12"] <- p.adjust(p_values_12$p_values_12, method = "BH")
    #  # 
    #  # all_adjusted_p_values <- cbind(unlist(adjusted_p_values[["p_values_12"]]),
    #  #                                unlist(adjusted_p_values[["p_values_13"]]),
    #  #                                unlist(adjusted_p_values[["p_values_14"]]),
    #  #                                unlist(adjusted_p_values[["p_values_15"]]))
    #  # colnames(all_adjusted_p_values) <- c("adj_pvalues_1_2","adj_pvalues_1_3","adj_pvalues_1_4","adj_pvalues_1_5")
    #  # 
    #  # 
    #  
    #  ##### p_values_for_all_ratio changed it later to adjusted!!!!
    #  df_volcano <- final_imputed_data %>% 
    #    select(!contains("mean_abundance")) %>%
    #    mutate(as.data.frame(p_values_for_all_ratio)) %>%
    #    #mutate(across(everything(), ~replace(., is.infinite(.), NA)))
    #    rowwise() %>%
    #    filter(across(where(is.numeric), ~!is.infinite(.)))
    #  
    #  
    #  pvalues_df <- merging_all_pvalues %>%
    #    pivot_longer(cols = starts_with("p_values"), names_to = "pvalues_col", values_to = "pvalues_value") %>%
    #    select(common_col_for_merging,Pool,isomericity,pvalues_col,pvalues_value)
    #  
    #  # Pivot the columns containing "Ratio"
    #  ratio_df <- merging_all_pvalues %>%
    #    pivot_longer(cols = ends_with("_Ratio"), names_to = "Ratio_col", values_to = "Ratio_value") %>%
    #    select(Ratio_col,Ratio_value)
    #  
    #  # Pivot the columns containing "exp_FC"
    #  exp_FC_df <- merging_all_pvalues %>%
    #    pivot_longer(cols = starts_with("exp_FC"), names_to = "exp_FC_col", values_to = "exp_FC_value") %>%
    #    select(exp_FC_col,exp_FC_value)
    #  
    #  volcano_final <- bind_cols(pvalues_df,ratio_df,exp_FC_df)
    # 
    #  actual_ratio <- c(2,10,20,100)
    #  actual_ratio <- data.frame(Ratio_col=unique(volcano_final$Ratio_col),log2(c(2,10,20,100)))
    #  
    #  log10_p_thresholds <- data.frame(ratio_col=unique(volcano_final$Ratio_col), -log10(p_thresholds))
    #  
    #  
    #  #### DETERMINE THRESHOLD OF EACH P VALUE COMPARISON BASED on 
    #  ####    WHERE WE CAN SEE THE CHANGE FROM TRUE to FALSE
    #  
    #  ## THIS WILL BE ADDED AS AN EXTRA LINE FOR EACH COMPARISION 
    #  ##     WHEN I MERGED ALL VOLCANO PLOTS INTO ONE
    #  
    #  volcano_final1 <- volcano_final %>% select(-Ratio_value) %>%
    #  mutate(isomericity = ifelse(is.na(isomericity), "False Positive", isomericity)) %>%
    #    unite(Pool_new, Pool, isomericity,sep = "_",remove = FALSE) %>%
    #    unite('new_col_coloring',Pool_new,Ratio_col,sep = "_",remove = FALSE) %>%
    #    mutate(new_col_coloring = if_else(grepl("ISO-REF", new_col_coloring), "ISO-REF", new_col_coloring)) %>%
    #    mutate(new_col_coloring = if_else(grepl("unexpected", new_col_coloring), "unexpected", new_col_coloring))
    #  
    #  #write.table(volcano_final1,file = paste0(file_path,"/","exp2wo_FAIMS_DDA_Proline_volcano_plot_data.tsv"),sep = "\t",row.names = F)
    #  
    #  ggplot(volcano_final1, aes(x = log2(exp_FC_value), y = -log10(pvalues_value))) +
    #    geom_point(aes(shape = Pool_new, color = new_col_coloring), size = 2.5) +
    #    #geom_line(aes(color = new_col_coloring), size = 1) +  # Add color aesthetic to geom_line()
    #    scale_color_manual(values = c("ISO-REF" = "#000000", "unexpected" = "#999999",
    #                                  "Others_multi_A1-A2_Ratio" = "#CC79A7",
    #                                  "Others_mono_A1-A2_Ratio" = "#CC79A7",
    #                                  "Others_multi_A1-A3_Ratio" = "#E69F00",
    #                                  "Others_mono_A1-A3_Ratio" = "#E69F00",
    #                                  "Others_multi_A1-A4_Ratio" = "#56B4E9",
    #                                  "Others_mono_A1-A4_Ratio" = "#56B4E9",
    #                                  "Others_multi_A1-A5_Ratio" = "#009E73",
    #                                  "Others_mono_A1-A5_Ratio" = "#009E73"),
    #                       labels = c('Non-variant', 'Variant non-isomeric A1 vs A2',
    #                                  'Variant non-isomeric A1 vs A3',
    #                                  'Variant non-isomeric A1 vs A4',
    #                                  'Variant non-isomeric A1 vs A5',
    #                                  'Variant isomeric A1 vs A2',
    #                                  'Variant isomeric A1 vs A3',
    #                                  'Variant isomeric A1 vs A4',
    #                                  'Variant isomeric A1 vs A5',
    #                                  'Unexpected')) +
    #    scale_shape_manual(values = c(16, 15, 12, 17),
    #                       labels = c('Non-variant', 'Variant isomeric', 'Variant non-isomeric', 'Unexpected')) +
    #    scale_y_continuous(limits = c(0, 7.2), breaks = seq(0, 7.2, by = 0.8)) +
    #    scale_x_continuous(limits = c(-7,7)) +
    #    #scale_y_continuous(breaks = seq(0, max(-log10(volcano_final1$pvalues_value)), length.out = 21)) +
    #    theme_bw() +
    #    theme(legend.text = element_text(size = 15),
    #          axis.title.x = element_text(size = 15),
    #          axis.title.y = element_text(size = 15),
    #          plot.title = element_text(size = 30),
    #          legend.title = element_text(size = 15),
    #          axis.text.x = element_text(size = 15),
    #          axis.title = element_text(size = 15),
    #          axis.text.y = element_text(size = 15)) +
    #    expand_limits(x = 0, y = 0) +
    #    geom_vline(data = actual_ratio, aes(xintercept = actual_ratio$log2.c.2..10..20..100..),color=c("#CC79A7","#E69F00","#56B4E9","#009E73"), size = 1, show.legend = FALSE) +
    #    geom_hline(data = log10_p_thresholds, aes(yintercept = log10_p_thresholds$X.log10.p_thresholds.),color=c("#CC79A7","#E69F00","#56B4E9","#009E73"), size = 1, linetype = 2, show.legend = FALSE)+ 
    #    labs(title = "Experiment 2 - DDA no FAIMS processed by Proline", color = "Classes", shape="Type")
    #  
    #    #coord_fixed()
    #  
    #  final_imputed_data %>% 
    #    pivot_longer(cols = starts_with("mean_"),
    #                 names_to = "Abundance_col",
    #                 values_to = "Abundance_value") %>%
    #  ggplot(aes(x = Abundance_col, y = common_col_for_merging, fill = Abundance_value)) +
    #    geom_tile(color = "white") +
    #    scale_fill_viridis() +
    #    labs(x = "Condition", y = "Peptide", fill = "Abundance") +
    #    ggtitle("Heatmap")
    #  
    #  
    #  
    #  ggplot(df_volcano, aes(x=df_volcano$`exp_FC_A1/A4`, y=-log10(pvalues_1_4),color=Pool)) + 
    #    geom_point(aes(Pool)) + geom_hline(yintercept = pvalue_threshold, size=1) + 
    #    geom_vline(xintercept = actual_ratio[3], size=1)
    #  #### Distribution of experimental quantitative ratios with two different pools
    #  
    # 
    #  
    #  library(ggplot2)
    #  
    #  ## This can be used as an object name
    #  df_for_figure <- final_imputed_data %>% 
    #    select(contains(col_sel_stat_analysis)) %>%
    #    pivot_longer(cols = col_sel_stat_analysis[3:11], names_to = "Theo_Exp", values_to = "Ratios",values_drop_na = T) %>%
    #    #select(!contains(c("common_col_for_merging"))) %>%
    #    filter(grepl("exp_FC",Theo_Exp)) %>%
    #    filter_at(vars(Pool), all_vars(!is.na(.))) %>%
    #    filter_at(vars(Ratios), all_vars(!is.infinite(.)))
    #  
    #  p1 <- gg_density(data_set = df_for_figure, 
    #                   x_df = df_for_figure$Ratios,
    #                   fill_df = df_for_figure$Theo_Exp,
    #                   color_df = df_for_figure$Pool,
    #                   header="Distribution of experimental quantitative ratios with two different pools",
    #                   facet_df = "Theo_Exp",
    #                   x_lab = "log10(Ratios)",
    #                   color_lab= "Pool",
    #                   fill_lab = "Ratios")
    #  
    #  
    #  #### Distribution of mean abundance of every sample with two different pools
    #  
    #  figure_with_mean_abundance <- final_imputed_data %>%
    #    select(starts_with("mean_") | contains(c("Pool", "common_col_for_merging"))) %>%
    #    pivot_longer(cols = starts_with("mean"),
    #                 names_to = "Theo_Exp",
    #                 values_to = "Mean_abundance",
    #                 values_drop_na = T) %>%
    #    filter_at(vars(Pool), all_vars(!is.na(.)))
    #  
    #  
    #  p2 <- gg_density(data_set = figure_with_mean_abundance, 
    #                   x_df = figure_with_mean_abundance$Mean_abundance,
    #                   fill_df = figure_with_mean_abundance$Theo_Exp,
    #                   color_df = figure_with_mean_abundance$Pool,
    #                   header="Distribution of mean abundance of every sample with two different pools",
    #                   facet_df = "Theo_Exp",
    #                   x_lab = "log10(Ratios)",
    #                   color_lab= "Pool",
    #                   fill_lab = "Ratios")
    #  
    #  # final_imputed_normalized_data %>%
    #  #   select(starts_with("mean_") | contains(c("Pool", "common_col_for_merging"))) %>%  
    #  #   filter_at(vars(Pool), all_vars(!is.na(.))) %>%
    #  
    #  #### MANUAL PLOTTING   
    #  #
    #  #   ggplot(aes(x =log10(mean_abundances_A1) , y = log10(mean_abundances_A5), color=Pool)) +
    #  #   geom_point()+
    #  #   #facet_wrap(vars(df_for_figure$common_col_for_merging))  + 
    #  #   geom_smooth(formula = y ~ x,method = "loess", colour = "green", fill = "green") +
    #  #   theme_minimal() +
    #  #   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
    #  #         axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
    #  #         plot.title = element_text(size=20),
    #  #         legend.title=element_text(size=15),
    #  #         axis.text=element_text(size=15),
    #  #         axis.title=element_text(size=15)
    #  #   )
    #  #### PLOTTING WITH FUNCTION
    #  #
    #  # gg_density_mean_abun <- function(data_set = final_imputed_normalized_data,
    #  #                                  x_df = final_imputed_normalized_data$mean_abundances_A1,
    #  #                                  y_df = final_imputed_normalized_data$mean_abundances_A5,
    #  #                                  color_df =final_imputed_normalized_data$Pool,
    #  #                                  #header,
    #  #                                  facet_df= "common_col_for_merging")
    #  
    #  ### BOX-PLOT: Experimental Quantity Ratio of Synthetic Peptides  
    #  
    #  p3 <- gg_boxplt_exp_ratio(data_set = df_for_figure, 
    #                            x_df = df_for_figure$Theo_Exp,
    #                            y_df = df_for_figure$Ratios,
    #                            fill_df = df_for_figure$Pool,
    #                            header="Experimental Quantity Ratio of Synthetic Peptides",
    #                            x_lab="Sample Names",
    #                            y_lab="Abundance Ratios",
    #                            fill_lab = "Pool")
    #  
    #  ### HALF-BOX-PLOT & HALF-SCATTER-PLOT: Experimental Quantity Ratio of Synthetic Peptides  
    #  library(gghalves)
    #  
    #  p4 <- gg_half_boxplt_exp_ratio(data_set = df_for_figure, 
    #                                 x_df = df_for_figure$Theo_Exp,
    #                                 y_df = df_for_figure$Ratios,
    #                                 fill_df = df_for_figure$Pool,
    #                                 header="Experimental Quantity Ratio of Synthetic Peptides",
    #                                 x_lab="Sample Names",
    #                                 y_lab="Abundance Ratios",
    #                                 fill_lab = "Pool")
    #  
    #  
    #  
    #  ### VIOLIN-PLOT: Experimental Quantity Ratio of Synthetic Peptides   
    #  
    #  ### TODO: fix y scaling without trimming 
    #  p5 <- gg_violin_exp_ratio(data_set = df_for_figure, 
    #                            x_df = df_for_figure$Theo_Exp,
    #                            y_df = df_for_figure$Ratios,
    #                            fill_df = df_for_figure$Pool,
    #                            header="Experimental Quantity Ratio of Synthetic Peptides",
    #                            x_lab="Sample Names",
    #                            y_lab="Abundance Ratios",
    #                            fill_lab = "Pool",
    #                            trim=TRUE)
    #  
    #  ### SAVE ALL PLOT AUTOMATICALLY (without giving a custom file name)
    #  sapply(1:5,function(x) ggsave(filename = paste0("p",x,".tiff"),
    #                                width = 50, height = 40, 
    #                                path = file_path,
    #                                units = "cm",
    #                                get(paste0("p",x)),
    #                                device = "tiff", #".svg"
    #  ))
    #  
    #  ### SAVE PLOTS SEPARATELY (with custom file name)
    #  # ggsave(filename = "Experimental_Quantity_Ratio_of_Synthetic_Peptides",
    #  #        path = "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/",
    #  #        plot = p5,
    #  #        device = "tiff")
    #  
    #  
    #  # library(limma)
    #  # sample_size <- 5
    #  # # Do t-test
    #  # # The code below does t-test for each row. Because of that, multiple test correction (like BH, Bonferoni)
    #  # 
    #  # ## 1 ## TIDY DATA FOR T-TEST
    #  # df_for_ttest<- final_imputed_normalized_data %>%
    #  #   select(starts_with("mean_") | contains(c("Pool", "common_col_for_merging"))) %>%  
    #  #   filter_at(vars(Pool), all_vars(!is.na(.)))
    #  # 
    #  # ## 2 ## FUNCTION FOR T-TEST ALL CONC. ACROSS A1
    #  # ttest_func <- function(A1,rest){
    #  #   t.test(A1,rest,alternative = "two.sided", var.equal = TRUE)}
    #  # 
    #  # p_values_for_all_ratio <- data.frame(1:dim(df_for_ttest)[1])
    #  # 
    #  # ## 3 ## AUTOMIZED T-TEST ALL CONC. ACROSS A1
    #  # for(k in 2:sample_size){
    #  #   
    #  #   assign(paste0("t_test_res_A1_to_A",k) , 
    #  #          lapply(df_for_ttest$mean_abundances_A1, ttest_func,
    #  #                 rest=df_for_ttest[,k]))
    #  #   
    #  #   assign(paste0("p_values_for_ratio",k),
    #  #          lapply(get(paste0("t_test_res_A1_to_A",k)), function (x) x[c('p.value')]))
    #  #   
    #  #   p_values_for_all_ratio <- cbind(p_values_for_all_ratio,
    #  #                                   as.data.frame(unlist(get(paste0("p_values_for_ratio",k)))))
    #  # }
    #  # 
    #  # 
    #  # 
    #  # #adjusted_p_values <- as.data.frame(p.adjust(p_values_A1_vs_A5, method = "BH", n = length(p_values_A1_vs_A5)))
    #  # adjusted_p_values <- lapply(p_values_for_all_ratio[,-1], function(x) p.adjust (x, method = "BH", n = length(x)))
    #  # 
    #  # 
    #  # 
    #  # 
    #  # 
    #  # #plot(correct_identifed_peps$mean_abundances_log10_A1,correct_identifed_peps$mean_abundances_log10_A5, pch = 16, col = "blue")
    #  # #abline(h = mean(na.omit(correct_identifed_peps$mean_abundances_log10_A1)) - mean(na.omit(correct_identifed_peps$mean_abundances_log10_A5)), col = "red")
    #  # 
    #  # boxplot(correct_identifed_peps$mean_abundances_log10_A1,correct_identifed_peps$mean_abundances_log10_A5)
    #  # ggplot(correct_identifed_peps, aes(x = mean_abundances_log10_A1, y = mean_abundances_log10_A5)) +
    #  #   geom_point(aes(color = "red", size = 5))# +
    #  # #scale_color_discrete(name = "") +
    #  # #scale_size_discrete(name = "Size")
    #  # # 
    #  # # Take -log10() of results
    #  # ggplot(log_10_filtered_abundances_rowMeans_A1_A5,aes(x=log_10_filtered_abundances_rowMeans_A1_A5$mean_abundances_log10_A1,
    #  #                                                      y =log_10_filtered_abundances_rowMeans_A1_A5$mean_abundances_log10_A5,
    #  # )) +
    #  #   geom_point()+
    #  #   #facet_wrap(vars(df_for_figure$common_col_for_merging))  + 
    #  #   #geom_smooth(method = "lm", colour = "green", fill = "green") +
    #  #   theme_light()
    #  # 
    #  # library(gginference)
    #  # ggttest(t.test(na.omit(correct_identifed_peps$mean_abundances_log10_A1),na.omit(correct_identifed_peps$mean_abundances_log10_A5), alternative = "two.sided", var.equal = TRUE))
    #  # 
    #  # 
    #  # 
    #  # df_for_figure_exp <- correct_identifed_peps[,c(11,74:77)]
    #  # df_for_figure <- correct_identifed_peps[,c(11:15,74:77)]
    #  # melt_df_for_figure_exp <- melt(df_for_figure_exp)
    #  # 
    #  # df_for_figure_theo <- correct_identifed_peps[,c(11:15)]
    #  # melt_df_for_figure_theo <- melt(df_for_figure_theo)
    #  # 
    #  # p1 <- ggplot(melt_df_for_figure_theo,aes(x =melt_df_for_figure_theo$variable , y =log2(melt_df_for_figure_theo$value),fill = melt_df_for_figure_theo$Pool)  ) +
    #  #   geom_boxplot() +
    #  #   theme_light() + #scale_y_continuous(limits = c(-60, 60),breaks = seq(-60, 60, by = 20)) +
    #  #   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
    #  #         axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
    #  #         plot.title = element_text(size=20),
    #  #         legend.title=element_text(size=15),
    #  #         axis.text=element_text(size=15),
    #  #         axis.title=element_text(size=15)
    #  #   ) +   stat_boxplot(geom = "errorbar") + ggtitle("Theoretical Quantity Ratio of Synthetic Peptides") +
    #  #   scale_x_discrete(labels=c("A1/A2","A1/A2","A1/A4","A1/A5")) +
    #  #   labs(x="Sample Names",y="Abundance Ratios", color="Pool Names")+
    #  #   scale_fill_brewer(palette="Set1")
    #  # 
    #  # p2 <-ggplot(melt_df_for_figure_exp,aes(x =melt_df_for_figure_exp$variable , y =log2(melt_df_for_figure_exp$value),fill = melt_df_for_figure_exp$Pool)  ) +
    #  #   geom_boxplot() +
    #  #   theme_light() +
    #  #   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
    #  #         axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
    #  #         plot.title = element_text(size=20),
    #  #         legend.title=element_text(size=15),
    #  #         axis.text=element_text(size=15),
    #  #         axis.title=element_text(size=15)
    #  #   ) +   stat_boxplot(geom = "errorbar") + ggtitle("Experimental  Quantity Ratio of Synthetic Peptides") +
    #  #   scale_x_discrete(labels=c("A1/A2","A1/A2","A1/A4","A1/A5")) +
    #  #   labs(x="Sample Names",y="Abundance Ratios", color="Pool Names") +
    #  #   scale_fill_brewer(palette="Set1")
    #  # 
    #  # p3 <-ggplot(melt_df_for_figure_exp,aes(x =melt_df_for_figure_exp$variable , y =melt_df_for_figure_exp$value,color = melt_df_for_figure_exp$Pool)  ) +
    #  #   geom_point() +
    #  #   theme_light() +
    #  #   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
    #  #         axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
    #  #         plot.title = element_text(size=20),
    #  #         legend.title=element_text(size=15),
    #  #         axis.text=element_text(size=15),
    #  #         axis.title=element_text(size=15)
    #  #   ) +   ggtitle("Experimental Quantity Ratio of Synthetic Peptides") +
    #  #   scale_x_discrete(labels=c("A1/A2","A1/A3","A1/A4","A1/A5")) +
    #  #   labs(x="Sample Names",y="Abundance Ratios", color="Pool Names") +
    #  #   scale_color_brewer(palette="Set1")
    #  # 
    #  # p4 <-ggplot(melt_df_for_figure_theo,aes(x =melt_df_for_figure_theo$variable , y =melt_df_for_figure_theo$value,color = melt_df_for_figure_theo$Pool)  ) +
    #  #   geom_point() +
    #  #   theme_light() +
    #  #   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
    #  #         axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
    #  #         plot.title = element_text(size=20),
    #  #         legend.title=element_text(size=15),
    #  #         axis.text=element_text(size=15),
    #  #         axis.title=element_text(size=15),
    #  #   ) +   ggtitle("Theoretical Quantity Ratio of Synthetic Peptides") +
    #  #   scale_x_discrete(labels=c("A1/A2","A1/A3","A1/A4","A1/A5")) +
    #  #   labs(x="Sample Names",y="Abundance Ratios", color="Pool Names") +
    #  #   scale_color_brewer(palette="Set1")
    #  # 
    #  # 
    #  # library(ggpubr)
    #  # ggarrange(p1, p2, common.legend = TRUE, legend="right")
    #  # ggarrange(p4, p3, common.legend = TRUE, legend="right")
    #  # 
    #  # 
    #  # iso_exp <- df_for_figure_exp %>% filter(grepl("ISO-REF",Pool))
    #  # melt_iso_theo <- melt(iso_exp)
    #  # ggplot(melt_iso_theo,aes(x =melt_iso_theo$variable , y =log2(melt_iso_theo$value),color = melt_iso_theo$Pool)  ) +
    #  #   geom_violin()  + #geom_point() +
    #  #   theme_light() +
    #  #   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
    #  #         axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
    #  #         plot.title = element_text(size=20),
    #  #         legend.title=element_text(size=15),
    #  #         axis.text=element_text(size=15),
    #  #         axis.title=element_text(size=15)
    #  #   ) +   ggtitle("Experimental Quantity Ratio of Synthetic Peptides") +
    #  #   scale_x_discrete(labels=c("A1/A2","A1/A3","A1/A4","A1/A5")) +
    #  #   labs(x="Sample Names",y="Abundance Ratios", color="Pool Names") +
    #  #   scale_color_brewer(palette="Set1")
    #  # 
    #  # 
    #  # other_exp <- df_for_figure_exp %>% filter(grepl("Others",Pool))
    #  # melt_other_exp <- melt(other_exp)
    #  # ggplot(melt_other_exp,aes(x =melt_other_exp$variable , y =log2(melt_other_exp$value),color = melt_other_exp$Pool)  ) +
    #  #   geom_violin() + #geom_point() +
    #  #   theme_light() +
    #  #   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
    #  #         axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
    #  #         plot.title = element_text(size=20),
    #  #         legend.title=element_text(size=15),
    #  #         axis.text=element_text(size=15),
    #  #         axis.title=element_text(size=15)
    #  #   ) +   ggtitle("Experimental Quantity Ratio of Synthetic Peptides") +
    #  #   scale_x_discrete(labels=c("A1/A2","A1/A3","A1/A4","A1/A5")) +
    #  #   labs(x="Sample Names",y="Abundance Ratios", color="Pool Names") +
    #  #   scale_color_brewer(palette="Set1")
    #  
    #  
    
  
  
  
}

