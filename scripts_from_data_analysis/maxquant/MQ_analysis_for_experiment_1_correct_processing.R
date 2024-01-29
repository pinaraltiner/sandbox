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
library(gtools)
###############################################
## MAXQUANT DATA ANALYSIS

### Source code was taken from here: https://rdrr.io/github/singjc/mstools/src/R/getModificationPosition.R
## The code was modified based on what I want and based on software input tyoe (DIANN-and MaxQuant)
#source("D:/dev/Desktop_copy/PHD/data_analysis/scripts/getModificationPositionMQ_func.R")
#source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/get_modification_func/getModificationPositionMQ_func_edit.R")
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/get_modification_func/getModificationPosition_general_change_condition_current_mod_sequence_MQ_Spectronaut.R")
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/ggplot/ggplot_functions.R")
#all_files <- list.files(paste0(file_path,all_dirs[i],"/"),pattern = ".xlsx")

#comb_result <- NULL
#comb_ecoli <- NULL
final_MQ_pep_quant_analysis_exp1 <- function(file_path,
                                             curr_dir,
                                             #file_name,
                                             theo_file_path,
                                             theo_file_name,
                                             sheet_theo_name,
                                             selected_species,
                                             mapping_file,
                                             exp_id,
                                             #exp_design,
                                             background_species,
                                             acquisiton_type,
                                             software_name){
  
  map_df <- read.table(file = mapping_file,sep = "\t",header = T)
  
  assign(paste0("tmp"), read.table(file=paste0(file_path,curr_dir,"/evidence.txt"),sep = "\t",header = T))
    
  if (any(curr_dir %in% c("exp1_withFAIMS_withEcoli_OXPAL230421"))) {
    tmp <- tmp %>%
      mutate(Experiment = ifelse(Raw.file == "OXPAL230421_08_-45","E1-M2-coli-inj1",Experiment)) 
    }else{}
  
  if(any(curr_dir %in% c("exp1_withFAIMS_withEcoli_OXPAL230421","exp1_withFAIMS_woEcoli"))){
    
    tmp <- tmp %>% separate(Raw.file,into=c("Raw.file","CV_type"),sep = '_-') %>%
      mutate(Experiment_with_CV=paste(Experiment,CV_type,sep = "_"))
  }
  else{}
  
  
  phospho_tmp <- tmp %>% 
    filter(grepl("Phospho",Modified.sequence)) %>%
    filter(grepl(selected_species,Proteins) & !grepl("CON__",Proteins)) %>%
    #drop_na(Modified.sequence) %>%
    mutate(curr_dir) %>%
    mutate(acq_type=acquisiton_type) %>%
    mutate(soft_name=software_name) 
  
  if (background_species %in% c("ECOLI","Escherichia coli")){
    
    ecoli_tmp <- tmp %>%
      filter(!grepl(selected_species,Proteins) & !grepl("CON_", Proteins)) %>%
      mutate(curr_dir) %>%
      mutate(acq_type=acquisiton_type) %>%
      mutate(soft_name=software_name) #%>%
      #group_by(Sequence,Experiment) %>%
      #slice(which.max(Intensity)) %>%
      #ungroup() %>%
      #distinct(Sequence,.keep_all = TRUE) %>%
      #separate(Experiment,into =c("Exp_id","Sample_id","tmp","inj_id"),sep = "-",remove = F) %>%
      #mutate(new_col=paste(Exp_id,Sample_id,sep = "_"))
    ## REMOVE DUPLICATE SEQUENCES in each sample
    
    comb_ecoli_dist <- ecoli_tmp %>% 
      rename(raw_file=Raw.file) %>%
      full_join(map_df,by="raw_file") %>%
      rename(pool_id_map_df=pool_id) %>%
      #separate(pool_id,into = c("expid","sample_id","Ecoli","inj"),sep = "_",remove = F) %>%
      #mutate(new_col=paste(expid,sample_id,Ecoli,sep = "_")) %>%
      group_by(Sequence,raw_file) %>% 
      distinct(Sequence,.keep_all = TRUE)%>% 
      ungroup() 
    
    plt1 <- gg_barplt_id_pep_count(data_set = comb_ecoli_dist,
                                   x_df = comb_ecoli_dist$raw_file,
                                   fill_df = comb_ecoli_dist$pool_id_map_df,
                                   ymax = 20000,
                                   header = paste("Total number of identified synthetic phospho-sequences across each pool",sep=" "),
                                   caption_lab = "Wrong Sequences were removed for the futher analysis.",
                                   x_lab = "Sample id",
                                   fill_lab =  "Sample id",
                                   y_lab = "Number of identified Sequence",
                                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name),
                                   size_num=7) + 
      scale_fill_brewer(palette = "Dark2") + theme(axis.text.x = element_text(size =7,angle = 90))
    
    write.table(comb_ecoli_dist,file=paste0(file_path,curr_dir,"/",
                                            "exp",
                                            exp_id,acquisiton_type,
                                            "merge_identified_ecoli_seq.tsv"),
                sep = "\t",col.names = T,row.names = F)
    }else{}
        
  #### PHOSPHO-POS EXTRACTION ####

    df2 <- apply(X = as.data.frame(phospho_tmp[,"Modified.sequence"]),1,function(x){getModificationPosition_general(mod_seq = x,software_name = "MQ_214")})
    
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
    
    final_results_with_common_col <- phospho_tmp %>%
      mutate(result_with_common_col) %>%
      mutate(extracted_values = sapply(str_extract_all(Phospho..STY..Probabilities,
                                                       "\\(\\d+(\\.\\d+)?\\)"), function(x) {
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
      )%>% relocate(extracted_values,.after = Sequence)
        #drop_na(extracted_values) %>%
    ## READ THEO LIST
    pep_list_w_theo <- read.xlsx(paste0(theo_file_path, theo_file_name), sheet = sheet_theo_name)
    pep_list_w_theo_quant <- pep_list_w_theo[,-1]
    common_col_theo_quant <- as.data.frame(paste(pep_list_w_theo_quant$Phosphopeptide.sequence,
                                                 pep_list_w_theo_quant$modified.position.in.peptide, sep = "_"))
    
    ## ADD COMMON COLUMN TO MERGE WITH EXP. DATA
    colnames(common_col_theo_quant) <- "pep_with_pos"
    pep_list_w_theo_quant_new <- cbind(common_col_theo_quant,pep_list_w_theo_quant)
    
    pep_list_w_theo_pep_comp <- pep_list_w_theo_quant_new %>% 
      select(!Phosphopeptide.sequence) %>%
      #rename(Sequence = Phosphopeptide.sequence) %>% 
      rename(pool_id_theo_list=pool_id) 
    
    pep_list_w_theo_seq_comp <- pep_list_w_theo_quant_new %>% 
      select(Phosphopeptide.sequence,pool_id) %>%
      #select(!pep_with_pos) %>%
      rename(Sequence = Phosphopeptide.sequence) %>% 
      mutate(pool_id_theo_list=pool_id) 
    
    # Extracted values
    comb_result_pos <- final_results_with_common_col %>% 
      select(Sequence, pep_with_pos,extracted_values,
             Phospho..STY..Probabilities,
             Experiment,Raw.file,
             `curr_dir`, MS.MS.scan.numbers,
             Retention.time,Intensity, #Marked.as
             Proteins) %>%
      rename(raw_file=Raw.file) %>%
      full_join(map_df,by="raw_file") #%>%
      #rename(pool_id_map_df=pool_id) 
      
    
    #### TOTAL NUM. OF PHOSPHO-SEQUENCES ####
    ### Correct mapping was done using "map_df".
    comb_result_seq <- comb_result_pos %>% 
      full_join(pep_list_w_theo_seq_comp,by=c("Sequence","pool_id")) %>%
      select(!c(pool_id,sample_name)) %>%
      ## RE-JOINING TO COUNT CORRECT, INCORRECT and MISSING SEQ. 
      full_join(map_df,by="raw_file") %>%
      rename(pool_id_map_df=pool_id) %>%
      mutate(Pool_for_seq_merge=ifelse(pool_id_theo_list==pool_id_map_df,"Correct Seq.","Wrong Seq.")) %>% #within theo list
      #mutate_at("Pool_for_seq_merge", ~replace_na(.,"Unexpected Seq.")) %>%
      mutate(Pool_for_seq_merge= ifelse(is.na(pool_id_map_df),"Missing",Pool_for_seq_merge)) %>%
      mutate(Pool_for_seq_merge=ifelse(is.na(pool_id_theo_list),"Wrong Seq.",Pool_for_seq_merge)) %>% # out of theo. list
      mutate(acq_type=acquisiton_type) %>%
      mutate(soft_name=software_name)
    
    ### Duplicate sequences were removed.
    comb_result_dist <- comb_result_seq %>%
      relocate(pool_id_map_df,pool_id_theo_list,Pool_for_seq_merge,.after = Sequence) %>% #
      group_by(raw_file,Sequence) %>%
      distinct(Sequence,.keep_all = T) %>%
      ungroup() 
    
    #### WRITE THE OBJECT AS TSV ####
    write.table(comb_result_seq,file=paste0(file_path,curr_dir,"/",
                                            "exp",
                                            exp_id,acquisiton_type,
                                            "merge_theo_list_with_id_phospho_seq.tsv"),
                sep = "\t",col.names = T,row.names = F)
    
    #### VISUALIZATION OF TOTAL NUM. OF PHOSPHO-SEQ ####
    plt3 <- gg_barplt_id_pep_count_stack(data_set = comb_result_dist,
                                   x_df = comb_result_dist$sample_name,
                                   fill_df = comb_result_dist$Pool_for_seq_merge,
                                   ymax = 250,
                                   size_num=10,
                                   header = paste("Total number of identified synthetic phospho-sequences across each pool",sep=" "),
                                   caption_lab = "Wrong Sequences were removed for the futher analysis.",
                                   x_lab = "Sample id",
                                   fill_lab =  "Sample id",
                                   y_lab = "Number of identified Sequence",
                                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) +
      
      scale_fill_brewer(palette = "Dark2") + theme(axis.text.x = element_text(angle = 90))
    
    #### EXTRACTION OF CORRECT SEQ. ####
    comb_result_seq_dist_cor <- comb_result_dist %>%
      filter(!grepl("Wrong",Pool_for_seq_merge) & !grepl("Missing",Pool_for_seq_merge))
    
    #### VISUALIZATION OF TOTAL NUM. OF CORRECTLY IDENTIFIED PHOSPHO-SEQ ####
    
    plt2 <- gg_barplt_id_pep_count_stack(data_set =comb_result_seq_dist_cor,
                                   x_df =comb_result_seq_dist_cor$sample_name,
                                   fill_df = comb_result_seq_dist_cor$Pool_for_seq_merge,
                                   ymax = 250,
                                   size_num=10,
                                   header = paste("Total number of identified synthetic phospho-sequences across each pool",sep=" "),
                                   caption_lab = "Duplicates were removed for the futher analysis.",
                                   x_lab = "Sample id",
                                   fill_lab =  "Sample id",
                                   y_lab = "Number of identified Sequence",
                                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) +
      scale_fill_brewer(palette = "Paired") + theme(axis.text.x = element_text(angle = 90))
    
    
    #### TOTAL NUM. OF PHOSPHO-PEPTIDES ####
    ### Only correct sequences were kept.
    ### Correct mapping was done using "pool_id_map_df" and "pep_with_pos".
    ### Duplicate peptides were removed.
    
    comb_result_pep <- comb_result_seq %>% 
      filter(!grepl("Wrong",Pool_for_seq_merge) & !grepl("Missing",Pool_for_seq_merge)) %>%
      
      select(!c(Pool_for_seq_merge,pool_id_theo_list,pool_id_map_df,sample_name)) %>%
      full_join(map_df,by="raw_file") %>%
      rename(pool_id_map_df=pool_id) %>%
      full_join(pep_list_w_theo_pep_comp,by="pep_with_pos") %>% ## IF FULL_JOIN IS USED,
      group_by(pep_with_pos,raw_file) %>%
      distinct(pep_with_pos,.keep_all = TRUE)%>%
      ungroup() %>%
      mutate(Pool_for_pep_merge=ifelse(pool_id_theo_list==pool_id_map_df,"Correct","Wrong Loc.within theo list")) %>%
      mutate(Pool_for_pep_merge= ifelse(is.na(pool_id_map_df),"Missing",Pool_for_pep_merge)) %>%
      mutate(Pool_for_pep_merge=ifelse(is.na(pool_id_theo_list),"Wrong Loc. out of theo. list",Pool_for_pep_merge)) %>%
      relocate(pool_id_map_df,pool_id_theo_list,Pool_for_pep_merge,.after = Sequence) %>% 
      mutate(acq_type=acquisiton_type) %>%
      mutate(soft_name=software_name) 
      # 
      # select(!c(Pool_for_seq_merge, pool_id_theo_list)) %>%
      # full_join(pep_list_w_theo_quant_new,by="pep_with_pos") %>% ## IF FULL_JOIN IS USED,
      # group_by(pep_with_pos,raw_file) %>%
      # distinct(pep_with_pos,.keep_all = TRUE)%>%
      # ungroup() %>%
      # mutate(Pool_for_pep_merge=ifelse(pool_id_theo_list==pool_id_map_df,"Correct","Wrong Loc.within theo list")) %>%
      # mutate(Pool_for_pep_merge= ifelse(is.na(pool_id_map_df),"Missing",Pool_for_pep_merge)) %>%
      # mutate(Pool_for_pep_merge=ifelse(is.na(pool_id_theo_list),"Wrong Loc. out of theo. list",Pool_for_pep_merge)) %>%
      # mutate(acq_type=acquisiton_type) %>%
      # mutate(soft_name=software_name) 
      # 
    
    write.table(comb_result_pep ,file=paste0(file_path,curr_dir,"/merge_theo_list_id_STY_max_int.tsv"),
                sep = "\t",col.names = T,row.names = F)
    
    
    #### VISUALIZATION OF TOTAL NUM. OF IDENTIFIED & LOCALIZED PHOSPHO-PEP ####
    plt4 <- gg_barplt_id_pep_count_stack(data_set =comb_result_pep,
                                   x_df =comb_result_pep$sample_name,
                                   fill_df = comb_result_pep$Pool_for_pep_merge,
                                   ymax = 250,
                                   size_num=10,
                                   header = paste("Number of Identified phospho-sites for each pool \n","Experiment",
                                                  exp_id, acquisiton_type),
                                   caption_lab = "Duplicates were removed for the futher analysis. \n Assessment was done by selecting only correct sequences",
                                   x_lab = "Sample id",
                                   fill_lab =  "Localization Accuracy",
                                   y_lab = "Number of Localized Phospho-peptides",
                                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) +
      scale_fill_brewer(palette = "Dark2") + theme(axis.text.x = element_text(angle = 90))
    
    comb_result_pep_filt <- comb_result_pep %>% filter(extracted_values >= 0.75)
    
    plt10 <- gg_barplt_id_pep_count_stack(data_set =comb_result_pep_filt,
                                         x_df =comb_result_pep_filt$sample_name,
                                         fill_df = comb_result_pep_filt$Pool_for_pep_merge,
                                         ymax = 250,
                                         size_num=10,
                                         header = paste("Number of Identified phospho-sites for each pool \n after applying localization threshold \n","Experiment",
                                                        exp_id, acquisiton_type),
                                         caption_lab = "Duplicates were removed for the futher analysis. \n Assessment was done by selecting only correct sequences",
                                         x_lab = "Sample id",
                                         fill_lab =  "Localization Accuracy",
                                         y_lab = "Number of Localized Phospho-peptides",
                                         subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) +
      scale_fill_brewer(palette = "Dark2") + theme(axis.text.x = element_text(angle = 90))
    
    
    
    #### EXTRACTION OF CORRECT PEP. ####
    comb_result_pep_cor <- comb_result_pep %>%
      filter(grepl("Correct",Pool_for_pep_merge))
    
    #### VISUALIZATION OF TOTAL NUM. OF CORRECTLY IDENTIFIED & LOCALIZED PHOSPHO-PEP ####
    plt5 <- gg_barplt_id_pep_count_stack(data_set =comb_result_pep_cor,
                                   x_df =comb_result_pep_cor$sample_name,
                                   fill_df = comb_result_pep_cor$Pool_for_pep_merge,
                                   ymax = 250,
                                   size_num=10,
                                   header = paste("Distribution of Correctly Identified phospho-sites for each pool \n","Experiment",
                                                  exp_id, acquisiton_type),
                                   caption_lab = "Duplicates were removed for the futher analysis. \n Assessment was done by selecting only correct sequences",
                                   x_lab = "Sample id",
                                   fill_lab =  "Localization Accuracy",
                                   y_lab = "Number of Localized Phospho-peptides",
                                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) +
      scale_fill_brewer(palette = "Dark2") + theme(axis.text.x = element_text(angle = 90))
    
    if(background_species == "ECOLI"){
      
      comb_result_pep_filtered <- comb_result_pep %>% 
        separate(sample_name,into = c("exp","samp","coli","inj"),sep = "_",remove = F) %>%
        filter(!grepl("Missing",Pool_for_pep_merge)) %>%
        mutate(samp_name_wo_inj=paste(exp,samp,coli,sep = "-")) %>%
        group_by(samp_name_wo_inj,pep_with_pos) %>%
        slice(which.max(extracted_values)) %>%
        ungroup()
      
    }else if (background_species ==""){
      
      comb_result_pep_filtered <- comb_result_pep %>% 
        separate(sample_name,into = c("exp","samp","inj"),sep = "_",remove = F) %>%
        filter(!grepl("Missing",Pool_for_pep_merge)) %>%
        mutate(samp_name_wo_inj=paste(exp,samp,sep = "-")) %>%
        group_by(samp_name_wo_inj,pep_with_pos) %>%
        slice(which.max(extracted_values)) %>%
        ungroup()
      
    }else{}
    
    #### LOCALIZATION ACCURACY ASSESSMENT (ROC LIKE PLOT GENERATION) #### 
    ### Column selection
    comb_result_pep_rmv_miss <- comb_result_pep_filtered %>% 
      select(pep_with_pos,Pool_for_pep_merge,extracted_values,pool_id_map_df)
    
    ### Custom localization threshold determination and
    ### Counting total number of correct and wrong localization at a given threshold
    threshold <- seq(0,1,length=100)
    final_df <- NULL
    
    for (i in 1:length(threshold)){
      
      df <- comb_result_pep_rmv_miss %>% subset(extracted_values > threshold[i]) %>% 
        count(Pool_for_pep_merge) %>% 
        mutate(threshold_val = threshold[i])
      
      final_df <-bind_rows(final_df,df)
      rm(df)
    }
    
    write.table(final_df, file=paste0(file_path,curr_dir,"/","Exp",exp_id,"_num_STY_w_scores.tsv"),sep = "\t",col.names = T,row.names = F)
    
    #wrong_df <- final_df %>% filter(grepl("Wrong",Pool_for_pep_merge))
    
    ### Filtering only the correct ones
    correct_df <- final_df %>% filter(grepl("Correct",Pool_for_pep_merge))
    
    ### Visualization
    plt6 <- ggplot(correct_df,aes(y=n, x=threshold_val)) + geom_line(size=2) +
      scale_fill_brewer(palette = "Dark2") +
      labs(title=paste("Number of Correctly Identified phospho-sites at various thresholds \n","Experiment",
                       exp_id, acquisiton_type),x="Sample id",y="Number of identified Sequence",
           subtitle = acquisiton_type) + 
      theme_minimal() +
      theme(legend.text = element_text(size=30), 
            axis.title.x = element_text(size=30),
            axis.title.y = element_text(size=30),
            plot.title = element_text(size=35),
            plot.subtitle = element_text(size = 25),
            plot.caption = element_text(size = 25),
            legend.title=element_text(size=30),
            axis.text.x = element_text(size=20,angle = 90),
            axis.text.y = element_text(size = 30),
            axis.title=element_text(size=30)) + 
      #scale_x_continuous(limits = c(0, 100),breaks = seq(from = 0, to = 100, by = 10)) +  # Set the ticks for the x-axis
      scale_y_continuous(limits = c(0, 180),breaks = seq(from = 0, to = 180, by = 10))
    
    ##################  ADJACENT SEQ ASSESSMENT ####################
    
    
    ### All possible combination of S,T,and Y were prepared.
    elements <- c("S", "T", "Y")
    ## 1 was removed because it is the same as "element object"
    for(i in 2:5){ ## 5 might be selected as a parameter of the main function
      assign(paste0("perm",i),
             as.data.frame(permutations(n = length(elements), r = i, v = elements, repeats.allowed = TRUE)))
      
    }
    
    perm2[,"comb"] <- as.data.frame(paste0(perm2$V1,perm2$V2))
    
    perm3[,"comb"] <- as.data.frame(paste0(perm3$V1,perm3$V2,perm3$V3))
    
    perm4[,"comb"] <- as.data.frame(paste0(perm4$V1,perm4$V2,perm4$V3,perm4$V4))
    
    perm5[,"comb"] <- as.data.frame(paste0(perm5$V1,perm5$V2,perm5$V3,perm5$V4,perm5$V5))
    
    ### Combination of STY possibilities 
    possibilites <-perm2 %>% bind_rows(perm3,perm4,perm5) %>% select(comb)
    possibilites <- possibilites$comb
    
    ## The function find matches to extract adjacent a.a 
    extract_pos_adj <- function(query,patt_text){
      tmp <- paste0(unlist(str_match_all(query,pattern=patt_text)))
      result <- tmp[!is.na(tmp)]
      return(result)
    }
    
    ## The function find matches to extract adjacent a.a 
    extract_pos_non_adj <- function(query,patt_text){
      tmp <- paste0(unlist(str_match_all(query, patt_text)))
      result <- tmp[!is.na(tmp)]
      return(result)
    }
    
    
    # comb_pep_adj_dist_df_neg <- comb_result_pep %>%
    #   select(Sequence) %>%
    #   mutate(
    #     Contains_STY = ifelse(sapply(Sequence, function(x) str_count(x, pattern = "(S(?=(S|T|Y)(S|T|Y))|T(?=(S|T|Y)(S|T|Y))|Y(?=(S|T|Y)(S|T|Y)))")) == 2,
    #                           sapply(Sequence, function(y) extract_pos(query = y,patt_text = perm4$comb)),""))
    # 
    # comb_pep_adj_dist_df_neg <- comb_result_pep %>%
    #   filter((sapply(Sequence, function(query) grepl(query, pattern = possibilites)) == TRUE)) %>%
    #   mutate(rowSums(sapply(Sequence,function(query) str_count(query,possibilites)))) %>%
    #   relocate(.after = Sequence)
    # # sum(sapply(Sequence, function(query) str_count(query, possibilites))),
    
    ### This regex was designed to check every S,T and Y in each sequence
    big_pattern ="(S(?=(S|T|Y)(S|T|Y))|T(?=(S|T|Y)(S|T|Y))|Y(?=(S|T|Y)(S|T|Y)))"
    
    ### ADJACENT COUNT FROM 3 to 5
    for(i in 1:3){
      assign(paste0("comb_result_pep_adj_",i+2),
             comb_result_pep_filtered %>% #comb_result_pep %>%
               #select(Sequence) %>%
               filter(sapply(Sequence, function(x) str_count(x, pattern = big_pattern)) == i) %>%
               mutate(STY_adj = sapply(Sequence, function(y) extract_pos_adj(query = y,patt_text = get(paste0("perm",i+2))$comb))) %>%
               mutate(STY_len=nchar(STY_adj)) %>%
               mutate(is_adj= "adjacent"))
    }
    
    ### ADJACENT COUNT FROM 2
    comb_result_pep_adj_2 <- comb_result_pep_filtered %>%#comb_result_pep
      #select(Sequence) %>%
      #filter(!(sapply(Sequence, function(x) str_count(x, pattern = big_pattern)) == 1) & 
      #!(sapply(Sequence, function(x) str_count(x, pattern = big_pattern)) == 2)) %>%
      
      filter(sapply(Sequence, function(x) str_count(x, pattern = big_pattern)) < 1 & 
               sapply(Sequence, function(x) str_count(x, pattern = "(S(?=S|T|Y)|T(?=S|T|Y)|Y(?=S|T))+")) == 1) %>%
      
      #filter(sapply(Sequence, function(x) str_count(x, pattern = "(S(?=S|T|Y)|T(?=S|T|Y)|Y(?=S|T))+")) == 2) %>%
      ### LIST MIGHT BE ADDED TO COLLECT ALL PAIRS 
      mutate(STY_adj = sapply(Sequence, function(x) extract_pos_adj(query = x,patt_text = perm2$comb))) %>%
      mutate(STY_len=nchar(STY_adj)) %>%
      mutate(is_adj= "adjacent")
    
    ### NON-ADJACENT COUNT
    comb_result_pep_non_adj <- comb_result_pep_filtered %>% #comb_result_pep
      filter(sapply(Sequence, function(x) str_count(x, pattern = big_pattern)) < 1 & 
               sapply(Sequence, 
                      function(x) str_count(x,
                                            pattern = "(S(?=S|T|Y)|T(?=S|T|Y)|Y(?=S|T))+")) == 0) %>%
      rowwise() %>%
      mutate(STY_adj = list(sapply(Sequence, function(x) extract_pos_non_adj(query = x, patt_text = elements)))) %>%
      #mutate(STY_adj = paste0(sapply(Sequence, function(x) extract_pos(query = x,patt_text = sing_STY))))
      mutate(STY_len=length(STY_adj)) %>%
      mutate(is_adj= "non_adjacent")
    
    comb_result_pep_non_adj$STY_adj_new <- sapply(comb_result_pep_non_adj$STY_adj, paste, collapse = "")
    
    ### COMBINATION OF ALL ADJACENT ONES
    comb_result_pep_adj_all <- comb_result_pep_non_adj %>%
      select(!c(STY_adj)) %>%
      rename(STY_adj = STY_adj_new) %>%
      bind_rows(comb_result_pep_adj_5,
                comb_result_pep_adj_4,
                comb_result_pep_adj_3,
                comb_result_pep_adj_2) %>%
      select(Sequence,
             STY_adj,
             STY_len,
             Pool_for_pep_merge,
             extracted_values,
             pep_with_pos,
             pool_id_map_df,
             sample_name,
             raw_file,is_adj) 
    
    ### TOTAL COUNT OF ALL SEQUENCES WITHOUT CONSIDERING ADJ_COUNT
    plt7 <- comb_result_pep_adj_all %>% tibble() %>%
      filter(grepl("Correct",Pool_for_pep_merge)) %>%
      count(is_adj) %>%
      #mutate(n_new = ifelse(Pool_for_seq_merge!="Correct",(-1*n),n)) %>%
      ggplot(aes(x=n,y=n, fill=is_adj)) +
      geom_col() + #geom_text(aes(label=after_stat(count)),stat = "count", position=position_dodge(width=0.9), vjust=-0.25,size=10)+
      geom_text(aes(label = n), position = position_stack(vjust = 0.5),size=10) +
      labs(x = "Adjacent STY count", y = "Phospho Peptide Count",
           fill="Match between \n pool and raw file",
           title = "Comparison of accuracy of phospho-peptides with number of holding adjacent amino acids",
           caption = paste(software_name,"Experiment",exp_id,acquisiton_type,sep = " ")) +
      theme_minimal() +
      scale_fill_manual(values = c("#377EB8", "#4DAF4A"))+ 
      theme(legend.text = element_text(size=30), 
            axis.title.x = element_text(size=30),
            axis.title.y = element_text(size=30),
            plot.title = element_text(size=35),
            plot.subtitle = element_text(size = 25),
            plot.caption = element_text(size = 25),
            legend.title=element_text(size=30),
            axis.text.x = element_text(size=20),
            axis.text.y = element_blank(),
            #axis.text.y = element_text(size = 30),
            axis.title=element_text(size=30)) #+
    #scale_x_continuous(limits = c(0,7),breaks = seq(from =0, to = 7, by = 1))
    
    ### TOTAL COUNT OF ALL SEQUENCES WITH CONSIDERING ADJ_COUNT AND ACCURACY
    plt8 <- comb_result_pep_adj_all %>%
      mutate(pep_class=Pool_for_pep_merge) %>%
      mutate(pep_class=ifelse(Pool_for_pep_merge!= "Correct","Wrong",pep_class)) %>%
      #filter(grepl("Correct",Pool_for_pep_merge)) %>%
      rowwise() %>%
      #distinct(Sequence,.keep_all = TRUE) %>%
      #group_by(pep_with_pos,new_col) %>%
      #count(pep_with_pos)
      #mutate(ptm_score = ifelse(Pool_for_pep_merge == "Wrong Localization",(ptm_score*-1),ptm_score)) %>%
      group_by(is_adj,pep_class) %>%
      count(STY_len) %>% #Pool_for_pep_merge
      ungroup() %>%
      mutate(n_label=n) %>%
      #mutate(n=ifelse(Pool_for_pep_merge =="Wrong Localization", (n* (-1)),n)) %>%
      mutate(n=ifelse(is_adj =="non_adjacent", (n* (-1)),n)) %>%
      mutate(STY_len= ifelse(STY_len < 0, (-1*STY_len),STY_len)) %>%
      #mutate(abs_n=c("141","79","82","26\n")) %>%
      ggplot( aes(x= STY_len,y=n,fill=interaction(is_adj,pep_class))) +
      geom_col() + labs(x = "Count of STY amino acids", y = "Total count",
                        title = "Distributions of the number of unique phospho-peptides holding adjacent and non-adjacent STY residues",
                        caption = paste(software_name,"Experiment",exp_id,acquisiton_type,sep = " "),fill="Is peptide adjacent?") +
      theme_minimal() +
      scale_fill_brewer(palette = "Paired",labels=c("Adjacent & Correct", "Non adjacent & Correct", "Adjacent & Wrong","Non adjacent & Wrong")) +
      #scale_fill_discrete() +
      theme(legend.text = element_text(size=30), 
            axis.title.x = element_text(size=30),
            axis.title.y = element_text(size=30),
            plot.title = element_text(size=35),
            plot.subtitle = element_text(size = 25),
            plot.caption = element_text(size = 25),
            legend.title=element_text(size=30),
            axis.text.x = element_text(size=20),
            axis.text.y = element_blank(),
            #axis.text.y = element_text(size = 30),
            axis.title=element_text(size=30)) +
      geom_text(aes(label =n_label), position = position_stack(vjust = 0.5),size=15)
    
    ### DISTRIBUTION OF ALL SEQUENCES WITH CONSIDERING ADJ_COUNT AND ACCURACY 
    plt9 <- comb_result_pep_adj_all %>%
      mutate(pep_class=Pool_for_pep_merge) %>%
      mutate(pep_class=ifelse(Pool_for_pep_merge!= "Correct","Wrong",pep_class)) %>%
      subset(!(STY_len == 1 & is_adj == "non_adjacent")) %>%
      
      ggplot( aes(x= extracted_values,fill=interaction(is_adj))) +
      geom_bar(stat="count") + labs(x = "Count of STY amino acids", y = "Total count",
                                     title = "Comparison of having an adjacent a.a effect of localization accuracy",
                                     caption = paste(software_name,"Experiment",exp_id,acquisiton_type,sep = " "),fill="Is peptide adjacent?") +
      facet_wrap(~pep_class) +
      theme_minimal() +
      #scale_fill_brewer(palette = "Paired")+
      scale_fill_manual(values = c("#377EB8", "#4DAF4A"))+ 
      theme(legend.text = element_text(size=30), 
            axis.title.x = element_text(size=30),
            axis.title.y = element_text(size=30),
            strip.text = element_text(size=30),
            plot.title = element_text(size=35),
            plot.subtitle = element_text(size = 25),
            plot.caption = element_text(size = 25),
            legend.title=element_text(size=30),
            axis.text.x = element_text(size=20),
            axis.text.y = element_text(size = 30),
            axis.title=element_text(size=30)) 
    
    #################################################
    plot_obj <- ls(pattern="plt")
    plot_obj <- plot_obj[!is.na(plot_obj)]
    sapply(1:length(plot_obj),function(x) ggsave(filename = paste0("p",x,".png"),
                                                 width = 60, height = 45, 
                                                 path = paste0(file_path,curr_dir,"/outputs_after_mapping_change/"),
                                                 units = "cm",
                                                 get(plot_obj[x]),
                                                 device = "png", #".svg"
    ))
    
    
    
    
    
} 
    
  