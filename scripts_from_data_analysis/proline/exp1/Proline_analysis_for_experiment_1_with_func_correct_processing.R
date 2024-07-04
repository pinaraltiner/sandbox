#library(PhosR)
library(stringr)
library(dplyr)
library(data.table)
library(openxlsx)
library(tidyr)
library(ggplot2)
library(gtools)
###############################################
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/ggplot/ggplot_functions.R")
#source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/roc_curve/new_roc_curve_generation_with_custom_threshold.R")

final_proline_pep_quant_analysis_syn <- function(file_path,
                                                 curr_dir,
                                                 sheet_name,
                                                 theo_file_path,
                                                 theo_file_name,
                                                 sheet_theo_name,
                                                 background_species,
                                                 selected_species,
                                                 exp_id,
                                                 #exp_design,
                                                 #fdr_threshold,
                                                 acquisiton_type,
                                                 software_name,
                                                 #test_type,
                                                 #num_reps,
                                                 #actual_ratio,
                                                 #subtitle,
                                                 mapping_file){
  
  
  map_df <- read.table(file = mapping_file,sep = "\t",header = T)
  
  # sample_size <- length(map_df$pool_id) / num_reps
  # sample_names <- paste0("A",1:sample_size)
  # comparisons <- NULL
  # for (i in 1:sample_size){
  #   tmp <- paste0(sample_names[1], "/",sample_names[i])
  #   comparisons[i] <- tmp
  #   rm(tmp)
  # }
  # comparisons <- comparisons[-1]
  
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
  
  
  all_files <- list.files(paste0(file_path,curr_dir,"/"),pattern = ".xlsx")

  comb_result <- NULL
  comb_ecoli <- NULL
  
  for(i in 1:length(all_files)){
    
    assign(paste0("tmp"), read.xlsx(paste0(file_path,curr_dir,"/",all_files[i]),sheet ="Best PSM from protein sets" ))
    
    phospho_tmp <- tmp %>% 
      filter(grepl("Phospho",modifications)) %>%
      filter(grepl(selected_species,accession) & !grepl("CON__",accession)) %>%
      mutate(all_files[i]) %>%
      mutate(phospho_pos = proline_phospho_pos_extraction(modifications)) %>%
      mutate(pep_with_pos = paste(sequence,phospho_pos,sep = "_"))
      #distinct(pep_with_pos,.keep_all = TRUE) ## The reason of putting this line
                                              ## is for eliminating the combine 
                                              ## version of phospho. Such as;
                                              ## N-term phospho, oxidation phospho, etc.
    
    #assign(paste0(all_files[i]), phospho_tmp)
    assign(paste0("comb_result"),bind_rows(comb_result, phospho_tmp))
    
    if(background_species == "ECOLI"){
      ecoli_tmp <- tmp %>% 
        filter(grepl(background_species,accession) & !grepl("CON__",accession)) %>%
        mutate(all_files[i]) %>%
        mutate(acq_type=acquisiton_type) %>%
        mutate(soft_name=software_name)
        #distinct(sequence,.keep_all = TRUE)
      
      assign(paste0("comb_ecoli"), bind_rows(comb_ecoli,ecoli_tmp))
    }else{
      
    }
    
    rm(tmp,phospho_tmp,ecoli_tmp)
    
  }

  if(background_species == "ECOLI"){
    ## REMOVE DUPLICATE SEQUENCES in each sample
    comb_ecoli_dist <- comb_ecoli %>% 
      separate(spectrum_title,into = c(paste0("tmp",1:6),"spec_tit","tmp7"),sep = ";",remove = FALSE) %>%
      select(!c(paste0("tmp",1:7))) %>%
      separate(spec_tit,into = c("tmp","raw_file"),sep = ":") %>%
      select(!tmp) %>%
      full_join(map_df,by="raw_file") %>%
      separate(pool_id,into = c("expid","sample_id","Ecoli","inj"),sep = "_",remove = F) %>%
      mutate(new_col=paste(expid,sample_id,Ecoli,sep = "_")) %>%
      group_by(sequence,raw_file) %>% 
      distinct(sequence,.keep_all = TRUE)%>% 
      ungroup() 
    
    plt1 <- gg_barplt_id_pep_count(data_set = comb_ecoli_dist,
                                    x_df = comb_ecoli_dist$raw_file,
                                    fill_df = comb_ecoli_dist$new_col,
                                    ymax = 20000,
                                    header = paste("Total number of identified synthetic phospho-sequences across each pool",sep=" "),
                                    caption_lab = "Wrong Sequences were removed for the futher analysis.",
                                    x_lab = "Sample id",
                                    fill_lab =  "Sample id",
                                    y_lab = "Number of identified Sequence",
                                    subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name),
                                    size_num=7) + 
      scale_fill_brewer(palette = "Dark2") + theme(axis.text.x = element_text(angle = 90))
    
    write.table(comb_ecoli_dist,file=paste0(file_path,curr_dir,"/", software_name,
                                      "experiment",
                                      exp_id,acquisiton_type,
                                      "merge_identified_ecoli_sequences.tsv"),
                sep = "\t",col.names = T,row.names = F)
    
    
  }else{
    
  }
  
  #pattern <- "=\\s*([0-9.]+)"
  #pattern <- "Phospho\\s*\\([^)]+\\)\\s*=\\s*([0-9.]+)"
  pattern= "Phospho\\s*\\([^)]+\\)\\s*=\\s*([0-9.]+)"
  # Use str_extract_all to extract all matches in each string
  
  if(curr_dir =="exp1_no_ecoli_with_FAIMS"){ # This is an exception found 
                                                # only in this type
    comb_result_pos <- comb_result %>%
      separate(spectrum_title,into = c("spec_tit","tmp1","tmp2"),sep = ";",remove = FALSE) %>%
      select(!c(paste0("tmp",1:2))) %>%
      separate(spec_tit,into = c("tmp","raw_file_path"),sep = ":") %>%
      select(!tmp) %>%
      separate(raw_file_path, into = c(paste0("tmp",1:7),"raw_file"),sep = "/") %>%
      select(!c(paste0("tmp",1:7))) %>%
      mutate(raw_file=str_remove(raw_file, "\"")) %>%
      mutate(raw_file=str_remove(raw_file, ".raw")) %>%
      full_join(map_df,by="raw_file") %>%
      #rename(pool_id_map_df=pool_id) %>% #paste(expid,sample_id,Ecoli,sep = "_")) %>% 
      mutate(extracted_values = str_extract_all(ptm_sites_confidence, pattern)) %>%
      mutate(ptm_val= lapply(extracted_values, function(matches) {
        matches1 <- gsub("Phospho\\s*\\([^)]+\\)\\s*=\\s*", "", matches)
        if(length(matches1 > 1)){
          #matches1 <- as.numeric(paste(matches1, collapse = "&"))
          matches1 <- as.numeric(max(matches1))
        }else{
          #matches1 <- as.numeric(matches1)
          as.numeric(matches1)
        }
      })) %>%
      relocate(c(extracted_values,ptm_val),.after = ptm_sites_confidence)
    
  }else{
    # Extracted values
    comb_result_pos <- comb_result %>% 
      separate(spectrum_title,into = c(paste0("tmp",1:6),"spec_tit","tmp7"),sep = ";",remove = FALSE) %>%
      select(!c(paste0("tmp",1:7))) %>%
      separate(spec_tit,into = c("tmp","raw_file"),sep = ":") %>%
      select(!tmp) %>%
      full_join(map_df,by="raw_file") %>%
      #separate(pool_id,into = c("expid","sample_id","Ecoli","inj"),sep = "_") %>%
      #rename(pool_id_map_df=pool_id) %>% #paste(expid,sample_id,Ecoli,sep = "_")) %>% 
      
      mutate(extracted_values = str_extract_all(ptm_sites_confidence, pattern)) %>%
      mutate(ptm_val= lapply(extracted_values, function(matches) {
        matches1 <- gsub("Phospho\\s*\\([^)]+\\)\\s*=\\s*", "", matches)
        if(length(matches1 > 1)){
          #matches1 <- as.numeric(paste(matches1, collapse = "&"))
          matches1 <- as.numeric(max(matches1))
        }else{
          #matches1 <- as.numeric(matches1)
          as.numeric(matches1)
        }
      })) %>%
      relocate(c(extracted_values,ptm_val),.after = ptm_sites_confidence)
  }

  # ## READ THEO LIST
  # pep_list_w_theo <- read.xlsx(paste0(theo_file_path, theo_file_name), sheet = sheet_theo_name)
  # pep_list_w_theo_quant <- pep_list_w_theo[,-1]
  # common_col_theo_quant <- as.data.frame(paste(pep_list_w_theo_quant$Phosphopeptide.sequence,
  #                                              pep_list_w_theo_quant$modified.position.in.peptide, sep = "_"))
  # 
  # ## ADD COMMON COLUMN TO MERGE WITH EXP. DATA
  # colnames(common_col_theo_quant) <- "pep_with_pos"
  # pep_list_w_theo_quant_new <- cbind(common_col_theo_quant,pep_list_w_theo_quant)
  # 
  # pep_list_w_theo_quant_new <- pep_list_w_theo_quant_new %>% 
  #   rename(pool_id_theo_list=pool_id) 
  # 
  # ## EXTRACTION OF UNIQUE SEQUENCES
  # pep_list_w_theo_unique <- pep_list_w_theo %>% 
  #   distinct(Phosphopeptide.sequence,.keep_all = TRUE) %>%
  #   rename(sequence = Phosphopeptide.sequence) %>% 
  #   select(sequence,pool_id) %>%
  #   rename(pool_id_theo_list=pool_id) #%>%
  #   #mutate(Pool_for_seq_merge="Correct")
  # 
  
  
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
    rename(sequence = Phosphopeptide.sequence) %>% 
    mutate(pool_id_theo_list=pool_id) 
  
  #### TOTAL NUM. OF PHOSPHO-SEQUENCES ####
    ### Correct mapping was done using "map_df".
  comb_result_seq <- comb_result_pos %>% 
    select(sequence,
           ptm_score,
           pep_with_pos,
           pool_id,
           #pool_id_map_df,
           sample_name,
           `all_files[i]`,
           accession,
           raw_file) %>% 
    full_join(pep_list_w_theo_seq_comp,by=c("sequence","pool_id")) %>%
    select(!c(pool_id,sample_name)) %>%
    ## RE-JOINING TO COUNT CORRECT, INCORRECT and MISSING SEQ. 
    full_join(map_df,by="raw_file") %>%
    rename(pool_id_map_df=pool_id) %>%
    mutate(Pool_for_seq_merge=ifelse(pool_id_theo_list==pool_id_map_df,"Correct Seq.","Wrong Seq")) %>% #. within theo list
    #mutate_at("Pool_for_seq_merge", ~replace_na(.,"Unexpected Seq.")) %>%
    mutate(Pool_for_seq_merge= ifelse(is.na(pool_id_map_df),"Missing",Pool_for_seq_merge)) %>%
    mutate(Pool_for_seq_merge=ifelse(is.na(pool_id_theo_list),"Wrong Seq.",Pool_for_seq_merge)) %>% # out of theo. list
    mutate(acq_type=acquisiton_type) %>%
    mutate(soft_name=software_name)
  

  
    ### Duplicate sequences were removed.
  comb_result_dist <- comb_result_seq %>%
    relocate(pool_id_map_df,pool_id_theo_list,Pool_for_seq_merge,.after = sequence) %>% #
    group_by(raw_file,sequence) %>%
    distinct(sequence,.keep_all = T) %>%
    ungroup() 
  
  #### WRITE THE OBJECT AS TSV ####
  write.table(comb_result_seq,file=paste0(file_path,curr_dir,"/", software_name,
                                                "experiment",
                                                exp_id,acquisiton_type,
                                                "merge_theo_list_with_identified_phospho_sequences.tsv"),
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
    filter(!grepl("Unexpected Seq.",Pool_for_seq_merge) & !grepl("Missing Seq.",Pool_for_seq_merge))
  
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
    relocate(pool_id_map_df,pool_id_theo_list,Pool_for_pep_merge,.after = sequence) %>% 
    mutate(acq_type=acquisiton_type) %>%
    mutate(soft_name=software_name) 
  
  write.table(comb_result_pep ,file=paste0(file_path,curr_dir,"/merge_theo_list_id_phospho_sites_max_int.tsv"),
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
    
    
  ## APPLYING LOCALIZATION THRESHOLD
  #comb_result_pep_thres  <- comb_result_pep_cor %>% 
    #filter(ptm_score >0.75) #%>%
   
  ## Because of that: Error in write.table(comb_result_pep_cor, file = paste0(file_path, curr_dir,: unimplemented type 'list' in 'EncodeElement'
  #comb_result_pep_cor_write <- apply(comb_result_pep_cor, 2,as.character)
  
  #write.table(comb_result_pep_cor_write,file=paste0(file_path,curr_dir,"/","merge_theo_list_id_phospho_sites_max_int.tsv"),
              #sep = "\t",col.names = T,row.names = F)
  #software_name,"experiment",exp_id,acquisiton_type,"merge_theo_list_id_phospho_sites_max_int.tsv"
  #comb_result_pep_max_int
  
  #### REMOVE DUPLICATE PEPTIDES COMING FROM DIFFERENT INJECTION #### 
    ### BY SELECTING BEST PTM_SCORED PEPTIDES
  if(background_species == "ECOLI"){
    
    comb_result_pep_filtered <- comb_result_pep %>% 
      filter(!grepl("Missing",Pool_for_pep_merge)) %>%
      separate(sample_name,into = c("exp","samp","coli","inj"),sep = "_",remove = F) %>%
      mutate(samp_name_wo_inj=paste(exp,samp,coli,sep = "-")) %>%
      group_by(samp_name_wo_inj,pep_with_pos) %>%
      slice(which.max(ptm_score)) %>%
      ungroup()
    
  }else if (background_species ==""){
    
    comb_result_pep_filtered <- comb_result_pep %>% 
      filter(!grepl("Missing",Pool_for_pep_merge)) %>%
      separate(sample_name,into = c("exp","samp","inj"),sep = "_",remove = F) %>%
      mutate(samp_name_wo_inj=paste(exp,samp,sep = "-")) %>%
      group_by(samp_name_wo_inj,pep_with_pos) %>%
      slice(which.max(ptm_score)) %>%
      ungroup()

  }else{}
  
  #### LOCALIZATION ACCURACY ASSESSMENT (ROC LIKE PLOT GENERATION) #### 
     ### Column selection
  comb_result_pep_rmv_miss <- comb_result_pep_filtered %>% 
    select(pep_with_pos,Pool_for_pep_merge,ptm_score,pool_id_map_df)
  
  ### Custom localization threshold determination and
  ### Counting total number of correct and wrong localization at a given threshold
  threshold <- seq(0,1,length=100)
  final_df <- NULL
  
  for (i in 1:length(threshold)){
    
    df <- comb_result_pep_rmv_miss %>% subset(ptm_score > threshold[i]) %>% 
      count(Pool_for_pep_merge) %>% 
      mutate(threshold_val = threshold[i])
    
    final_df <-bind_rows(final_df,df)
    rm(df)
  }
  
  write.table(final_df, file=paste0(file_path,curr_dir,"/","Experiment",exp_id,software_name,"_number_of_sites_with_scores.tsv"),sep = "\t",col.names = T,row.names = F)
  
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
             filter(sapply(sequence, function(x) str_count(x, pattern = big_pattern)) == i) %>%
             mutate(STY_adj = sapply(sequence, function(y) extract_pos_adj(query = y,patt_text = get(paste0("perm",i+2))$comb))) %>%
             mutate(STY_len=nchar(STY_adj)) %>%
             mutate(is_adj= "adjacent"))
  }
  
  ### ADJACENT COUNT FROM 2
  comb_result_pep_adj_2 <- comb_result_pep_filtered %>%#comb_result_pep
    #select(Sequence) %>%
    #filter(!(sapply(Sequence, function(x) str_count(x, pattern = big_pattern)) == 1) & 
    #!(sapply(Sequence, function(x) str_count(x, pattern = big_pattern)) == 2)) %>%
    
    filter(sapply(sequence, function(x) str_count(x, pattern = big_pattern)) < 1 & 
             sapply(sequence, function(x) str_count(x, pattern = "(S(?=S|T|Y)|T(?=S|T|Y)|Y(?=S|T))+")) == 1) %>%
    
    #filter(sapply(Sequence, function(x) str_count(x, pattern = "(S(?=S|T|Y)|T(?=S|T|Y)|Y(?=S|T))+")) == 2) %>%
    ### LIST MIGHT BE ADDED TO COLLECT ALL PAIRS 
    mutate(STY_adj = sapply(sequence, function(x) extract_pos_adj(query = x,patt_text = perm2$comb))) %>%
    mutate(STY_len=nchar(STY_adj)) %>%
    mutate(is_adj= "adjacent")

  ### NON-ADJACENT COUNT
  comb_result_pep_non_adj <- comb_result_pep_filtered %>% #comb_result_pep
    filter(sapply(sequence, function(x) str_count(x, pattern = big_pattern)) < 1 & 
             sapply(sequence, 
                    function(x) str_count(x,
                                          pattern = "(S(?=S|T|Y)|T(?=S|T|Y)|Y(?=S|T))+")) == 0) %>%
    rowwise() %>%
    mutate(STY_adj = list(sapply(sequence, function(x) extract_pos_non_adj(query = x, patt_text = elements)))) %>%
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
    select(sequence,
           STY_adj,
           STY_len,
           Pool_for_pep_merge,
           ptm_score,
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
    
    ggplot( aes(x= ptm_score,fill=interaction(is_adj))) +
    geom_bar(stat="count")  + labs(x = "Count of STY amino acids", y = "Total count",
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




#   
#   quant_peptides <- read.xlsx(paste0(file_path,file_name), sheet = sheet_name)
#   
#   pep_list_w_theo <- read.xlsx(paste0(theo_file_path, theo_file_name), sheet = sheet_theo_name)
#   pep_list_w_theo_quant <- pep_list_w_theo[,-1]
# 
#   
#   pep_list_w_theo_unique <- pep_list_w_theo %>% 
#     distinct(Phospopeptide.sequence,.keep_all = TRUE) %>%
#     rename(sequence = Phospopeptide.sequence) %>% 
#     select(sequence,Pool) %>%
#     rename(Pool_for_seq_merge=Pool)
#   
#   if(background_species == ""){
#     
#   }else{
#     ecoli_seq_dist <- quant_peptides_cor_abun %>% 
#       select(sequence, modifications, accession) %>%
#       filter(grepl(background_species,accession)) %>%
#       distinct(sequence,.keep_all = T) %>%
#       mutate(species=background_species) 
#     
#     ecoli_seq <- quant_peptides_cor_abun %>% 
#       select(sequence, modifications, accession) %>%
#       filter(grepl(background_species,accession)) %>%
#       mutate(species=background_species) 
#   }
#   
#  
#   
#   all_seq <- quant_peptides_cor_abun %>% 
#     filter(grepl("Phospho",modifications) & grepl(selected_spcies,accession)) %>%
#     #distinct(sequence, .keep_all = T) %>%
#     separate(accession, into = c("protein","species"),remove = F,sep="_") %>%
#     full_join(pep_list_w_theo_unique,by="sequence") %>%
#     mutate_at("Pool_for_seq_merge", ~replace_na(.,"Unexpected")) %>%
#     mutate(Pool_for_seq_merge= ifelse(is.na(species),"missing",Pool_for_seq_merge)) %>%
#     filter(!grepl("Unexpected",Pool_for_seq_merge)) %>%
#     bind_rows(ecoli_seq) %>% 
#     mutate(Pool_for_seq_merge= ifelse(is.na(Pool_for_seq_merge),background_species,Pool_for_seq_merge))
#   #mutate(acq_type=acquisiton_type) %>%
#   #mutate(soft_name=software_name)
#   
#   all_seq_syn <- all_seq %>%
#     select(sequence, modifications, Pool_for_seq_merge) %>%
#     distinct(sequence, .keep_all = T) %>%
#     bind_rows(ecoli_seq_dist) %>%
#     mutate(Pool_for_seq_merge= ifelse(is.na(Pool_for_seq_merge),background_species,Pool_for_seq_merge)) %>%
#     mutate(acq_type=acquisiton_type) %>%
#     mutate(soft_name=software_name)
#   
#   p13 <- gg_barplt_id_pep_count(data_set = all_seq_syn,
#                                 x_df = all_seq_syn$Pool_for_seq_merge,
#                                 fill_df = all_seq_syn$Pool_for_seq_merge,
#                                 ymax = 20000,
#                                 header = paste("Total number of identified phosphorylated", selected_spcies,"and", background_species,"across each sample",sep=" "),
#                                 caption_lab = "NA values are removed.",
#                                 x_lab = "Sample id",
#                                 fill_lab =  "Sample id",
#                                 y_lab = "Number of identified peptides",
#                                 subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
#   
#   write.table(all_seq_syn, file=paste0(file_path,"Experiment2",software_name,"_number_of_unique_sequence_for_each_species.txt"),sep = "\t",col.names = T,row.names = F)
#   #################################################
#   
#   
#   
# }  
#   
#   
# # Creation of common column merging peptide sequence and phospho positions -> experimental data
# common_col_exp <- as.data.frame(paste(id_syn_phospho_pep$sequence, phospho_ptm_pos_df, sep = "_"))
# colnames(common_col_exp) <- "common_col_for_merging"
# 
# # Bind it to the quant data
# syn_phospho_pep_new <- cbind(common_col_exp,id_syn_phospho_pep)
# 
# 
# #### THEORETICAL SEQUENCE DATA ####
# sheet_theo_name <- "final_pep_list_with_pos"
# theo_path = "D:/dev/Desktop_copy/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/"
# theo_file = "pep_list_ordered_in_pool_id_with_pos.xlsx"
# pep_list_w_theo <- read.xlsx(paste0(theo_path, theo_file), sheet = sheet_theo_name)
# #pep_list_w_theo <- pep_list_w_theo[,-1]
# 
# 
# common_col_theo <- as.data.frame(paste(pep_list_w_theo$Sequence,
#                                        pep_list_w_theo$modified.position.in.peptide, sep = "_"))
# colnames(common_col_theo) <- "common_col_merging"
# pep_list_w_theo_new <- cbind(common_col_theo,pep_list_w_theo)
# colnames(pep_list_w_theo_new)[3] <- "sequence"
# 
# #####
# #
# final_results_with_common_col <- mutate(result_with_common_col,id_syn_phospho_pep)
# #
# # Combine data frame (we will continue with this for further step)
# df_merge_phosphosite <- final_results_with_common_col %>% 
#   left_join(pep_list_w_theo_new,by="common_col_merging")
# #####
# 
# mapping_path <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_1/MaxQuant/"
# mapping_file <- "exp1_wo_FAIMS_wo_Ecoli_mapping.txt"
# ## MAPPING_PATH
# raw_files_order <- as.data.frame(read_tsv(paste0(mapping_path,mapping_file)))
# 
# 
# ### EXPERIMENTAL DATA SORTING BASED ON POOL ID
# for (i in 1:dim(raw_files_order)[1]){
#   
#   query <- paste0("_", strsplit(raw_files_order[i,"raw_file"], "_")[[1]][2])
#   assign(paste0(raw_files_order[i,"pool_id"],query), syn_phospho_pep_new %>% filter(grepl(query,spectrum_title)))
#   
#   rm(query)
#   
# }
# pool_size <- 8
# ## TO MERGE RAWS FILES ARE ASSOCIATED TO THE SAME POOL_ID
# for (j in 1:pool_size){
#   
#   object_list <- mget(ls(pattern=paste0("pool",j)))
#   assign(paste0("pool",j,"_all"), do.call(rbind, object_list))
#   
#   rm(object_list)
#   
#   # new_col_name <- "common_col_merging"
#   # pool_j_all <- get(paste0("pool", j, "_all"))
#   # pool_j_all[[new_col_name]] <- paste(pool_j_all$Stripped.Sequence, pool_j_all$value, sep = "_")
#   #assign(paste0("pool", j, "_all"), pool_j_all)
#   
#   ## WITHOUT INSIDE FOR LOOP
#   #pool1_all[,"common_col_merging"] <- paste(pool1_all$Stripped.Sequence,
#   #pool1_all$value, sep="_")
#   
#   tmp <- pep_list_w_theo_new %>%
#     filter(grepl(paste0("pool",j), pool_id)) %>%
#     ## DON'T WORRY ABOUT WARNING MESSAGE,
#     ## DUPLICATES WILL BE ELIMINATED IN THE NEXT STEP
#     full_join(get(paste0("pool",j,"_all")),by="sequence") ## TODO:TRY left_join() without filtering duplicates
#   
#   # ## REMOVE MULTIPLE PHOSPHOSITE THAT ARE RELATED TO THE SAME PEPTIDE
#   assign(paste0("pool",j,"merged"),tmp %>%
#            filter(duplicated(sequence) == FALSE))
#   
#   ## UNEXPECTEDLY IDENTIFIED - FALSE POSITIVES
#   assign(paste0("unexp_id_pool",j),get(paste0("pool",j,"merged")) %>% 
#            filter(is.na(spectrum_title) == FALSE & is.na(Well.position) == TRUE) %>%
#            mutate(type=paste0("unexpected_pool",j)))
#   
#   ## CORRECTLY IDENTIFIED - TRUE POSTIVIES
#   assign(paste0("correct_id_pool",j),get(paste0("pool",j,"merged")) %>% 
#            filter(is.na(spectrum_title) == FALSE & is.na(Well.position) == FALSE)%>%
#            mutate(type=paste0("correct_pool",j)))
#   
#   ### MISSED IDENTIFIED - FALSE NEGATIVES
#   assign(paste0("missed_id_pool",j),get(paste0("pool",j,"merged")) %>% 
#            filter(is.na(spectrum_title) == TRUE & is.na(Well.position) == FALSE)%>%
#            mutate(type=paste0("missed_pool",j)))
#   
# }
# 
# type_of_id <- c("correct_id_pool","missed_id_pool","unexp_id_pool")  
# 
# for (k in 1:length(type_of_id)){
#   
#   assign(paste0("final_",type_of_id[k],"list"),mget(ls(pattern=paste0(type_of_id[k]))))
#   assign(paste0("ffinal_",type_of_id[k]), do.call(rbind, get(paste0("final_",type_of_id[k],"list"))))
#   #assign(paste0("ffinal_",type_of_id[k]),cbind(get(paste0("ffinal_",type_of_id[k])),paste0(type_of_id[k])))
#   
# }
# 
# 
# final_all_id_list <- mget(ls(pattern=paste0("ffinal")))
# final_all_id <- do.call(rbind, final_all_id_list)
# 
# 
# output_path <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_1/Proline/exp1_wo_Ecoli_noFAIMS/"
# output_file <- "modified_processable_version_proline_output.txt"
# write_tsv(final_all_id,paste0(output_path,output_file))
# 
# library(viridis)
# library(ggplot2)
# final_all_id  %>% 
#   group_by(type) %>% 
#   summarise(n()) %>%
#   separate(type,c('type_name','pool_id')) %>%
#   ggplot( aes(fill=type_name,x=pool_id,y=`n()`)) + geom_bar(stat="identity",position = "dodge")+
#   geom_text(aes(label=`n()`),colour = "gray", size = 4,position = position_dodge(0.9),vjust=1.6) +
#   theme_minimal() +
#   labs(title="Number of identified peptides using Proline") + 
#   theme(legend.text = element_text(size=15), 
#         axis.title.x = element_text(size = 15),
#         axis.title.y = element_text(size = 15),
#         plot.title = element_text(size=20),
#         legend.title=element_text(size=15),
#         axis.text=element_text(size=15),
#         axis.title=element_text(size=15)) + 
#   scale_fill_viridis(discrete = TRUE,option = "viridis")
# 
# 
# 
# 
# 
# 
# 
# intensity_match <- final_results_with_common_col %>%
#   select(common_col_merging,Intensity,Experiment) %>% 
#   tibble() %>%
#   mutate(row = row_number()) %>%
#   pivot_wider(names_from = "Experiment",
#               values_from = "Intensity")
# 
# pep_list_with_pool_id <- read.xlsx(paste0("D:/dev/Desktop_copy/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/",
#                                           "pep_list_ordered_in_pool_id_with_pos.xlsx"), sheet = "final_pep_list_with_pos")
# 
# pep_list_with_pool_id[,"common_col_merging"] <- paste(pep_list_with_pool_id$sequence,
#                                                       pep_list_with_pool_id$modified.position.in.peptide,
#                                                       sep="_")
# 
# 
# 
# 
# 
# 
# 
# 
