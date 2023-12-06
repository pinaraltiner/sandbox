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
      filter(grepl("Phospho",modifications) & grepl(selected_species,accession)) %>%
      mutate(all_files[i]) %>%
      mutate(phospho_pos = proline_phospho_pos_extraction(modifications)) %>%
      mutate(pep_with_pos = paste(sequence,phospho_pos,sep = "_")) %>%
      distinct(pep_with_pos,.keep_all = TRUE) ## The reason of putting this line
                                              ## is for eliminating the combine 
                                              ## version of phospho. Such as;
                                              ## N-term phospho, oxidation phospho, etc.
    
    #assign(paste0(all_files[i]), phospho_tmp)
    assign(paste0("comb_result"),bind_rows(comb_result, phospho_tmp))
    
    if(background_species == "ECOLI"){
      ecoli_tmp <- tmp %>% 
        filter(grepl(background_species,accession)) %>%
        mutate(all_files[i]) %>%
        mutate(acq_type=acquisiton_type) %>%
        mutate(soft_name=software_name) %>%
        distinct(sequence,.keep_all = TRUE)
      
      assign(paste0("comb_ecoli"), bind_rows(comb_ecoli,ecoli_tmp))
    }else{
      
    }
    
    rm(tmp,phospho_tmp,ecoli_tmp)
    
  }

  if(background_species == "ECOLI"){
    
    comb_ecoli <- comb_ecoli %>% 
      separate(spectrum_title,into = c(paste0("tmp",1:6),"spec_tit","tmp7"),sep = ";",remove = FALSE) %>%
      select(!c(paste0("tmp",1:7))) %>%
      separate(spec_tit,into = c("tmp","raw_file"),sep = ":") %>%
      select(!tmp) %>%
      full_join(map_df,by="raw_file") %>%
      separate(pool_id,into = c("expid","sample_id","Ecoli","inj"),sep = "_") %>%
      mutate(new_col=paste(expid,sample_id,Ecoli,sep = "_")) %>%
      group_by(Sequence,new_col) %>% 
      distinct(Sequence,.keep_all = TRUE)%>% 
      ungroup() 
    
    plot6 <- gg_barplt_id_pep_count(data_set = comb_ecoli,
                                    x_df = comb_ecoli$raw_file,
                                    fill_df = comb_ecoli$new_col,
                                    ymax = 20000,
                                    header = paste("Total number of identified synthetic phospho-sequences across each pool",sep=" "),
                                    caption_lab = "Wrong Sequences were removed for the futher analysis.",
                                    x_lab = "Sample id",
                                    fill_lab =  "Sample id",
                                    y_lab = "Number of identified Sequence",
                                    subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) + 
      scale_fill_brewer(palette = "Dark2")
    
    write.table(comb_ecoli,file=paste0(file_path,curr_dir,"/", software_name,
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
    comb_result <- comb_result %>% 
      separate(spectrum_title,into = c("spec_tit","tmp1","tmp2"),sep = ";",remove = FALSE) %>%
      select(!c(paste0("tmp",1:2))) %>%
      separate(spec_tit,into = c("tmp","raw_file_path"),sep = ":") %>%
      select(!tmp) %>%
      separate(raw_file_path, into = c(paste0("tmp",1:7),"raw_file"),sep = "/") %>%
      select(!c(paste0("tmp",1:7))) %>%
      mutate(raw_file=str_remove(raw_file, "\"")) %>%
      full_join(map_df,by="raw_file") %>%
      separate(pool_id,into = c("expid","sample_id","Ecoli","inj"),sep = "_") %>%
      mutate(new_col=paste(expid,sample_id,Ecoli,sep = "_")) %>% 
      
      mutate(extracted_values = str_extract_all(ptm_sites_confidence, pattern)) %>%
      mutate(ptm_val= lapply(extracted_values, function(matches) {
        matches1 <- gsub("Phospho\\s*\\([^)]+\\)\\s*=\\s*", "", matches)
        if(length(matches1 > 1)){
          matches1 <- paste(matches1, collapse = "&")
          #matches1 <- as.numeric(max(matches1))
        }else{
          #matches1 <- as.numeric(matches1)
          matches1
        }
      })) %>%
      relocate(c(extracted_values,ptm_val),.after = ptm_sites_confidence)
    
  }else{
    # Extracted values
    comb_result <- comb_result %>% 
      separate(spectrum_title,into = c(paste0("tmp",1:6),"spec_tit","tmp7"),sep = ";",remove = FALSE) %>%
      select(!c(paste0("tmp",1:7))) %>%
      separate(spec_tit,into = c("tmp","raw_file"),sep = ":") %>%
      select(!tmp) %>%
      full_join(map_df,by="raw_file") %>%
      separate(pool_id,into = c("expid","sample_id","Ecoli","inj"),sep = "_") %>%
      mutate(new_col=paste(expid,sample_id,Ecoli,sep = "_")) %>% 
      
      mutate(extracted_values = str_extract_all(ptm_sites_confidence, pattern)) %>%
      mutate(ptm_val= lapply(extracted_values, function(matches) {
        matches1 <- gsub("Phospho\\s*\\([^)]+\\)\\s*=\\s*", "", matches)
        if(length(matches1 > 1)){
          matches1 <- as.numeric(paste(matches1, collapse = "&"))
          #matches1 <- as.numeric(max(matches1))
        }else{
          #matches1 <- as.numeric(matches1)
          as.numeric(matches1)
        }
      })) %>%
      relocate(c(extracted_values,ptm_val),.after = ptm_sites_confidence)
    
    
  }
 
  
  ## READ THEO LIST
  pep_list_w_theo <- read.xlsx(paste0(theo_file_path, theo_file_name), sheet = sheet_theo_name)
  pep_list_w_theo_quant <- pep_list_w_theo[,-1]
  common_col_theo_quant <- as.data.frame(paste(pep_list_w_theo_quant$Phosphopeptide.sequence,
                                               pep_list_w_theo_quant$modified.position.in.peptide, sep = "_"))
  
  ## ADD COMMON COLUMN TO MERGE WITH EXP. DATA
  colnames(common_col_theo_quant) <- "pep_with_pos"
  pep_list_w_theo_quant_new <- cbind(common_col_theo_quant,pep_list_w_theo_quant)
  
  
  ## EXTRACTION OF UNIQUE SEQUENCES
  pep_list_w_theo_unique <- pep_list_w_theo %>% 
    group_by(pool_id) %>%
    summarise(unique(Phosphopeptide.sequence)) %>%
    rename(Sequence = 	
             `unique(Phosphopeptide.sequence)`) %>% 
    select(Sequence,pool_id) %>%
    mutate(Pool_for_seq_merge="Correct")
   #rename(Pool_for_seq_merge=pool_id)
  
  
  if(background_species == "ECOLI"){
    
    comb_ecoli_result_seq <- comb_result %>% 
      select(sequence,ptm_score,new_col,inj,`all_files[i]`,
                                             accession,raw_file) %>% #Intensity, #Marked.as,
      rename(Sequence =sequence) %>%
      group_by(raw_file,Sequence) %>%
      summarise(unique(Sequence),
                #pep_with_pos,
                ptm_score,
                new_col,inj,
                `all_files[i]`,
                accession) %>%
      mutate(acq_type=acquisiton_type) %>%
      mutate(soft_name=software_name) %>%
      mutate(new_col_inj = paste(new_col,inj,sep = "_")) %>%
      ## THIS PART IS NEW!
      group_by(Sequence,raw_file) %>% 
      distinct(Sequence,.keep_all = TRUE)%>% 
      ungroup() 
    
    write.table(comb_ecoli_result_seq,file=paste0(file_path,curr_dir,"/", software_name,
                                           "experiment",
                                           exp_id,acquisiton_type,
                                           "merge_theo_list_with_identified_phospho_sequences.tsv"),
                sep = "\t",col.names = T,row.names = F)
    
    plot1 <- gg_barplt_id_pep_count(data_set = comb_ecoli_result_seq,
                                    x_df = comb_ecoli_result_seq$new_col ,
                                    fill_df = comb_ecoli_result_seq$Pool_for_seq_merge,
                                    ymax = 250,
                                    header = paste("Total number of identified synthetic phospho-sequences across each pool",sep=" "),
                                    caption_lab = "Wrong Sequences were removed for the futher analysis.",
                                    x_lab = "Sample id",
                                    fill_lab =  "Sample id",
                                    y_lab = "Number of identified Sequence",
                                    subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
    
  }else{
    
  }
  
    comb_result_seq <- comb_result %>%
      select(sequence, pep_with_pos,ptm_score,new_col,inj,`all_files[i]`,
             accession,raw_file) %>% #Intensity, #Marked.as,
      rename(Sequence =sequence) %>%
      group_by(raw_file,Sequence) %>%
      summarise(unique(Sequence),
                pep_with_pos,
                ptm_score,
                new_col,inj,
                `all_files[i]`,
                accession
      ) %>%
      full_join(pep_list_w_theo_unique,by="Sequence") %>%
      mutate_at("Pool_for_seq_merge", ~replace_na(.,"Unexpected")) %>%
      mutate(Pool_for_seq_merge= ifelse(is.na(raw_file),"missing",Pool_for_seq_merge)) %>%
      #filter(!grepl("Unexpected",Pool_for_seq_merge))
      mutate(acq_type=acquisiton_type) %>%
      mutate(soft_name=software_name)
    
    # comb_result_seq_dist <- comb_result_seq %>%
    #   group_by(Sequence,raw_file) %>% ##pool_id
    #   distinct(Sequence,.keep_all = TRUE)%>%
    #   ungroup()
    # 
    # #comb_result_seq_dist <- comb_result_seq %>% distinct(Sequence,.keep_all = TRUE)
    # plot1 <- gg_barplt_id_pep_count(data_set =comb_result_seq_dist,
    #                                 x_df =comb_result_seq_dist$raw_file,
    #                                 fill_df = comb_result_seq_dist$Pool_for_seq_merge ,
    #                                 ymax = 250,
    #                                 header = paste("Total number of identified synthetic phospho-sequences across each pool",sep=" "),
    #                                 caption_lab = "Wrong Sequences were removed for the futher analysis.",
    #                                 x_lab = "Sample id",
    #                                 fill_lab =  "Sample id",
    #                                 y_lab = "Number of identified Sequence",
    #                                 subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) +
    #   scale_fill_brewer(palette = "Paired") + theme(axis.text.x = element_text(angle = 90))
    
    comb_result_seq_dist_cor <- comb_result_seq %>%
      filter(grepl("Correct",Pool_for_seq_merge)) %>%
      group_by(Sequence,raw_file) %>%
      distinct(Sequence,.keep_all = TRUE)%>%
      ungroup()
    
    plot1 <-   gg_barplt_id_pep_count(data_set =comb_result_seq_dist_cor,#comb_result_seq_dist,
                                      x_df =comb_result_seq_dist_cor$raw_file,#comb_result_seq_dist$raw_file,
                                      fill_df = comb_result_seq_dist_cor$Pool_for_seq_merge,#comb_result_seq_dist$Pool_for_seq_merge ,
                                      ymax = 250,
                                      header = paste("Total number of identified synthetic phospho-sequences across each pool",sep=" "),
                                      caption_lab = "Duplicates were removed for the futher analysis.",
                                      x_lab = "Sample id",
                                      fill_lab =  "Sample id",
                                      y_lab = "Number of identified Sequence",
                                      subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) +
      scale_fill_brewer(palette = "Paired") + theme(axis.text.x = element_text(angle = 90))
    
    
    plot23 <- ggplot(comb_result_seq_dist_cor, aes(x=comb_result_seq_dist_cor$raw_file,
                                              fill=comb_result_seq_dist_cor$pool_id )) +
      geom_bar() + theme_bw() +
      geom_text(color="black",aes(label=after_stat(count)),
                show.legend = F,stat = "count",size=8,
                position = position_stack(vjust = 0.5)) +
      scale_fill_brewer(palette = "Dark2") +
      labs(title=paste("Number of Identified phospho-sequence for each pool \n","Experiment",
                       exp_id, acquisiton_type),x="Sample id",y="Number of identified Sequence",
           subtitle = "\n Contamination assessment 2",fill="Pool id") +
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
            axis.title=element_text(size=30))
    
    
    
    comb_result_pep <- comb_result %>% #comb_result_seq %>%
      #filter(!grepl("Unexpected",Pool_for_seq_merge) & !grepl("missing",Pool_for_seq_merge)) %>%
      #select(!c(Pool_for_seq_merge, pool_id)) %>%
      #group_by(`all_files[i]`,pep_with_pos) %>%
      #summarise(unique(pep_with_pos),
      #Sequence,
      #ptm_score,new_col,inj,
      #`all_files[i]`,accession) %>%
      #ungroup() %>%
      full_join(pep_list_w_theo_quant_new,by="pep_with_pos") %>% ## IF FULL_JOIN IS USED,
      group_by(pep_with_pos,raw_file) %>%
      distinct(pep_with_pos,.keep_all = TRUE)%>%
      ungroup() %>%
      mutate(Pool_for_pep_merge=NA) %>%
      mutate(Pool_for_pep_merge=ifelse(is.na(Pool),"Wrong Localization","Correct")) %>%  ## MISSING PEPTIDE INFO COULD BE OBTAINED
      mutate(Pool_for_pep_merge= ifelse(is.na(raw_file),"missing",Pool_for_pep_merge)) %>%
      mutate(Pool=ifelse(is.na(Pool),"Unexpected",Pool)) %>%
      mutate(Pool= ifelse(is.na(raw_file),"missing",Pool)) %>%
      mutate(isomericity=ifelse(is.na(isomericity),"Unexpected",isomericity)) %>%
      mutate(isomericity=ifelse(is.na(raw_file),"missing",isomericity)) %>%
      #drop_na(pool_id) ## ELIMINATION OF UNEXPECTED PEPTIDES
      mutate(acq_type=acquisiton_type) %>%
      mutate(soft_name=software_name) %>%
      rename(Sequence = sequence)
    
    plot21 <- ggplot(comb_result_pep, aes(x=raw_file,fill=Pool_for_pep_merge)) +
      geom_bar() + theme_bw() +
      geom_text(color="black",aes(label=after_stat(count)),
                show.legend = F,stat = "count",size=8,
                position = position_stack(vjust = 0.5)) +
      scale_fill_brewer(palette = "Dark2") +
      labs(title=paste("Number of Identified phospho-sites for each pool \n","Experiment",
                       exp_id, acquisiton_type),x="Sample id",y="Number of identified Sequence",
           subtitle = "Contaminant assessment 1",fill="Pool id") +
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
            axis.title=element_text(size=30))
    
    comb_result_pep_cor <- comb_result_pep %>%
      filter(grepl("Correct",Pool_for_pep_merge))
    
    plot23 <- ggplot(comb_result_pep_cor, aes(x=comb_result_pep_cor$raw_file,fill=pool_id)) +
      geom_bar() + theme_bw() +
      geom_text(color="black",aes(label=after_stat(count)),
                show.legend = F,stat = "count",size=8,
                position = position_stack(vjust = 0.5)) +
      scale_fill_brewer(palette = "Dark2") +
      labs(title=paste("Distribution of Correctly Identified phospho-sites for each pool \n","Experiment",
                       exp_id, acquisiton_type),x="Sample id",y="Number of identified Sequence",
           subtitle = "Contaminant assessment 3") +
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
            axis.title=element_text(size=30))
    
  ## NAs are removed
  comb_result_pep_thres  <- comb_result_pep_cor %>% 
    #group_by(pep_with_pos,raw_file) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    filter(ptm_score >0.75) %>%
    #slice(which.max(ptm_score)) %>% ## ELIMINATE MULTIPLE CHARGES
    ungroup() 
  ## Because of that: Error in write.table(comb_result_pep_cor, file = paste0(file_path, curr_dir,: unimplemented type 'list' in 'EncodeElement'
  comb_result_pep_cor_write <- apply(comb_result_pep_cor, 2,as.character)
  
  write.table(comb_result_pep_cor_write,file=paste0(file_path,curr_dir,"/","merge_theo_list_id_phospho_sites_max_int.tsv"),
              sep = "\t",col.names = T,row.names = F)
  #software_name,"experiment",exp_id,acquisiton_type,"merge_theo_list_id_phospho_sites_max_int.tsv"
  #comb_result_pep_max_int
  
  
  plot5 <- comb_result_pep_thres %>% #distinct(pep_with_pos,.keep_all = T) %>%
    ggplot( aes(x=raw_file,fill=pool_id)) +
    geom_bar() + theme_bw() +
    geom_text(color="black",aes(label=after_stat(count)),
              show.legend = F,stat = "count",size=8,
              position = position_stack(vjust = 0.5)) +
    scale_fill_brewer(palette = "Dark2") +
    labs(title=paste("Distribution of Correctly Identified phospho-sequence for each pool \n after applying localization threshold","Experiment",
                     exp_id, acquisiton_type),x="Sample id",y="Number of identified Sequence",subtitle = "Contaminant assessment 4") +
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
          axis.title=element_text(size=30))
  # 
  # plot5 <- comb_result_pep_max_int %>% #distinct(pep_with_pos,.keep_all = T) %>%
  #   ggplot( aes(x=Pool,fill=isomericity)) +
  #   geom_bar() + theme_bw() +
  #   geom_text(color="black",aes(label=after_stat(count)),
  #             show.legend = F,stat = "count",size=8,
  #             position = position_stack(vjust = 0.5)) + 
  #   scale_fill_brewer(palette = "Dark2") +
  #   labs(title=paste("Number of Identified phospho-sites based on pools of Experiment 2",
  #                    acquisiton_type),x="Sample id",y="Number of identified Sequence") +
  #   theme_minimal() +
  #   theme(legend.text = element_text(size=30), 
  #         axis.title.x = element_text(size=30),
  #         axis.title.y = element_text(size=30),
  #         plot.title = element_text(size=35),
  #         plot.subtitle = element_text(size = 25),
  #         plot.caption = element_text(size = 25),
  #         legend.title=element_text(size=30),
  #         axis.text.x = element_text(size=20,angle = 90),
  #         axis.text.y = element_text(size = 30),
  #         axis.title=element_text(size=30))
  # 
  comb_result_pep_rmv_miss <- comb_result_pep %>% filter(!grepl("missing",Pool_for_pep_merge)) #%>%
    #select(pep_with_pos,Sequence,ptm_score,raw_file,pool_id,Pool_for_pep_merge)
  
  pool_list <- c(paste0("pool",1:8))
  raw_list <- unique(comb_result_pep_rmv_miss$raw_file)
  correct_match_poolid_rawfile <- NULL
  
  for (i in 1:length(unique(comb_result_pep_rmv_miss$raw_file))){
    
    assign(paste0("pool",i), comb_result_pep_rmv_miss %>% 
             select(pep_with_pos,
                    Sequence,
                    pool_id,
                    raw_file,
                    new_col,
                    ptm_score,
                    Pool_for_pep_merge) %>%
             rowwise() %>%
             filter(grepl(pool_list[i],pool_id)) %>%
             mutate(cor_pool=ifelse(raw_file == raw_list[i] && pool_id == pool_list[i],"Correct_match","Wrong match"))) #%>%
             #mutate(cor_pool=ifelse(is.na(pool_id[i]), "Unexpected",cor_pool)) %>%
             #filter(grepl("Correct_match",cor_pool) | grepl( "Unexpected",cor_pool)))
             #filter(grepl("Correct_match",cor_pool)))
    
    correct_match_poolid_rawfile <- rbind(correct_match_poolid_rawfile, get(paste0("pool",i)))
  }

  
  threshold <- seq(0,1,length=100)
  final_df <- NULL
  
  for (i in 1:length(threshold)){
    # assign(paste0("tmp",i), comb_result_pep %>% 
    #   filter(max_value_column > threshold[i]) %>%
    #   count(Pool_for_pep_merge,pool_id.x))
    df <- correct_match_poolid_rawfile %>% subset(ptm_score > threshold[i]) %>% 
      
      count(cor_pool) %>% 
      
      mutate(threshold_val = threshold[i])
    
    final_df <-bind_rows(final_df,df)
    rm(df)
  }
  
  write.table(final_df, file=paste0(file_path,curr_dir,"/","Experiment",exp_id,software_name,"_number_of_sites_with_scores.tsv"),sep = "\t",col.names = T,row.names = F)
  
  #wrong_df <- final_df %>% filter(grepl("Wrong",Pool_for_pep_merge))
  correct_df <- final_df %>% filter(grepl("Correct_match",cor_pool))
  
  plot4 <- ggplot(correct_df,aes(y=n, x=threshold_val)) + geom_line(size=2) +
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
  
  
  ############################################
  elements <- c("S", "T", "Y")
                 ## 1 was removed because it is the same as "element object"
  for(i in 2:5){ ## 5 might be selected as a parameter of the main function
    assign(paste0("perm",i),
           as.data.frame(permutations(n = length(elements), r = i, v = elements, repeats.allowed = TRUE)))
    
    #perm <- as.data.frame(perm)
    #assign(get(paste0("perm",i)[,"comb"]),as.data.frame(get(paste0("perm",i))$V1))
    
  }
  #perm <- permutations(n = length(elements), r = 1, v = elements, repeats.allowed = TRUE)
  #perm <- as.data.frame(perm)
  #perm[,"comb"] <- as.data.frame(paste0(perm$V1))
  
  #perm2 <- permutations(n = length(elements), r = 2, v = elements, repeats.allowed = TRUE)
  #perm2 <- as.data.frame(perm2)
  perm2[,"comb"] <- as.data.frame(paste0(perm2$V1,perm2$V2))
  #colnames(perm2)[3] <- "comb"
  
  #perm3 <- permutations(n = length(elements), r = 3, v = elements, repeats.allowed = TRUE)
  #perm3 <- as.data.frame(perm3)
  perm3[,"comb"] <- as.data.frame(paste0(perm3$V1,perm3$V2,perm3$V3))
  #colnames(perm3)[4] <- "comb"
  
  #perm4 <- permutations(n = length(elements), r = 4, v = elements, repeats.allowed = TRUE)
  #perm4 <- as.data.frame(perm4)
  perm4[,"comb"] <- as.data.frame(paste0(perm4$V1,perm4$V2,perm4$V3,perm4$V4))
  
  #perm5 <- permutations(n = length(elements), r = 5, v = elements, repeats.allowed = TRUE)
  #perm5 <- as.data.frame(perm5)
  perm5[,"comb"] <- as.data.frame(paste0(perm5$V1,perm5$V2,perm5$V3,perm5$V4,perm5$V5))
  
  possibilites <-perm2 %>% bind_rows(perm3,perm4,perm5) %>% select(comb)
  possibilites <- possibilites$comb
  
  #sing_STY <-perm$comb
  
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
  # 
  big_pattern ="(S(?=(S|T|Y)(S|T|Y))|T(?=(S|T|Y)(S|T|Y))|Y(?=(S|T|Y)(S|T|Y)))"
  
  for(i in 1:3){
    assign(paste0("comb_result_pep_adj_",i+2),
          comb_result_pep %>%
             #select(Sequence) %>%
             filter(sapply(Sequence, function(x) str_count(x, pattern = big_pattern)) == i) %>%
             mutate(STY_adj = sapply(Sequence, function(y) extract_pos_adj(query = y,patt_text = get(paste0("perm",i+2))$comb))) %>%
             mutate(STY_len=nchar(STY_adj)) %>%
             mutate(is_adj= "adjacent"))
      
  }
  
  comb_result_pep_adj_2 <- comb_result_pep %>%#comb_result_pep
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

  
  comb_result_pep_non_adj <- comb_result_pep %>% 
    filter(sapply(Sequence, function(x) str_count(x, pattern = big_pattern)) < 1 & 
             sapply(Sequence, 
                    function(x) str_count(x,
                                          pattern = "(S(?=S|T|Y)|T(?=S|T|Y)|Y(?=S|T))+")) == 0) %>%
    rowwise() %>%
    mutate(STY_adj = list(sapply(Sequence, function(x) extract_pos_non_adj(query = x, patt_text = elements)))) %>%
    #mutate(STY_adj = paste0(sapply(Sequence, function(x) extract_pos(query = x,patt_text = sing_STY))))
    mutate(STY_len=length(STY_adj)) %>%
    mutate(is_adj= "non_adjacent")
    
  
  # comb_result_pep_adj_5 <- comb_result_seq %>% #comb_result_pep
  #   #select(Sequence) %>%
  #   filter(sapply(Sequence, function(x) str_count(x, pattern = big_pattern)) == 3) %>%
  #   mutate(STY_adj = sapply(Sequence, function(y) extract_pos(query = y,patt_text = perm5$comb))) %>%
  #   mutate(STY_len=nchar(STY_adj))
  # 
  # 
  # comb_result_pep_adj_4 <- comb_result_seq %>% #comb_result_pep
  #   #select(Sequence) %>%
  #   filter(sapply(Sequence, function(x) str_count(x, pattern = big_pattern)) == 2) %>%
  #   mutate(STY_adj = sapply(Sequence, function(y) extract_pos(query = y,patt_text = perm4$comb))) %>%
  #   mutate(STY_len=nchar(STY_adj))
  # 
  # comb_result_pep_adj_3 <- comb_result_seq %>%#comb_result_pep
  #   #select(Sequence) %>%
  #   filter(sapply(Sequence, function(x) str_count(x, pattern = big_pattern)) == 1) %>%
  #   mutate(STY_adj =sapply(Sequence, function(x) extract_pos(query = x,patt_text = perm3$comb))) %>%
  #   mutate(STY_len=nchar(STY_adj))
  #  
  #mutate(STY_adj =sapply(Sequence, function(x) extract_pos(query = x,patt_text = perm2$com0b)[1])) %>%
  
  ### FOR THE SAKE OF COMBINIG ALL ADJACENT OBJECTS ABOVE, ONE A.A WAS SELECT, IF THERE IS NO ADJACENT A.A
  # mutate(STY_adj_new =ifelse(sapply(STY_adj, function(x) any(length(x) == 0)), sapply(Sequence,function(x) extract_pos(query = x,
  #                                                                                                                   patt_text = sing_STY)[1]),NA)) %>%
  # #mutate(STY_pair =ifelse(sapply(STY_adj,function(x) length(x) > 1), sapply(STY_adj, function (x) length(x)),1)) %>%
  # ### TODO: ELIMINATE REPETATIVE CODE
  # mutate(STY_adj_new = ifelse(is.na(STY_adj_new),sapply(Sequence, function(x) extract_pos(query = x,patt_text = perm2$comb)[1]), STY_adj_new)) %>%
  
  comb_result_pep_non_adj$STY_adj_new <- sapply(comb_result_pep_non_adj$STY_adj, paste, collapse = "")
  
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
           ptm_score,
           pep_with_pos,
           pool_id,
           new_col,
           raw_file,is_adj) #%>%
    #mutate(is_adj= ifelse(STY_len == 1, "non-adjacent","adjacent")) %>%
    #tibble()
  
  # 
  correct_match_poolid_rawfile_adj_seq <- NULL

  for (i in 1:length(unique(comb_result_pep_adj_all$raw_file))){

    assign(paste0("pool",i), comb_result_pep_adj_all %>%
             rowwise() %>%
             filter(grepl(pool_list[i],pool_id)) %>%
             mutate(cor_pool=ifelse(raw_file == raw_list[i] && pool_id == pool_list[i],"Correct_match","Wrong match"))) #%>%
             #mutate(cor_pool=ifelse(is.na(pool_id[i]), "Unexpected",cor_pool)) %>%
             #filter(grepl("Correct_match",cor_pool) | grepl( "Unexpected",cor_pool)))
             #filter(grepl("Correct_match",cor_pool)))

    correct_match_poolid_rawfile_adj_seq <- rbind(correct_match_poolid_rawfile_adj_seq, get(paste0("pool",i)))
  }
  
  plot8 <- correct_match_poolid_rawfile_adj_seq %>% tibble() %>%
    group_by(cor_pool) %>%
    count(is_adj) %>%
    #mutate(n_new = ifelse(Pool_for_seq_merge!="Correct",(-1*n),n)) %>%
    ggplot(aes(x=is_adj,y=n, fill=cor_pool)) +
    geom_col() + #geom_text(aes(label=after_stat(count)),stat = "count", position=position_dodge(width=0.9), vjust=-0.25,size=10)+
    geom_text(aes(label = n), position = position_stack(vjust = 0.5),size=10) +
  #geom_text(aes(label = n), vjust = -0.5,size=10) +
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
  
  ###################################################
  
  # comb_pep_adj_dist_df_neg <- comb_result_pep %>% mutate(
  #   STY_adj = rowSums(sapply(possibilites, function(query) grepl(query, Sequence))) > 0,
  #   STY_len = rowSums(sapply(possibilites, function(query) str_count(Sequence, query)))) %>% 
  #   mutate(STY_len = ifelse(STY_adj == FALSE,(rowSums(sapply(sing_STY, function(query) str_count(Sequence, query)))*-1),STY_len))
  # #mutate(stat=ifelse(STY_count < 0,"negative","positive"))
  # 
  # comb_pep_adj_dist_df_pos <- comb_result_pep %>% mutate(
  #   STY_adj = rowSums(sapply(possibilites, function(query) grepl(query, Sequence))) > 0,
  #   STY_len = rowSums(sapply(possibilites, function(query) str_count(Sequence, query)))) %>% 
  #   mutate(STY_len = ifelse(STY_adj == FALSE,rowSums(sapply(sing_STY, function(query) str_count(Sequence, query))),STY_len))
  # 
  # comb_result_pep_adj <- comb_result_pep_adj_all %>%
  #   filter(grepl("non-adjacent",is_adj))
  # 
  # comb_result_pep_adj_and_nonadj <- comb_pep_adj_dist_df_neg %>% 
  #   filter(STY_adj==FALSE) %>%
  #   select(Sequence,
  #          STY_adj,
  #          STY_len,
  #          Pool_for_pep_merge,
  #          ptm_score,
  #          pep_with_pos,
  #          pool_id,
  #          new_col,
  #          raw_file) %>%
  #   mutate(STY_adj=ifelse(STY_adj == TRUE,"adjacent","non-adjacent")) %>%
  #   rename(is_adj = STY_adj) %>%
  #   mutate(STY_adj=NA) %>%
  #   bind_rows(comb_result_pep_adj) %>%
  #   filter(!grepl("missing",Pool_for_pep_merge))
  #   
  
  
  #mutate(stat=ifelse(STY_count < 0,"negative","positive"))
  ###############################################################
  
  #comb_result_pep_adj_all_corr_match <- NULL
  
  # for (i in 1:length(raw_list)){
  #   
  #   assign(paste0("pool",i), comb_result_pep_adj_all %>%
  #            rowwise() %>%
  #            mutate(cor_pool=ifelse(raw_file == raw_list[i] && pool_id == pool_list[i],"Correct_match","Wrong match")))#%>%
  #            #mutate(cor_pool=ifelse(is.na(pool_id[i]), "Unexpected",cor_pool)) %>%
  #            #filter(grepl("Correct_match",cor_pool) | grepl( "Unexpected",cor_pool)))
  #            #filter(grepl("Correct_match",cor_pool)))
  #   
  #   comb_result_pep_adj_all_corr_match <- rbind(comb_result_pep_adj_all_corr_match, get(paste0("pool",i)))
  #   
  # }
  #     
  # perm <-comb_result_pep_adj_all_corr_match %>%
  #   filter(grepl("Wrong match",cor_pool)) %>%
  #   group_by(pep_with_pos,raw_file) %>%
  #   distinct(pep_with_pos,.keep_all = TRUE)
                                                            
  
  ############################ ###########################
  #sum(sapply(possibilites, function(query) grepl(query,"SFNGSLTTKNVAVDELSR",fixed = F)))
  
  # comb_pep_adj_dist_df_pos <- comb_result_pep %>% mutate(
  #   Contains_STY = rowSums(sapply(possibilites, function(query) grepl(query, Sequence))) > 0,
  #   STY_count = rowSums(sapply(possibilites, function(query) str_count(Sequence, query)))) %>% 
  #   mutate(STY_count = ifelse(Contains_STY == FALSE,rowSums(sapply(sing_STY, function(query) str_count(Sequence, query))),STY_count))
  # #mutate(stat=ifelse(STY_count < 0,"negative","positive"))
  # 
  # 
  # gg_comb_pep_adj_dist_df <- comb_pep_adj_dist_df_neg %>% 
  #   group_by(Contains_STY) %>% 
  #   count(STY_count) %>%
  #   mutate(ifelse(Contains_STY==FALSE,(-1*(n)),n)) %>%
  #   mutate(ifelse(STY_count < 0, abs(STY_count),STY_count))
  # y=`ifelse(Contains_STY == FALSE, (-1 * (n)), n)`,
  #x=ifelse(STY_count < 0, abs(STY_count), STY_count)`,
  #fill=Contains_STY
  
  # plot7 <- ggplot(correct_match_poolid_rawfile_adj_seq, aes(x= STY_len,fill=is_adj)) +
  #   geom_bar() +
  #   geom_text(aes(label=after_stat(count)),stat = "count", position=position_dodge(width=0.9), vjust=-0.25,size=10)+
  #   labs(x = "Contains_STY", y = "STY_count",
  #        title = "Total number of adjacent and non-adjacent STY residues",caption = paste(software_name,"Experiment",exp_id,acquisiton_type,sep = " ")) +
  #   theme_minimal() +
  #   scale_fill_brewer(palette = "Set1")+
  #   theme(legend.text = element_text(size=30),
  #         axis.title.x = element_text(size=30),
  #         axis.title.y = element_text(size=30),
  #         plot.title = element_text(size=35),
  #         plot.subtitle = element_text(size = 25),
  #         plot.caption = element_text(size = 25),
  #         legend.title=element_text(size=30),
  #         axis.text.x = element_text(size=20),
  #         axis.text.y = element_blank(),
  #         #axis.text.y = element_text(size = 30),
  #         axis.title=element_text(size=30)) +
  #   scale_x_continuous(limits = c(-6,7),breaks = seq(from =-6, to = 7, by = 1))

  ##### ptmRS SCORE DISTRIBUTION BASED ON ADJACENT DISTANCES ##### 
  
  # comb_pep_adj_dist_prop <- comb_pep_adj_dist_df_pos %>%
  #   mutate(STY_prop = as.numeric(sapply(ptm_score, function(x) {
  #     str_extract(x, "\\d+(?:\\.\\d+)?(?=;|$)") #%>%
  #     #ifelse((x<10), str_extract(x, "\\d+\\.\\d{2,3}|\\d+"), .)
  #   })))
  #
  
  plot7 <- correct_match_poolid_rawfile_adj_seq %>%
    rowwise() %>%
    #distinct(Sequence,.keep_all = TRUE) %>%
    #group_by(pep_with_pos,new_col) %>%
    #count(pep_with_pos)
    #mutate(ptm_score = ifelse(Pool_for_pep_merge == "Wrong Localization",(ptm_score*-1),ptm_score)) %>%
    group_by(is_adj,cor_pool) %>%
    count(STY_len) %>% #Pool_for_pep_merge
    ungroup() %>%
    mutate(n_label=n) %>%
    #mutate(n=ifelse(Pool_for_pep_merge =="Wrong Localization", (n* (-1)),n)) %>%
    mutate(n=ifelse(is_adj =="non_adjacent", (n* (-1)),n)) %>%
    mutate(STY_len= ifelse(STY_len < 0, (-1*STY_len),STY_len)) %>%
    #mutate(abs_n=c("141","79","82","26\n")) %>%
    ggplot( aes(x= STY_len,y=n,fill=interaction(is_adj,cor_pool))) +
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

  
  plot9 <- correct_match_poolid_rawfile_adj_seq %>% 
    subset(!(STY_len == 1 & is_adj == "non_adjacent")) %>%
    
    ggplot( aes(x= ptm_score,fill=interaction(is_adj))) +
    geom_density(alpha=0.8) + labs(x = "Count of STY amino acids", y = "Total count",
                      title = "Comparison of having an adjacent a.a effect of localization accuracy",
                      caption = paste(software_name,"Experiment",exp_id,acquisiton_type,sep = " "),fill="Is peptide adjacent?") +
    facet_wrap(~cor_pool) +
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
  plot_obj <- ls(pattern="plot")
  plot_obj <- plot_obj[!is.na(plot_obj)]
  sapply(1:length(plot_obj),function(x) ggsave(filename = paste0("p",x,".tiff"),
                                               width = 60, height = 45, 
                                               path = paste0(file_path,curr_dir,"/outputs_with_new_script/"),
                                               units = "cm",
                                               get(plot_obj[x]),
                                               device = "tiff", #".svg"
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
