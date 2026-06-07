#library(PhosR)
library(stringr)
library(dplyr)
#library(data.table) # REMOVED: unused
library(openxlsx)
library(tidyr)
library(ggplot2)
library(ggpattern)
library(gtools)
###############################################

source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/ggplot/ggplot_functions.R")
#source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/roc_curve/new_roc_curve_generation_with_custom_threshold.R")

final_pd_pep_quant_analysis_exp1 <- function(file_path,
                                             #file_name,
                                             theo_file_path,
                                             theo_file_name,
                                             sheet_theo_name,
                                             selected_species,
                                             exp_id,
                                             #exp_design,
                                             background_species,
                                             mapping_file,
                                             acquisiton_type,
                                             software_name){
  
  intended_dir <-paste0(file_path,"/output_final_aft_mapp_func_mod_col_sel")#output_final")
  
  if(dir.exists(intended_dir)){
    new_path <- intended_dir
    
  }else{
    dir.create(intended_dir)
    new_path <- list.dirs(intended_dir)
    
  }
  
  map_df <- read.table(file = mapping_file,sep = "\t",header = T)
  
  # proline_phospho_pos_extraction <- function(df){
  #   # Extraction of phospho positions from quant peptides object
  #   phospho_ptm_pos <- lapply(df, function(each_ptm_protein_positions) {
  #     
  #     ptm_list <- as.list(strsplit(each_ptm_protein_positions,"; ", fixed=TRUE)[[1]]) # Split ptm_protein_position depending on ";"
  #     ptm_list <- ptm_list[grepl("Phospho", ptm_list, fixed = TRUE)] # Extract only which contains "Phospho"
  #     phospho_positions <- lapply(ptm_list, function(ptm) { 
  #       #sub('Phospho \\(([A-Z]\\d+)\\)', "\\d+", ptm) #then, remove "Phospho" and remain only positions
  #       as.character(str_extract(ptm, "\\d+"))
  #       
  #     })
  #     
  #     phospho_positions_as_str <- paste(phospho_positions, collapse="&") #combine each position with "|"
  #     
  #   })
  #   # Data conversion 
  #   phospho_ptm_pos_df <- t(as.data.frame(phospho_ptm_pos))
  #   rownames(phospho_ptm_pos_df) <- 1:length(phospho_ptm_pos_df)
  #   
  #   return(phospho_ptm_pos_df)
  #   
  # }
  ## EXTRACTION OF PTM SITE PROBS
  get_max_value <- function(input_string) {
    numeric_values <- max(as.numeric(str_extract(unlist(strsplit(unlist(strsplit(input_string, "; ")),":")),"\\d+")))
    max(as.numeric(numeric_values))
  }
  ###################
  
  extract_phospho_numbers <- function(input_string) {
    phospho_part <- regmatches(input_string, gregexpr("Phospho \\[[^]]+\\]", input_string))
    
    if (length(phospho_part) == 0 || length(phospho_part[[1]]) == 0) {
      return(list(NA, NA))
    }
    
    phospho_text <- phospho_part[[1]]
    letter_numbers <- gregexpr("[A-Z](\\d+)", phospho_text)
    extracted_letters <- regmatches(phospho_text, letter_numbers)[[1]]
    extracted_values <- gregexpr("\\((\\d+(?:\\.\\d+)?)\\)", phospho_text)
    extracted_values <- regmatches(phospho_text, extracted_values)[[1]]
    
    # Remove unwanted characters
    extracted_letters <- gsub("[A-Z]", "", extracted_letters)  # Remove letters
    extracted_values <- gsub("\\(|\\)", "", extracted_values)  # Remove parentheses
    
    # Handle empty cases correctly
    if (length(extracted_letters) == 0) {
      extracted_letters <- NA
    } else {
      extracted_letters <- paste(extracted_letters, collapse = "&")
    }
    
    if (length(extracted_values) == 0) {
      extracted_values <- NA
    } else {
      extracted_values <- paste(max(as.numeric(extracted_values)), collapse = "&")
    }
    
    return(list(extracted_letters, extracted_values))
  }
  
  
  # 
  # # Define the improved extraction function that handles decimal phospho scores
  # extract_phospho_info <- function(input_string) {
  #   
  #   # Step 1: Split the input string by semicolon ";"
  #   split_entries <- strsplit(input_string, ";")[[1]]
  #   
  #   # Initialize vectors to store positions and phospho scores
  #   positions <- c()
  #   phospho_scores <- c()
  #   
  #   # Step 2: Loop through each entry and extract information
  #   for (i in seq_along(split_entries)) {
  #     
  #     # Clean the entry by trimming any leading/trailing spaces
  #     entry <- trimws(split_entries[i])
  #     
  #     # Split the string by white space
  #     split_string <- strsplit(entry, " ")[[1]]
  #     
  #     # Extract the modification part for further analysis
  #     modification_part <- split_string[3]
  #     
  #     # Extract the position (number after 'S', 'T', or 'Y') and the phospho score
  #     position <- as.numeric(sub(".*[STY](\\d+).*", "\\1", modification_part))
  #     
  #     # Modify the regex to capture decimal numbers inside parentheses
  #     phospho_score <- as.numeric(sub(".*\\(([0-9]+\\.?[0-9]*)\\).*", "\\1", modification_part))
  #     
  #     # Append the extracted position and phospho score to respective vectors
  #     positions <- c(positions, position)
  #     phospho_scores <- c(phospho_scores, phospho_score)
  #   }
  #   
  #   # Combine positions using "&"
  #   combined_positions <- paste(positions, collapse = "&")
  #   
  #   # Get the highest phospho score
  #   max_phospho_score <- max(phospho_scores)
  #   
  #   # Return the combined positions and the highest phospho score
  #   return(list(
  #     positions = combined_positions,
  #     phospho_score = max_phospho_score
  #   ))
  # }
  
  all_dirs <- list.files(file_path,pattern = '.txt')
  
  #all_files <- list.files(paste0(file_path,all_dirs[i],"/"),pattern = ".xlsx")
  comb_result <- NULL
  comb_ecoli <- NULL
  
  for(i in 1:length(all_dirs)){
    assign(paste0("tmp"), read.table(file=paste0(file_path,all_dirs[i]),sep = "\t",header = T))
    
    if(grepl("PeptideIsoforms", all_dirs[i]) |grepl("PeptideGroups", all_dirs[i]) ){
      
      phospho_tmp <- tmp %>% 
        filter(grepl("Phospho",Modifications) & grepl(selected_species,Master.Protein.Descriptions)) %>%
        mutate(raw_files= all_dirs[i]) %>%
        separate(raw_files,into = c(paste0("tmp",1:4)),sep = "-") %>%
        select(!c(paste0("tmp",2:4))) %>%
        rename(raw_file=tmp1) %>%
        #filter(!grepl("positions not distinguishable",Modification.Pattern)) %>%
        rowwise() %>% # Ensure that the function is applied row by row
        mutate(results = list(extract_phospho_numbers(Modifications)),
               phospho_pos = results[[1]],
               phospho_score =results[[2]]) %>%
        select(!results) %>%
        
        ungroup() %>% # After rowwise, return to normal dataframe operations
        mutate(pep_with_pos = paste(Sequence,phospho_pos,sep = "_")) %>%
        drop_na(phospho_score)
      
    }else{
      
      phospho_tmp <- tmp %>% 
        filter(grepl("Phospho",Modifications) & grepl(selected_species,Master.Protein.Descriptions)) %>%
        mutate(raw_files=all_dirs[i]) %>%
        separate(raw_files,into = c(paste0("tmp",1:3)),sep = "-") %>%
        select(!c(paste0("tmp",2:4))) %>% 
        rename(raw_file=tmp1) %>%
        mutate(phospho_pos = proline_phospho_pos_extraction(Modifications)) %>%
        mutate(pep_with_pos = paste(Sequence,phospho_pos,sep = "_")) #%>%
      #distinct(pep_with_pos,.keep_all = TRUE)
    }
    
    assign(paste0(all_dirs[i]), phospho_tmp)
    assign(paste0("comb_result"),bind_rows(comb_result, phospho_tmp))
    
    if(background_species == "Escherichia coli"){
      ecoli_tmp <- tmp %>% 
        filter(grepl(background_species,Master.Protein.Descriptions)) %>%
        mutate(raw_files=all_dirs[i]) %>%
        separate(raw_files,into = c(paste0("tmp",1:3)),sep = "-") %>%
        rename(raw_file=tmp1) %>%
        mutate(acq_type=acquisiton_type) %>%
        mutate(soft_name=software_name) %>%
        distinct(Sequence,.keep_all = TRUE)
      
      assign(paste0("comb_ecoli"), bind_rows(comb_ecoli,ecoli_tmp))
    }else{
      
    }
    
    rm(tmp,phospho_tmp,ecoli_tmp)
    
  }
  
  
  if(background_species == "Escherichia coli"){
    
    comb_ecoli_dist <- comb_ecoli %>% 
      #mutate(raw_file=str_remove_all(Spectrum.File,".raw")) %>%
      full_join(map_df,by="raw_file") %>%
      #separate(`all_dirs[i]`,into = c("tmp","injs","tmp2"),sep = "_") %>%
      #select(!c(tmp,tmp2)) %>%
      #separate(injs,into = c("tmp","expid","sample_id","Ecoli","inj"),sep = "-") %>%
      #mutate(new_col=paste(expid,sample_id,Ecoli,sep = "_")) %>%
      ## THIS PART IS NEW!
      group_by(Sequence,raw_file) %>% 
      distinct(Sequence,.keep_all = TRUE)%>% 
      ungroup() 
    
    plot1 <- gg_barplt_id_pep_count(data_set = comb_ecoli_dist,
                                    x_df = comb_ecoli_dist$raw_file,
                                    fill_df = comb_ecoli_dist$pool_id,
                                    ymax = 20000,
                                    size_num = 8,
                                    header = paste("Total number of identified synthetic phospho-sequences across each pool",sep=" "),
                                    caption_lab = "Wrong Sequences were removed for the futher analysis.",
                                    x_lab = "Sample id",
                                    fill_lab =  "Sample id",
                                    y_lab = "Number of identified Sequence",
                                    subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) + 
      scale_fill_brewer(palette = "Dark2") + theme(axis.text.x = element_text(angle = 90))
    
    write.table(comb_ecoli_dist,file=paste0(new_path, "/",software_name, ##file_path
                                            "experiment",
                                            exp_id,
                                            "merge_id_ecoli_seq.tsv"),
                sep = "\t",col.names = T,row.names = F)
    
    
  }else{
    
  }
  
  ## READ THEO LIST
  pep_list_w_theo <- read.xlsx(paste0(theo_file_path, theo_file_name), sheet = sheet_theo_name)
  #pep_list_w_theo_quant <- pep_list_w_theo[,-1]
  #common_col_theo_quant <- as.data.frame(paste(pep_list_w_theo_quant$Phosphopeptide.sequence,
  #pep_list_w_theo_quant$modified.position.in.peptide, sep = "_"))
  
  ## ADD COMMON COLUMN TO MERGE WITH EXP. DATA
  #colnames(common_col_theo_quant) <- "pep_with_pos"
  #pep_list_w_theo_quant_new <- cbind(common_col_theo_quant,pep_list_w_theo_quant)
  
  # pep_list_w_theo_pep_comp <- pep_list_w_theo %>% 
  #   select(!c(Sequence,Phosphopeptide.sequence)) %>%
  #   #rename(Sequence = Phosphopeptide.sequence) %>% 
  #   rename(pool_id_theo_list=pool_id) 
  # 
  # pep_list_w_theo_seq_comp <- pep_list_w_theo %>% 
  #   select(Sequence,pool_id) %>%
  #   #select(!pep_with_pos) %>%
  #   #rename(sequence = Sequence) %>% 
  #   mutate(pool_id_theo_list=pool_id) 
  
  pep_list_w_theo_seq_map <- pep_list_w_theo %>%
    select(Sequence,Neutral.mass.SH.Cys) %>%
    rename(Neutral_mass=Neutral.mass.SH.Cys) %>%
    mutate(Neutral_mass = round(Neutral_mass, 0)) %>%
    mutate(Neutral_mass_theo=Neutral_mass) #%>%
  #mutate(pool_id_theo_list_seq=pool_id) 
  
  
  pep_list_w_theo_seq_map_pool <-  pep_list_w_theo %>%
    select(Sequence,pool_id) %>%
    mutate(pool_id_theo_list_seq=pool_id) 
  
  pep_list_w_theo_pep_map_pool <-  pep_list_w_theo %>%
    select(pep_with_pos,pool_id) %>%
    mutate(pool_id_theo_list_pep=pool_id) 
  
  if(grepl("PeptideIsoforms", all_dirs[i]) |grepl("PeptideGroups", all_dirs[i]) ){
    
    # Extracted values
    comb_result_pos <- comb_result %>% 
      select(Sequence, pep_with_pos,phospho_pos,phospho_score,raw_file,
             #First.Scan,Confidence,
             #Intensity, #Marked.as
             Protein.Accessions) %>%
      separate(pep_with_pos,into = c("seq","Positions"),sep = "_",remove = F) %>%
      #mutate(raw_file=str_remove_all(Spectrum.File,".raw")) %>%
      full_join(map_df,by="raw_file") %>%
      rename(Experiment=sample_name) %>%
      #separate(pool_id,into = c("expid","sample_id","Ecoli","inj"),sep = "_") %>%
      mutate(pool_id_map_df=pool_id)
  }else{
    # Extracted values
    comb_result_pos <- comb_result %>% 
      separate(pep_with_pos,into = c("seq","Positions"),sep = "_",remove = F) %>%
      select(Sequence, pep_with_pos,phospho_pos,ptmRS.Best.Site.Probabilities,
             First.Scan,raw_file,
             Confidence,Intensity, #Marked.as
             Protein.Accessions) %>%
      #mutate(raw_file=str_remove_all(Spectrum.File,".raw")) %>%
      full_join(map_df,by="raw_file") %>%
      rename(Experiment=sample_name) %>%
      #separate(pool_id,into = c("expid","sample_id","Ecoli","inj"),sep = "_") %>%
      mutate(pool_id_map_df=pool_id)
  }
  
  ##################################################################################################################
  ##################################################################################################################
  ### Mapping by precursor mass + phospho  and sequence ####
  
  amino_acid_table <- read.delim("D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/amino_acid_table.txt")
  mapping_df <- as.data.frame(cbind(amino_acid_table$X1.letter.code,amino_acid_table$Monoisotopic.Mass))
  colnames(mapping_df) <- c("letters","mono_isotopic")
  
  mapping_df$mono_isotopic <- as.numeric(mapping_df$mono_isotopic)
  # Apply the function to each sequence in the sequences data frame
  
  comb_result_pos$Sum <- mapply(calculate_sum,
                                comb_result_pos$Sequence,
                                comb_result_pos$Positions, MoreArgs = list(mapping_df = mapping_df))
  
  comb_result_mz_seq_map <- comb_result_pos %>%
    mutate(Neutral_mass = round(Sum, 0)) %>%
    mutate(Neutral_mass_res=Neutral_mass) %>%
    full_join(pep_list_w_theo_seq_map,by=c("Neutral_mass","Sequence"),relationship="many-to-many") %>% # #pep_with_pos
    mutate(map_seq = ifelse(is.na(Neutral_mass_res), "missing", NA)) %>%
    mutate(map_seq = ifelse(is.na(Neutral_mass_theo), "Subset or different mz", map_seq)) %>%
    mutate(map_seq = ifelse(!is.na(Neutral_mass_res) & !is.na(Neutral_mass_theo) & Neutral_mass_res == Neutral_mass_theo , # pool_id_map_df==pool_id_theo_list_seq
                            "Correct seq. & Correct mz",map_seq)) %>%
    select(!c(Positions,Sum)) #%>%
  ####### THESE LAST THREE ROWS ARE COMMENTED BECAUSE THERE IS NO INTENSITY INFO
  #group_by(pep_with_pos,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
  #slice(which.max(Intensity)) %>%
  # ungroup() 
  
  
  comb_result_seq_pool <- comb_result_mz_seq_map %>% 
    full_join(pep_list_w_theo_seq_map_pool,by=c("Sequence","pool_id")) %>% # #pep_with_pos
    mutate(map_seq_pool=ifelse(map_seq!='Subset or different mz' & pool_id_theo_list_seq==pool_id_map_df, "Correct Seq. & mz & pool",NA)) %>%
    mutate(map_seq_pool=ifelse(map_seq!='Subset or different mz' & is.na(pool_id_theo_list_seq), "Correct Seq. & mz but wrong pool",map_seq_pool)) %>%
    
    mutate(map_seq_pool=ifelse(map_seq== "Subset or different mz",map_seq,map_seq_pool)) %>%
    mutate(map_seq_pool = ifelse(is.na(Experiment), "missing", map_seq_pool))
  
  
  # #### TOTAL NUM. OF PHOSPHO-SEQUENCES ####
  # ### Correct mapping was done using "map_df".
  # 
  #   comb_result_seq <- comb_result_pos %>%
  #   full_join(pep_list_w_theo_seq_comp,by=c("Sequence","pool_id")) %>%
  #   select(!c(pool_id,sample_name)) %>%
  #   ## RE-JOINING TO COUNT CORRECT, INCORRECT and MISSING SEQ. 
  #   full_join(map_df,by="raw_file") %>%
  #   rename(pool_id_map_df=pool_id) %>%
  #   mutate(Pool_for_seq_merge=ifelse(pool_id_theo_list==pool_id_map_df,"Correct Seq.","Wrong Seq.")) %>% #within theo list
  #   #mutate_at("Pool_for_seq_merge", ~replace_na(.,"Unexpected Seq.")) %>%
  #   mutate(Pool_for_seq_merge= ifelse(is.na(pool_id_map_df),"Missing",Pool_for_seq_merge)) %>%
  #   mutate(Pool_for_seq_merge=ifelse(is.na(pool_id_theo_list),"Wrong Seq.",Pool_for_seq_merge)) %>% # out of theo. list
  #   mutate(acq_type=acquisiton_type) %>%
  #   mutate(soft_name=software_name)
  # 
  #### WRITE THE OBJECT AS TSV ####
  #### WRITE THE OBJECT AS TSV ####
  
  
  ### Duplicate sequences were removed.
  comb_result_dist <- comb_result_seq_pool %>%
    #relocate(pool_id_map_df,pool_id_theo_list,Pool_for_seq_merge,.after = Sequence) %>% #
    group_by(raw_file,Sequence,Neutral_mass_res) %>%
    distinct(Sequence,.keep_all = T) %>%
    ungroup() 
  
  write.table(comb_result_seq_pool,file=paste0(new_path,"/", #file_path,curr_dir
                                               "exp",
                                               exp_id,#acquisiton_type,
                                               "merge_theo_list_with_identified_phospho_seq.tsv"),
              sep = "\t",col.names = T,row.names = F)
  
  write.table(comb_result_dist,file=paste0(new_path,"/", software_name, #file_path,curr_dir
                                           "experiment",
                                           exp_id,#acquisiton_type,
                                           "merge_theo_list_with_identified_phospho_seq_unique_ones.tsv"),
              sep = "\t",col.names = T,row.names = F)
  
  #### VISUALIZATION OF TOTAL NUM. OF PHOSPHO-SEQ ####
  plot2 <-   gg_barplt_id_pep_count_stack(data_set =comb_result_dist,#comb_result_seq_dist,
                                          x_df =comb_result_dist$Experiment,#comb_result_seq_dist$raw_file,
                                          fill_df = comb_result_dist$map_seq_pool,#comb_result_seq_dist$Pool_for_seq_merge ,
                                          ymax = 250,
                                          size_num = 10,
                                          header = paste("Total number of identified synthetic phospho-sequences across each pool",sep=" "),
                                          caption_lab = "Duplicates were removed for the futher analysis.",
                                          x_lab = "Sample id",
                                          fill_lab =  "Sample id",
                                          y_lab = "Number of identified Sequence",
                                          subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) +
    scale_fill_brewer(palette = "Paired") + theme(axis.text.x = element_text(angle = 90))
  
  #### EXTRACTION OF CORRECT SEQ. ####
  comb_result_seq_dist_cor <- comb_result_dist %>%
    #filter(!grepl("Wrong",Pool_for_seq_merge) & !grepl("Missing",Pool_for_seq_merge))
    filter(!grepl("Subset or different mz",map_seq) & !grepl("Missing",map_seq))
  
  #### VISUALIZATION OF TOTAL NUM. OF CORRECTLY IDENTIFIED PHOSPHO-SEQ ####
  plot3 <- gg_barplt_id_pep_count_stack(data_set = comb_result_seq_dist_cor,
                                        x_df = comb_result_seq_dist_cor$Experiment,
                                        fill_df = comb_result_seq_dist_cor$map_seq_pool,
                                        ymax = 250,
                                        size_num=10,
                                        header = paste("Total number of identified synthetic phospho-sequences across each pool",sep=" "),
                                        caption_lab = "Wrong Sequences were removed for the futher analysis.",
                                        x_lab = "Sample id",
                                        fill_lab =  "Sample id",
                                        y_lab = "Number of identified Sequence",
                                        subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) +
    
    scale_fill_brewer(palette = "Dark2") + theme(axis.text.x = element_text(angle = 90))
  
  
  #### TOTAL NUM. OF PHOSPHO-PEPTIDES ####
  ### Only correct sequences were kept.
  ### Correct mapping was done using "pool_id_map_df" and "pep_with_pos".
  ### Duplicate peptides were removed.
  if(grepl("PeptideIsoforms", all_dirs[i]) |grepl("PeptideGroups", all_dirs[i]) ){
    
    
    comb_result_pep_poolwise <- comb_result_seq_pool %>% 
      filter(!grepl("Subset or different mz",map_seq) & !grepl("Missing",map_seq))%>%
      filter(!grepl('Correct Seq. & mz but wrong pool',map_seq_pool)) %>%
      #filter(!grepl("Wrong",Pool_for_seq_merge) & !grepl("Missing",Pool_for_seq_merge)) %>%
      left_join(pep_list_w_theo_pep_map_pool ,by=c("pep_with_pos","pool_id")) %>%
      
      mutate(map_pep=ifelse(pool_id_map_df==pool_id_theo_list_pep,"Correct Seq. & Correct Loc.",NA)) %>%
      mutate(map_pep=ifelse(map_seq_pool=="Correct Seq. & mz & pool"& is.na(pool_id_theo_list_pep),"Correct Seq. & Wrong Loc.",map_pep)) %>%
      mutate(map_pep=ifelse(map_seq_pool=="missing","missing",map_pep)) %>%
      rename(max_phospho=phospho_score) %>%
      group_by(pep_with_pos,raw_file) %>%
      distinct(pep_with_pos,.keep_all = TRUE)%>%
      ungroup()
    
    #mutate(map_pep=ifelse(is.na(pool_id_theo_list_pep),"Correct Seq. & Wrong Loc.",map_seq_pool)) %>%
    #select(!c(Positions,Sum)) %>%
    
    ####### THESE LAST THREE ROWS ARE COMMENTED BECAUSE THERE IS NO INTENSITY INFO
    # group_by(pep_with_pos,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    # #slice(which.max(Intensity)) %>%
    #ungroup() 
    
    # comb_result_pep <- comb_result_seq %>% 
    #   rename(max_phospho=phospho_score) %>%
    #   #mutate(max_phospho=sapply(`ptmRS.Best.Site.Probabilities`,function(row) get_max_value(row))) %>%
    #   #relocate(max_phospho,.after = ptmRS.Best.Site.Probabilities) %>%
    #   filter(!grepl("Wrong",Pool_for_seq_merge) & !grepl("Missing",Pool_for_seq_merge)) %>%
    #   select(!c(Pool_for_seq_merge,pool_id_theo_list,pool_id_map_df,sample_name)) %>%
    #   full_join(map_df,by="raw_file") %>%
    #   rename(pool_id_map_df=pool_id) %>%
    #   full_join(pep_list_w_theo_pep_comp,by="pep_with_pos") %>% ## IF FULL_JOIN IS USED,
    #   group_by(pep_with_pos,raw_file) %>%
    #   distinct(pep_with_pos,.keep_all = TRUE)%>%
    #   ungroup() %>%
    #   mutate(map_pep=ifelse(pool_id_theo_list==pool_id_map_df,"Correct","Wrong Loc.within theo list")) %>%
    #   mutate(map_pep= ifelse(is.na(pool_id_map_df),"Missing",map_pep)) %>%
    #   mutate(map_pep=ifelse(is.na(pool_id_theo_list),"Wrong Loc. out of theo. list",map_pep)) %>%
    #   relocate(pool_id_map_df,pool_id_theo_list,map_pep,.after = Sequence) %>% 
    #   mutate(acq_type=acquisiton_type) %>%
    #   mutate(soft_name=software_name)
    
  }else{
    comb_result_pep_poolwise <- comb_result_seq_pool %>% 
      mutate(max_phospho=sapply(`ptmRS.Best.Site.Probabilities`,function(row) get_max_value(row))) %>%
      relocate(max_phospho,.after = ptmRS.Best.Site.Probabilities) %>%
      # filter(!grepl("Wrong",Pool_for_seq_merge) & !grepl("Missing",Pool_for_seq_merge)) %>%
      # select(!c(Pool_for_seq_merge,pool_id_theo_list,pool_id_map_df,sample_name)) %>%
      # full_join(map_df,by="raw_file") %>%
      # rename(pool_id_map_df=pool_id) %>%
      # full_join(pep_list_w_theo_pep_comp,by="pep_with_pos") %>% ## IF FULL_JOIN IS USED,
      # group_by(pep_with_pos,raw_file) %>%
      # distinct(pep_with_pos,.keep_all = TRUE)%>%
      # ungroup() %>%
      # mutate(map_pep=ifelse(pool_id_theo_list==pool_id_map_df,"Correct","Wrong Loc.within theo list")) %>%
      # mutate(map_pep= ifelse(is.na(pool_id_map_df),"Missing",map_pep)) %>%
      # mutate(map_pep=ifelse(is.na(pool_id_theo_list),"Wrong Loc. out of theo. list",map_pep)) %>%
      # relocate(pool_id_map_df,pool_id_theo_list,map_pep,.after = Sequence) %>% 
      # mutate(acq_type=acquisiton_type) %>%
      # mutate(soft_name=software_name)
      filter(!grepl("Subset or different mz",map_seq) & !grepl("Missing",map_seq))%>%
      filter(!grepl('Correct Seq. & mz but wrong pool',map_seq_pool)) %>%
      #filter(!grepl("Wrong",Pool_for_seq_merge) & !grepl("Missing",Pool_for_seq_merge)) %>%
      left_join(pep_list_w_theo_pep_map_pool ,by=c("pep_with_pos","pool_id")) %>%
      mutate(map_pep=ifelse(pool_id_map_df==pool_id_theo_list_pep,"Correct Seq. & Correct Loc.",NA)) %>%
      mutate(map_pep=ifelse(map_seq_pool=="Correct Seq. & mz & pool"& is.na(pool_id_theo_list_pep),"Correct Seq. & Wrong Loc.",map_pep)) %>%
      mutate(map_pep=ifelse(map_seq_pool=="missing","missing",map_pep)) %>%
      group_by(pep_with_pos,raw_file) %>%
      distinct(pep_with_pos,.keep_all = TRUE)%>%
      ungroup()
  }
  
  
  
  # filter(!grepl("Unexpected Seq.",Pool_for_seq_merge) & !grepl("Missing Seq.",Pool_for_seq_merge)) %>%
  # select(!c(Pool_for_seq_merge, pool_id_theo_list)) %>%
  # full_join(pep_list_w_theo_quant_new,by="pep_with_pos") %>% ## IF FULL_JOIN IS USED,
  
  # mutate(map_pep=ifelse(pool_id_theo_list==pool_id_map_df,"Correct","Wrong Loc.within theo list")) %>%
  # mutate(map_pep= ifelse(is.na(pool_id_map_df),"Missing",map_pep)) %>%
  # mutate(map_pep=ifelse(is.na(pool_id_theo_list),"Wrong Loc. out of theo. list",map_pep)) %>%
  # mutate(acq_type=acquisiton_type) %>%
  # mutate(soft_name=software_name) 
  # 
  write.table(comb_result_pep_poolwise ,file=paste0(new_path,"/merge_theo_list_id_phospho_sites_only_unique_ones.tsv"), #file_path,curr_dir
              sep = "\t",col.names = T,row.names = F)
  
  #### VISUALIZATION OF TOTAL NUM. OF IDENTIFIED & LOCALIZED PHOSPHO-PEP ####
  plot4 <- gg_barplt_id_pep_count_stack(data_set =comb_result_pep_poolwise,
                                        x_df =comb_result_pep_poolwise$Experiment,
                                        fill_df = comb_result_pep_poolwise$map_pep,
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
  comb_result_pep_cor <- comb_result_pep_poolwise %>%
    #filter(grepl("Correct",map_pep))
    filter(grepl("Correct Seq. & mz & pool",map_seq_pool))
  
  #### VISUALIZATION OF TOTAL NUM. OF CORRECTLY IDENTIFIED & LOCALIZED PHOSPHO-PEP ####
  plot5 <- gg_barplt_id_pep_count_stack(data_set =comb_result_pep_cor,
                                        x_df =comb_result_pep_cor$Experiment,
                                        fill_df = comb_result_pep_cor$map_pep,
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
  
  
  if(background_species == "Escherichia coli"){
    
    comb_result_pep_filtered <- comb_result_pep_poolwise %>% 
      filter(!grepl("missing",map_pep)) %>%
      separate(Experiment,into = c("exp","samp","coli","inj"),sep = "_",remove = F) %>%
      mutate(samp_name_wo_inj=paste(exp,samp,coli,sep = "-")) %>%
      group_by(samp_name_wo_inj,pep_with_pos) %>%
      slice(which.max(max_phospho)) %>%
      ungroup()
    
  }else if (background_species ==""){
    
    comb_result_pep_filtered <- comb_result_pep_poolwise %>% 
      filter(!grepl("missing",map_pep)) %>%
      separate(Experiment,into = c("exp","samp","inj"),sep = "_",remove = F) %>%
      mutate(samp_name_wo_inj=paste(exp,samp,sep = "-")) %>%
      group_by(samp_name_wo_inj,pep_with_pos) %>%
      slice(which.max(max_phospho)) %>%
      ungroup()
    
  }else{}
  
  write.table(comb_result_pep_filtered ,file=paste0(new_path,"/merge_theo_list_id_phospho_sites_max_ptm_score.tsv"), #file_path,curr_dir
              sep = "\t",col.names = T,row.names = F)
  
  #### LOCALIZATION ACCURACY ASSESSMENT (ROC LIKE PLOT GENERATION) #### 
  ### Column selection
  comb_result_pep_rmv_miss <- comb_result_pep_filtered %>% 
    select(pep_with_pos,map_pep,max_phospho,pool_id_map_df)
  
  
  ### Custom localization threshold determination and
  ### Counting total number of correct and wrong localization at a given threshold
  threshold <- 1:100
  final_df <- NULL
  
  for (i in 1:length(threshold)){
    
    df <- comb_result_pep_rmv_miss %>% subset(max_phospho > threshold[i]) %>% 
      count(map_pep) %>% 
      mutate(threshold_val = threshold[i])
    
    final_df <-bind_rows(final_df,df)
    rm(df)
  }
  
  write.table(final_df, file=paste0(new_path,"/","Experiment",exp_id,software_name,"_num_sites_with_scores.tsv"),sep = "\t",col.names = T,row.names = F) #file_path,curr_dir
  
  
  #wrong_df <- final_df %>% filter(!grepl("Correct",map_pep))
  ## Threshold values were divided by 100 to make all thresholds the same range
  correct_df <- final_df %>% filter(grepl("Correct Seq. & Correct Loc",map_pep)) %>% #filter(grepl("Correct",map_pep)) %>% 
    mutate(new_threshold_val=threshold_val/100)
  
  ### Visualization
  plot6 <- ggplot(correct_df,aes(y=n, x=new_threshold_val)) + geom_line(size=2) +
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
  pep_list_w_theo_adj_map <- pep_list_w_theo %>% select(pep_with_pos,S_count,T_count,Y_count,row_sum)
  ### COMBINATION OF ALL ADJACENT ONES
  comb_result_pep_adj_all <- comb_result_pep_non_adj %>%
    select(!c(STY_adj)) %>%
    rename(STY_adj = STY_adj_new) %>%
    bind_rows(comb_result_pep_adj_5,
              comb_result_pep_adj_4,
              comb_result_pep_adj_3,
              comb_result_pep_adj_2) %>%
    left_join(pep_list_w_theo_adj_map,by = "pep_with_pos") %>%
    select(Sequence,
           S_count,
           T_count,
           Y_count,
           row_sum,
           STY_adj,
           STY_len,
           map_pep,
           map_seq_pool,
           max_phospho,
           pep_with_pos,
           pool_id_map_df,
           #sample_name,
           raw_file,is_adj) 
  
  write.table(comb_result_pep_adj_all,file = paste0(new_path,"/all_adj&nonadj_corr_wrong_loc_nolocthreshold.txt"),sep = "\t",col.names = T,row.names = F)
  
  ### TOTAL COUNT OF ALL SEQUENCES WITHOUT CONSIDERING ADJ_COUNT
  plot7 <- comb_result_pep_adj_all %>% tibble() %>%
    filter(grepl("Correct Seq. & Correct Loc",map_pep)) %>% ##filter(grepl("Correct",map_pep)) %>%
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
    scale_fill_manual(values = c("#a1d76a", "#e9a3c9"))+ 
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
  
  #scale_x_continuous(limits = c(0,7),breaks = seq(from =0, to = 7, by = 1))
  
  ### TOTAL COUNT OF ALL SEQUENCES WITH CONSIDERING ADJ_COUNT AND ACCURACY
  plot8 <- comb_result_pep_adj_all %>%
    mutate(pep_class=map_pep) %>%
    #mutate(pep_class=ifelse(map_pep!= "Correct","Wrong",pep_class)) %>%
    #filter(grepl("Correct",map_pep)) %>%
    rowwise() %>%
    #distinct(Sequence,.keep_all = TRUE) %>%
    #group_by(pep_with_pos,new_col) %>%
    #count(pep_with_pos)
    #mutate(ptm_score = ifelse(map_pep == "Wrong Localization",(ptm_score*-1),ptm_score)) %>%
    group_by(is_adj,pep_class) %>%
    count(STY_len) %>% #map_pep
    ungroup() %>%
    mutate(n_label=n) %>%
    #mutate(n=ifelse(map_pep =="Wrong Localization", (n* (-1)),n)) %>%
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
  
  plot11 <- comb_result_pep_adj_all %>%
    mutate(pep_class=map_pep) %>%
    # mutate(map_pep=ifelse(map_pep!="Correct","Wrong",map_pep)) %>%
    #filter(!grepl("Wrong Loc. out of theo. list",map_pep)) %>%
    #filter(ptm_score > 0.75) %>%
    count(row_sum,is_adj,map_pep) %>% mutate(n_label=n) %>%
    mutate(n=ifelse(is_adj =="non_adjacent", (n* (-1)),n)) %>%
    mutate(STY_len= ifelse(row_sum < 0, (-1*row_sum),row_sum)) %>%
    ggplot(aes(x=STY_len,y=n,fill=map_pep,pattern=is_adj)) +
    geom_col_pattern(alpha=0.8,
                     color = "white", 
                     pattern_fill = "white",
                     pattern_angle = 45,
                     pattern_density = 0.1,
                     pattern_spacing = 0.025,
                     pattern_key_scale_factor = 0.6) +
    scale_pattern_manual(values = c(adjacent = "stripe", non_adjacent = "none")) +
    geom_text(aes(label =n_label), position = position_stack(vjust = 0.5),size=5) +
    scale_fill_manual(values = c("#a1d76a", "#e9a3c9")) + 
    theme_minimal()+
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
          axis.title=element_text(size=30)) +
    ggtitle("Count of peptides having adjacent STYs or not as a function of total number of STY amino acids") +
    labs(subtitle = "No Localization filter was applied.",x="Length of STY a.a in the sequence",y="Count",fill="Is correctly localized?")
  
  plot12 <- comb_result_pep_adj_all %>% 
    filter(!grepl("non_adjacent",is_adj)) %>%
    mutate(pep_class=map_pep) %>%
    #mutate(map_pep=ifelse(map_pep!="Correct","Wrong",map_pep)) %>%
    #filter(!grepl("Wrong Loc. out of theo. list",map_pep)) %>%
    #filter(ptm_score > 0.75) %>%
    count(STY_len,is_adj,map_pep) %>% mutate(n_label=n) %>%
    mutate(n=ifelse(map_pep =="Wrong", (n* (-1)),n)) %>%
    mutate(STY_len= ifelse(STY_len < 0, (-1*STY_len),STY_len)) %>%
    ggplot(aes(x=STY_len,y=n,fill=map_pep)) +
    geom_col() +
    geom_text(aes(label =n_label), position = position_stack(vjust = 0.5),size=5) +
    scale_fill_manual(values = c("#a1d76a", "#e9a3c9")) + 
    theme_minimal()+
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
          axis.title=element_text(size=30)) + 
    ggtitle("Count of peptides having adjacent STYs as a function of length of the STY pattern") +
    labs(subtitle = "No Localization filter was applied.",x="Length of STY pattern",y="Count",fill="Is correctly localized?")
  
  ### DISTRIBUTION OF ALL SEQUENCES WITH CONSIDERING ADJ_COUNT AND ACCURACY 
  plot9 <- comb_result_pep_adj_all %>%
    mutate(pep_class=map_pep) %>%
    #mutate(pep_class=ifelse(map_pep!= "Correct","Wrong",pep_class)) %>%
    subset(!(STY_len == 1 & is_adj == "non_adjacent")) %>%
    
    ggplot( aes(x= max_phospho,fill=interaction(is_adj))) +
    geom_bar(stat = "count") + labs(x = "Count of STY amino acids", y = "Total count",
                                    title = "Comparison of having an adjacent a.a effect of localization accuracy",
                                    caption = paste(software_name,"Experiment",exp_id,acquisiton_type,sep = " "),fill="Is peptide adjacent?") +
    facet_wrap(~pep_class) +
    theme_minimal() +
    #scale_fill_brewer(palette = "Paired")+
    scale_fill_manual(values = c("#a1d76a", "#e9a3c9"))+ 
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
  plt_obj <- ls(pattern="plot")
  plt_obj <- plt_obj[!is.na(plt_obj)]
  sapply(1:length(plt_obj),function(x) ggsave(filename = paste0("p",x,".png"),
                                              width = 60, height = 45, 
                                              path = paste0(file_path,"/output_final_aft_mapp_func_mod_col_sel/"),#output_final #outputs_after_mapping_change
                                              units = "cm",
                                              get(plt_obj[x]),
                                              device = "png", #".svg"
  ))
  
  
  
}
