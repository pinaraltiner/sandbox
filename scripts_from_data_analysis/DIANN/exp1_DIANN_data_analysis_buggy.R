#library(PhosR)
library(stringr)
library(dplyr)
library(tibble)
library(tibble)
#library(data.table) # REMOVED: unused
library(openxlsx)
library(tidyr)
library(ggplot2)
#library(tidyverse) # REMOVED: replaced by individual package imports already present
library(purrr)
library(gtools)
library(ggpattern)
library(nanoparquet)
###############################################
# CONTAINER_DISABLED: "D:/dev/Pinar/PHD/wet_lab_experiments/DIA_analysis/Exp_1/Exp1_with_Ecoli/DIANN_v1.9/noMBR/"
# CONTAINER_DISABLED: 
source("/home/claude/post_processing_parquet/additional_functions/ggplot_functions.R")
source("/home/claude/post_processing_parquet/additional_functions/get_modification_func/getModificationPosition_func_for_all_mods_fixed.R")
# CONTAINER_DISABLED: 
# CONTAINER_DISABLED: file_paths <- paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DIA_analysis/Exp_1/Exp1_with_Ecoli/DIANN_v1.9/noMBR/",
# CONTAINER_DISABLED:                c("using_PD_Exp1_DDA_speclib_MaxVarMod4/",
# CONTAINER_DISABLED:                  "misclv2_MaxVarMod4_charge14_FDR1/",
# CONTAINER_DISABLED:                  "maxVarMod3_wrong_speclib_naming/misclv2_MaxVarMod3_charge14_FDR100/",
# CONTAINER_DISABLED:                  "misclv2_MaxVarMod4_charge14_FDR100/",
# CONTAINER_DISABLED:                  "misclv2_MaxVarMod5_charge14_FDR100/","run_dir_comb_fasta_msclv2_MaxVarMod3Charge14/"))
# CONTAINER_DISABLED: #"D:/dev/Pinar/PHD/wet_lab_experiments/DIA_analysis/Exp_1/DIANN/Exp1_with_Ecoli/DIANN_v1.9/noMBR/misclv2_MaxVarMod5_charge14/")
# CONTAINER_DISABLED: #"D:/dev/Pinar/PHD/wet_lab_experiments/DIA_analysis/Exp_1/DIANN/Exp1_with_Ecoli/DIANN_v1.9/withMBR/") #noMBR
# CONTAINER_DISABLED: file_names <- c("report.parquet","misclv2_MaxVarMod4_charge14_FDR1_report.parquet",
# CONTAINER_DISABLED:                 "misclv2_MaxVarMod3_charge14_FDR100report.parquet",
# CONTAINER_DISABLED:                 "report.parquet",
# CONTAINER_DISABLED:                 "Exp1_with_ecoli_misclv2_MaxVarMod5_charge14_FDR100.parquet",
# CONTAINER_DISABLED:                 "comb_fasta_msclv2_MaxVarMod3Charge14_report.parquet")#) "Exp1_with_ecoli_misclv2_MaxVarMod5_charge14_FDR100.parquet"
# CONTAINER_DISABLED: #"misclv2_MaxVarMod5_charge14_Exploris_report.tsv")#"same_params_and_lib_exp2_redesigned_MaxVarMods4_FDR1_Misclgvg2_MBR.tsv")
# CONTAINER_DISABLED: #"same_params_and_lib_exp2_redesigned_MaxVarMods4_FDR1_Misclgvg2_noMBR.tsv")
# CONTAINER_DISABLED: mapping <- "D:/dev/Pinar/PHD/data_analysis/DIA_mapping/exp1_batch2_mapping_btw_rawfile_pool_id.txt"
# CONTAINER_DISABLED: 
# CONTAINER_DISABLED: #"D:\dev\Pinar\PHD\wet_lab_experiments\DIA_data_analysis\experiment_1\DIANN\Exp1_with_Ecoli"
# CONTAINER_DISABLED: 
# CONTAINER_DISABLED: selected_spcies="HUMAN"
# CONTAINER_DISABLED: background_species= "ECOLI"
# CONTAINER_DISABLED: theo_file_path="D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/"
# CONTAINER_DISABLED: theo_file_name="Synthetic peptides list_theo_conc_corrected_isomericity_new_with_plates.xlsx"
# CONTAINER_DISABLED: sheet_theo_name ="ISOREF_REF2_Others" #"ISO-refOTHER with FC_correct"
# CONTAINER_DISABLED: acquisiton_type="DIA no FAIMS Exploris"
# CONTAINER_DISABLED: subtitle = #"version 1.9 & no MBR charge 1 & 4 MaxVarMod 5 "
# CONTAINER_DISABLED: #fdr_threshold = 0.05
# CONTAINER_DISABLED: 
# CONTAINER_DISABLED: #exp_design=experiment_name
# CONTAINER_DISABLED: exp_id=1
# CONTAINER_DISABLED: software_name="DIA-NN"
# CONTAINER_DISABLED: num_reps=3
# CONTAINER_DISABLED: #test_type="limma"
# CONTAINER_DISABLED: file_path <- file_paths[6]
# CONTAINER_DISABLED: file_name <- file_names[6]
# CONTAINER_DISABLED: 
# CONTAINER_DISABLED: final_diann_pep_quant_analysis_syn_exp1(file_path = file_path,
# CONTAINER_DISABLED:                                         file_name = file_name,
# CONTAINER_DISABLED:                                         sheet_name="",
# CONTAINER_DISABLED:                                         theo_file_path=theo_file_path,
# CONTAINER_DISABLED:                                         theo_file_name=theo_file_name,
# CONTAINER_DISABLED:                                         sheet_theo_name=sheet_theo_name,
# CONTAINER_DISABLED:                                         background_species="ECOLI",
# CONTAINER_DISABLED:                                         selected_spcies="HUMAN",
# CONTAINER_DISABLED:                                         exp_id=1,
# CONTAINER_DISABLED:                                         mapping=mapping,
# CONTAINER_DISABLED:                                         #exp_design,
# CONTAINER_DISABLED:                                         #fdr_threshold,
# CONTAINER_DISABLED:                                         acquisiton_type=acquisiton_type,
# CONTAINER_DISABLED:                                         software_name="DIA-NN",
# CONTAINER_DISABLED:                                         #test_type,
# CONTAINER_DISABLED:                                         num_reps=3,
# CONTAINER_DISABLED:                                         #actual_ratio,
# CONTAINER_DISABLED:                                         subtitle="version 1.9 & no MBR charge 1 & 4 MaxVarMod 3 FDR 1% using comb. FASTA")


final_diann_pep_quant_analysis_syn_exp1_fixed <- function(file_path,
                                                          file_name,
                                                          sheet_name,
                                                          theo_file_path,
                                                          theo_file_name,
                                                          sheet_theo_name,
                                                          background_species,
                                                          selected_spcies,
                                                          exp_id,
                                                          mapping,
                                                          #exp_design,
                                                          #fdr_threshold,
                                                          acquisiton_type,
                                                          software_name,
                                                          #test_type,
                                                          num_reps,
                                                          #actual_ratio,
                                                          subtitle){
  
  intended_dir <-paste0(file_path,"output_final_aft_mapp_func/")
  
  if(dir.exists(intended_dir)){
    new_path <- intended_dir
    
  }else{
    dir.create(intended_dir)
    new_path <- list.dirs(intended_dir)
    
  }
  
  mapping_file <- read.table(mapping,sep = "\t",header = T)
  exp_design <- mapping_file$sample_name
  
  sample_size <- length(exp_design) / num_reps
  sample_names <- paste0("M",1:sample_size)
  # comparisons <- NULL
  # for (i in 1:sample_size){
  #   tmp <- paste0(sample_names[1], "/",sample_names[i])
  #   comparisons[i] <- tmp
  #   rm(tmp)
  # }
  # comparisons <- comparisons[-1]
  # 
  quant_peptides <- as.data.frame(read_parquet(paste0(file_path,file_name)))
  
  ## Version-aware intensity column: v1.9 uses PG.Normalised, v2.3.2 uses Precursor.Normalised
  if (!"PG.Normalised" %in% colnames(quant_peptides)) {
    quant_peptides$PG.Normalised <- quant_peptides$Precursor.Normalised
  }
  
  quant_peptides_with_cond <- quant_peptides %>% 
    rename(raw_file=Run) %>%
    left_join(mapping_file,by="raw_file") %>%
    mutate(pool_id_map_df=pool_id) %>%
    rename("Sequence"="Stripped.Sequence") %>%
    rename("Intensity"= "Precursor.Normalised")# %>%
  #filter(Global.Q.Value < 0.05)
  
  ## PHOSPHO-FILTERING
  quant_phospho <- quant_peptides_with_cond %>% 
    filter(grepl("HUMAN", Protein.Names) & !grepl("CON__",Protein.Names)) %>%
    filter(grepl("UniMod:21",Modified.Sequence)) #%>%
  #filter(PTM.Q.Value < 0.05 )
  
  ## THEORETICAL PEPTIDE LIST
  pep_list_w_theo <- read.xlsx(paste0(theo_file_path, theo_file_name), sheet = sheet_theo_name)
  
  # ##### OLD AMPPING SYTLE #####
  
  # pep_list_w_theo_pep_comp <- pep_list_w_theo %>% 
  #   select(!c(Sequence,Phosphopeptide.sequence)) %>%
  #   #rename(Sequence = Phosphopeptide.sequence) %>% 
  #   rename(pool_id_theo_list=pool_id) 
  
  # ##### OLD AMPPING SYTLE #####
  
  # pep_list_w_theo_seq_comp <- pep_list_w_theo %>% 
  #   select(Sequence,pool_id) %>%
  #   #select(!pep_with_pos) %>%
  #   #rename(sequence = Sequence) %>% 
  #   mutate(pool_id_theo_list=pool_id) 
  #   
  # 
  # ##### OLD AMPPING SYTLE #####
  
  # pep_list_w_theo_unique <- pep_list_w_theo %>% select(Phosphopeptide.sequence, Pool) %>%
  #   distinct(Phosphopeptide.sequence,.keep_all = TRUE) %>%
  #   rename(Sequence = Phosphopeptide.sequence) %>%
  #   rename(Pool_for_seq_merge=Pool) %>%
  #   mutate(situation="Correct")
  
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
  
  # pep_list_w_theo_pep_map <- pep_list_w_theo %>%
  #   select(pep_with_pos,pool_id) %>%
  #   rename(pool_id_theo_list_pep=pool_id) 
  # 
  ## input sequence will be like this: (UniMod:1)AGGKPS(UniMod:21)QS(UniMod:21)PSQEAAGEAVLGAK
  
  df2 <- apply(quant_phospho[,"Modified.Sequence",drop=FALSE],1,getModificationPosition_)
  
  #apply(X = as.data.frame(quant_phospho[,"Modified.Sequence"]),1,function(x){getModificationPosition_general(mod_seq = x,software_name = )})
  
  test <- map_dfr(df2, enframe)
  
  ## OUTPUT FORMAT DOES NOT SUITABLE FOR DISTINGUSING BTW MODS 
  results1 <- map_dfr(df2, ~ enframe(.x)) %>%
    filter(grepl("modification",name) | grepl("pep_seq", name)) %>%
    mutate(value = map_chr(value, str_c, collapse="&")) %>%
    mutate(mods=case_when(grepl("(UniMod:21",fixed = T,name) ~ "phospho",
                          grepl("(UniMod:35",fixed = T,name) ~ "Oxidation",
                          grepl("(UniMod:1",fixed = T,name) ~ "N-term_Acetyl",
                          TRUE ~ ""))
  
  # map_dfr(df2, ~ enframe(.x) %>%
  #           mutate(value = map_chr(value, str_c, collapse="&")) %>% 
  #           filter(grepl("modification_\\(UniMod\\:21", name)) %>%
  #           add_column(modifications="Phospho") )
  
  
  ## RE-SHAPING DATA TO HAVE ONE SEQ with ALL MODS in ONE ROW
  
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
    full_join(filter(results_with_index, name == "pep_seq"), by = "id") %>% 
    ungroup() %>%
    ## SELECTING USEFUL COLUMNS
    select(c(name.x,value.x,mods.x,value.y))
  
  result_with_common_col <- reshaped_results %>%
    filter(grepl("modification_(UniMod:21)",fixed = T,name.x)) %>%
    mutate(pep_with_pos = paste(value.y,value.x, sep = "_")) %>%
    select(!value.y) %>%
    rename("mods_unimod_id" = "name.x",
           "phospho_positions" ="value.x",
           "all_mods_with_mod_type" = "mods.x")
  
  
  final_results_with_common_col <- mutate(result_with_common_col,quant_phospho)
  
  
  # Extracted values
  comb_result_pos <- final_results_with_common_col %>% 
    separate(pep_with_pos,into = c("seq","Positions"),sep = "_",remove = F) %>%
    select(Sequence, pep_with_pos,PTM.Site.Confidence,Positions,
           sample_name,raw_file,Q.Value, #### !!!  Q.Value is only necessary for 100% FDR RUNS !!! ####
           RT,Intensity, #Marked.as
           Protein.Names,pool_id_map_df,pool_id) %>%
    rename(Experiment=sample_name)
  #rename(ptm_score=extracted_values) %>%
  
  ##################################################################################################################
  ##################################################################################################################
  ### Mapping by precursor mass + phospho  and sequence ####
  
  amino_acid_table <- read.delim("/home/claude/post_processing_parquet/amino_acid_table.txt")
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
    select(!c(Positions,Sum)) %>%
    group_by(pep_with_pos,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    slice(which.max(Intensity)) %>%
    ungroup() 
  
  
  comb_result_seq_pool <- comb_result_mz_seq_map %>% 
    full_join(pep_list_w_theo_seq_map_pool,by=c("Sequence","pool_id")) %>% # #pep_with_pos
    mutate(map_seq_pool=ifelse(map_seq!='Subset or different mz' & pool_id_theo_list_seq==pool_id_map_df, "Correct Seq. & mz & pool",NA)) %>%
    mutate(map_seq_pool=ifelse(map_seq!='Subset or different mz' & is.na(pool_id_theo_list_seq), "Correct Seq. & mz but wrong pool",map_seq_pool)) %>%
    
    mutate(map_seq_pool=ifelse(map_seq== "Subset or different mz",map_seq,map_seq_pool)) %>%
    mutate(map_seq_pool = ifelse(is.na(Experiment), "missing", map_seq_pool))
  
  
  #
  ##################################################################################################################
  ##################################################################################################################
  
  
  #### TOTAL NUM. OF PHOSPHO-SEQUENCES ####
  ### Correct mapping was done using "map_df".
  # comb_result_seq <- comb_result_pos %>% 
  #   full_join(pep_list_w_theo_seq_comp,by=c("Sequence","pool_id")) %>%
  #   select(!c(pool_id,Experiment)) %>%
  #   ## RE-JOINING TO COUNT CORRECT, INCORRECT and MISSING SEQ. 
  #   full_join(mapping_file,by="raw_file") %>%
  #   rename(pool_id_map_df=pool_id) %>%
  #   mutate(Pool_for_seq_merge=ifelse(pool_id_theo_list==pool_id_map_df,"Correct Seq.","Wrong Seq.")) %>% #within theo list
  #   #mutate_at("Pool_for_seq_merge", ~replace_na(.,"Unexpected Seq.")) %>%
  #   mutate(Pool_for_seq_merge= ifelse(is.na(pool_id_map_df),"Missing",Pool_for_seq_merge)) %>%
  #   mutate(Pool_for_seq_merge=ifelse(is.na(pool_id_theo_list),"Wrong Seq.",Pool_for_seq_merge)) %>% # out of theo. list
  #   mutate(acq_type=acquisiton_type) %>%
  #   mutate(soft_name=software_name)
  # 
  ### Duplicate sequences were removed.
  ### Duplicate sequences were removed.
  comb_result_dist <- comb_result_seq_pool %>%
    # relocate(pool_id_map_df,pool_id_theo_list,Pool_for_seq_merge,.after = Sequence) %>% #
    group_by(raw_file,Sequence,Neutral_mass_res) %>%
    
    distinct(Sequence,.keep_all = T) %>%
    ungroup() 
  
  #### WRITE THE OBJECT AS TSV ####
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
  plot3 <- gg_barplt_id_pep_count_stack(data_set = comb_result_dist,
                                        x_df = comb_result_dist$Experiment,
                                        fill_df = comb_result_dist$map_seq_pool,
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
    filter(!grepl("Subset or different mz",map_seq) & !grepl("Missing",map_seq))
  
  #### VISUALIZATION OF TOTAL NUM. OF CORRECTLY IDENTIFIED PHOSPHO-SEQ ####
  
  plot2 <- gg_barplt_id_pep_count_stack(data_set =comb_result_seq_dist_cor,
                                        x_df =comb_result_seq_dist_cor$Experiment,
                                        fill_df = comb_result_seq_dist_cor$map_seq_pool,
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
  
  
  comb_result_pep_poolwise <- comb_result_seq_pool %>% 
    filter(!grepl("Subset or different mz",map_seq)) %>%
    filter(!grepl('Correct Seq. & mz but wrong pool',map_seq_pool)) %>%
    left_join(pep_list_w_theo_pep_map_pool ,by=c("pep_with_pos","pool_id")) %>%
    
    mutate(map_pep=ifelse(pool_id_map_df==pool_id_theo_list_pep,"Correct Seq. & Correct Loc.",NA)) %>%
    mutate(map_pep=ifelse(map_seq_pool=="Correct Seq. & mz & pool"& is.na(pool_id_theo_list_pep),"Correct Seq. & Wrong Loc.",map_pep)) %>%
    mutate(map_pep=ifelse(map_seq_pool=="missing","missing",map_pep)) %>%
    #mutate(map_pep=ifelse(is.na(pool_id_theo_list_pep),"Correct Seq. & Wrong Loc.",map_seq_pool)) %>%
    #select(!c(Positions,Sum)) %>%
    group_by(pep_with_pos,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    slice(which.max(Intensity)) %>%
    ungroup() %>%
    rename(ptm_score=PTM.Site.Confidence)
  
  #distinct(pep_with_pos,.keep_all = T) %>%
  #count(map_loc)
  # comb_result_pep <- comb_result_seq %>% 
  #   filter(!grepl("Wrong",Pool_for_seq_merge) & !grepl("Missing",Pool_for_seq_merge)) %>%
  #   
  #   select(!c(Pool_for_seq_merge,pool_id_theo_list,pool_id_map_df,sample_name)) %>%
  #   full_join(mapping_file,by="raw_file") %>%
  #   rename(pool_id_map_df=pool_id) %>%
  #   full_join(pep_list_w_theo_pep_comp,by="pep_with_pos") %>% ## IF FULL_JOIN IS USED,
  #   group_by(pep_with_pos,raw_file) %>%
  #   distinct(pep_with_pos,.keep_all = TRUE)%>%
  #   ungroup() %>%
  #   mutate(map_seq_pool=ifelse(pool_id_theo_list==pool_id_map_df,"Correct","Wrong Loc.within theo list")) %>%
  #   mutate(map_seq_pool= ifelse(is.na(pool_id_map_df),"Missing",map_seq_pool)) %>%
  #   mutate(map_seq_pool=ifelse(is.na(pool_id_theo_list),"Wrong Loc. out of theo. list",map_seq_pool)) %>%
  #   relocate(pool_id_map_df,pool_id_theo_list,map_seq_pool,.after = Sequence) %>% 
  #   mutate(acq_type=acquisiton_type) %>%
  #   mutate(soft_name=software_name) %>%
  #   rename(ptm_score=PTM.Site.Confidence)
  # 
  write.table(comb_result_pep_poolwise ,file=paste0(new_path,"/merge_theo_list_id_phospho_sites_only_unique_ones.tsv"), #file_path,curr_dir
              sep = "\t",col.names = T,row.names = F)
  
  #write.table(comb_result_pep ,file=paste0(new_path,"/merge_theo_list_id_STY_max_int.tsv"), #file_path,curr_dir
  #sep = "\t",col.names = T,row.names = F)
  
  
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
  
  comb_result_pep_filt <- comb_result_pep_poolwise  %>%
    filter(ptm_score >= 0.75)
  
  plot10 <- gg_barplt_id_pep_count_stack(data_set =comb_result_pep_filt,
                                         x_df =comb_result_pep_filt$Experiment,
                                         fill_df = comb_result_pep_filt$map_pep,
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
  comb_result_pep_cor <- comb_result_pep_poolwise %>%
    filter(grepl("Correct Seq. & mz & pool",map_seq_pool))
  
  #### VISUALIZATION OF TOTAL NUM. OF CORRECTLY IDENTIFIED & LOCALIZED PHOSPHO-PEP ####
  plot5 <- gg_barplt_id_pep_count_stack(data_set =comb_result_pep_poolwise,
                                        x_df =comb_result_pep_poolwise$Experiment,
                                        fill_df = comb_result_pep_poolwise$map_pep,
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
  
  
  comb_result_pep_filtered <- comb_result_pep_poolwise %>% 
    separate(Experiment,into = c("exp","samp","inj"),sep = "_",remove = F) %>%
    #filter(!grepl("Missing",map_seq_pool)) %>%
    mutate(samp_name_wo_inj=paste(exp,samp,sep = "-")) %>%
    group_by(samp_name_wo_inj,pep_with_pos) %>%
    slice(which.max(ptm_score)) %>%
    ungroup()
  
  
  
  write.table(comb_result_pep_filtered ,file=paste0(new_path,"/merge_theo_list_id_phospho_sites_max_ptm_score.tsv"), #file_path,curr_dir
              sep = "\t",col.names = T,row.names = F)
  
  #### LOCALIZATION ACCURACY ASSESSMENT (ROC LIKE PLOT GENERATION) #### 
  ### Column selection
  comb_result_pep_rmv_miss <- comb_result_pep_filtered %>% 
    select(pep_with_pos,map_pep,ptm_score,pool_id_map_df)# %>%
  #select(pep_with_pos,map_seq_pool,ptm_score,pool_id_map_df)# %>%
  #rename(ptm_score=extracted_values)
  
  ### Custom localization threshold determination and
  ### Counting total number of correct and wrong localization at a given threshold
  threshold <- seq(0,1,length=100)
  final_df <- NULL
  
  for (i in 1:length(threshold)){
    
    df <- comb_result_pep_rmv_miss %>% subset(ptm_score > threshold[i]) %>% 
      count(map_pep) %>% 
      mutate(threshold_val = threshold[i])
    
    final_df <-bind_rows(final_df,df)
    rm(df)
  }
  
  write.table(final_df, file=paste0(new_path,"/","Experiment",exp_id,software_name,"_num_sites_with_scores.tsv"),sep = "\t",col.names = T,row.names = F) #file_path,curr_dir
  ### Filtering only the correct ones
  correct_df <- final_df %>% filter(grepl("Correct Seq. & Correct Loc.",map_pep))#"Correct Seq. & mz & pool",map_pep))
  
  ### Visualization
  plot6 <- ggplot(correct_df,aes(y=n, x=threshold_val)) + geom_line(size=2) +
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
          axis.title=element_text(size=30)) #+ 
  #scale_x_continuous(limits = c(0, 100),breaks = seq(from = 0, to = 100, by = 10)) +  # Set the ticks for the x-axis
  #scale_y_continuous(limits = c(0, 420),breaks = seq(from = 0, to = 420, by = 50))
  
  
  barplt_df <- comb_result_pep_filtered  %>%
    select(Sequence,Experiment,Intensity,pep_with_pos,Protein.Names,ptm_score,map_seq_pool) %>%
    #filter(Global.Q.Value < 0.08) %>%
    separate(Experiment, into = c("Exp_id","Sample_id","Rep_id"),sep = "_",remove = F) %>%
    #mutate(sample_rep_id_seq = paste(pep_with_pos, Sample_id,Rep_id, sep = "_")) %>%
    group_by(pep_with_pos,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    slice(which.max(Intensity)) %>%
    ungroup()
  
  ####### ADDITIONAL PLOT TO DISPLAY MISSING and UNEXPECTED PEPTIDES ########
  
  plot12 <- gg_barplt_id_pep_count(data_set = barplt_df,
                                   x_df = barplt_df$map_seq_pool,
                                   fill_df = barplt_df$map_seq_pool,
                                   ymax = 20000,
                                   size_num = 10,
                                   header = "Total number of phospho-site across each sample before filtering",
                                   caption_lab = "NA values are removed.",
                                   x_lab = "Sample id",
                                   fill_lab =  "Sample id",
                                   y_lab = "Number of identified peptides",
                                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
  barplt_df_filt <-barplt_df  %>% filter(ptm_score >= 0.75)
  plot13 <- gg_barplt_id_pep_count(data_set = barplt_df_filt,
                                   x_df = barplt_df_filt$map_seq_pool,
                                   fill_df = barplt_df_filt$map_seq_pool,
                                   ymax = 20000,
                                   size_num = 10,
                                   header = "Total number of phospho-site across each sample after loc. filtering",
                                   caption_lab = "NA values are removed.",
                                   x_lab = "Sample id",
                                   fill_lab =  "Sample id",
                                   y_lab = "Number of identified peptides",
                                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
  
  
  barplt_df_ecoli <- quant_peptides_with_cond %>% 
    filter(grepl(background_species, Protein.Names)) %>%
    select(Sequence,sample_name,Intensity,Protein.Names) %>%
    separate(sample_name, into = c("Exp_id","Sample_id", "Rep_id"), sep = "_",remove = F) %>%
    #mutate(sample_rep_id_seq = paste(Sequence,Experiment,sep = "_")) %>%
    group_by(Sequence,sample_name) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    slice(which.max(Intensity)) %>%
    ungroup()
  
  plot1 <- gg_barplt_id_pep_count(data_set = barplt_df_ecoli,
                                  x_df = barplt_df_ecoli$Sample_id,
                                  fill_df = barplt_df_ecoli$Rep_id,
                                  ymax = 20000,
                                  size_num = 10,
                                  header = "Total number of quantified phospho-site across each sample",
                                  caption_lab = "NA values are removed.",
                                  x_lab = "Sample id",
                                  fill_lab =  "Sample id",
                                  y_lab = "Number of identified peptides",
                                  subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
  
  ### This regex was designed to check every S,T and Y in each sequence
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
    bind_rows(#comb_result_pep_adj_5 disabled for withMBR
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
           map_seq_pool,
           map_pep,
           ptm_score,
           pep_with_pos,
           pool_id_map_df,
           Experiment,
           raw_file,is_adj) 
  
  write.table(comb_result_pep_adj_all,file = paste0(new_path,"/all_adj&nonadj_corr_wrong_loc_nolocthreshold.txt"),sep = "\t",col.names = T,row.names = F)
  
  ### TOTAL COUNT OF ALL SEQUENCES WITHOUT CONSIDERING ADJ_COUNT
  plot7 <- comb_result_pep_adj_all %>% tibble() %>%
    filter(grepl("Correct Seq. & Correct Loc",map_pep)) %>%
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
  
  ### TOTAL COUNT OF ALL SEQUENCES WITH CONSIDERING ADJ_COUNT AND ACCURACY
  plot8 <- comb_result_pep_adj_all %>%
    mutate(pep_class=map_pep) %>%
    #mutate(pep_class=ifelse(map_seq_pool!= "Correct Seq. & mz & pool","Wrong",pep_class)) %>%
    #filter(grepl("Correct",map_seq_pool)) %>%
    rowwise() %>%
    #distinct(Sequence,.keep_all = TRUE) %>%
    #group_by(pep_with_pos,new_col) %>%
    #count(pep_with_pos)
    #mutate(ptm_score = ifelse(map_seq_pool == "Wrong Localization",(ptm_score*-1),ptm_score)) %>%
    group_by(is_adj,pep_class) %>%
    count(STY_len) %>% #map_seq_pool
    ungroup() %>%
    mutate(n_label=n) %>%
    #mutate(n=ifelse(map_seq_pool =="Wrong Localization", (n* (-1)),n)) %>%
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
    #mutate(map_seq_pool=ifelse(map_seq_pool!="Correct Seq. & mz & pool","Wrong",map_seq_pool)) %>%
    #filter(!grepl("Wrong Loc. out of theo. list",map_seq_pool)) %>%
    #filter(ptm_score > 0.75) %>%
    count(row_sum,is_adj,map_seq_pool) %>% mutate(n_label=n) %>%
    mutate(n=ifelse(is_adj =="non_adjacent", (n* (-1)),n)) %>%
    mutate(STY_len= ifelse(row_sum < 0, (-1*row_sum),row_sum)) %>%
    ggplot(aes(x=STY_len,y=n,fill=map_seq_pool,pattern=is_adj)) +
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
  
  plot14 <- comb_result_pep_adj_all %>% filter(!grepl("non_adjacent",is_adj)) %>%
    mutate(pep_class=map_pep) %>%
    #mutate(map_seq_pool=ifelse(map_seq_pool!="Correct Seq. & mz & pool","Wrong",map_seq_pool)) %>%
    #filter(!grepl("Wrong Loc. out of theo. list",map_seq_pool)) %>%
    #filter(ptm_score > 0.75) %>%
    count(STY_len,is_adj,map_seq_pool) %>% mutate(n_label=n) %>%
    mutate(n=ifelse(map_seq_pool =="Wrong", (n* (-1)),n)) %>%
    mutate(STY_len= ifelse(STY_len < 0, (-1*STY_len),STY_len)) %>%
    ggplot(aes(x=STY_len,y=n,fill=map_seq_pool)) +
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
    #mutate(pep_class=ifelse(map_seq_pool!= "Correct Seq. & mz & pool","Wrong",pep_class)) %>%
    subset(!(STY_len == 1 & is_adj == "non_adjacent")) %>%
    
    ggplot( aes(x= ptm_score,fill=interaction(is_adj))) +
    geom_bar(stat="count") + labs(x = "Count of STY amino acids", y = "Total count",
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
  sapply(1:length(plt_obj),function(x) tryCatch(
    ggsave(filename = paste0("p",x,".png"),
           width = 60, height = 45, 
           path = paste0(new_path),
           units = "cm",
           get(plt_obj[x]),
           device = "png", dpi = 150,
    ),
    error = function(e) message("Warning: plot ", plt_obj[x], " failed to render: ", e$message)
  ))
  
}
