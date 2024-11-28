#library(PhosR)
library(stringr)
library(dplyr)
library(data.table)
library(openxlsx)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(purrr)
library(gtools)
library(ggpattern)
###############################################
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/ggplot/ggplot_functions.R")
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/get_modification_func/getModificationPosition_func_for_all_mods.R")

file_path <- c("D:/dev/Pinar/PHD/wet_lab_experiments/DIA_analysis/Exp_1/DIANN/Exp1_with_Ecoli/DIANN_v1.9/noMBR/misclv2_MaxVarMod5_charge14_FDR100/")
#"D:/dev/Pinar/PHD/wet_lab_experiments/DIA_analysis/Exp_1/DIANN/Exp1_with_Ecoli/DIANN_v1.9/noMBR/misclv2_MaxVarMod5_charge14/")
#"D:/dev/Pinar/PHD/wet_lab_experiments/DIA_analysis/Exp_1/DIANN/Exp1_with_Ecoli/DIANN_v1.9/withMBR/") #noMBR
file_name <- c("Exp1_with_ecoli_misclv2_MaxVarMod5_charge14_FDR100.tsv")#"misclv2_MaxVarMod5_charge14_FDR100report.tsv")
#"misclv2_MaxVarMod5_charge14_Exploris_report.tsv")#"same_params_and_lib_exp2_redesigned_MaxVarMods4_FDR1_Misclgvg2_MBR.tsv")
#"same_params_and_lib_exp2_redesigned_MaxVarMods4_FDR1_Misclgvg2_noMBR.tsv")
mapping <- "D:/dev/Pinar/PHD/data_analysis/DIA_mapping/exp1_batch2_mapping_btw_rawfile_pool_id.txt"

#"D:\dev\Pinar\PHD\wet_lab_experiments\DIA_data_analysis\experiment_1\DIANN\Exp1_with_Ecoli"

selected_spcies="HUMAN"
background_species= "ECOLI"
theo_file_path="D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/"
theo_file_name="Synthetic peptides list_theo_conc_corrected_isomericity_new_with_plates.xlsx"
sheet_theo_name ="ISOREF_REF2_Others" #"ISO-refOTHER with FC_correct"
acquisiton_type="DIA no FAIMS Exploris"
subtitle = #"version 1.9 & no MBR charge 1 & 4 MaxVarMod 5 "
#fdr_threshold = 0.05

#exp_design=experiment_name
exp_id=1
software_name="DIA-NN"
num_reps=3
#test_type="limma"

final_diann_pep_quant_analysis_syn_exp1(file_path = file_path,
                                        file_name = file_name,
                                        sheet_name="",
                                        theo_file_path=theo_file_path,
                                        theo_file_name=theo_file_name,
                                        sheet_theo_name=sheet_theo_name,
                                        background_species="ECOLI",
                                        selected_spcies="HUMAN",
                                        exp_id=1,
                                        mapping=mapping,
                                        #exp_design,
                                        #fdr_threshold,
                                        acquisiton_type=acquisiton_type,
                                        software_name="DIA-NN",
                                        #test_type,
                                        num_reps=3,
                                        #actual_ratio,
                                        subtitle="version 1.9 & no MBR charge 1 & 4 MaxVarMod 4 FDR100%")


final_diann_pep_quant_analysis_syn_exp1 <- function(file_path,
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
  
  intended_dir <-paste0(file_path,"output_final/")
  
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
  quant_peptides <- read_tsv(paste0(file_path,file_name),col_names = T)
  
  quant_peptides_with_cond <- quant_peptides %>% 
    rename(raw_file=Run) %>%
    left_join(mapping_file,by="raw_file") %>%
    rename("Sequence"="Stripped.Sequence") %>%
    rename("Intensity"= "PG.Normalised")# %>%
  #filter(Global.Q.Value < 0.05)
  
  ## PHOSPHO-FILTERING
  quant_phospho <- quant_peptides_with_cond %>% 
    filter(grepl("HUMAN", Protein.Names) & !grepl("CON__",Protein.Names)) %>%
    filter(grepl("UniMod:21",Modified.Sequence)) #%>%
  #filter(PTM.Q.Value < 0.05 )
  
  ## THEORETICAL PEPTIDE LIST
  pep_list_w_theo <- read.xlsx(paste0(theo_file_path, theo_file_name), sheet = sheet_theo_name)
  
  pep_list_w_theo_pep_comp <- pep_list_w_theo %>% 
    select(!c(Sequence,Phosphopeptide.sequence)) %>%
    #rename(Sequence = Phosphopeptide.sequence) %>% 
    rename(pool_id_theo_list=pool_id) 
  
  pep_list_w_theo_seq_comp <- pep_list_w_theo %>% 
    select(Sequence,pool_id) %>%
    #select(!pep_with_pos) %>%
    #rename(sequence = Sequence) %>% 
    mutate(pool_id_theo_list=pool_id) 
  
  
  ## input sequence will be like this: (UniMod:1)AGGKPS(UniMod:21)QS(UniMod:21)PSQEAAGEAVLGAK
  
  df2 <- apply(quant_phospho[,"Modified.Sequence"],1,getModificationPosition_)
  
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
    select(Sequence, pep_with_pos,PTM.Site.Confidence,
           sample_name,raw_file,Q.Value, #### !!!  Q.Value is only necessary for 100% FDR RUNS !!! ####
           RT,Intensity, #Marked.as
           Protein.Names,pool_id) %>%
    rename(Experiment=sample_name)
    #rename(ptm_score=extracted_values) %>%

  #### TOTAL NUM. OF PHOSPHO-SEQUENCES ####
  ### Correct mapping was done using "map_df".
  comb_result_seq <- comb_result_pos %>% 
    full_join(pep_list_w_theo_seq_comp,by=c("Sequence","pool_id")) %>%
    select(!c(pool_id,Experiment)) %>%
    ## RE-JOINING TO COUNT CORRECT, INCORRECT and MISSING SEQ. 
    full_join(mapping_file,by="raw_file") %>%
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
  write.table(comb_result_seq,file=paste0(new_path,"/", #file_path,curr_dir
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
  
  plot2 <- gg_barplt_id_pep_count_stack(data_set =comb_result_seq_dist_cor,
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
    full_join(mapping_file,by="raw_file") %>%
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
    mutate(soft_name=software_name) %>%
    rename(ptm_score=PTM.Site.Confidence)
  
  write.table(comb_result_pep ,file=paste0(new_path,"/merge_theo_list_id_phospho_sites_only_unique_ones.tsv"), #file_path,curr_dir
              sep = "\t",col.names = T,row.names = F)
  
  #write.table(comb_result_pep ,file=paste0(new_path,"/merge_theo_list_id_STY_max_int.tsv"), #file_path,curr_dir
  #sep = "\t",col.names = T,row.names = F)
  
  
  #### VISUALIZATION OF TOTAL NUM. OF IDENTIFIED & LOCALIZED PHOSPHO-PEP ####
  plot4 <- gg_barplt_id_pep_count_stack(data_set =comb_result_pep,
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

   comb_result_pep_filt <- comb_result_pep %>% filter(ptm_score >= 0.75)
  
  plot10 <- gg_barplt_id_pep_count_stack(data_set =comb_result_pep_filt,
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
  plot5 <- gg_barplt_id_pep_count_stack(data_set =comb_result_pep_cor,
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
  

    comb_result_pep_filtered <- comb_result_pep %>% 
      separate(sample_name,into = c("exp","samp","inj"),sep = "_",remove = F) %>%
      filter(!grepl("Missing",Pool_for_pep_merge)) %>%
      mutate(samp_name_wo_inj=paste(exp,samp,sep = "-")) %>%
      group_by(samp_name_wo_inj,pep_with_pos) %>%
      slice(which.max(ptm_score)) %>%
      ungroup()
    
 
  
  write.table(comb_result_pep_filtered ,file=paste0(new_path,"/merge_theo_list_id_phospho_sites_max_ptm_score.tsv"), #file_path,curr_dir
              sep = "\t",col.names = T,row.names = F)
  
  #### LOCALIZATION ACCURACY ASSESSMENT (ROC LIKE PLOT GENERATION) #### 
  ### Column selection
  comb_result_pep_rmv_miss <- comb_result_pep_filtered %>% 
    select(pep_with_pos,Pool_for_pep_merge,ptm_score,pool_id_map_df)# %>%
  #rename(ptm_score=extracted_values)
  
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
  
  write.table(final_df, file=paste0(new_path,"/","Experiment",exp_id,software_name,"_num_sites_with_scores.tsv"),sep = "\t",col.names = T,row.names = F) #file_path,curr_dir
  ### Filtering only the correct ones
  correct_df <- final_df %>% filter(grepl("Correct",Pool_for_pep_merge))
  
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
          axis.title=element_text(size=30)) + 
    #scale_x_continuous(limits = c(0, 100),breaks = seq(from = 0, to = 100, by = 10)) +  # Set the ticks for the x-axis
    scale_y_continuous(limits = c(0, 180),breaks = seq(from = 0, to = 180, by = 10))
  
  
  barplt_df <- comb_result_pep_filtered  %>%
    select(Sequence,sample_name,Intensity,pep_with_pos,Protein.Names,ptm_score,Pool_for_pep_merge) %>%
    rename(Experiment=sample_name) %>%
    #filter(Global.Q.Value < 0.08) %>%
    separate(Experiment, into = c("Exp_id","Sample_id","Rep_id"),sep = "_",remove = F) %>%
    #mutate(sample_rep_id_seq = paste(pep_with_pos, Sample_id,Rep_id, sep = "_")) %>%
    group_by(pep_with_pos,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    slice(which.max(Intensity)) %>%
    ungroup()
  
  ####### ADDITIONAL PLOT TO DISPLAY MISSING and UNEXPECTED PEPTIDES ########
  
  plot12 <- gg_barplt_id_pep_count(data_set = barplt_df,
                                x_df = barplt_df$Pool_for_pep_merge,
                                fill_df = barplt_df$Pool_for_pep_merge,
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
                                x_df = barplt_df_filt$Pool_for_pep_merge,
                                fill_df = barplt_df_filt$Pool_for_pep_merge,
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
  
  ### COMBINATION OF ALL ADJACENT ONES
  comb_result_pep_adj_all <- comb_result_pep_non_adj %>%
    select(!c(STY_adj)) %>%
    rename(STY_adj = STY_adj_new) %>%
    bind_rows(#comb_result_pep_adj_5 disabled for withMBR
              comb_result_pep_adj_4,
              comb_result_pep_adj_3,
              comb_result_pep_adj_2) %>%
    select(Sequence,
           S_count,
           T_count,
           Y_count,
           row_sum,
           STY_adj,
           STY_len,
           Pool_for_pep_merge,
           ptm_score,
           pep_with_pos,
           pool_id_map_df,
           sample_name,
           raw_file,is_adj) 
  
  write.table(comb_result_pep_adj_all,file = paste0(new_path,"/all_adj&nonadj_corr_wrong_loc_nolocthreshold.txt"),sep = "\t",col.names = T,row.names = F)
  
  ### TOTAL COUNT OF ALL SEQUENCES WITHOUT CONSIDERING ADJ_COUNT
  plot7 <- comb_result_pep_adj_all %>% tibble() %>%
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
  
  plot11 <- comb_result_pep_adj_all %>% mutate(Pool_for_pep_merge=ifelse(Pool_for_pep_merge!="Correct","Wrong",Pool_for_pep_merge)) %>%
    #filter(!grepl("Wrong Loc. out of theo. list",Pool_for_pep_merge)) %>%
    #filter(ptm_score > 0.75) %>%
    count(row_sum,is_adj,Pool_for_pep_merge) %>% mutate(n_label=n) %>%
    mutate(n=ifelse(is_adj =="non_adjacent", (n* (-1)),n)) %>%
    mutate(STY_len= ifelse(row_sum < 0, (-1*row_sum),row_sum)) %>%
    ggplot(aes(x=STY_len,y=n,fill=Pool_for_pep_merge,pattern=is_adj)) +
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
    mutate(Pool_for_pep_merge=ifelse(Pool_for_pep_merge!="Correct","Wrong",Pool_for_pep_merge)) %>%
    #filter(!grepl("Wrong Loc. out of theo. list",Pool_for_pep_merge)) %>%
    #filter(ptm_score > 0.75) %>%
    count(STY_len,is_adj,Pool_for_pep_merge) %>% mutate(n_label=n) %>%
    mutate(n=ifelse(Pool_for_pep_merge =="Wrong", (n* (-1)),n)) %>%
    mutate(STY_len= ifelse(STY_len < 0, (-1*STY_len),STY_len)) %>%
    ggplot(aes(x=STY_len,y=n,fill=Pool_for_pep_merge)) +
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
    mutate(pep_class=Pool_for_pep_merge) %>%
    mutate(pep_class=ifelse(Pool_for_pep_merge!= "Correct","Wrong",pep_class)) %>%
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
  sapply(1:length(plt_obj),function(x) ggsave(filename = paste0("p",x,".png"),
                                              width = 60, height = 45, 
                                              path = paste0(new_path),
                                              units = "cm",
                                              get(plt_obj[x]),
                                              device = "png", #".svg"
  ))
  
}
  