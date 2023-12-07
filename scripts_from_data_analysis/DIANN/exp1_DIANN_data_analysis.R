#library(PhosR)
library(stringr)
library(dplyr)
library(data.table)
library(openxlsx)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(purrr)
###############################################
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/ggplot/ggplot_functions.R")
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/get_modification_func/getModificationPosition_func_for_all_mods.R")

file_path <- "D:/dev/Pinar/PHD/wet_lab_experiments/DIA_data_analysis/experiment_1/DIANN/Exp1_with_Ecoli/"
file_name <- "report.tsv" #"exp2_unimod_report.tsv"
mapping <- "D:/dev/Pinar/PHD/data_analysis/DIA_data_processing/mapping_files/exp1_batch2_mapping_btw_rawfile_pool_id.txt"

"D:\dev\Pinar\PHD\wet_lab_experiments\DIA_data_analysis\experiment_1\DIANN\Exp1_with_Ecoli"

selected_spcies="HUMAN"
background_species= "ECOLI"
theo_file_path= "D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/"
theo_file_name="Synthetic peptides list_theo_conc_corrected_isomericity_new.xlsx" #"Synthetic peptides list_theo_conc_corrected_pool_id_iso_count_final.xlsx"
sheet_theo_name ="ISO-refOTHER with FC_correct" #"ISO-ref and OTHER with FC"
acquisiton_type="DIA no FAIMS Exploris"
subtitle = ""
#fdr_threshold = 0.05
actual_ratio = c(2,10,20,100)
#exp_design=experiment_name
exp_id=1
software_name="DIANN"
num_reps=3
#test_type="limma"


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
  
  
  mapping_file <- read.table(mapping,sep = "\t",header = T)
  exp_design <- mapping_file$Experiment
  
  sample_size <- length(exp_design) / num_reps
  sample_names <- paste0("A",1:sample_size)
  comparisons <- NULL
  for (i in 1:sample_size){
    tmp <- paste0(sample_names[1], "/",sample_names[i])
    comparisons[i] <- tmp
    rm(tmp)
  }
  comparisons <- comparisons[-1]
  
  quant_peptides <- read_tsv(paste0(file_path,file_name),col_names = T)
  
  quant_peptides_with_cond <- quant_peptides %>% 
    left_join(mapping_file,by="Run") %>%
    rename("Sequence"="Stripped.Sequence") %>%
    rename("Intensity"= "PG.Normalised")# %>%
  #filter(Global.Q.Value < 0.05)
  
  ## PHOSPHO-FILTERING
  quant_phospho <- quant_peptides_with_cond %>% 
    filter(grepl("HUMAN", Protein.Names) & !grepl("CON__",Protein.Names)) %>%
    filter(grepl("UniMod:21",Modified.Sequence)) #%>%
  #filter(PTM.Q.Value < 0.05 )
  
  ## THEORETICAL PEPTIDE LIST
  pep_list_w_theo_quant <- read.xlsx(paste0(theo_file_path, theo_file_name), sheet = sheet_theo_name)
  pep_list_w_theo_quant <- pep_list_w_theo_quant[,-1]
  
  common_col_theo_quant <- as.data.frame(paste(pep_list_w_theo_quant$Phosphopeptide.sequence,
                                               pep_list_w_theo_quant$modified.position.in.peptide, sep = "_"))
  colnames(common_col_theo_quant) <- "pep_with_pos"
  pep_list_w_theo_quant_new <- cbind(common_col_theo_quant,pep_list_w_theo_quant)
  
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
  
  barplt_df <- final_results_with_common_col %>%
    select(Sequence,Experiment,Intensity,pep_with_pos,Protein.Names,Global.Q.Value) %>%
    #filter(Global.Q.Value < 0.08) %>%
    separate(Experiment, into = c("Exp_id","Sample_id","Rep_id"),sep = "_",remove = F) %>%
    #mutate(sample_rep_id_seq = paste(pep_with_pos, Sample_id,Rep_id, sep = "_")) %>%
    group_by(pep_with_pos,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    slice(which.max(Intensity)) %>%
    ungroup()
  
  ####### ADDITIONAL PLOT TO DISPLAY MISSING and UNEXPECTED PEPTIDES ########
  df_merge_syn <- barplt_df %>%
    select(pep_with_pos,Experiment,Intensity, Protein.Names) %>% 
    pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
    full_join(pep_list_w_theo_quant_new,by="pep_with_pos") %>% 
    mutate_at("Pool", ~replace_na(.,"Unexpected")) %>%
    mutate(Pool= ifelse(is.na(Protein.Names),"missing",Pool)) %>%
    select(pep_with_pos,starts_with(exp_design),Pool) %>%
    mutate(soft_name=software_name, ion_mobility=acquisiton_type)
    
    
  
  p11 <- gg_barplt_id_pep_count(data_set = df_merge_syn,
                                x_df = df_merge_syn$Pool,
                                fill_df = df_merge_syn$Pool,
                                ymax = 20000,
                                header = "Total number of quantified phospho-site across each sample",
                                caption_lab = "NA values are removed.",
                                x_lab = "Sample id",
                                fill_lab =  "Sample id",
                                y_lab = "Number of identified peptides",
                                subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
  
  
  barplt_df_ecoli <- quant_peptides_with_cond %>% 
    filter(grepl(background_species, Protein.Names)) %>%
    select(Sequence,Experiment,Intensity,Protein.Names) %>%
    separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "_",remove = F) %>%
    #mutate(sample_rep_id_seq = paste(Sequence,Experiment,sep = "_")) %>%
    group_by(Sequence,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    slice(which.max(Intensity)) %>%
    ungroup()
  
  p2 <- gg_barplt_id_pep_count(data_set = barplt_df_ecoli,
                               x_df = barplt_df_ecoli$Sample_id,
                               fill_df = barplt_df_ecoli$Rep_id,
                               ymax = 20000,
                               header = "Total number of quantified phospho-site across each sample",
                               caption_lab = "NA values are removed.",
                               x_lab = "Sample id",
                               fill_lab =  "Sample id",
                               y_lab = "Number of identified peptides",
                               subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
  
  p1 <- gg_barplt_id_pep_count(data_set = barplt_df,
                               x_df = barplt_df$Sample_id,
                               fill_df = barplt_df$Rep_id,
                               ymax = 200,
                               header = "Total number of quantified phospho-site across each sample",
                               caption_lab = "NA values are removed.",
                               x_lab = "Sample id",
                               fill_lab =  "Sample id",
                               y_lab = "Number of identified peptides",
                               subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
  
}
  