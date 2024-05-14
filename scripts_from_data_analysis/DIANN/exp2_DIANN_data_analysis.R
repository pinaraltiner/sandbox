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
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/roc_curve/new_roc_curve_generation_with_custom_threshold.R")

final_diann_pep_quant_analysis_syn <- function(file_path,
                                                     file_name,
                                                     sheet_name,
                                                     theo_file_path,
                                                     theo_file_name,
                                                     sheet_theo_name,
                                                     background_species,
                                                     selected_spcies,
                                                     exp_id,
                                                     loc_filter,
                                                     loc_filter_opt,
                                                     numerator, ## Which sample id was selected as numerator for quant process
                                                     mapping,
                                                     #exp_design,
                                                     fdr_threshold,
                                                     acquisiton_type,
                                                     software_name,
                                                     test_type,
                                                     num_reps,
                                                     actual_ratio,
                                                     subtitle,
                                               size_variying_pep){
  
  
  mapping_file <- read.table(mapping,sep = "\t",header = T)
  exp_design <- mapping_file$Experiment
  
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

  exp_design <- str_sort(exp_design)
  
  quant_peptides <- read_tsv(paste0(file_path,file_name),col_names = T)
  
  quant_peptides_with_cond <- quant_peptides %>% 
    filter(!grepl("CON__",Protein.Names)) %>%
    left_join(mapping_file,by="Run") %>%
    rename("Sequence"="Stripped.Sequence") %>%
    rename("Intensity"= "PG.Normalised")
  
  imputed_values <- quant_peptides_with_cond  %>%
    group_by(Experiment) %>% 
    summarise(first_quantile=quantile(Intensity,probs=0.01,na.rm=TRUE))
  
  imputed_values_vec <- as.vector(imputed_values$first_quantile)
  
  ## PHOSPHO-FILTERING
  quant_phospho <- quant_peptides_with_cond %>% 
    filter(grepl("HUMAN", Protein.Names) & !grepl("CON__",Protein.Names)) %>%
    filter(grepl("UniMod:21",Modified.Sequence))
  
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
  
  pep_list_w_theo_unique <- pep_list_w_theo %>% select(Phosphopeptide.sequence, Pool) %>%
    distinct(Phosphopeptide.sequence,.keep_all = TRUE) %>%
    rename(Sequence = Phosphopeptide.sequence) %>%
    rename(Pool_for_seq_merge=Pool) %>%
    mutate(situation="Correct")

  ############################
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
    select(Sequence,Experiment,Intensity,pep_with_pos,Protein.Names,Run,PTM.Site.Confidence) %>%
    separate(Experiment, into = c("Exp_id","Sample_id","Rep_id"),sep = "-",remove = F) %>%
    #mutate(sample_rep_id_seq = paste(pep_with_pos, Sample_id,Rep_id, sep = "_")) %>%
    group_by(pep_with_pos,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    slice(which.max(Intensity)) %>%
    ungroup()
  
  ######################################################################
  
  site_prob <- barplt_df %>% 
    select(PTM.Site.Confidence, pep_with_pos,Experiment) %>% 
    pivot_wider(names_from = "Experiment",
                values_from = "PTM.Site.Confidence") %>% 
    relocate(exp_design) %>%
    rowwise() %>%
    mutate(max_value = max(c_across(all_of(exp_design)), na.rm = TRUE)) %>%
    select(!exp_design)
  #mutate(max_value=ifelse(!is.numeric(max_value),NA,max_value))
  #filter(!grepl(-Inf,max_value))
  
  ecoli_seq_dist <- quant_peptides_with_cond %>% 
    select(Sequence,Modified.Sequence, Protein.Names) %>%
    filter(grepl(background_species,Protein.Names)) %>%
    mutate(species=background_species) %>%
    distinct(Sequence,.keep_all = T)
  
  
  ecoli_seq <- quant_peptides_with_cond %>% 
    select(Sequence,Modified.Sequence, Protein.Names) %>%
    filter(grepl(background_species,Protein.Names)) %>%
    mutate(species=background_species)
  

  all_seq <- quant_peptides_with_cond %>%
    filter(grepl("UniMod:21",Modified.Sequence)) %>%
    filter(grepl(selected_spcies,Protein.Names) & !grepl("CON__",Protein.Names)) %>%
    separate(Protein.Names, into = c("protein","species"),remove = F,sep="_") %>%
    full_join(pep_list_w_theo_unique,by="Sequence") %>%
    mutate_at("Pool_for_seq_merge", ~replace_na(.,"Unexpected")) %>%
    mutate(Pool_for_seq_merge= ifelse(is.na(species),"missing",Pool_for_seq_merge)) %>%
    filter(!grepl("Unexpected",Pool_for_seq_merge)) %>%
    bind_rows(ecoli_seq) %>% 
    mutate(Pool_for_seq_merge= ifelse(is.na(Pool_for_seq_merge),background_species,Pool_for_seq_merge))
    
  all_seq_syn <- all_seq %>%
    select(Sequence, Modified.Sequence, Pool_for_seq_merge,PTM.Site.Confidence) %>%
    distinct(Sequence, .keep_all = T) %>%
    bind_rows(ecoli_seq_dist) %>%
    mutate(Pool_for_seq_merge= ifelse(is.na(Pool_for_seq_merge),background_species,Pool_for_seq_merge)) %>%
    mutate(acq_type=acquisiton_type) %>%
    mutate(soft_name=software_name)
  
  
  ########## ########## ########## ########## ########## ########## ########## ##########
  
  all_corr_seq <- quant_peptides_with_cond %>% 
    filter(grepl("21",Modified.Sequence)) %>%
    filter(grepl(selected_spcies,Protein.Names) & !grepl("CON__",Protein.Names)) %>%
    mutate(species=selected_spcies) %>% 
    group_by(Experiment) %>%
    distinct(Sequence,.keep_all = T) %>%
    full_join(pep_list_w_theo_unique,by="Sequence") %>%
    mutate_at("situation", ~replace_na(.,"Unexpected")) %>%
    mutate(situation= ifelse(is.na(species),"missing",situation)) %>%
    filter(!grepl("Unexpected",situation) & !grepl("missing",situation)) %>%
    select(Sequence,Modified.Sequence,PTM.Site.Confidence,Experiment,situation,Pool_for_seq_merge) %>%
    separate(Experiment, into = c("exp_id","samp_id","rep_id"),sep = "-")
  
  plot15 <- gg_barplt_id_pep_count(data_set = all_corr_seq,
                                   x_df = all_corr_seq$samp_id,
                                   fill_df = all_corr_seq$rep_id,
                                   ymax = nrow(all_corr_seq),
                                   size_num=10,
                                   header = paste("Total number of correctly identified ",selected_spcies,"phospho-sequence","across each sample",sep=" "),
                                   caption_lab = "Mapping was done without considering phospho-positions.",
                                   x_lab = "Sample id",
                                   fill_lab =  "Sample id",
                                   y_lab = "Number of identified peptides",
                                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) +
    theme(axis.text.x = element_text(angle = 90))
  
  ########## ########## ########## ########## ########## ########## ########## ##########
  
  
   ## This is the same as p13 in DDA data analysis
  plot12 <- gg_barplt_id_pep_count(data_set = all_seq_syn,
                                x_df = all_seq_syn$Pool_for_seq_merge,
                                fill_df = all_seq_syn$Pool_for_seq_merge,
                                ymax = 20000,
                                size_num = 10,
                                header = paste("Total number of identified phosphorylated", selected_spcies,"and", background_species,"sequence across each sample",sep=" "),
                                caption_lab = "NA values are removed. (p13)",
                                x_lab = "Sample id",
                                fill_lab =  "Sample id",
                                y_lab = "Number of identified sequence",
                                subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
  write.table(all_seq_syn, file=paste0(file_path,"Experiment2",software_name,"_number_of_unique_sequence_for_each_species.txt"),sep = "\t",col.names = T,row.names = F)

  
  barplt_df_ecoli <- quant_peptides_with_cond %>% 
    filter(grepl(background_species, Protein.Names)) %>%
    select(Sequence,Experiment,Intensity,Protein.Names) %>%
    separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F) %>%
    mutate(sample_rep_id_seq = paste(Sequence, Sample_id,Rep_id, sep = "_")) %>%
    group_by(sample_rep_id_seq,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    slice(which.max(Intensity)) %>%
    ungroup()
  
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
                               subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
  # barplt_prot <- quant_phospho %>%
  #   select(Sequence,Experiment,Intensity,Protein.Names,Run) %>%
  #   separate(Experiment, into = c("Exp_id","Sample_id","Rep_id"),sep = "-",remove = F) %>%
  #   #mutate(sample_rep_id_seq = paste(pep_with_pos, Sample_id,Rep_id, sep = "_")) %>%
  #   group_by(Protein.Names,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
  #   slice(which.max(Intensity)) %>%
  #   ungroup()
  # 
  # plotx <- gg_barplt_id_pep_count(data_set = barplt_prot,
  #                                 x_df = barplt_prot$Sample_id,
  #                                 fill_df = barplt_prot$Rep_id,
  #                                 ymax = 20000,
  #                                 size_num = 10,
  #                                 header = "Total number of quantified phospho-site across each sample",
  #                                 caption_lab = "NA values are removed.",
  #                                 x_lab = "Sample id",
  #                                 fill_lab =  "Sample id",
  #                                 y_lab = "Number of identified peptides",
  #                                 subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
  barplt_prot_ecoli <- quant_peptides_with_cond %>% 
    filter(grepl(background_species, Protein.Names)) %>%
    select(Sequence,Experiment,Intensity,Protein.Names) %>%
    separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F) %>%
    #mutate(sample_rep_id_seq = paste(Sequence, Sample_id,Rep_id, sep = "_")) %>%
    group_by(Protein.Names,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    slice(which.max(Intensity)) %>%
    ungroup() %>%
    mutate(type=subtitle)
  
  plot2 <- gg_barplt_id_pep_count(data_set = barplt_prot_ecoli,
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
  
  
  write.table(barplt_prot_ecoli, file=paste0(file_path,"Exp2",software_name,"_number_of_unique_",background_species,"proteins_",".txt"),sep = "\t",col.names = T,row.names = F)
  
  
  
  ## THIS RESHAPING IS ONLY FOR ELIMINATION OF MULTIPLE PHOSPHO-SITES and ECOLI PEPTIDES
  ## ELIMINATION STEP IS NOT NECESSARY FOR ECOLI, 1st STRATEGY can be used only (this will decrease lines of code)
  
  barplt_df_wide <- barplt_df %>%  ## If you select "charge" column, it will bring multiple rows for one seq
    select(Sequence,Experiment,Intensity,pep_with_pos,Protein.Names) %>%
    pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
    mutate(species=selected_spcies) %>%
    relocate(exp_design,.after = "pep_with_pos")
  
  barplt_df_ecoli_wide <- barplt_df_ecoli %>% 
    select(Sequence,Experiment, Protein.Names,Intensity) %>%
    pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
    mutate(species=background_species) %>% 
    relocate(exp_design,.after = c("Protein.Names","Sequence"))
  
  #################################################
  
  if(loc_filter_opt == TRUE){
    
    barplt_df_wide <- barplt_df %>%  ## If you select "charge" column, it will bring multiple rows for one seq
      select(Sequence,Experiment,Intensity,pep_with_pos,Protein.Names) %>%
      pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
      relocate(exp_design,.after = where(is.character)) %>%
      left_join(site_prob, by="pep_with_pos") %>%
      mutate(species=selected_spcies) %>%
      filter(max_value >= loc_filter) %>%
      relocate(pep_with_pos, .after = Sequence) 
    #relocate(exp_design,.after = "Proteins")
  }else{
    barplt_df_wide <- barplt_df %>%  ## If you select "charge" column, it will bring multiple rows for one seq
      select(Sequence,Experiment,Intensity,pep_with_pos,Protein.Names) %>%
      pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
      left_join(site_prob, by="pep_with_pos") %>%
      mutate(species=selected_spcies) %>%
      #filter(max_value >= loc_filter) %>%
      relocate(pep_with_pos, .after = Sequence) 
  }
  
  #################################################
  
  ####### ADDITIONAL PLOT TO DISPLAY MISSING and UNEXPECTED PEPTIDES ########
  df_merge_syn <- barplt_df_wide %>%
    full_join(pep_list_w_theo,by="pep_with_pos") %>% 
    mutate_at("Pool", ~replace_na(.,"Unexpected")) %>%
    mutate(Pool= ifelse(is.na(Protein.Names),"missing",Pool)) %>%
    select(pep_with_pos,starts_with(exp_design),Pool) %>%
    mutate(soft_name=software_name, ion_mobility=acquisiton_type)
  
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
  
  write.table(df_merge_syn,file = paste0(file_path,"Count_of_missing_unexpected_correct_phospho-sites_",
                                         software_name,"_Experiment",exp_id,".txt"),
              sep = "\t",row.names = F)
  
  ########## ########## ########## ########## ########## ########## ########## ##########
  
  # Nothing is changed
  filtered_abundances<-barplt_df_wide[rowSums(!is.na(select(barplt_df_wide,starts_with(exp_design))))>0,]
  filtered_abundances_ecoli <-barplt_df_ecoli_wide[rowSums(!is.na(select(barplt_df_ecoli_wide,starts_with(exp_design))))>0,]
  library(kableExtra)
  
  na_phospho_selected <- apply(X = is.na(filtered_abundances %>% select(Sequence, Protein.Names,starts_with("E2"))), MARGIN = 2, FUN = sum)
  na_ecoli <- apply(X = is.na(filtered_abundances_ecoli %>% select(Sequence,Protein.Names,starts_with("E2"))), MARGIN = 2, FUN = sum)
  
  na_table <- bind_rows(na_phospho_selected,na_ecoli)
  na_table$species <- c(selected_spcies,background_species)
  na_table$total <- c(dim(filtered_abundances)[1],dim(filtered_abundances_ecoli)[1])
  
  na_table %>% select(-Sequence) %>%
    kbl(caption = paste("Number of NA values across all samples \n Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) %>%
    kable_material(c("striped", "hover")) %>%
    kable_styling(bootstrap_options = "striped", full_width = F, position = "left", font_size = 12) %>%
    kable_minimal(full_width = F) %>%
    footnote(general =  paste("This table was created after elimination of multiple charges by selecting either phospho-sites of", selected_spcies, "and", background_species, "sequences \n that has the highest abundace."),
             # number = c("Footnote 1; ", "Footnote 2; "),
             # alphabet = c("Footnote A; ", "Footnote B; "),
             # symbol = c("Footnote Symbol 1; ", "Footnote Symbol 2")
             footnote_as_chunk = T, title_format = c("italic", "underline")) %>%
    #as_image(width = 8) %>%
    save_kable(paste0(file_path,"outputs_with_new_script/table1.png"))
  
  ## REMOVE SEQUENCE COLUMN AFTER NA TABLE
  filtered_abundances <- filtered_abundances %>% select(!Sequence)
  
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
    mutate(sample_id_seq = paste(Sequence, sample_ids, sep = "_"))
  
  ecoli_density_plt<- quant_peptides_ECOLI_density_plt %>%
    select(contains(c("sample_ids","Intensity","species"))) 
  
  
  colnames(abundances_rowMeans) <- paste0("mean abundance",1:sample_size)
  
  
  ### RATIO SUPPRESION ASSESSMENT ### 
  ### For this calculation, wrong localizations were removed. 
  
  pep_list_w_theo_sel <- pep_list_w_theo %>%
    select(pep_with_pos,isomericity,Pool,pool_id) 
  
  ratio_supp_fixed <- filtered_abundances %>%
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
  
  ratio_supp_spiked <- filtered_abundances %>%
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
  
  
  
  
  
  #### MEAN ABUNDANCE RATIO WITH  DENSITY PLOT ####
  ### BEFORE IMPUTATION ###
  quant_phospho_density_plt <- filtered_abundances %>%
    select(!starts_with("E")) %>%
    bind_cols(abundances_rowMeans) %>% 
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
                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
  
  # Impute missing values
  for (j in 1:length(imputed_values_vec)){
    # Number NA
    #num_NA <- length(abundances_for_impute_all[,j+2][is.na(abundances_for_impute_all[,j+2])])
    
    filtered_abundances[,j+2][is.na(filtered_abundances[,j+2])] <- imputed_values_vec[j]
    filtered_abundances_ecoli[,j+2][is.na(filtered_abundances_ecoli[,j+2])] <- imputed_values_vec[j]
    
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
  abundances_all_aft_imputation <- filtered_abundances_ecoli %>%
    rename_with(~ paste0("pep_with_pos"), matches("^Seq")) %>%
    bind_rows(filtered_abundances) 
  
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
  
  final_imputed_data_syn <- final_imputed_data %>% filter(grepl(selected_spcies,species))
  
  final_imputed_data_ecoli <- final_imputed_data  %>% filter(!grepl(selected_spcies, species))
  
  write.table(final_imputed_data, 
              file =paste0(file_path,"/outputs_with_new_script/final_imputed_data_DIANN_",
                           exp_id, 
                           acquisiton_type,".txt"),sep = "\t",row.names = F)
  
  df_merge <- final_imputed_data_syn %>%
    left_join(pep_list_w_theo,by="pep_with_pos") %>% 
    mutate_at("Pool", ~replace_na(.,"Unexpected")) %>%
    bind_rows(final_imputed_data_ecoli) %>%
    mutate_at("Pool", ~replace_na(.,background_species)) 
  
  #write.table(final_imputed_data, file = "final_imputed_normalized_data_PAL _T_cell_Exp3_( 5 conc 3reps)_NoFAIMS_DDA_with_cont_230206_2023-02-07_0947.txt",sep = "\t",row.names = F)
  
  df_mean_ab_after_impt <- df_merge %>% 
    select(contains("aft_imp") | contains("species"),Pool) %>%
    tibble() %>% 
    #separate(accession, into = c("uniprot_id", "species", "position"), remove = F) %>%
    pivot_longer(cols = contains("aft_imp"),
                 names_to = "Mean_abundance",
                 values_to = "values")# %>%
  #mutate_at("Pool", ~replace_na(.,background_species))
  
  
  df_FC_ratio_after_impt <- df_merge %>% 
    select(starts_with("exp_") | contains("species"),Pool) %>%
    #separate(accession, into = c("uniprot_id", "species", "position"), remove = F) %>%
    tibble() %>% 
    pivot_longer(cols = starts_with("exp_"),
                 names_to = "exp_FC",
                 values_to = "values")# %>%
  #mutate_at("Pool", ~replace_na(.,background_species))
  
  plot4 <- gg_density(data_set = df_mean_ab_after_impt, 
                   x_df = df_mean_ab_after_impt$values,
                   fill_df = df_mean_ab_after_impt$Mean_abundance,
                   color_df = df_mean_ab_after_impt$Pool,
                   header="Distribution of mean abundance of every sample after imputation",
                   facet_df = "Mean_abundance",
                   x_lab = "log10(values)",
                   color_lab= "",
                   fill_lab = "Sample Names",
                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
  ############################## ############################

  ### BOX-PLOT: Experimental Quantity Ratio of Phospho Peptides  
  
  df_FC_ratio_absErr_after_impt <- df_merge %>% 
    select(pep_with_pos,starts_with("exp_") | contains("species"),Pool,isomericity) %>%
    #separate(accession, into = c("uniprot_id", "species", "position"), remove = F) %>%
    tibble() %>% 
    pivot_longer(cols = starts_with("exp_"),
                 names_to = "exp_FC",
                 values_to = "values") %>%
    mutate_at("Pool", ~replace_na(.,background_species)) %>%
    mutate(actual_ratio_val= case_when(grepl(comparisons[1],exp_FC) ~ as.numeric(actual_ratio[1]),
                                       grepl(comparisons[2],exp_FC) ~as.numeric(actual_ratio[2]),
                                       grepl(comparisons[3],exp_FC) ~as.numeric(actual_ratio[3]),
                                       grepl(comparisons[4],exp_FC) ~as.numeric(actual_ratio[4]))) %>%
    drop_na(actual_ratio_val) %>%
    mutate(actual_ratio_val=ifelse(Pool=="ECOLI",1,actual_ratio_val)) %>%
    mutate(actual_ratio_val=ifelse((Pool=="Fixed_isomeric" | Pool=="Fixed_non-isomeric"),1,actual_ratio_val)) %>%
    mutate(log2_fc_values =log2(values)) %>%
    filter(!grepl("Unexpected",Pool)) %>%
    #filter(!grepl("ECOLI",species)) %>%
    mutate(Pool=ifelse(Pool=="Diluted",paste0(Pool,"_",isomericity),Pool)) %>%
    #mutate(lower_bound = quantile(log2_fc_values, 0.25) - 1.5 * IQR(log2_fc_values),
    # upper_bound = quantile(log2_fc_values, 0.75) + 1.5 * IQR(log2_fc_values)) %>%
    mutate(acq_type=acquisiton_type) %>%
    mutate(soft_name=software_name)
  
  ranges <- df_FC_ratio_absErr_after_impt %>% 
    group_by(exp_FC,Pool) %>%
    summarise(firstQ=(quantile(log2_fc_values,probs = 0.25)),
              iqr_val=IQR(log2_fc_values),
              thirdQ=(quantile(log2_fc_values,probs = 0.75))) %>%
    mutate(lower_bound=firstQ - 1.5*iqr_val) %>%
    mutate(upper_bound=thirdQ + 1.5*iqr_val) %>%
    ungroup() %>%
    mutate(Pool_FC=paste0(Pool,"_",exp_FC)) %>%
    select(Pool_FC,lower_bound,upper_bound)
  
  df_FC_ratio_absErr_after_impt_filt <- df_FC_ratio_absErr_after_impt %>% #filter(log2_fc_values >= lower_bound & log2_fc_values <= upper_bound) %>%
    mutate(Pool_FC = paste0(Pool,"_",exp_FC)) %>%
    left_join(ranges,by = "Pool_FC") %>%
    group_by(Pool_FC) %>%
    filter(log2_fc_values >= lower_bound & log2_fc_values <= upper_bound) %>%
    mutate(abs_err = abs(log2(actual_ratio_val)-log2_fc_values)) %>%
    mutate(rel_err=(abs_err/abs(actual_ratio_val))*100) 
  
  
  
  
  write.table(df_FC_ratio_absErr_after_impt_filt,file = paste0(file_path,"AbsError_FC",software_name,"_",acquisiton_type,".txt"),sep = 
                "\t",col.names = T,row.names = F)
  
  
  library(gghalves)
  plot18 <- gg_half_boxplt_exp_ratio_nolog(data_set = df_FC_ratio_absErr_after_impt_filt , 
                                           x_df = df_FC_ratio_absErr_after_impt_filt$exp_FC,
                                           y_df = df_FC_ratio_absErr_after_impt_filt$rel_err,
                                           fill_df = df_FC_ratio_absErr_after_impt_filt$Pool,
                                           header="Relative Absolute Error (%) using Experimental Quantity Ratio of Phospho Peptides",
                                           x_lab="Sample Names",
                                           y_lab="Relative Absolute Error (%)",
                                           fill_lab = "Pool Names",
                                           
                                           subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) + 
    #scale_y_continuous(limits = c(0,200)) + #breaks = seq(from =0, to=100,by=10)) +
    
    scale_fill_manual(values = c("#6A51A3","#6A51A3","grey68","#E6550D","#E6550D"))+ #"grey68"
    scale_color_manual(values =  c("#6A51A3","#6A51A3","grey68","#E6550D","#E6550D")) #"grey68"
  #geom_half_point(alpha = 1, show.legend = TRUE, aes(color=df_FC_ratio_absErr_after_impt_filt$Pool,shape=df_FC_ratio_absErr_after_impt_filt$Pool))#+ scale_shape_manual(values = c(15,24,8,15,24)) 
  
  plot19 <- gg_violin_exp_ratio_nolog(data_set = df_FC_ratio_absErr_after_impt_filt, 
                                      x_df = df_FC_ratio_absErr_after_impt_filt$exp_FC,
                                      y_df = df_FC_ratio_absErr_after_impt_filt$log2_fc_values,
                                      fill_df = df_FC_ratio_absErr_after_impt_filt$Pool,
                                      trim=TRUE,
                                      header="Fold change Ratio of every sample based on pool names after imputation",
                                      x_lab="Sample Names",
                                      y_lab="Fold change values",
                                      fill_lab = "Pool Names",
                                      subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name,"\n after filtering outliers")
  ) + scale_fill_manual(values = c("#6A51A3","#6A51A3","grey68","#E6550D","#E6550D"))+ #"grey68"
    scale_color_manual(values =  c("#6A51A3","#6A51A3","grey68","#E6550D","#E6550D")) 
  #scale_y_continuous(limits = c(0,200)) + #breaks = seq(from =0, to=100,by=10)) 
  
  
  ############################## ############################
  
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
  
  # p15 <- gg_raincloud(data_set = df_mean_ab_after_impt,
  #                    x_df = df_mean_ab_after_impt$Mean_abundance,
  #                    y_df = df_mean_ab_after_impt$values,
  #                    fill_df = df_mean_ab_after_impt$Mean_abundance,
  #                    header = "Distribution of mean abundance of every sample after imputation",
  #                    x_lab = "Sample Names",
  #                    y_lab = " Density of log10(Mean Abundance)",
  #                    fill_lab = "Sample Names",
  #                    caption_lab = "",
  #                    subtitle_txt = "")
  # 
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
  
  plot6 <- gg_boxplt_exp_ratio(data_set = df_FC_ratio_after_impt, 
                            x_df = df_FC_ratio_after_impt$exp_FC,
                            y_df = df_FC_ratio_after_impt$values,
                            fill_df = df_FC_ratio_after_impt$Pool,
                            header="Experimental Quantity Ratio of Phospho Peptides",
                            x_lab="Sample Names",
                            y_lab="Abundance Ratios",
                            fill_lab = "Sample Names",
                            subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
  
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
          ttest_func(select(stat_analysis, contains("A1-") & contains("log10_"))[j,],
                     select(stat_analysis, contains(paste0("A", i, "-")) & contains("log10_"))[j,])
        } else if (test_type == "wilcoxon") {
          wilcox.test(select(stat_analysis, contains("A1-") & contains("log10_"))[j,],
                      select(stat_analysis, contains(paste0("A", i, "-")) & contains("log10_"))[j,])
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
    merge_stat_df_final <- merge_stat_df
    
    
  }else if(test_type=="limma"){
    library(limma)
    design_matrix <- model.matrix(~factor(c(rep(2,num_reps),rep(1,num_reps))))
    merge_stat_df <-NULL
    for ( i in 1:sample_size){
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
      unite('Pool_new',Pool.x, isomericity,sep = "_",remove = FALSE) %>%
      mutate(Pool_new=ifelse(Pool_new=="Fixed_isomeric_isomeric","Fixed_isomeric",Pool_new)) %>%
      mutate(Pool_new=ifelse(Pool_new=="Fixed_non-isomeric_nonisomeric","Fixed_non-isomeric",Pool_new)) %>%
      unite('new_col_coloring',Pool_new,ratio,sep = "_",remove = FALSE) %>%
      #mutate(new_col_coloring = if_else(grepl("Fixed_isomeric", new_col_coloring), "Fixed_isomeric", new_col_coloring)) %>%
      #mutate(new_col_coloring = if_else(grepl("Fixed_non-isomeric_nonisomeric", new_col_coloring), "Fixed_nonisomeric", new_col_coloring)) %>%
      mutate(new_col_coloring = if_else(grepl("Unexpected", new_col_coloring), "Unexpected", new_col_coloring)) %>%
      rename(pep_with_pos=pep_with_pos.x) %>% rename(Pool=Pool.x) %>%
      #unite(Pool_new, Pool.x, isomericity,sep = "_",remove = FALSE) %>%
      #unite('new_col_coloring',Pool_new,ratio,sep = "_",remove = FALSE) %>%
      #mutate(new_col_coloring = if_else(grepl("Fixed_multi", new_col_coloring), "Fixed_multi", new_col_coloring)) %>%
      #mutate(new_col_coloring = if_else(grepl("Fixed_mono", new_col_coloring), "Fixed_mono", new_col_coloring)) %>%
      #mutate(new_col_coloring = if_else(grepl("Unexpected", new_col_coloring), "Unexpected", new_col_coloring)) %>%
      #rename(pep_with_pos=pep_with_pos.x) %>% rename(Pool=Pool.x) #%>%
      filter(!grepl("A1/A6",fold_change_comp))
    
  }else{
    print("Statistical test could not be assessed. Check the input files!")
  }
  
  ##TODO: Find more logical way to map these values!!!
  
  ranges <- merge_stat_df_final %>% 
    mutate(log2_fc_values =log2(fold_change_values)) %>%
    group_by(fold_change_comp,Pool_new) %>%
    summarise(firstQ=(quantile(log2_fc_values,probs = 0.25)),
              iqr_val=IQR(log2_fc_values),
              thirdQ=(quantile(log2_fc_values,probs = 0.75))) %>%
    mutate(lower_bound=firstQ - 1.5*iqr_val) %>%
    mutate(upper_bound=thirdQ + 1.5*iqr_val) %>%
    ungroup() %>%
    mutate(Pool_FC=paste0(Pool_new,"_",fold_change_comp)) %>%
    select(Pool_FC,lower_bound,upper_bound)
  
  merge_stat_df_final_filt <- merge_stat_df_final %>% 
    mutate(Pool_FC=paste0(Pool_new,"_",fold_change_comp)) %>%
    mutate(log2_fc_values =log2(fold_change_values)) %>%
    left_join(ranges,by = "Pool_FC") %>%
    group_by(Pool_FC) %>%
    filter(log2_fc_values >= lower_bound & log2_fc_values <= upper_bound) %>%
    ungroup()
  
  ##TODO: Find more logical way to map these values!!!
  comparisons <- comparisons[-5]
  
  if(length(comparisons) == 4){
    actual_ratio_col <- merge_stat_df_final %>%
      select(fold_change_comp) %>% distinct() %>%
      mutate(actual_ratio_val= case_when(grepl(comparisons[1],fold_change_comp) ~ actual_ratio[1],
                                         grepl(comparisons[2],fold_change_comp) ~actual_ratio[2],
                                         grepl(comparisons[3],fold_change_comp) ~actual_ratio[3],
                                         grepl(comparisons[4],fold_change_comp) ~actual_ratio[4]))
    actual_ratio_col_filt <- merge_stat_df_final_filt %>%
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
    actual_ratio_col_filt <- merge_stat_df_final_filt %>%
      select(fold_change_comp) %>% distinct() %>%
      mutate(actual_ratio_val= case_when(grepl(comparisons[1],fold_change_comp) ~ actual_ratio[1],
                                         grepl(comparisons[2],fold_change_comp) ~actual_ratio[2],
                                         grepl(comparisons[3],fold_change_comp) ~actual_ratio[3],
                                         grepl(comparisons[4],fold_change_comp) ~actual_ratio[4]))#,
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
    group_by(fold_change_comp, new_col_coloring) %>%
    filter(P.Value < 0.05) %>% 
    count(new_col_coloring) %>% left_join(actual_ratio_col)
  
  point_count_y_axis_filt <- merge_stat_df_final_filt %>%
    group_by(fold_change_comp, new_col_coloring) %>%
    filter(P.Value < 0.05) %>% 
    count(new_col_coloring) %>% left_join(actual_ratio_col_filt)
  
  
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
  
  # Apply the function to calculate y_pos within each group
  point_count_y_axis_filt$y_pos <- unlist(by(point_count_y_axis_filt$fold_change_comp, 
                                             point_count_y_axis_filt$fold_change_comp, calculate_y_pos))
  
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
  
  actual_ratio_col <- actual_ratio_col %>% bind_cols(col)
  
  actual_ratio_col_filt <- actual_ratio_col_filt %>% bind_cols(col)
  
  mapped_coloring_dat <- as.data.frame(mapped_coloring)
  
  point_count_y_axis <- mapped_coloring_dat %>% 
    bind_cols(labels) %>%
    rename(colors=1,new_col_coloring=2) %>%
    right_join(point_count_y_axis,by="new_col_coloring")
  
  point_count_y_axis_filt <-mapped_coloring_dat %>% 
    bind_cols(labels) %>%
    rename(colors=1,new_col_coloring=2) %>%
    right_join(point_count_y_axis_filt,by="new_col_coloring")
  
  plot20 <- ggplot(merge_stat_df_final_filt,aes(x =log2(merge_stat_df_final_filt$fold_change_values), y = -log10(merge_stat_df_final_filt$P.Value))) +
    geom_point(size = 4, aes(color = new_col_coloring,shape=Pool)) + # Pool_new might be use to 
    #scale_shape_identity() +                                        # differentiate peptides are found as Unexpected and reference 
    # in their associated concentration 
    # Pool can be used to show only difference btw pools same as coloring 'less complex visualization)
    #geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "red") +
    #scale_fill_manual(values=setNames(mapped_coloring, labels)) + 
    #geom_line(aes(color =setNames(mapped_coloring, labels)), size = 1) +  # Add color aesthetic to geom_line()
    scale_color_manual(values =setNames(mapped_coloring, labels)) +
    #scale_y_continuous(limits = c(0, max(-log10(merge_stat_df_final$adj.P.Val))), breaks = seq(0, max(-log10(merge_stat_df_final$adj.P.Val)), by = 0.8)) +
    #scale_x_continuous(limits = c(min(log2(merge_stat_df_final$fold_change_values)),max(log2(merge_stat_df_final$fold_change_values)))) +#facet_wrap(~ratio) +
    
    #scale_x_continuous(breaks = seq(from =round(min(log2(merge_stat_df_final$fold_change_values))), to=(round(max(log2(merge_stat_df_final$fold_change_values)))+2),by=2)) +
    #scale_y_continuous(breaks = seq(from =round(min(-log10(merge_stat_df_final$P.Value))), to=(round(max(-log10(merge_stat_df_final$P.Value)))+2),by=2)) +
    scale_y_continuous(
      limits = c(0,12), 
      breaks = seq(0, 12,2)
    ) +
    scale_x_continuous(
      limits = c(-9,11), 
      breaks = seq(-9, 11,2)
    ) +
    
    #scale_y_continuous(breaks = seq(0, max(-log10(volcano_final1$pvalues_value)), length.out = 21)) +
    theme_bw() +
    theme(legend.text = element_text(size = 45),
          axis.title.x = element_text(size = 45),
          axis.title.y = element_text(size = 45),
          plot.title = element_text(size = 55),
          legend.title = element_text(size = 45),
          axis.text.x = element_text(size = 45),
          axis.title = element_text(size = 45),
          axis.text.y = element_text(size = 45),
          plot.subtitle = element_text(size = 45)) +
    labs( y= "-log10(p values)", x="log2(fold change)",title = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name), subtitle = paste("Limma was used \n", subtitle,"\n after filtering outliers")) +
    geom_vline(data = actual_ratio_col_filt, aes(xintercept = log2(actual_ratio_val), show.legend = FALSE),color=col,size=1.5) +
    geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "red",size=1.5) + 
    geom_label(data = point_count_y_axis_filt, aes(x = log2(actual_ratio_val), y = y_pos, fill=new_col_coloring,label = n),color="white",size=14,show.legend = FALSE) +
    scale_fill_manual(values =setNames(mapped_coloring, labels))
  
  
  plot9 <- ggplot(merge_stat_df_final,aes(x =log2(merge_stat_df_final$fold_change_values), y = -log10(merge_stat_df_final$P.Value))) +
    geom_point(size = 4, aes(color = new_col_coloring,shape=Pool)) + # Pool_new might be use to 
    #scale_shape_identity() +                                        # differentiate peptides are found as Unexpected and reference 
    # in their associated concentration 
    # Pool can be used to show only difference btw pools same as coloring 'less complex visualization)
    #geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "red") +
    #scale_fill_manual(values=setNames(mapped_coloring, labels)) + 
    #geom_line(aes(color =setNames(mapped_coloring, labels)), size = 1) +  # Add color aesthetic to geom_line()
    scale_color_manual(values =setNames(mapped_coloring, labels)) +
    scale_y_continuous(
      limits = c(0,12), 
      breaks = seq(0, 12,2)
    ) +
    scale_x_continuous(
      limits = c(-9,11), 
      breaks = seq(-9, 11,2)
    ) +
    
    #scale_y_continuous(limits = c(0, max(-log10(merge_stat_df_final$adj.P.Val))), breaks = seq(0, max(-log10(merge_stat_df_final$adj.P.Val)), by = 0.8)) +
    #scale_x_continuous(limits = c(min(log2(merge_stat_df_final$fold_change_values)),max(log2(merge_stat_df_final$fold_change_values)))) +#facet_wrap(~ratio) +
    #scale_x_continuous(breaks = seq(from =round(min(log2(merge_stat_df_final$fold_change_values))), to=(round(max(log2(merge_stat_df_final$fold_change_values)))+2),by=2)) +
    #scale_y_continuous(breaks = seq(from =round(min(-log10(merge_stat_df_final$P.Value))), to=(round(max(-log10(merge_stat_df_final$P.Value)))+2),by=2)) +
    #scale_y_continuous(limits = c(0, 11), breaks = seq(0,11, by = 2)) +
    #scale_x_continuous(limits = c(-7,11),breaks = seq(-7,11, by = 2))+
    #scale_y_continuous(breaks = seq(0, max(-log10(volcano_final1$pvalues_value)), length.out = 21)) +
    theme_bw() +
    theme(legend.text = element_text(size = 45),
          axis.title.x = element_text(size = 45),
          axis.title.y = element_text(size = 45),
          plot.title = element_text(size = 55),
          legend.title = element_text(size = 45),
          axis.text.x = element_text(size = 45),
          axis.title = element_text(size = 45),
          axis.text.y = element_text(size = 45),
          plot.subtitle = element_text(size = 45)) +
    labs( y= "-log10(p values)", x="log2(fold change)",title = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name), subtitle = paste("Limma was used \n", subtitle)) +
    geom_vline(data = actual_ratio_col, aes(xintercept = log2(actual_ratio_val), show.legend = FALSE),color=col,size=1.5) +
    geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "red",size=1.5) + 
    geom_label(data = point_count_y_axis, aes(x = log2(actual_ratio_val), y = y_pos, fill=new_col_coloring,label = n),color="white",size=14,show.legend = FALSE) +
    scale_fill_manual(values =setNames(mapped_coloring, labels))
  
  merge_stat_df_final_text <- merge_stat_df_final %>% mutate(soft_name=paste0(software_name)) %>% mutate(acq_type=paste0(acquisiton_type))
  
  write.table(merge_stat_df_final_text,file = paste0(file_path,"volcano_plot_",software_name,"_",acquisiton_type,".txt"),sep = 
                "\t",col.names = T,row.names = F)
  
  merge_stat_df_final_text <- merge_stat_df_final_filt %>% mutate(soft_name=paste0(software_name)) %>% mutate(acq_type=paste0(acquisiton_type))
  write.table(merge_stat_df_final_text,file = paste0(file_path,"volcano_plot_",software_name,"_",acquisiton_type,"filtered.txt"),sep = 
                "\t",col.names = T,row.names = F)
  
  
  df_roc_filt <- merge_stat_df_final_filt %>%
    select(pep_with_pos, Pool,P.Value) %>%
    mutate(Pool = if_else(grepl("Diluted_isomeric", Pool), "Diluted", Pool)) %>%
    mutate(Pool = if_else(grepl("Diluted_nonisomeric", Pool), "Diluted", Pool))
  #filter(!grepl("Unexpected",Pool))
  ### ROC analysis custom func
  df_roc_order_filt <- df_roc_filt[order(df_roc_filt$P.Value),]
  
  df_roc_func_filt <- compute_roc_curve(df=df_roc_order_filt, flag = "Diluted",expected = length(comparisons)*size_variying_pep_size)
  
  df_roc <- merge_stat_df_final %>%
    select(pep_with_pos, Pool,P.Value) %>%
    select(pep_with_pos, Pool,P.Value) %>%
    mutate(Pool = if_else(grepl("Diluted_isomeric", Pool), "Diluted", Pool)) %>%
    mutate(Pool = if_else(grepl("Diluted_nonisomeric", Pool), "Diluted", Pool))
  #filter(!grepl("Unexpected",Pool))
  
  ### ROC analysis custom func
  df_roc_order <- df_roc[order(df_roc$P.Value),]
  
  df_roc_func <- compute_roc_curve(df=df_roc_order, flag = "Diluted",expected = (length(comparisons)*size_variying_pep_size))
  
  
  intended_dir <-paste0(file_path,"outputs_new")
  
  if(dir.exists(intended_dir)){
    new_path <- intended_dir
    
  }else{
    dir.create(intended_dir)
    new_path <- list.dirs(intended_dir)
    
  }
  
  write.table(df_roc_func, file = paste0(new_path,"/new_custom_Roc_analysis_",exp_id,"_",software_name,"_",".txt"),sep = "\t",row.names = F)
  write.table(df_roc_func_filt, file = paste0(new_path,"/new_custom_Roc_analysis_",exp_id,"_",software_name,"_","filtered.txt"),sep = "\t",row.names = F)
  
  
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
  
  
  #roc_raw_non_var <- roc(df_roc$non_var, df_roc$P.Value)
  #tpr_and_fpr_non_var  <- cbind(roc_raw_non_var$sensitivities,
  #roc_raw_non_var$specificities,
  #"Non-variant Pool")
  
  roc_plt_df <- as.data.frame(tpr_and_fpr_variant) %>% 
    #bind_rows(as.data.frame(tpr_and_fpr_non_var)) 
    bind_cols(software_name)
  
  colnames(roc_plt_df) <- c("sensitivity", "fpr","Pool_type","Software_name")
  
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
  
  write.table(roc_plt_df, file = paste0(file_path,"pRoc_analysis_",exp_id,"_",software_name,"_",".txt"),sep = "\t",row.names = F)
  
 
  # #### ROC Analysis by using pROC
  # df_roc <- merge_stat_df_final %>%
  #   select(pep_with_pos, Pool,P.Value)
  # 
  # df_roc$variant <- ifelse(df_roc$Pool == "Others", TRUE, FALSE)
  # df_roc$non_var <- ifelse(df_roc$Pool == "ISO-REF", TRUE, FALSE)
  # 
  # library(pROC)
  # # Calculate ROC curve for raw p-values
  # roc_raw_variant <- roc(df_roc$variant, df_roc$P.Value)
  # tpr_and_fpr_variant  <- cbind(roc_raw_variant$sensitivities,
  #                               roc_raw_variant$specificities,
  #                               "Variant Pool")
  # 
  # 
  # roc_raw_non_var <- roc(df_roc$non_var, df_roc$P.Value)
  # tpr_and_fpr_non_var  <- cbind(roc_raw_non_var$sensitivities,
  #                               roc_raw_non_var$specificities,
  #                               "Non-variant Pool")
  # 
  # roc_plot_df <- as.data.frame(tpr_and_fpr_variant) %>% 
  #   bind_rows(as.data.frame(tpr_and_fpr_non_var)) %>% bind_cols(software_name)
  # 
  # colnames(roc_plot_df) <-    c("sensitivity", "specificity","Pool_type","Software_name")
  # 
  # 
  # p10 <- roc_plot_df %>% group_by(Pool_type) %>% 
  #   ggplot( aes(y=as.numeric(sensitivity), x = as.numeric(specificity), color=Pool_type)) +
  #   geom_path(size=1.5) +  scale_x_reverse() + theme_bw() +
  #   theme(legend.text = element_text(size = 20),
  #         axis.title.x = element_text(size = 20),
  #         axis.title.y = element_text(size = 20),
  #         plot.title = element_text(size = 25),
  #         legend.title = element_text(size = 20),
  #         axis.text.x = element_text(size = 20),
  #         axis.text.y = element_text(size = 20)) +
  #   scale_color_brewer(palette = "Dark2") +
  #   labs(y="True Positive Rate \n (Sensitivity)", x="False Positive Rate \n (Specificity)",
  #        title =  paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name), 
  #        subtitle = paste(subtitle), color="Pool Type")
  # 
  # write.table(roc_plot_df, file = paste0(file_path,"Roc_analysis_",exp_id,"_",software_name,"_",".txt"),sep = "\t",row.names = F) #acquisiton_type ## IT WAS TOO LONG-> GIVES AN ERROR
  # 
  
  
  plt_obj <- ls(pattern="plot")
  plt_obj <- plt_obj[!is.na(plt_obj)]
  sapply(1:length(plt_obj),function(x) ggsave(filename = paste0("p",x,".png"),
                                               width = 60, height = 45, 
                                               path = paste0(file_path,"/outputs_new/"),
                                               units = "cm",
                                               get(plt_obj[x]),
                                               device = "png", #".svg"
  ))
  
  
  
  
}
  