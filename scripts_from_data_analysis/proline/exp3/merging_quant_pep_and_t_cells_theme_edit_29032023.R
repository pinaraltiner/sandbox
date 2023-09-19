#
final_proline_pep_quant_analysis_bio <- function(file_path,
                                             file_name,
                                             sheet_name,
                                             theo_file_path,
                                             theo_file_name,
                                             sheet_theo_name,
                                             selected_spcies,
                                             acquisiton_type,
                                             exp_id,
                                             software_name,
                                             test_type,
                                             exp_design,
                                             num_reps){
  #library(PhosR)
  library(stringr)
  library(dplyr)
  library(data.table)
  library(openxlsx)
  library(ggplot2)
  library(tidyr)
  
  source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/ggplot/ggplot_functions.R")
  
  sample_size <- length(exp_design) / num_reps
  sample_names <- paste0("A",1:sample_size)
  comparisons <- NULL
  for (i in 1:sample_size){
    tmp <- paste0(sample_names[1], "_vs_",sample_names[i])
    comparisons[i] <- tmp
    rm(tmp)
  }
  
  quant_peptides <- read.xlsx(paste0(file_path,file_name), sheet = sheet_name)
  
  abundances_for_impute <- quant_peptides %>% 
    select(starts_with("abundance_")) %>% ## spectrum_title remove it because it was not make it as rownames (has duplicates)
    rename_with(~exp_design,matches("abundance"))
  
  quant_peptides_cor_abun <- quant_peptides %>% 
    rename_with(~exp_design,matches("^abundance"))
  
  quant_phospho_peptides <- quant_peptides_cor_abun %>% 
    filter(grepl(selected_spcies,accession)) %>% 
    filter(grepl("Phospho",modifications)) 
  
  quant_peptides_ECOLI <- quant_peptides_cor_abun %>% 
    filter(grepl("ECOLI",accession)) %>% 
    select(sequence,modifications,accession,spectrum_title,starts_with(exp_design))
    
  
  phospho_ptm_pos <- lapply(quant_phospho_peptides$modifications, function(each_ptm_protein_positions) {
    
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
  
  quant_phospho_peptides$pep_with_pos <- paste(quant_phospho_peptides$sequence, phospho_ptm_pos_df, sep = "_")

  
  # Nothing is changed
  filtered_abundances<-quant_phospho_peptides[rowSums(!is.na(select(quant_phospho_peptides,starts_with(exp_design))))>0,]
  filtered_abundances_ecoli <-quant_peptides_ECOLI[rowSums(!is.na(select(quant_peptides_ECOLI,starts_with(exp_design))))>0,]

  df_id_pep <- filtered_abundances %>% 
    select(sequence,modifications, pep_with_pos, starts_with(exp_design),accession) %>%
    tibble() %>%
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
  pivot_longer(cols = starts_with("E3"), 
               values_to = "intensity",
               names_to = "sample_ids",
               values_drop_na = T) %>%
    separate(sample_ids, into = c("Exp_id","Sample_id", "Rep_id"), sep = "_",remove = F)# %>% 
  #group_by(Sample_id) %>%
  #count()
  
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
  barplt_df <- df_id_pep %>%
    mutate(sample_rep_id_seq = paste(pep_with_pos, Sample_id,Rep_id, sep = "_"))%>%
    group_by(sample_rep_id_seq) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    slice(which.max(intensity)) %>%
    ungroup()
  
  ## SAME STRATEGIES ABOVE (3rd) WAS APPLIED TO BACKGROUND AS WELL
  barplt_df_ecoli <- filtered_abundances_ecoli %>% 
    select(sequence,modifications, sequence, starts_with(exp_design),accession) %>%
    pivot_longer(cols = starts_with("E3"), 
                 values_to = "intensity",
                 names_to = "sample_ids",
                 values_drop_na = T) %>%
    separate(sample_ids, into = c("Exp_id","Sample_id", "Rep_id"), sep = "_",remove = F) %>%
    mutate(sample_rep_id_seq = paste(sequence, Sample_id,Rep_id, sep = "_")) %>%
    group_by(sample_rep_id_seq,sample_ids) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    slice(which.max(intensity)) %>%
    ungroup()
  
  # barplt_df_ecoli <- filtered_abundances_ecoli %>% 
  #   pivot_longer(cols = starts_with("E3"), 
  #                values_to = "intensity",
  #                names_to = "sample_ids",
  #                values_drop_na = T) %>%
  #   separate(sample_ids, into = c("Exp_id","Sample_id", "Rep_id"), sep = "_") %>%
  #   mutate(sample_rep_id_seq = paste(sequence, Sample_id,Rep_id, sep = "_")) %>%
  #   filter(duplicated(sample_rep_id_seq)==FALSE)
  
  

  p1 <- gg_barplt_id_pep_count(data_set = barplt_df,
                         x_df = barplt_df$Sample_id,
                         fill_df = barplt_df$Rep_id,
                         ymax = 20000,
                         header = "Total number of quantified phospho-site across each sample",
                         caption_lab = "NA values are removed.",
                         x_lab = "Sample id",
                         fill_lab =  "Sample id",
                         y_lab = "Number of identified peptides",
                         subtitle_txt = "")
  
  
  p2 <- gg_barplt_id_pep_count(data_set = barplt_df_ecoli,
                         x_df = barplt_df_ecoli$Sample_id,
                         fill_df = barplt_df_ecoli$Rep_id,
                         ymax = 20000,
                         header = "Total number of quantified Ecoli across each sample",
                         caption_lab = "NA values are removed.",
                         x_lab = "Sample id",
                         fill_lab =  "Sample id",
                         y_lab = "Number of identified peptides",
                         subtitle_txt = "")
  
  ## THIS RESHAPING IS ONLY FOR ELIMINATION OF MULTIPLE PHOSPHO-SITES and ECOLI PEPTIDES
   ## ELIMINATION STEP IS NOT NECESSARY FOR ECOLI, 1st STRATEGY can be used only (this will decrease lines of code)
  barplt_df_wide <- barplt_df %>% 
    select(pep_with_pos, accession, sample_ids,intensity) %>%
    pivot_wider(names_from = "sample_ids",values_from = "intensity")

  barplt_df_ecoli_wide <- barplt_df_ecoli %>% 
    select(sequence, accession, sample_ids,intensity) %>%
    pivot_wider(names_from = "sample_ids",values_from = "intensity")
  
  ## DENSITY PLOT OF BEFORE IMPUTATION 
  
  sample_size <- 5
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
  

  quant_peptides_ECOLI_density_plot <- barplt_df_ecoli_wide %>%
    select(!starts_with("E")) %>%
    bind_cols(abundances_ecoli_rowMeans) %>%
    tibble() %>%
    rename_with(~ paste0("mean_abun",1:5), matches("^row")) %>%
    pivot_longer(cols = starts_with("mean"), 
                 values_to = "intensity",
                 names_to = "sample_ids",
                 values_drop_na = T) %>%
    mutate(sample_id_seq = paste(sequence, sample_ids, sep = "_"))
  
  ecoli_density_plot<- quant_peptides_ECOLI_density_plot %>%
    select(contains(c("sample_ids","intensity","accession"))) 
  
  
  colnames(abundances_rowMeans) <- paste0("mean_abun",1:sample_size)
  
  #### MEAN ABUNDANCE RATIO WITH  DENSITY PLOT ####
      ### BEFORE IMPUTATION ###
  quant_phospho_density_plot <- barplt_df_wide %>%
    select(!starts_with("E")) %>%
    bind_cols(abundances_rowMeans) %>% 
    #rename_with(~ paste0("mean_abun",1:5), matches("^row")) %>%
    tibble() %>% #mutate(pep_with_pos = sequence) %>% ###  At this stage, no need for phospho-position#   
    pivot_longer(cols = starts_with("mean"),
                 names_to = "sample_ids",
                 values_to = "intensity",
                 values_drop_na = T) %>%
    mutate(sample_id_seq = paste(pep_with_pos, sample_ids, sep = "_"))
  
  density_df <-quant_phospho_density_plot %>%
    select(c(sample_ids,intensity,accession)) %>%
    bind_rows(ecoli_density_plot) %>%
    separate(accession, into = c("prot_id","species"),remove = F)
  
  
  p3 <- gg_density(data_set = density_df, 
             x_df = density_df$intensity,
             fill_df = density_df$species,
             color_df = NULL,
             header="Distribution of mean abundance of every sample before imputation",
             facet_df = "sample_ids",
             x_lab = "log10(intensities)",
             color_lab= "",
             fill_lab = "species",
             subtitle_txt = "")
  
  
  library(kableExtra)
  na_phospho_mouse <- apply(X = is.na(filtered_abundances), MARGIN = 2, FUN = sum)
  na_ecoli <- apply(X = is.na(filtered_abundances_ecoli), MARGIN = 2, FUN = sum)
  
  na_table <- bind_rows(na_phospho_mouse,na_ecoli)
  na_table$species <- c("MOUSE","ECOLI")
  na_table$total <- c(dim(filtered_abundances)[1],dim(filtered_abundances_ecoli)[1])
  
  na_table %>% select(-common_col_for_merging) %>%
    kbl(caption = paste("Number of NA values across all samples \n Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) %>%
    kable_material(c("striped", "hover")) %>%
    kable_styling(bootstrap_options = "striped", full_width = F, position = "left", font_size = 12) %>%
    kable_minimal(full_width = F) %>%
    footnote(general = "This table was created after elimination of multiple charges by selecting either phospho-sites of mouse or background sequences \n that has the highest abundace.",
             # number = c("Footnote 1; ", "Footnote 2; "),
             # alphabet = c("Footnote A; ", "Footnote B; "),
             # symbol = c("Footnote Symbol 1; ", "Footnote Symbol 2")
             footnote_as_chunk = T, title_format = c("italic", "underline")) %>%
    #as_image(width = 8) %>%
    save_kable(paste0(file_path,"table1.png"))
  
  
  abundances_for_before_impt <- barplt_df_wide %>%
    bind_rows(barplt_df_ecoli_wide) %>% select(exp_design)
  
  # Calculate 1 percent quantile of each sample
  impute_values <- apply(abundances_for_impute, 2 , quantile , probs = 0.01 , na.rm = TRUE )
  
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
    rename_with(~ paste0("pep_with_pos"), matches("^seq")) %>%
    bind_rows(barplt_df_wide) 
    #select(sequence, accession) %>% #spectrum_title
    #bind_cols(abundances_for_before_impt)
  
  # Take log10 
  #log_10_abundances_only_for_impute <- log10(abundances_only_for_impute)
  #colnames(log_10_abundances_only_for_impute) <- paste0(colnames(abundances_only_for_impute),"_log10")
  
  ## No need is right now. I will cont with "log_10_abundances_only_for_impute" for rowMeans and FC
  #log_10_filtered_abundances <- cbind(filtered_abundances,log_10_abundances_only_for_impute)
  
  # Take mean of triplicates of each sample 
  # Ask sample_size additional parameter
  
  
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
  for(An in 2:cols){
    filtered_abundances_rowMeans[,paste0("exp_FC_A1/A",An)] <- filtered_abundances_rowMeans[,1]/filtered_abundances_rowMeans[,An]
    
  }
  # To calculate all binary combination in the data frame
  #mat <- do.call(cbind, lapply(cols, function(xj) 
  #  sapply(cols, function(xi) (filtered_abundances_rowMeans[, xj]/(filtered_abundances_rowMeans[, xj])))))
  #colnames(mat) <-  outer(names(filtered_abundances_rowMeans), names(filtered_abundances_rowMeans), paste0)
  
  final_imputed_data <- cbind(abundances_all_aft_imputation, filtered_abundances_rowMeans,filtered_abundances_log10,filtered_abundances_log10_rowMeans) #filtered_abundances
  
  #write.table(final_imputed_data, file = "final_imputed_normalized_data_PAL _T_cell_Exp3_( 5 conc 3reps)_NoFAIMS_DDA_with_cont_230206_2023-02-07_0947.txt",sep = "\t",row.names = F)
  
  df_mean_ab_after_impt <- final_imputed_data %>% 
    select(contains("aft_imp") | contains("accession")) %>%
    tibble() %>% 
    separate(accession, into = c("uniprot_id", "species"), remove = F) %>%
    pivot_longer(cols = contains("aft_imp"),
                 names_to = "Mean_abundance",
                 values_to = "values")
    
  df_FC_ratio_after_impt <- final_imputed_data %>% 
    select(starts_with("exp_")| contains("accession")) %>%
    separate(accession, into = c("uniprot_id", "species"), remove = F) %>%
    tibble() %>% 
    pivot_longer(cols = starts_with("exp_"),
                 names_to = "exp_FC",
                 values_to = "values")
  
  p4 <- gg_density(data_set = df_mean_ab_after_impt, 
                   x_df = df_mean_ab_after_impt$values,
                   fill_df = df_mean_ab_after_impt$Mean_abundance,
                   color_df = df_mean_ab_after_impt$species,
                   header="Distribution of mean abundance of every sample after imputation",
                   facet_df = "Mean_abundance",
                   x_lab = "log10(values)",
                   color_lab= "",
                   fill_lab = "Sample Names",
                   subtitle_txt = "")
  
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
  
  
  
  p5 <- gg_density(data_set = df_FC_ratio_after_impt,
                   x_df = df_FC_ratio_after_impt$values,
                   fill_df = df_FC_ratio_after_impt$exp_FC,
                   color_df = df_FC_ratio_after_impt$species,
                   header="Distribution of Fold change Ratio of every sample after imputation",
                   facet_df = "exp_FC",
                   x_lab = "log10(values)",
                   color_lab= "",
                   fill_lab = "Sample Names",
                   subtitle_txt = "")
  # 

  ### BOX-PLOT: Experimental Quantity Ratio of Phospho Peptides  
  
  p6 <- gg_boxplt_exp_ratio(data_set = df_FC_ratio_after_impt, 
                            x_df = df_FC_ratio_after_impt$exp_FC,
                            y_df = df_FC_ratio_after_impt$values,
                            fill_df = df_FC_ratio_after_impt$species,
                            header="Experimental Quantity Ratio of T-cell Phospho Peptides",
                            x_lab="Sample Names",
                            y_lab="Abundance Ratios",
                            fill_lab = "Sample Names",
                            subtitle_txt = "")
  
  
  ### HALF-BOX-PLOT & HALF-SCATTER-PLOT: Experimental Quantity Ratio of Synthetic Peptides  
  library(gghalves)
  
  p7 <- gg_half_boxplt_exp_ratio(data_set = df_FC_ratio_after_impt, 
                                 x_df = df_FC_ratio_after_impt$exp_FC,
                                 y_df = df_FC_ratio_after_impt$values,
                                 fill_df = df_FC_ratio_after_impt$species,
                                 header="Experimental Quantity Ratio of T-cell Phospho Peptides with Background",
                                 x_lab="Sample Names",
                                 y_lab="Abundance Ratios",
                                 fill_lab = "Sample Names",
                                 subtitle_txt = "")
  
  ### VIOLIN-PLOT: Experimental Quantity Ratio of Synthetic Peptides   
  
  ### TODO: fix y scaling without trimming 
  p8 <- gg_violin_exp_ratio(data_set = df_FC_ratio_after_impt, 
                            x_df = df_FC_ratio_after_impt$exp_FC,
                            y_df = df_FC_ratio_after_impt$values,
                            fill_df = df_FC_ratio_after_impt$species,
                            header="Experimental Quantity Ratio of T-cell Phospho Peptides with Background",
                            x_lab="Sample Names",
                            y_lab="Abundance Ratios",
                            fill_lab = "Sample Names",
                            trim=TRUE,
                            subtitle_txt = "")
  
  
  
  # ############################ # Do t-test # ############################
  # 
  # # The code below does t-test for each row. Because of that, multiple test correction (like BH, Bonferoni)
  # p_values_for_all_ratio <- NULL
  # for(k in 2:sample_size){
  #   
  #   assign(paste0("t_test_res_A1_to_A",k) , lapply(na.omit(correct_identifed_peps$mean_abundances_A1),
  #                                                  na.omit(correct_identifed_peps[,paste0("mean_abundances_A",k)], 
  #                                                          function(x,y) t.test(x,y,alternative = "two.sided", var.equal = TRUE))))
  #   assign(paste0("p_values_A1_vs_A",k), cbind(p_values_for_all_ratio, get(paste0("t_test_res_A1_to_A",k))$p.value))
  #   
  # }
  # 
  # #adjusted_p_values <- as.data.frame(p.adjust(p_values_A1_vs_A5, method = "BH", n = length(p_values_A1_vs_A5)))
  # adjusted_p_values <- lapply(p_values_for_all_ratio, function(x) p.adjust (x, method = "BH", n = length(x)))
  # 
  # #plot(correct_identifed_peps$mean_abundances_log10_A1,correct_identifed_peps$mean_abundances_log10_A5, pch = 16, col = "blue")
  # #abline(h = mean(na.omit(correct_identifed_peps$mean_abundances_log10_A1)) - mean(na.omit(correct_identifed_peps$mean_abundances_log10_A5)), col = "red")
  # 
  # boxplot(correct_identifed_peps$mean_abundances_log10_A1,correct_identifed_peps$mean_abundances_log10_A5)
  # ggplot(correct_identifed_peps, aes(x = mean_abundances_log10_A1, y = mean_abundances_log10_A5)) +
  #   geom_point(aes(color = "red", size = 5))# +
  # #scale_color_discrete(name = "") +
  # #scale_size_discrete(name = "Size")
  # # 
  # # Take -log10() of results
  # ggplot(log_10_filtered_abundances_rowMeans_A1_A5,aes(x=log_10_filtered_abundances_rowMeans_A1_A5$mean_abundances_log10_A1,
  #                                                      y =log_10_filtered_abundances_rowMeans_A1_A5$mean_abundances_log10_A5,
  # )) +
  #   geom_point()+
  #   #facet_wrap(vars(df_for_figure$common_col_for_merging))  + 
  #   #geom_smooth(method = "lm", colour = "green", fill = "green") +
  #   theme_light()
  # 
  # library(gginference)
  # ggttest(t.test(na.omit(correct_identifed_peps$mean_abundances_log10_A1),na.omit(correct_identifed_peps$mean_abundances_log10_A5), alternative = "two.sided", var.equal = TRUE))
  # 
  # 
  # 
  # 
  
  ttest_func <- function(x, y) {
    # if (sum(!is.na(x)) < 2 | sum(!is.na(y)) < 2) {
    #   return(NA)
    # }else{
    #   
    # }
    t.test(x, y,var.equal=TRUE,alternative = c("two.sided"))$p.value
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
  stat_analysis <- final_imputed_data %>%
    select(pep_with_pos,accession, starts_with("log10_") | starts_with("mean_log10_") | starts_with("exp_FC")) #spectrum_title
  
  rownames(stat_analysis) <- paste0(stat_analysis$common_col_for_merging,"@",stat_analysis$Pool,"@",(1:nrow(stat_analysis))) #stat_analysis$spectrum_title,"@"
  # if(test_type== "t.test" | test_type== "wilcoxon"){
  #     library(multtest)
  #   all_pvalues <- NULL
  #   all_adjust_pval <- NULL
    # for (i in 2:sample_size){
    #   p_values_tmp <- NULL
    #   
    #   for(j in 1:dim(stat_analysis)[1]){
    #     
    #     if(test_type=="t.test"){
    #       p_values_tmp[j] <- ttest_func(select(stat_analysis,contains("A1_") & contains("log10_"))[j,],     ### FOR DIFFERENT KIND OF EXP SETUP, 
    #                                     select(stat_analysis,contains(paste0("A",i,"_")) & contains("log10_"))[j,] ) ## It should be defined as an input.
    #       
    #     }else if(test_type=="wilcoxon"){
    #       p_values_tmp[j] <- wilcox.test(select(stat_analysis,contains("A1_") & contains("log10_"))[j,],     ### FOR DIFFERENT KIND OF EXP SETUP, 
    #                                      select(stat_analysis,contains(paste0("A",i,"_")) & contains("log10_"))[j,])
    #       ## It should be defined as an input.
    #     }
    #     
    #   }
    #   p_values_tmp <- as.data.frame(p_values_tmp)
    #   colnames(p_values_tmp) <- paste0("pvalues_A1","/","A",i)
    #   all_pvalues <- bind_cols(all_pvalues,p_values_tmp)
    #   rownames(all_pvalues) <- row.names(stat_analysis)
    #   rm(p_values_tmp)
    #   
    #   #### p-value adjust at a given FDR ####
    #   #procs<-c("Bonferroni","Holm","Hochberg","SidakSS","SidakSD","BH","BY","ABH","TSBH")
    #   adjust_pval_tmp <-  mt.rawp2adjp(all_pvalues[,paste0("pvalues_A1","/","A",i)],
    #                                    proc="BH", alpha=0.05)
    #   qval <- data.frame(adjust_pval_tmp$adjp,adjust_pval_tmp$index)[order(adjust_pval_tmp$index),2]
    #   qval <- as.data.frame(qval)
    #   colnames(qval) <- paste0("adjust_pval_A1","/","A",i)
    #   all_adjust_pval <- bind_cols(all_adjust_pval,qval)
    #   rownames(all_adjust_pval) <- row.names(stat_analysis)
    #   rm(adjust_pval_tmp,qval)
    #   
    # }
    if (test_type == "t.test" || test_type == "wilcoxon") {
        all_pvalues <- matrix(NA, nrow = nrow(stat_analysis), ncol = sample_size - 1)
        all_adjust_pval <- matrix(NA, nrow = nrow(stat_analysis), ncol = sample_size - 1)
        
        for (i in 2:sample_size) {
            col_A1 <- select(stat_analysis, contains("A1-") & contains("log10_"))
            col_Ai <- select(stat_analysis, contains(paste0("A", i, "-")) & contains("log10_"))
            
            p_values_tmp <- vector("numeric", length = nrow(stat_analysis))
            
            for (j in seq_along(p_values_tmp)) {
                if (test_type == "t.test") {
                    p_values_tmp[j] <- ttest_func(col_A1[j,], col_Ai[j,])
                } else if (test_type == "wilcoxon") {
                    p_values_tmp[j] <- wilcox.test(col_A1[j,], col_Ai[j,])
                }
            }
            
            all_pvalues[, i - 1] <- p_values_tmp
            
            adjust_pval_tmp <- mt.rawp2adjp(all_pvalues[, i - 1], proc = "BH", alpha = 0.05)
            qval <- data.frame(adjust_pval_tmp$adjp, adjust_pval_tmp$index)[order(adjust_pval_tmp$index), 2]
            all_adjust_pval[, i - 1] <- qval
        }
        
        all_pvalues <- as.data.frame(all_pvalues)
        all_adjust_pval <- as.data.frame(all_adjust_pval)
        
        colnames(all_pvalues) <- paste0("pvalues_A1/A", 2:sample_size)
        colnames(all_adjust_pval) <- paste0("adjust_pval_A1/A", 2:sample_size)
        
        rownames(all_pvalues) <- row.names(stat_analysis)
        rownames(all_adjust_pval) <- row.names(stat_analysis)
    
    all_pvalues_common_col <- stat_analysis %>% 
      select(pep_with_pos, accession) %>% #spectrum_title
      bind_cols(all_pvalues) %>% 
      pivot_longer(cols = starts_with("pvalues_"), values_to = "pvalues", names_to ="p_ratios") %>%
      separate(p_ratios, into = c("tmp","ratio"),sep = "_") %>%
      select(!tmp) %>%
      mutate(common_col = paste(pep_with_pos,accession,ratio,1:((sample_size-1)*nrow(stat_analysis)),sep="@")) #spectrum_title
    
    ## THE BEST WAY TO DO is this:
    merge_stat_df <- final_imputed_data %>%
      select(pep_with_pos, accession, starts_with("exp_FC")) %>%
      pivot_longer(cols = starts_with("exp_FC"), values_to = "fold_change_values", names_to ="fold_change_ratios") %>%
      separate(fold_change_ratios, into = c("tmp","tmp1","ratio"),sep = "_") %>%
      select(!c(tmp,tmp1)) %>%
      mutate(common_col = paste(pep_with_pos,accession,1:((sample_size-1)*nrow(final_imputed_data)),sep="@")) %>%
      bind_cols(all_pvalues_common_col$ratio,all_pvalues_common_col$pvalues) %>%
      rename_with(.col =6 , ~"ratio1") %>%
      rename_with(.col=7, ~ "pvalues") %>%
      separate(accession, into = c("prot_id","species"),sep = "_")
    
    
    
    p9_t_test <- ggplot(merge_stat_df ,aes(x =log2(fold_change_values), y = -log10(merge_stat_df$pvalues), color=ratio)) +
      geom_point(size = 2,aes(shape=Pool)) + #, aes(shape=merge_stat_df_final$species)
      #facet_wrap(~ratio) +
      #scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 0.8)) +
      #scale_x_continuous(limits = c(-3,3),breaks = seq(-3, 3, by = 0.8)) +
      scale_color_brewer(palette = "Set1") +
      #scale_y_continuous(breaks = seq(0, max(-log10(volcano_final1$pvalues_value)), length.out = 21)) +
      theme_bw() +
      theme(legend.text = element_text(size = 15),
            axis.title.x = element_text(size = 15),
            axis.title.y = element_text(size = 15),
            plot.title = element_text(size = 30),
            legend.title = element_text(size = 15),
            axis.text.x = element_text(size = 15),
            axis.title = element_text(size = 15),
            axis.text.y = element_text(size = 15)) +
      labs(title =  paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name), subtitle = "T-test was used")
    #expand_limits(x = 0, y = 0) +
    #geom_vline(data = actual_ratio, aes(xintercept = actual_ratio$X.1....log2.c.2..10..20..100..., size = 1, show.legend = FALSE)) + #color=c("#CC79A7","#E69F00","#56B4E9","#009E73")
    #geom_hline(data = log10_p_thresholds, aes(yintercept = log10_p_thresholds$X.log10.p_thresholds.),color=c("#CC79A7","#E69F00","#56B4E9","#009E73"), size = 1, linetype = 2, show.legend = FALSE)+ 
    #
    
    
    
    
    ### MERGING I: All used columns are merged and used to combine pvalues and ratios
          ## Before changing the shape of data 
          ## This eliminates any duplication and mismatching btw dataframes
          ## We are 100% sured that all pvalues are associated with its ratio.
      ## Shape of volcano plot looks quite weird esspecially in the A1vsA5.
    
  ## THIS IS CORRECT BUT UNNECESSARILY LONG WAY
    # test <- stat_analysis %>%
    #   select(pep_with_pos, accession, starts_with("exp_FC")) %>% #spectrum_title
    #   bind_cols(all_pvalues) %>%
    #   pivot_longer(cols = starts_with("exp_FC"), values_to = "fold_change_values", names_to ="fold_change_ratios") %>%
    #   separate(fold_change_ratios, into = c("tmp","tmp1","ratio"),sep = "_") %>%
    #   select(!c(tmp,tmp1)) %>%
    #   mutate(common_col = paste(pep_with_pos,accession,ratio,1:((sample_size-1)*nrow(stat_analysis)),sep="@")) %>%
    #   left_join(all_pvalues_common_col, by="common_col") %>%
    #   separate(accession.x, into = c("uniprot_id","species"),sep = "_")
    # 
    # 
    ### MERGING II: Previous version of combining pvalues and ratios
        ## No errors were appeared when I double pivot_longer()
        ## But there is no common column between those values so that
        ## We cannot know which pvalues correspond to which ratio (MAYBE WE DON'T NEED IT)
    ## Shape of volcano plot looks better when I used this way.
    
    ## IT's NOT CORRECT ##
    
      # merge_stat_df1 <- stat_analysis %>%
      #   select(pep_with_pos, accession, starts_with("exp_FC")) %>%
      #   bind_cols(all_pvalues) %>%
      #   pivot_longer(cols = starts_with("exp_FC"), values_to = "fold_change_values", names_to ="fold_change_ratios") %>%
      #   pivot_longer(cols = starts_with("pvalues_"), values_to = "pvalues", names_to ="p_ratios") %>%
      #   separate(accession, into = c("uniprot_id","species"),sep = "_")
      
      
  }else if(test_type=="limma"){
    library(limma)
    ## DESIGN SETUP FOR LIMMA (TODO: supply this matrix as an imput makes it faster)
    ## MANUALLY DESIGN GENERATION this is unnecessary to have that but still I'll keep it
    # for (i in 0:num_reps){
    #   design[(1:num_reps+(num_reps*i)),paste0("A",(i+2))] <- 1
    # }
    # 
    # design[is.na(design)] <- 0
    # num_reps <- 3
    # sample_size <- length(exp_design) / num_reps
    # num_comparison <- sample_size -1 
    # 
    # design <- matrix(nrow = num_comparison*num_reps, ncol = sample_size)
    # colnames(design) <- paste0("A",1:sample_size)
    # design[,"A1"] <- 1
    # 
    
    ## DESIGN SETUP FOR LIMMA (TODO: supply this matrix as an imput makes it faster)
    # Create the design matrix
    
      ################################################################################
      #### ALL COMPARISON INTO ONE DESIGN MATRIX (TODO: FIX IT DOES NOT WORK NOW!####
      
      ###  adjust.method argument in the topTable function does not directly control 
      # the target false discovery rate (FDR) threshold.
      
      ## Since it is complicated design, it needs additionaly constrast matrix! ##
      
      #contrast_mat <- makeContrasts(A1-A2,A1-A3,A1-A4,A1-A5,levels=c("A1","A2","A3","A4","A5"))
      #fit <- lmFit(stat_analysis[,3:18], design_matrix)
      #fit <- contrasts.fit(fit, contrast_mat)
      #fit <- eBayes(fit)
      #test <- topTable(fit,number = Inf, adjust.method = "BH",p.value = 1)
      
      ## FURTHER DATA MANUPLATIONS - MIGHT BE USEFUL LATER ON ##
      
      # for (j in 2:length(comparisons)){
      #   assign(paste0("results_",comparisons[j]),
      #          topTable(data_limma, coef = comparisons[j],
      #                   number = Inf, adjust.method = "BH",p.value = 1))
      #   
      #   assign(paste0("results_",comparisons[j],"_",j),cbind(get(paste0("results_",comparisons[j])), as.data.frame(rownames(get(paste0("results_",comparisons[j]))))))
      #   
      # }
      # colnames(results_A1_vs_A2_2)[1:(length(colnames(results_A1_vs_A2_2))-1)] <- paste0(colnames(results_A1_vs_A2_2), "_12")
      # colnames(results_A1_vs_A3_3)[1:(length(colnames(results_A1_vs_A3_3))-1)] <- paste0(colnames(results_A1_vs_A3_3), "_13")
      # colnames(results_A1_vs_A4_4)[1:(length(colnames(results_A1_vs_A4_4))-1)] <- paste0(colnames(results_A1_vs_A4_4), "_14")
      # colnames(results_A1_vs_A5_5)[1:(length(colnames(results_A1_vs_A5_5))-1)] <- paste0(colnames(results_A1_vs_A5_5), "_15")
      # 
      # exclude_column <- "rownames(get(paste0(\"results_\", comparisons[j])))"
      # 
      # stat_df <- results_A1_vs_A2_2 %>% 
      #   left_join(results_A1_vs_A3_3,by=exclude_column) %>%
      #   left_join(results_A1_vs_A4_4,by=exclude_column) %>%
      #   left_join(results_A1_vs_A5_5, by=exclude_column) %>%
      #   select(!contains("accession_") & contains(c("P.value_","logFC_",exclude_column))) %>%
      #   separate(exclude_column, into = c("pep_with_pos", "spectrum_title", "accession","index"), sep = "@") %>%
      #   pivot_longer(cols = starts_with("P.value"), values_to = "pvalues", names_to ="p_ratios") %>%
      #   pivot_longer(cols = starts_with("logFC_"), values_to = "fold_change_values", names_to ="fold_change_ratios")
      # 
      
      ## THIS IS FIRST TIME THAT I TRIED THE RESULTS WERE OVERESTIMATED ##
      
      #all_pvalues <- topTable(data_limma,number = Inf, adjust.method = "BH",p.value = 1)
      #design_matrix <- model.matrix(~ 0 + factor(rep(sample_names, each = 3)))
      #colnames(design_matrix) <-comparisons
      #data_limma <- limma::eBayes(limma::lmFit(stat_analysis[,3:18], design_matrix))
      ##############################################################################
      
      #### ONE COMPARISION AT A TIME WITH SEPARATE DESIGN MATRIX ####
      
      design_matrix <- model.matrix(~factor(c(rep(2,num_reps),rep(1,num_reps))))
      merge_stat_df <-NULL
      for ( i in 2:sample_size){
        # Change only the colname iteratively makes fit to every comparison
        colnames(design_matrix) <- c("Intercept", paste0("A1-A",i))
        #print(colnames(design_matrix))
        # Col selection for each comparison
        assign(paste0("df_A1vsA",i),stat_analysis %>% select(1:2 | contains("log10_E3_A1") | contains(paste0("log10_E3_A",i))))
        # First, linear model was built
        assign(paste0("fit",i) ,lmFit(get(paste0("df_A1vsA",i))[,3:8], design_matrix))
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
      
      ## MULTIPLE PIVOTING IN ONE DATAFRAME DOES NOT WORK - RATHER THAN THAT, BIND_ROWS() WAS USED ABOVE.
      
      merge_stat_df_final <- merge_stat_df %>%  #
        rename_with(.col =1 , ~"mult_col") %>%
        rename_with(.col=8, ~ "ratios") %>%
        separate(mult_col, into = c("pep_with_pos","accession","id"),sep = "@") %>% #"spectrum_title"
        separate(accession, into = c("uniprot_id","species"),sep = "_")
      
      p9_limma <- gg_volcano(data_set = merge_stat_df_final,
                   x_df = merge_stat_df_final$logFC,  
                   y_df = -log10(merge_stat_df_final$P.Value),
                   color_df = merge_stat_df_final$species,
                   facet_df = "ratios",
                   header = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name), ## Specifications of the header will be asked as an input.
                   color_lab = "Classes",
                   subtitle_txt = "Limma was used",
                   x_lab = "log2(fold_change_values)",
                   y_lab = "-log10(p_values)")   
      
      #select(accession, spectrum_title, contains("P.value") | contains("logFC")) %>%
      #pivot_longer(cols = contains("P.value"), values_to = "pvalues", names_to ="p_ratios") %>%
      #pivot_longer(cols = contains("logFC"), values_to = "fold_change_values", names_to ="fold_change_ratios")
      
      # Filter the results based on the FDR threshold
      # Add column that defined values either above or lower the threshold as Boolean 
      
      #desired_fdr_threshold <- 0.005
      #merge_stat_df$sig <- merge_stat_df$adj.P.Val < desired_fdr_threshold
      #significant_results_A1_vs_A2 <- subset(results_A1_vs_A2, adj.P.Val <= desired_fdr_threshold)
      
  }else{
      print("Statistical test could not be assessed. Check the input files!")
    }
    
    #necessary_cols <- as.data.frame(rownames(all_pvalues))
    ## DON'T USE ALL COMPARISON AT ONCE
    # stat_df <- necessary_cols %>%
    #   rename(mult_col = 1) %>%
    #   separate(mult_col, c("pep_with_pos", "spectrum_title","accession","id"),sep = "@") %>%
    #   select(-accession) %>%
    #   bind_cols(all_pvalues) %>% tibble() ## tibble() will delete rownames! 
    #                                       ## Be careful that merging correct rownames to pvalues
    # 
    #results_A1_vs_A2$sig <- results_A1_vs_A2$adj.P.Val < desired_fdr_threshold
    ## IGNORE FOR NOW!!
    # for (j in 1:nrow(results_A1_vs_A2)){
    #   #results_A1_vs_A2[j,"test12"] <- (results_A1_vs_A2[j,"rank12"]/nrow(results_A1_vs_A2))*0.05
    #   results_A1_vs_A3[j,"test_bool"] <- results_A1_vs_A3[j,"adj.P.Val"] > results_A1_vs_A3[j,"P.Value"]
    # }
    
    # stat_df <- necessary_cols %>%
    #   rename(mult_col = 1) %>%
    #   separate(mult_col, c("pep_with_pos", "spectrum_title","accession","id"),sep = "@") %>%
    #   select(-accession) %>%
    #   bind_cols(all_pvalues) %>% tibble() ## tibble() will delete rownames! 
    # ## Be careful that merging correct rownames to pvalues
    # 
    
    
    # ggplot(data=merge_stat_df1,aes(x=log2(merge_stat_df1$fold_change_values),y=-log10(merge_stat_df1$pvalues),color=species)) +
    #   geom_point() + facet_wrap(~fold_change_ratios)
    # 
    sapply(1:8,function(x) ggsave(filename = paste0("p",x,".tiff"),
                                  width = 50, height = 40, 
                                  path = paste0(file_path,"/outputs_with_new_script/"),
                                  units = "cm",
                                  get(paste0("p",x)),
                                  device = "tiff", #".svg"
    ))
    
    ggsave(filename = paste0("p9_t_test.tiff"),
           width = 50, height = 40, 
           path = paste0(file_path,"/outputs_with_new_script/"),
           units = "cm",
           p9_t_test,
           device = "tiff")
    
    ggsave(filename = paste0("p9_limma.tiff"),
           width = 50, height = 40, 
           path = paste0(file_path,"/outputs_with_new_script/"),
           units = "cm",
           p9_limma,
           device = "tiff")
    
   
  }

  
  
  
  



