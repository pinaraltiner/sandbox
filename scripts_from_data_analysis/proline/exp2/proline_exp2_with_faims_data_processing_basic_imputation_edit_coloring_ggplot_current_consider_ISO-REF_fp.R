#library(PhosR)
library(stringr)
library(dplyr)
library(data.table)
library(openxlsx)
library(tidyr)
library(ggplot2)
###############################################
source("D:/dev/Desktop_copy/PHD/data_analysis/scripts/ggplot_functions.R")
source("D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/comparision_DDA_noFAIMS_mq_proline_pd/roc_curve_generation_proline_edit.R")
# Experiment 2
 file_path <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/exp2_re_injection/"
 file_name <- "PAL _Phosphopeptides exp2 ( 5 conc 3reps) DDA_230117 with Design_2023-06-07_1003.xlsx"


#file_path <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/with_FAIMS/"
#file_name <- "PAL _Exp2_( 5 conc 3reps)_withFAIMS_DDA_25052023_noDesign - correct_2023-06-13_1514.xlsx"
selected_spcies = "_HUMAN"

# Correct syn. peptide list for Experiment 2
theo_file_path <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/"
theo_file_name <- "Synthetic peptides list_theo_conc_added_pool_id_iso_count.xlsx"   ### TODO CHANGE me
sheet_theo_name <- "ISO-ref and OTHER with FC"


# Common constant objects
sample_size <- 5
sheet_name <- "Best PSM from protein sets"
#sheet_name <- "Quantified peptide ions"







final_proline_pep_quant_analysis <- function(file_path,
                                             file_name,
                                             sheet_name,
                                             theo_file_path,
                                             theo_file_name,
                                             sheet_theo_name){
  exp2_FAIMS_DDA_proline <- read.xlsx(paste0(file_path,file_name), sheet = sheet_name)
  
  experiment_name <- c("E2_A1_R1",
                       "E2_A1_R2",
                       "E2_A1_R3",
                       "E2_A2_R1",
                       "E2_A2_R2",
                       "E2_A2_R3",
                       "E2_A3_R1",
                       "E2_A3_R2",
                       "E2_A3_R3",
                       "E2_A4_R1",
                       "E2_A4_R2",
                       "E2_A4_R3",
                       "E2_A5_R1",
                       "E2_A5_R2",
                       "E2_A5_R3")
  
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
  abundances_for_impute <- exp2_FAIMS_DDA_proline %>% 
    select(starts_with("abundance_"))
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
  # 
  # Calculate 1 percent quantile of each sample
  impute_values <- apply(abundances_for_impute, 2 , quantile , probs = 0.01 , na.rm = TRUE )
  
  # Impute missing values
  for (j in 1:length(impute_values)){
    # Number NA
    num_NA <- length(abundances_for_impute[,j][is.na(abundances_for_impute[,j])])
    
    abundances_for_impute[,j][is.na(abundances_for_impute[,j])] <- impute_values[j]
    
    # After imputation number of imputed values
    num_imp <-length(abundances_for_impute[,j][(abundances_for_impute[,j]==impute_values[j])])
    
    # This is verification of imputation is done successfully
    # Because we expect to see that number of imputed values should be the same amount as number of NA
    print(setequal(num_NA,num_imp))
    #print(num_NA)
    #print(num_imp)
  }
  
  
  selected_cols <- colnames(exp2_FAIMS_DDA_proline[c(1:4,6,28)])
  #selected_cols <- colnames(exp2_noFAIMS_DDA_proline[c(1:4,15)])
  
  
  syn_phospho_pep_proline <- exp2_FAIMS_DDA_proline %>% 
    select(-starts_with("abundance_")) %>% 
    mutate(abundances_for_impute) %>%
    filter(grepl("Phospho", ptm_protein_positions)) %>% 
    filter(grepl("_HUMAN", accession)) %>%
    select(contains(selected_cols) | starts_with("abundance"))
  
  
  pep_list_w_theo_quant <- read.xlsx(paste0(theo_file_path, theo_file_name), sheet = sheet_theo_name)
  pep_list_w_theo_quant <- pep_list_w_theo_quant[,-1]
  
  
  # Extraction of phospho positions from quant peptides object
  phospho_ptm_pos <- lapply(syn_phospho_pep_proline$ptm_protein_positions, function(each_ptm_protein_positions) {
    
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
  
  
  # Creation of common column merging peptide sequence and phospho positions -> experimental data
  common_col_exp_quant <- as.data.frame(paste(syn_phospho_pep_proline$sequence, phospho_ptm_pos_df, sep = "_"))
  colnames(common_col_exp_quant) <- "common_col_for_merging"
  
  # Bind it to the quant data
  syn_phospho_pep_proline_new <- cbind(common_col_exp_quant,syn_phospho_pep_proline)
  
  # Creation of common column merging peptide sequence and phospho positions -> theoretical data
  common_col_theo_quant <- as.data.frame(paste(pep_list_w_theo_quant$Phosphopeptide.sequence,
                                               pep_list_w_theo_quant$modified.position.in.peptide, sep = "_"))
  colnames(common_col_theo_quant) <- "common_col_for_merging"
  pep_list_w_theo_quant_new <- cbind(common_col_theo_quant,pep_list_w_theo_quant)
  
  # Sorting columns A1 to A5 and R1 to R3
  sorted_colnames <- c("psm_count","abundance") # Abundance contains both w/wo raw_abundance
  for (i in 1:length(sorted_colnames)){
    assign(paste0("tmp",i), sort(colnames(select(syn_phospho_pep_proline_new,
                                                 contains(sorted_colnames[i])))))
  }
  
  sort_col <- c(tmp1,tmp2)
  
  sort_syn_phospho_pep_proline_new <- syn_phospho_pep_proline_new %>%
    relocate(all_of(sort_col), .after = 7) %>% rename_with(~ experiment_name, all_of(sort_col)) 
  
  
  ##### REMOVE REDUNDANCY OF PHOSPHO SEQUENCES DUE TO DIFFERENT CHARGES
  syn_phospho_pep_proline_unique_max_absum <- sort_syn_phospho_pep_proline_new %>% 
    #mutate(across(all_of(experiment_name), ~ ifelse(is.na(.), 0, .))) %>%
    rowwise() %>%
    mutate(row_sum = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
    group_by(common_col_for_merging) %>%
    #summarise(max_row_sum = max(row_sum, na.rm = TRUE)) %>%
    slice(which.max(row_sum)) %>%
    ungroup()
  
  

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
  df_merge <- syn_phospho_pep_proline_unique_max_absum %>% left_join(pep_list_w_theo_quant_new,by="common_col_for_merging")
  
  
  ## ALWAYS KEEP IT COMMENTED TO PREVENT OVERWRITE
  #write.xlsx(df_merge,file = "D:/dev/Desktop_copy/PHD/wet_lab_experiments/experiment2_quantification_peptide_level/Correct_imputation_exp2proline_and_peplist_theo_before_stast_analysis.xlsx")
  
  
  
  #### PROVEMENT OF SAMPLE PREPARATION WAS SUCCESSFULLY DONE!!
  # Take mean of triplicates of each sample 
  # Ask sample_size additional parameter
  sample_size <- 5
  filtered_abundances_rowMeans <- NULL
  filtered_abundances_log10 <- NULL
  for (k in 1:sample_size){
    # If separate version of row means is not needed, it can be commented later.
    # Separate row Means can be collected in temp object to merge in "log_10_filtered_abundances_rowMeans"
    
    assign(paste0("row_means_abundances_A",k),as.data.frame(rowMeans(df_merge %>% select(contains(paste0("E2_A",k))))))

    
    assign(paste0("log10_abundances_A",k),as.data.frame(log10(df_merge  %>% select(contains(paste0("E2_A",k)))))) 
                                                            
    
    filtered_abundances_rowMeans<- bind_cols(filtered_abundances_rowMeans,get(paste0("row_means_abundances_A",k)))
    
    filtered_abundances_log10 <- bind_cols(filtered_abundances_log10, get(paste0("log10_abundances_A",k)))
    
    
  }
  
  colnames(filtered_abundances_rowMeans) <- paste0("mean_abundances_A",1:sample_size)
  
  filtered_abundances_log10 <- filtered_abundances_log10 %>%
    rename_with(~ paste0("log10_", .x), everything())# %>%
    #mutate_all(function(x) ifelse(is.na(x), 0, x)) %>%
    #mutate_all(function(x) ifelse(is.infinite(x), 0, x))


  
  # Calculate Fold Change by keeping A1 constant (mean(S1)/mean(S2), etc.)
  cols <- ncol(filtered_abundances_rowMeans)
  for(An in 2:cols){
    filtered_abundances_rowMeans[,paste0("exp_FC_A1/A",An)] <- filtered_abundances_rowMeans[,1]/filtered_abundances_rowMeans[,An]
    
  }
  # To calculate all binary combination in the data frame
  #mat <- do.call(cbind, lapply(cols, function(xj) 
  #  sapply(cols, function(xi) (filtered_abundances_rowMeans[, xj]/(filtered_abundances_rowMeans[, xj])))))
  #colnames(mat) <-  outer(names(filtered_abundances_rowMeans), names(filtered_abundances_rowMeans), paste0)
  
  final_imputed_data <- cbind(df_merge,filtered_abundances_rowMeans,filtered_abundances_log10)
 # write.table(final_imputed_data,file="D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/with_FAIMS/exp2_with_FAIMS_Proline_after_merging_psm_count.tsv",sep = "\t")
  
  #test <- final_imputed_data
  
  unexpectedly_id_peps <- final_imputed_data %>% 
    filter_at(vars(Pool), all_vars(is.na(.)))
  
  
  output_path <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/comparision_DDA_noFAIMS_mq_proline_pd/"
  #write.xlsx(unexpectedly_id_peps,file = paste0(output_path,"Corrected_imputation_proline_unexpectedly_identified_phosphopeptides_exp2_noFAIMS_DDA.xlsx"))
  
  
  stat_analysis <- final_imputed_data %>%
    select(common_col_for_merging | starts_with("log10_E2"))
  
  
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
  
  
  # iterate through each column of the data set and calculate t-tests
  triplicate_indx <- c(1,3,4,6,7,9,10,12,13,15) +1
  
  
  p_values_12 <- NULL
  p_values_13 <- NULL
  p_values_14 <- NULL
  p_values_15 <- NULL
    
  ##TODO put this into another for loop to reduce redundancy of the code  
    for(j in 1:dim(stat_analysis)[1]){
      
      #if ((sum(stat_analysis[j,2:4])/3) !=  stat_analysis[j,2:4][1] | (sum(stat_analysis[j,5:7])/3) !=  stat_analysis[j,5:7][1]){
        
      p_values_12[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,5:7])
      #p_values_13[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,8:10])
      #p_values_14[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,11:13])
      #p_values_15[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,14:16])
      #}else{
        
      #}
      #if ((sum(stat_analysis[j,2:4])/3) !=  stat_analysis[j,2:4][1] | (sum(stat_analysis[j,8:10])/3) !=  stat_analysis[j,8:10][1]){
        
        p_values_13[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,8:10])

      #}else{
        
      #} 
      #if ((sum(stat_analysis[j,2:4])/3) !=  stat_analysis[j,2:4][1] | (sum(stat_analysis[j,11:13])/3) !=  stat_analysis[j,11:13][1]){
        
        
        p_values_14[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,11:13])
       
      #}else{
        
     # }
      #if ((sum(stat_analysis[j,2:4])/3) !=  stat_analysis[j,2:4][1] | (sum(stat_analysis[j,14:16])/3) !=  stat_analysis[j,14:16][1]){
        
        p_values_15[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,14:16])
     # }else{
        
     # } 
      }
  
  
  p_values_12 <- as.data.frame(p_values_12)
  p_values_13 <- as.data.frame(p_values_13)
  p_values_14 <- as.data.frame(p_values_14)
  p_values_15 <- as.data.frame(p_values_15)

  
    for (i in 2:sample_size){
      
      assign(paste0("p_values_1",i), 
             cbind(stat_analysis$common_col_for_merging,
                   as.data.frame(get(paste0("p_values_1",i)))))
      
    }
  
  
  test_num <- nrow(stat_analysis)
  
  for (j in 1:nrow(stat_analysis)){
    
    p_values_12 <- p_values_12[order(p_values_12$p_values_12),]
    colnames(p_values_12)[1] <- "common_col_for_merging"
    p_values_12["rank12"] <- 1:test_num
    p_values_12[j,"test12"] <- (p_values_12[j,"rank12"]/test_num)*0.05
    p_values_12[j,"test_bool"] <- p_values_12[j,"test12"] > p_values_12[j,"p_values_12"]
  
    
    p_values_13 <- p_values_13[order(p_values_13$p_values_13),]
    colnames(p_values_13)[1] <- "common_col_for_merging"
    p_values_13["rank13"] <- 1:test_num
    p_values_13[j,"test13"] <- (p_values_13[j,"rank13"]/test_num)*0.05
    p_values_13[j,"test_bool"] <- p_values_13[j,"test13"] > p_values_13[j,"p_values_13"]
    
    
    p_values_14 <- p_values_14[order(p_values_14$p_values_14),]
    colnames(p_values_14)[1] <- "common_col_for_merging"
    p_values_14["rank14"] <- 1:test_num
    p_values_14[j,"test14"] <- (p_values_14[j,"rank14"]/test_num)*0.05
    p_values_14[j,"test_bool"] <- p_values_14[j,"test14"] > p_values_14[j,"p_values_14"]
    
    
    p_values_15 <- p_values_15[order(p_values_15$p_values_15),]
    colnames(p_values_15)[1] <- "common_col_for_merging"
    p_values_15["rank15"] <- 1:test_num
    p_values_15[j,"test15"] <- (p_values_15[j,"rank15"]/test_num)*0.05
    p_values_15[j,"test_bool"] <- p_values_15[j,"test15"] > p_values_15[j,"p_values_15"]
    
  }
 ### TODO: Extract where the first FALSE was generated to use for threshold of each comparison.
  #p_thresholds <- c(2.100824e-02,2.381169e-02,3.119312e-02,2.567998e-02)
  p_thresholds<- c(2.045695e-02,2.359024e-02,2.751729e-02,2.703404e-02)
  col_sel_stat_analysis <- c("common_col_for_merging",
                             #"modifications",
                             #"accession",
                             "isomericity",
                             "isomeric_count",
                             #"pool_id",
                             "Pool",
                             "A1-A2_Ratio",
                             "A1-A3_Ratio",
                             "A1-A3_Ratio",
                             "A1-A4_Ratio",
                             "A1-A5_Ratio",
                             "exp_FC_A1/A2",
                             "exp_FC_A1/A3",
                             "exp_FC_A1/A4",
                             "exp_FC_A1/A5")
    
  

  for (i in 12:15){
    assign(paste0("complete_p_values_",i), 
            final_imputed_data %>% 
              select(col_sel_stat_analysis) %>%
              left_join(get(paste0("p_values_",i)),
                        by="common_col_for_merging") %>%
              mutate(ls_correctness = ifelse(test_bool == FALSE, 0, 1)))}
    
    assign(paste0("roc_data",i),compute_roc_curve(get(paste0("complete_p_values_",i)),flag = "Others",expected = 145 ))
    #assign(paste0("roc_data_unexpected",i),compute_roc_curve(get(paste0("complete_p_values_",i)),flag = "unexpected",expected = 100))
    #assign(paste0("roc_data_isoref",i),compute_roc_curve(get(paste0("complete_p_values_",i)),flag = "ISO-REF",expected = 37 ))
    
  common_colnames <-c("common_col_for_merging", "p_values","rank","test","test_bool")
  
  
  
  colnames(p_values_12) <- common_colnames
  colnames(p_values_13) <- common_colnames
  colnames(p_values_14) <- common_colnames
  colnames(p_values_15) <- common_colnames
  complete_pvalues_all <- rbind(p_values_12,p_values_13,p_values_14,p_values_15)
  roc_data_all <- compute_roc_curve(complete_pvalues_all,flag = "Others",expected = 145 )
   
  
  ## ROC CURVE GENERATION ### 
  roc_data121 <- cbind(roc_data12, "A1 vs A2","Others", "Proline")
  roc_data131 <- cbind(roc_data13, "A1 vs A3","Others", "Proline")
  roc_data141 <- cbind(roc_data14, "A1 vs A4","Others", "Proline")
  roc_data151 <- cbind(roc_data15, "A1 vs A5","Others", "Proline")
  
  colnames(roc_data121)[c(4:6)] <- c("Comparison","Pool","Software")
  colnames(roc_data131)[c(4:6)] <- c("Comparison","Pool","Software")
  colnames(roc_data141)[c(4:6)] <- c("Comparison","Pool","Software")
  colnames(roc_data151)[c(4:6)] <- c("Comparison","Pool","Software")
  
  roc_final <- rbind(roc_data121,roc_data131,roc_data141,roc_data151)
  
  #roc_curve_data_exp2_noFAIMS_MQ <- read.delim("D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/MQ_data_analysis/exp2_wo_FAIMS/exp2_wo_FAIMSwith_MBR/roc_curve_data_exp2_noFAIMS_MQ.txt")
  write.table(roc_final,file = "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/exp2_re_injection/exp2_no_FAIMS_proline_roc_curve_data_consider_isoref_fp.tsv",sep = "\t",row.names = F)
  #roc_curve_data_pro_faims <- read.delim("D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/exp2_re_injection/exp2wo_FAIMS_DDA_Proline_volcano_plot_data.tsv",sep = "\t")
  
 
  roc_proline_mq <- rbind(roc_curve_data_exp2_noFAIMS_MQ,roc_final)
  
  
  roc_data_unexpected121 <- cbind(roc_data_unexpected12, "A1 vs A2","unexcepted")
  roc_data_unexpected131 <- cbind(roc_data_unexpected13, "A1 vs A3","unexcepted")
  roc_data_unexpected141 <- cbind(roc_data_unexpected14, "A1 vs A4","unexcepted")
  roc_data_unexpected151 <- cbind(roc_data_unexpected15, "A1 vs A5","unexcepted")
  
  colnames(roc_data_unexpected121)[c(4,5)] <- c("Comparison","Pool")
  colnames(roc_data_unexpected131)[c(4,5)] <- c("Comparison","Pool")
  colnames(roc_data_unexpected141)[c(4,5)] <- c("Comparison","Pool")
  colnames(roc_data_unexpected151)[c(4,5)] <- c("Comparison","Pool")
  
  
  
  roc_final2 <- rbind(roc_data_unexpected121,roc_data_unexpected131,
                      roc_data_unexpected141,roc_data_unexpected151)
  
  roc_curve_combine_proline_mq_score_40 <- read.delim("D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/comparision_DDA_noFAIMS_mq_proline_pd/comparison_wrong_localized_pep_MQ_Proline_with_out_FAIMS/final_figures/roc_curve_combine_proline_mq_score_40.txt")
  
  line_types <- c("solid","dashed" )

  ##TODO: make it more professional
  ggplot(roc_final,aes(x=fdp,y=tpr)) + 
    geom_line(size=1.2,aes(color=Comparison,linetype=Software)) +
    scale_color_manual(values = c("#CC79A7", "#E69F00", "#56B4E9", "#009E73"),
                       labels=c('A1 vs A2','A1 vs A3','A1 vs A4','A1 vs A5')) +
    scale_linetype_manual(values = line_types)+
    theme_bw() +
    theme(legend.text = element_text(size=15), 
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          plot.title = element_text(size=30),
          legend.title=element_text(size=15),
          axis.text.x=element_text(size=15),
          axis.title=element_text(size=15),
          axis.text.y = element_text(size = 15)) + 
    expand_limits(x = 0, y = 0) +
    scale_y_continuous(limits = c(0,100)) +
    #facet_wrap(~Ratio_col,scales = "free_x") +
    labs(title = "Experiment 2 - DDA without FAIMS - \n Comparision of Softwares")
  
  ##TODO: make it more professional
  ggplot(roc_data12,aes(x=fdp,y=tpr)) + geom_line()
  
  
  num_wrong_local_pep_Proline_wo_FAIMS <- c("ARSRTPPSAPSQSR_3&11",
  "ARSRTPPSAPSQSR_3&13",
  "ATSLPSLDTPGELR_2",
  "FSDQAGPAIPTSNSYSK_14",
  "HTDDEMTGYVATR_12",
  "HTDDEMTGYVATR_7",
  "LMTGDTYTAHAGAK_8",
  "LPLTRSHNNFVAILDLPEGEHQYK_4",
  "NGSLKPGSSHR_9",
  "RLSSTSLASGHSVR_3&5",
  "RLSSTSLASGHSVR_4&5",
  "RLSSTSLASGHSVR_4&6",
  "RLSSTSLASGHSVR_5",
  "RLSSTSLASGHSVR_5&6",
  "RLSSTSLASGHSVR_6",
  "SFGSPNRAYTHQVVTR_1&10",
  "SNSTSSMSSGLPEQDR_1",
  "STGDPQGVIR_1",
  "STVASMMHR_2",
  "TAGTSFMMTPYVVTR_14",
  "TGMGSGSAGKEGGPFK_1",
  "TVSTSSQPEENVDR_4",
  "VIEDNEYTAR_8",
  "VSPSPTTYR_6",
  "VSPSPTTYR_7",
  "YATPQVIQAPGPR_3")
  
  num_wrong_local_pept_MQ_score_40 <- c("ARSRTPPSAPSQSR_3&13",
    "DIYSTDYYR_8",
    "ESKSSPRPTAEK_2",
    "FSDQAGPAIPTSNSYSK_14",
    "IQPAGNTSPR_7",
    "LPLTRSHNNFVAILDLPEGEHQYK_4",
    "LSYYEYDFER_3",
    "RLSSFVTK_7",
    "RLSSTSLASGHSVR_3&5",
    "RLSSTSLASGHSVR_3&6",
    "RLSSTSLASGHSVR_5",
    "RSMSPFRGPK_2",
    "SFGSPNRAYTHQVVTR_1&10",
    "STVASMMHR_2",
    "TAGTSFMMTPYVVTR_14",
    "TVSTSSQPEENVDR_4&5")
  
  png_width <- 1800
  png_height <- 1500
  
  venn.diagram(
    x = list(num_wrong_local_pep_Proline_wo_FAIMS, num_wrong_local_pept_MQ_score_40),
    #category.names = c("Proline_wo_FAIMS", "MQ_wo_FAIMS"),
    #fill = c("#00AFBB", "#D16103"),  # Specify custom colors for the sets
    alpha = 0.5,  # Set transparency level
    fontfamily = "sans",  # Specify font family
    fontface = "bold",  # Specify font face
    fontcolor = "black",  # Specify font color
    fontsize = 12,  # Specify font size
    cat.fontfamily = "sans",  # Specify category label font family
    cat.fontface = "bold",  # Specify category label font face
    cat.fontsize = 12,  # Specify category label font size
    cat.cex = 1.2,  # Specify category label expansion factor
    cex = 1,  # Specify overall expansion factor for the diagram
    filename = "venn.png",  # Optional: specify a filename to save the plot as an image
    width = png_width,
    height = png_height
  )
  
  
  
  
  merging_all_pvalues <- final_imputed_data %>% 
    select(col_sel_stat_analysis) %>%
    left_join(p_values_12,by="common_col_for_merging") %>%
    left_join(p_values_13,by="common_col_for_merging") %>%
    left_join(p_values_14,by="common_col_for_merging") %>%
    left_join(p_values_15,by="common_col_for_merging")
  
  merging_all_pvalues <- merging_all_pvalues %>%
    mutate(Pool=replace_na(Pool,"unexpected"))
  
  
  
  
  
    
  # 
  # 
  # 
  # 
  # 
  # #adjusted_p_values <- as.data.frame(p.adjust(p_values_A1_vs_A5, method = "BH", n = length(p_values_A1_vs_A5)))
  # adjusted_p_values <- lapply(stat_analysis, function(x) p.adjust (x, method = "BH", n = length(x)))
  # p_values_12[,"adj_pvalue12"] <- p.adjust(p_values_12$p_values_12, method = "BH")
  # 
  # all_adjusted_p_values <- cbind(unlist(adjusted_p_values[["p_values_12"]]),
  #                                unlist(adjusted_p_values[["p_values_13"]]),
  #                                unlist(adjusted_p_values[["p_values_14"]]),
  #                                unlist(adjusted_p_values[["p_values_15"]]))
  # colnames(all_adjusted_p_values) <- c("adj_pvalues_1_2","adj_pvalues_1_3","adj_pvalues_1_4","adj_pvalues_1_5")
  # 
  # 
  
  ##### p_values_for_all_ratio changed it later to adjusted!!!!
  df_volcano <- final_imputed_data %>% 
    select(!contains("mean_abundance")) %>%
    mutate(as.data.frame(p_values_for_all_ratio)) %>%
    #mutate(across(everything(), ~replace(., is.infinite(.), NA)))
    rowwise() %>%
    filter(across(where(is.numeric), ~!is.infinite(.)))
  
  
  pvalues_df <- merging_all_pvalues %>%
    pivot_longer(cols = starts_with("p_values"), names_to = "pvalues_col", values_to = "pvalues_value") %>%
    select(common_col_for_merging,Pool,isomericity,pvalues_col,pvalues_value)
  
  # Pivot the columns containing "Ratio"
  ratio_df <- merging_all_pvalues %>%
    pivot_longer(cols = ends_with("_Ratio"), names_to = "Ratio_col", values_to = "Ratio_value") %>%
    select(Ratio_col,Ratio_value)
  
  # Pivot the columns containing "exp_FC"
  exp_FC_df <- merging_all_pvalues %>%
    pivot_longer(cols = starts_with("exp_FC"), names_to = "exp_FC_col", values_to = "exp_FC_value") %>%
    select(exp_FC_col,exp_FC_value)
  
  volcano_final <- bind_cols(pvalues_df,ratio_df,exp_FC_df)

  actual_ratio <- c(2,10,20,100)
  actual_ratio <- data.frame(Ratio_col=unique(volcano_final$Ratio_col),log2(c(2,10,20,100)))
  
  log10_p_thresholds <- data.frame(ratio_col=unique(volcano_final$Ratio_col), -log10(p_thresholds))
  
  
  #### DETERMINE THRESHOLD OF EACH P VALUE COMPARISON BASED on 
  ####    WHERE WE CAN SEE THE CHANGE FROM TRUE to FALSE
  
  ## THIS WILL BE ADDED AS AN EXTRA LINE FOR EACH COMPARISION 
  ##     WHEN I MERGED ALL VOLCANO PLOTS INTO ONE
  
  volcano_final1 <- volcano_final %>% select(-Ratio_value) %>%
  mutate(isomericity = ifelse(is.na(isomericity), "False Positive", isomericity)) %>%
    unite(Pool_new, Pool, isomericity,sep = "_",remove = FALSE) %>%
    unite('new_col_coloring',Pool_new,Ratio_col,sep = "_",remove = FALSE) %>%
    mutate(new_col_coloring = if_else(grepl("ISO-REF", new_col_coloring), "ISO-REF", new_col_coloring)) %>%
    mutate(new_col_coloring = if_else(grepl("unexpected", new_col_coloring), "unexpected", new_col_coloring))
  
  #write.table(volcano_final1,file = paste0(file_path,"/","exp2wo_FAIMS_DDA_Proline_volcano_plot_data.tsv"),sep = "\t",row.names = F)
  
  ggplot(volcano_final1, aes(x = log2(exp_FC_value), y = -log10(pvalues_value))) +
    geom_point(aes(shape = Pool_new, color = new_col_coloring), size = 2.5) +
    #geom_line(aes(color = new_col_coloring), size = 1) +  # Add color aesthetic to geom_line()
    scale_color_manual(values = c("ISO-REF" = "#000000", "unexpected" = "#999999",
                                  "Others_multi_A1-A2_Ratio" = "#CC79A7",
                                  "Others_mono_A1-A2_Ratio" = "#CC79A7",
                                  "Others_multi_A1-A3_Ratio" = "#E69F00",
                                  "Others_mono_A1-A3_Ratio" = "#E69F00",
                                  "Others_multi_A1-A4_Ratio" = "#56B4E9",
                                  "Others_mono_A1-A4_Ratio" = "#56B4E9",
                                  "Others_multi_A1-A5_Ratio" = "#009E73",
                                  "Others_mono_A1-A5_Ratio" = "#009E73"),
                       labels = c('Non-variant', 'Variant non-isomeric A1 vs A2',
                                  'Variant non-isomeric A1 vs A3',
                                  'Variant non-isomeric A1 vs A4',
                                  'Variant non-isomeric A1 vs A5',
                                  'Variant isomeric A1 vs A2',
                                  'Variant isomeric A1 vs A3',
                                  'Variant isomeric A1 vs A4',
                                  'Variant isomeric A1 vs A5',
                                  'Unexpected')) +
    scale_shape_manual(values = c(16, 15, 12, 17),
                       labels = c('Non-variant', 'Variant isomeric', 'Variant non-isomeric', 'Unexpected')) +
    scale_y_continuous(limits = c(0, 7.2), breaks = seq(0, 7.2, by = 0.8)) +
    scale_x_continuous(limits = c(-7,7)) +
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
    expand_limits(x = 0, y = 0) +
    geom_vline(data = actual_ratio, aes(xintercept = actual_ratio$log2.c.2..10..20..100..),color=c("#CC79A7","#E69F00","#56B4E9","#009E73"), size = 1, show.legend = FALSE) +
    geom_hline(data = log10_p_thresholds, aes(yintercept = log10_p_thresholds$X.log10.p_thresholds.),color=c("#CC79A7","#E69F00","#56B4E9","#009E73"), size = 1, linetype = 2, show.legend = FALSE)+ 
    labs(title = "Experiment 2 - DDA no FAIMS processed by Proline", color = "Classes", shape="Type")
  
    #coord_fixed()
  
  final_imputed_data %>% 
    pivot_longer(cols = starts_with("mean_"),
                 names_to = "Abundance_col",
                 values_to = "Abundance_value") %>%
  ggplot(aes(x = Abundance_col, y = common_col_for_merging, fill = Abundance_value)) +
    geom_tile(color = "white") +
    scale_fill_viridis() +
    labs(x = "Condition", y = "Peptide", fill = "Abundance") +
    ggtitle("Heatmap")
  
  
  
  ggplot(df_volcano, aes(x=df_volcano$`exp_FC_A1/A4`, y=-log10(pvalues_1_4),color=Pool)) + 
    geom_point(aes(Pool)) + geom_hline(yintercept = pvalue_threshold, size=1) + 
    geom_vline(xintercept = actual_ratio[3], size=1)
  #### Distribution of experimental quantitative ratios with two different pools
  
 
  
  library(ggplot2)
  
  ## This can be used as an object name
  df_for_figure <- final_imputed_data %>% 
    select(contains(col_sel_stat_analysis)) %>%
    pivot_longer(cols = col_sel_stat_analysis[3:11], names_to = "Theo_Exp", values_to = "Ratios",values_drop_na = T) %>%
    #select(!contains(c("common_col_for_merging"))) %>%
    filter(grepl("exp_FC",Theo_Exp)) %>%
    filter_at(vars(Pool), all_vars(!is.na(.))) %>%
    filter_at(vars(Ratios), all_vars(!is.infinite(.)))
  
  p1 <- gg_density(data_set = df_for_figure, 
                   x_df = df_for_figure$Ratios,
                   fill_df = df_for_figure$Theo_Exp,
                   color_df = df_for_figure$Pool,
                   header="Distribution of experimental quantitative ratios with two different pools",
                   facet_df = "Theo_Exp",
                   x_lab = "log10(Ratios)",
                   color_lab= "Pool",
                   fill_lab = "Ratios")
  
  
  #### Distribution of mean abundance of every sample with two different pools
  
  figure_with_mean_abundance <- final_imputed_data %>%
    select(starts_with("mean_") | contains(c("Pool", "common_col_for_merging"))) %>%
    pivot_longer(cols = starts_with("mean"),
                 names_to = "Theo_Exp",
                 values_to = "Mean_abundance",
                 values_drop_na = T) %>%
    filter_at(vars(Pool), all_vars(!is.na(.)))
  
  
  p2 <- gg_density(data_set = figure_with_mean_abundance, 
                   x_df = figure_with_mean_abundance$Mean_abundance,
                   fill_df = figure_with_mean_abundance$Theo_Exp,
                   color_df = figure_with_mean_abundance$Pool,
                   header="Distribution of mean abundance of every sample with two different pools",
                   facet_df = "Theo_Exp",
                   x_lab = "log10(Ratios)",
                   color_lab= "Pool",
                   fill_lab = "Ratios")
  
  # final_imputed_normalized_data %>%
  #   select(starts_with("mean_") | contains(c("Pool", "common_col_for_merging"))) %>%  
  #   filter_at(vars(Pool), all_vars(!is.na(.))) %>%
  
  #### MANUAL PLOTTING   
  #
  #   ggplot(aes(x =log10(mean_abundances_A1) , y = log10(mean_abundances_A5), color=Pool)) +
  #   geom_point()+
  #   #facet_wrap(vars(df_for_figure$common_col_for_merging))  + 
  #   geom_smooth(formula = y ~ x,method = "loess", colour = "green", fill = "green") +
  #   theme_minimal() +
  #   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
  #         axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
  #         plot.title = element_text(size=20),
  #         legend.title=element_text(size=15),
  #         axis.text=element_text(size=15),
  #         axis.title=element_text(size=15)
  #   )
  #### PLOTTING WITH FUNCTION
  #
  # gg_density_mean_abun <- function(data_set = final_imputed_normalized_data,
  #                                  x_df = final_imputed_normalized_data$mean_abundances_A1,
  #                                  y_df = final_imputed_normalized_data$mean_abundances_A5,
  #                                  color_df =final_imputed_normalized_data$Pool,
  #                                  #header,
  #                                  facet_df= "common_col_for_merging")
  
  ### BOX-PLOT: Experimental Quantity Ratio of Synthetic Peptides  
  
  p3 <- gg_boxplt_exp_ratio(data_set = df_for_figure, 
                            x_df = df_for_figure$Theo_Exp,
                            y_df = df_for_figure$Ratios,
                            fill_df = df_for_figure$Pool,
                            header="Experimental Quantity Ratio of Synthetic Peptides",
                            x_lab="Sample Names",
                            y_lab="Abundance Ratios",
                            fill_lab = "Pool")
  
  ### HALF-BOX-PLOT & HALF-SCATTER-PLOT: Experimental Quantity Ratio of Synthetic Peptides  
  library(gghalves)
  
  p4 <- gg_half_boxplt_exp_ratio(data_set = df_for_figure, 
                                 x_df = df_for_figure$Theo_Exp,
                                 y_df = df_for_figure$Ratios,
                                 fill_df = df_for_figure$Pool,
                                 header="Experimental Quantity Ratio of Synthetic Peptides",
                                 x_lab="Sample Names",
                                 y_lab="Abundance Ratios",
                                 fill_lab = "Pool")
  
  
  
  ### VIOLIN-PLOT: Experimental Quantity Ratio of Synthetic Peptides   
  
  ### TODO: fix y scaling without trimming 
  p5 <- gg_violin_exp_ratio(data_set = df_for_figure, 
                            x_df = df_for_figure$Theo_Exp,
                            y_df = df_for_figure$Ratios,
                            fill_df = df_for_figure$Pool,
                            header="Experimental Quantity Ratio of Synthetic Peptides",
                            x_lab="Sample Names",
                            y_lab="Abundance Ratios",
                            fill_lab = "Pool",
                            trim=TRUE)
  
  ### SAVE ALL PLOT AUTOMATICALLY (without giving a custom file name)
  sapply(1:5,function(x) ggsave(filename = paste0("p",x,".tiff"),
                                width = 50, height = 40, 
                                path = file_path,
                                units = "cm",
                                get(paste0("p",x)),
                                device = "tiff", #".svg"
  ))
  
  ### SAVE PLOTS SEPARATELY (with custom file name)
  # ggsave(filename = "Experimental_Quantity_Ratio_of_Synthetic_Peptides",
  #        path = "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/",
  #        plot = p5,
  #        device = "tiff")
  
  
  # library(limma)
  # sample_size <- 5
  # # Do t-test
  # # The code below does t-test for each row. Because of that, multiple test correction (like BH, Bonferoni)
  # 
  # ## 1 ## TIDY DATA FOR T-TEST
  # df_for_ttest<- final_imputed_normalized_data %>%
  #   select(starts_with("mean_") | contains(c("Pool", "common_col_for_merging"))) %>%  
  #   filter_at(vars(Pool), all_vars(!is.na(.)))
  # 
  # ## 2 ## FUNCTION FOR T-TEST ALL CONC. ACROSS A1
  # ttest_func <- function(A1,rest){
  #   t.test(A1,rest,alternative = "two.sided", var.equal = TRUE)}
  # 
  # p_values_for_all_ratio <- data.frame(1:dim(df_for_ttest)[1])
  # 
  # ## 3 ## AUTOMIZED T-TEST ALL CONC. ACROSS A1
  # for(k in 2:sample_size){
  #   
  #   assign(paste0("t_test_res_A1_to_A",k) , 
  #          lapply(df_for_ttest$mean_abundances_A1, ttest_func,
  #                 rest=df_for_ttest[,k]))
  #   
  #   assign(paste0("p_values_for_ratio",k),
  #          lapply(get(paste0("t_test_res_A1_to_A",k)), function (x) x[c('p.value')]))
  #   
  #   p_values_for_all_ratio <- cbind(p_values_for_all_ratio,
  #                                   as.data.frame(unlist(get(paste0("p_values_for_ratio",k)))))
  # }
  # 
  # 
  # 
  # #adjusted_p_values <- as.data.frame(p.adjust(p_values_A1_vs_A5, method = "BH", n = length(p_values_A1_vs_A5)))
  # adjusted_p_values <- lapply(p_values_for_all_ratio[,-1], function(x) p.adjust (x, method = "BH", n = length(x)))
  # 
  # 
  # 
  # 
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
  # df_for_figure_exp <- correct_identifed_peps[,c(11,74:77)]
  # df_for_figure <- correct_identifed_peps[,c(11:15,74:77)]
  # melt_df_for_figure_exp <- melt(df_for_figure_exp)
  # 
  # df_for_figure_theo <- correct_identifed_peps[,c(11:15)]
  # melt_df_for_figure_theo <- melt(df_for_figure_theo)
  # 
  # p1 <- ggplot(melt_df_for_figure_theo,aes(x =melt_df_for_figure_theo$variable , y =log2(melt_df_for_figure_theo$value),fill = melt_df_for_figure_theo$Pool)  ) +
  #   geom_boxplot() +
  #   theme_light() + #scale_y_continuous(limits = c(-60, 60),breaks = seq(-60, 60, by = 20)) +
  #   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
  #         axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
  #         plot.title = element_text(size=20),
  #         legend.title=element_text(size=15),
  #         axis.text=element_text(size=15),
  #         axis.title=element_text(size=15)
  #   ) +   stat_boxplot(geom = "errorbar") + ggtitle("Theoretical Quantity Ratio of Synthetic Peptides") +
  #   scale_x_discrete(labels=c("A1/A2","A1/A2","A1/A4","A1/A5")) +
  #   labs(x="Sample Names",y="Abundance Ratios", color="Pool Names")+
  #   scale_fill_brewer(palette="Set1")
  # 
  # p2 <-ggplot(melt_df_for_figure_exp,aes(x =melt_df_for_figure_exp$variable , y =log2(melt_df_for_figure_exp$value),fill = melt_df_for_figure_exp$Pool)  ) +
  #   geom_boxplot() +
  #   theme_light() +
  #   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
  #         axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
  #         plot.title = element_text(size=20),
  #         legend.title=element_text(size=15),
  #         axis.text=element_text(size=15),
  #         axis.title=element_text(size=15)
  #   ) +   stat_boxplot(geom = "errorbar") + ggtitle("Experimental  Quantity Ratio of Synthetic Peptides") +
  #   scale_x_discrete(labels=c("A1/A2","A1/A2","A1/A4","A1/A5")) +
  #   labs(x="Sample Names",y="Abundance Ratios", color="Pool Names") +
  #   scale_fill_brewer(palette="Set1")
  # 
  # p3 <-ggplot(melt_df_for_figure_exp,aes(x =melt_df_for_figure_exp$variable , y =melt_df_for_figure_exp$value,color = melt_df_for_figure_exp$Pool)  ) +
  #   geom_point() +
  #   theme_light() +
  #   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
  #         axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
  #         plot.title = element_text(size=20),
  #         legend.title=element_text(size=15),
  #         axis.text=element_text(size=15),
  #         axis.title=element_text(size=15)
  #   ) +   ggtitle("Experimental Quantity Ratio of Synthetic Peptides") +
  #   scale_x_discrete(labels=c("A1/A2","A1/A3","A1/A4","A1/A5")) +
  #   labs(x="Sample Names",y="Abundance Ratios", color="Pool Names") +
  #   scale_color_brewer(palette="Set1")
  # 
  # p4 <-ggplot(melt_df_for_figure_theo,aes(x =melt_df_for_figure_theo$variable , y =melt_df_for_figure_theo$value,color = melt_df_for_figure_theo$Pool)  ) +
  #   geom_point() +
  #   theme_light() +
  #   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
  #         axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
  #         plot.title = element_text(size=20),
  #         legend.title=element_text(size=15),
  #         axis.text=element_text(size=15),
  #         axis.title=element_text(size=15),
  #   ) +   ggtitle("Theoretical Quantity Ratio of Synthetic Peptides") +
  #   scale_x_discrete(labels=c("A1/A2","A1/A3","A1/A4","A1/A5")) +
  #   labs(x="Sample Names",y="Abundance Ratios", color="Pool Names") +
  #   scale_color_brewer(palette="Set1")
  # 
  # 
  # library(ggpubr)
  # ggarrange(p1, p2, common.legend = TRUE, legend="right")
  # ggarrange(p4, p3, common.legend = TRUE, legend="right")
  # 
  # 
  # iso_exp <- df_for_figure_exp %>% filter(grepl("ISO-REF",Pool))
  # melt_iso_theo <- melt(iso_exp)
  # ggplot(melt_iso_theo,aes(x =melt_iso_theo$variable , y =log2(melt_iso_theo$value),color = melt_iso_theo$Pool)  ) +
  #   geom_violin()  + #geom_point() +
  #   theme_light() +
  #   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
  #         axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
  #         plot.title = element_text(size=20),
  #         legend.title=element_text(size=15),
  #         axis.text=element_text(size=15),
  #         axis.title=element_text(size=15)
  #   ) +   ggtitle("Experimental Quantity Ratio of Synthetic Peptides") +
  #   scale_x_discrete(labels=c("A1/A2","A1/A3","A1/A4","A1/A5")) +
  #   labs(x="Sample Names",y="Abundance Ratios", color="Pool Names") +
  #   scale_color_brewer(palette="Set1")
  # 
  # 
  # other_exp <- df_for_figure_exp %>% filter(grepl("Others",Pool))
  # melt_other_exp <- melt(other_exp)
  # ggplot(melt_other_exp,aes(x =melt_other_exp$variable , y =log2(melt_other_exp$value),color = melt_other_exp$Pool)  ) +
  #   geom_violin() + #geom_point() +
  #   theme_light() +
  #   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
  #         axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
  #         plot.title = element_text(size=20),
  #         legend.title=element_text(size=15),
  #         axis.text=element_text(size=15),
  #         axis.title=element_text(size=15)
  #   ) +   ggtitle("Experimental Quantity Ratio of Synthetic Peptides") +
  #   scale_x_discrete(labels=c("A1/A2","A1/A3","A1/A4","A1/A5")) +
  #   labs(x="Sample Names",y="Abundance Ratios", color="Pool Names") +
  #   scale_color_brewer(palette="Set1")
  
  
  
  
}

