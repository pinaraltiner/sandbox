#library(PhosR)
library(stringr)
library(dplyr)
library(data.table)
library(openxlsx)
library(ggplot2)
library(tidyr)
###############################################

source("D:/dev/Desktop_copy/PHD/data_analysis/scripts/ggplot_functions.R")

file_path_exp3 <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_3/Proline_data_analysis/quant_pep_06022023/"
file_name_exp3 <- "PAL _Tcell_Exp3_( 5 conc 3reps)_NoFAIMS_DDA_with_cont_230206.xlsx"
selected_spcies = "_MOUSE"
sheet_name <- "Quantified peptide ions"


final_proline_pep_quant_analysis_bio <- function(file_path,
                                             file_name,
                                             sheet_name,
                                             theo_file_path,
                                             theo_file_name,
                                             sheet_theo_name,
                                             selected_spcies,
                                             sample_size){
  quant_peptides <- read.xlsx(paste0(file_path_exp3,file_name_exp3), sheet = sheet_name)
  
  quant_peptides_with_all <- quant_peptides %>% 
    filter(grepl("_MOUSE",accession)) %>% 
    filter(grepl("Phospho",modifications)) 
  
  quant_peptides_ECOLI <- quant_peptides %>% 
    filter(grepl("ECOLI",accession)) %>% 
    select(sequence,ptm_protein_positions, starts_with("abundance"))
    
  
  phospho_ptm_pos <- lapply(quant_peptides_with_all$ptm_protein_positions, function(each_ptm_protein_positions) {
    
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
  
  quant_peptides_with_all$pep_with_pos <- paste(quant_peptides_with_all$sequence, phospho_ptm_pos_df, sep = "_")

  # Nothing is changed
  filtered_abundances<-quant_peptides_with_all[rowSums(!is.na(select(quant_peptides_with_all,starts_with("abundance_"))))>0,]
  
  df_id_pep <- filtered_abundances %>% 
    select(sequence,ptm_protein_positions, pep_with_pos, starts_with("abundance")) %>%
    tibble() %>%
    rename(A1_R1= 4, # Using column index to rename the colnames
           A1_R2= 5,
           A1_R3= 6,
           A2_R1= 7,
           A2_R2= 8,
           A2_R3= 9,
           A3_R1= 10,
           A3_R2= 11,
           A3_R3= 12,
           A4_R1= 13,
           A4_R2= 14,
           A4_R3= 15,
           A5_R1= 16,
           A5_R2= 17,
           A5_R3= 18) %>%
    pivot_longer(cols = starts_with("A"), 
                 values_to = "intensity",
                 names_to = "sample_ids",
                 values_drop_na = T) %>%
    separate(sample_ids, into = c("Sample_id", "Rep_id"), sep = "_")# %>% 
    #group_by(Sample_id) %>%
    #count()
  
  barplt_df <- df_id_pep %>%
    mutate(sample_rep_id_seq = paste(pep_with_pos, Sample_id,Rep_id, sep = "_"))%>%
    filter(duplicated(sample_rep_id_seq)==FALSE)
  
  p1 <- gg_barplt_id_pep_count(data_set = barplt_df,
                         x_df = barplt_df$Sample_id,
                         fill_df = barplt_df$Rep_id,
                         ymax = 20000,
                         header = "Total number of quantified phospho-site across each sample",
                         caption_lab = "NA values are removed.",
                         x_lab = "Sample id",
                         fill_lab =  "Sample id",
                         y_lab = "Number of identified peptides")
  
  barplt_df_ecoli <- quant_peptides_ECOLI %>% 
    mutate(sample_rep_id_seq = paste(sequence, Sample_id,Rep_id, sep = "_")) %>%
    filter(duplicated(sample_rep_id_seq)==FALSE)
  
  p2 <- gg_barplt_id_pep_count(data_set = barplt_df_ecoli,
                         x_df = barplt_df_ecoli$Sample_id,
                         fill_df = barplt_df_ecoli$Rep_id,
                         ymax = 20000,
                         header = "Total number of quantified Ecoli across each sample",
                         caption_lab = "NA values are removed.",
                         x_lab = "Sample id",
                         fill_lab =  "Sample id",
                         y_lab = "Number of identified peptides")

  abundances_only_for_impute <- filtered_abundances %>% select(starts_with("abundance_"))
  
  
  ## DENSITY PLOT OF BEFORE IMPUTATION 
  
  sample_size <- 5
  abundances_rowMeans <- NULL
  abundances_ecoli_rowMeans <- NULL
  for (k in 1:sample_size){
    # If separate version of row means is not needed, it can be commented later.
    # Separate row Means can be collected in temp object to merge in "log_10_filtered_abundances_rowMeans"
    assign(paste0("abundances_A",k),as.data.frame(rowMeans(abundances_only_for_impute  %>%
                                                             select(contains(paste0("A",k))) %>%
                                                             select(starts_with("abundance_")))))
    abundances_rowMeans<- bind_cols(abundances_rowMeans,get(paste0("abundances_A",k)))
    
    
    assign(paste0("ecoli_abundances_A",k),as.data.frame(rowMeans(quant_peptides_ECOLI  %>%
                                                             select(contains(paste0("A",k))) %>%
                                                             select(starts_with("abundance_")))))  
    
    abundances_ecoli_rowMeans<- bind_cols(abundances_ecoli_rowMeans,get(paste0("ecoli_abundances_A",k)))
    
  }
  
  
  quant_peptides_ECOLI_density_plot <- quant_peptides_ECOLI %>%
    select(!starts_with("abundance")) %>%
    bind_cols(abundances_ecoli_rowMeans) %>%
    tibble() %>%
    rename(mean_abundances_A1= 3, # Using column index to rename the colnames
           mean_abundances_A2= 4,
           mean_abundances_A3= 5,
           mean_abundances_A4= 6,
           mean_abundances_A5= 7) %>%
    pivot_longer(cols = starts_with("mean"), 
                 values_to = "intensity",
                 names_to = "sample_ids",
                 values_drop_na = T) %>%
    mutate(sample_id_seq = paste(sequence, sample_ids, sep = "_")) %>%
    filter(duplicated(sample_id_seq)==FALSE) %>%
    mutate(species="ECOLI") 
  
  ecoli_density_plot<- quant_peptides_ECOLI_density_plot %>%
    select(contains(c("sample_ids","intensity","species"))) 
  
  
  
  
  colnames(abundances_rowMeans) <- paste0("mean_abundances_A",1:sample_size)
  
  #### MEAN ABUNDANCE RATIO WITH  DENSITY PLOT ####
  figure_with_mean_abundance <- abundances_rowMeans %>% 
    tibble() %>% mutate(pep_with_pos = filtered_abundances$pep_with_pos) %>%
    mutate(species="MOUSE") %>%
    pivot_longer(cols = starts_with("mean"), 
                 names_to = "sample_ids",
                 values_to = "intensity",
                 values_drop_na = T)
  
  density_df <-figure_with_mean_abundance %>%
    select(c(sample_ids,intensity,species)) %>%
    bind_rows(ecoli_density_plot)
  
  
  p3 <- gg_density(data_set = density_df, 
             x_df = density_df$intensity,
             fill_df = density_df$species,
             color_df = NULL,
             header="Distribution of mean abundance of every sample before imputation",
             facet_df = "sample_ids",
             x_lab = "log10(intensities)",
             color_lab= "",
             fill_lab = "species")

  # Calculate 1 percent quantile of each sample
  impute_values <- apply(abundances_only_for_impute, 2 , quantile , probs = 0.01 , na.rm = TRUE )
  
  # Impute missing values
  for (j in 1:length(impute_values)){
    # Number NA
    num_NA <- length(abundances_only_for_impute[,j][is.na(abundances_only_for_impute[,j])])
    
    abundances_only_for_impute[,j][is.na(abundances_only_for_impute[,j])] <- impute_values[j]
    
    # After imputation number of imputed values
    num_imp <-length(abundances_only_for_impute[,j][(abundances_only_for_impute[,j]==impute_values[j])])
    
    # This is verification of imputation is done successfully
    # Because we expect to see that number of imputed values should be the same amount as number of NA
    # However there is also possibility that data has already had the same value just by chance.
    
    print(setequal(num_NA,num_imp))
    #print(num_NA)
    #print(num_imp)
    #print(impute_values[j])
  }
  
  
  # Take log10 
  #log_10_abundances_only_for_impute <- log10(abundances_only_for_impute)
  #colnames(log_10_abundances_only_for_impute) <- paste0(colnames(abundances_only_for_impute),"_log10")
  
  ## No need is right now. I will cont with "log_10_abundances_only_for_impute" for rowMeans and FC
  #log_10_filtered_abundances <- cbind(filtered_abundances,log_10_abundances_only_for_impute)
  
  # Take mean of triplicates of each sample 
  # Ask sample_size additional parameter
  
  filtered_abundances_rowMeans <- NULL
  for (k in 1:sample_size){
    # If separate version of row means is not needed, it can be commented later.
    # Separate row Means can be collected in temp object to merge in "log_10_filtered_abundances_rowMeans"
    assign(paste0("aft_imp_abundances_A",k),as.data.frame(rowMeans(abundances_only_for_impute  %>% select(contains(paste0("A",k))) %>%
                                                             select(starts_with("abundance_")))))
    filtered_abundances_rowMeans<- bind_cols(filtered_abundances_rowMeans,get(paste0("aft_imp_abundances_A",k)))
    
  }
  
  colnames(filtered_abundances_rowMeans) <- paste0("mean_abundances_aft_imp_A",1:sample_size)
  
  # Calculate Fold Change by keeping A1 constant (mean(S1)/mean(S2), etc.)
  cols <- ncol(filtered_abundances_rowMeans)
  for(An in 2:cols){
    filtered_abundances_rowMeans[,paste0("exp_FC_A1/A",An)] <- filtered_abundances_rowMeans[,1]/filtered_abundances_rowMeans[,An]
    
  }
  # To calculate all binary combination in the data frame
  #mat <- do.call(cbind, lapply(cols, function(xj) 
  #  sapply(cols, function(xi) (filtered_abundances_rowMeans[, xj]/(filtered_abundances_rowMeans[, xj])))))
  #colnames(mat) <-  outer(names(filtered_abundances_rowMeans), names(filtered_abundances_rowMeans), paste0)
  
  final_imputed_normalized_data <- cbind(filtered_abundances,filtered_abundances_rowMeans)
  
  #write.table(final_imputed_normalized_data, file = "final_imputed_normalized_data_PAL _T_cell_Exp3_( 5 conc 3reps)_NoFAIMS_DDA_with_cont_230206_2023-02-07_0947.txt",sep = "\t",row.names = F)
  
  df_mean_ab_after_impt <- final_imputed_normalized_data %>% 
    select(starts_with("mean_")| contains("common_col_for_merging")) %>%
    tibble() %>% 
    pivot_longer(cols = starts_with("mean"),
                 names_to = "Mean_abundance",
                 values_to = "values")
    
  df_FC_ratio_after_impt <- final_imputed_normalized_data %>% 
    select(starts_with("exp_")| contains("common_col_for_merging")) %>%
    tibble() %>% 
    pivot_longer(cols = starts_with("exp_"),
                 names_to = "exp_FC",
                 values_to = "values")
  
  p4 <- gg_density(data_set = df_mean_ab_after_impt, 
                   x_df = df_mean_ab_after_impt$values,
                   fill_df = df_mean_ab_after_impt$Mean_abundance,
                   color_df = NULL,
                   header="Distribution of mean abundance of every sample after imputation",
                   facet_df = "Mean_abundance",
                   x_lab = "log10(values)",
                   color_lab= "",
                   fill_lab = "Sample Names")
  
  p5 <- gg_raincloud(data_set = df_mean_ab_after_impt,
                     x_df = df_mean_ab_after_impt$Mean_abundance,
                     y_df = df_mean_ab_after_impt$values,
                     fill_df = df_mean_ab_after_impt$Mean_abundance,
                     header = "Distribution of mean abundance of every sample after imputation",
                     x_lab = "Sample Names",
                     y_lab = " Density of log10(Mean Abundance)",
                     fill_lab = "Sample Names",
                     caption_lab = "")
  
  
  
  # p4 <- gg_density(data_set = df_FC_ratio_after_impt, 
  #                  x_df = df_FC_ratio_after_impt$values,
  #                  fill_df = df_FC_ratio_after_impt$exp_FC,
  #                  color_df = NULL,
  #                  header="Distribution of Fold change Ratio of every sample after imputation",
  #                  facet_df = "exp_FC",
  #                  x_lab = "log10(values)",
  #                  color_lab= "",
  #                  fill_lab = "Sample Names")
  # 

  ### BOX-PLOT: Experimental Quantity Ratio of Phospho Peptides  
  
  p6 <- gg_boxplt_exp_ratio(data_set = df_FC_ratio_after_impt, 
                            x_df = df_FC_ratio_after_impt$exp_FC,
                            y_df = df_FC_ratio_after_impt$values,
                            fill_df = df_FC_ratio_after_impt$exp_FC,
                            header="Experimental Quantity Ratio of T-cell Phospho Peptides",
                            x_lab="Sample Names",
                            y_lab="Abundance Ratios",
                            fill_lab = "Sample Names")
  
  
  ### HALF-BOX-PLOT & HALF-SCATTER-PLOT: Experimental Quantity Ratio of Synthetic Peptides  
  library(gghalves)
  
  p7 <- gg_half_boxplt_exp_ratio(data_set = df_FC_ratio_after_impt, 
                                 x_df = df_FC_ratio_after_impt$exp_FC,
                                 y_df = df_FC_ratio_after_impt$values,
                                 fill_df = df_FC_ratio_after_impt$exp_FC,
                                 header="Experimental Quantity Ratio of T-cell Phospho Peptides",
                                 x_lab="Sample Names",
                                 y_lab="Abundance Ratios",
                                 fill_lab = "Sample Names")
  
  ### VIOLIN-PLOT: Experimental Quantity Ratio of Synthetic Peptides   
  
  ### TODO: fix y scaling without trimming 
  p8 <- gg_violin_exp_ratio(data_set = df_FC_ratio_after_impt, 
                            x_df = df_FC_ratio_after_impt$exp_FC,
                            y_df = df_FC_ratio_after_impt$values,
                            fill_df = df_FC_ratio_after_impt$exp_FC,
                            header="Experimental Quantity Ratio of T-cell Phospho Peptides",
                            x_lab="Sample Names",
                            y_lab="Abundance Ratios",
                            fill_lab = "Sample Names",
                            trim=TRUE)
  
  sapply(1:8,function(x) ggsave(filename = paste0("p",x,".tiff"),
                                width = 50, height = 40, 
                                path = file_path,
                                units = "cm",
                                get(paste0("p",x)),
                                device = "tiff", #".svg"
  ))
  
  
  
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
}

