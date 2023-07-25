#library(PhosR)
library(stringr)
library(dplyr)
library(data.table)
library(openxlsx)
library(ggplot2)
library(tidyr)
###############################################

source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/ggplot/ggplot_functions.R")
# 
# file_path_exp3 <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_3/Proline_data_analysis/quant_pep_06022023/"
# file_name_exp3 <- "PAL _Tcell_Exp3_( 5 conc 3reps)_NoFAIMS_DDA_with_cont_230206.xlsx"
# selected_spcies = "_MOUSE"
# sheet_name <- "Quantified peptide ions"
# 


file_path_exp3_tims <- "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_3/Proline_timsdata/"
file_name_exp3_tims <- "PAL TimsTOF_data_conversion_all_exp_06062023 E3 corrected_2023-06-28_1522.xlsx"
# Common constant objects
sample_size <- 5
#sheet_name <- "Quantified peptide ions"
sheet_name <- "Best PSM from protein sets"
file_path<- file_path_exp3_tims 
file_name <- file_name_exp3_tims
selected_spcies <-  "_MOUSE"
test_type <- c("t.test","wilcoxon","limma") ## one selection at a time

exp_design <- c("E3_A1_R1","E3_A1_R2",
                "E3_A1_R3",
                "E3_A2_R1",
                "E3_A2_R2",
                "E3_A2_R3",
                "E3_A3_R1",
                "E3_A3_R2",
                "E3_A3_R3",
                "E3_A4_R1",
                "E3_A4_R2",
                "E3_A4_R3",
                "E3_A5_R1",
                "E3_A5_R2",
                "E3_A5_R3")
final_proline_pep_quant_analysis_bio <- function(file_path,
                                             file_name,
                                             sheet_name,
                                             theo_file_path,
                                             theo_file_name,
                                             sheet_theo_name,
                                             selected_spcies,
                                             sample_size){
  quant_peptides <- read.xlsx(paste0(file_path,file_name), sheet = sheet_name)
  
  abundances_for_impute <- quant_peptides %>% 
    select(starts_with("abundance_")) %>% ## spectrum_title remove it because it was not make it as rownames (has duplicates)
    rename_with(~exp_design,matches("abundance"))
  
  quant_peptides_cor_abun <- quant_peptides %>% 
    rename_with(~exp_design,matches("^abundance"))
  
  quant_peptides_with_all <- quant_peptides_cor_abun %>% 
    filter(grepl("_MOUSE",accession)) %>% 
    filter(grepl("Phospho",modifications)) 
  
  quant_peptides_ECOLI <- quant_peptides_cor_abun %>% 
    filter(grepl("ECOLI",accession)) %>% 
    select(sequence,modifications, starts_with(exp_design))
    
  
  phospho_ptm_pos <- lapply(quant_peptides_with_all$modifications, function(each_ptm_protein_positions) {
    
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
  filtered_abundances<-quant_peptides_with_all[rowSums(!is.na(select(quant_peptides_with_all,starts_with(exp_design))))>0,]
  filtered_abundances_ecoli <-quant_peptides_ECOLI[rowSums(!is.na(select(quant_peptides_ECOLI,starts_with(exp_design))))>0,]
  
  df_id_pep <- filtered_abundances %>% 
    select(sequence,modifications, pep_with_pos, starts_with(exp_design)) %>%
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
    separate(sample_ids, into = c("Exp_id","Sample_id", "Rep_id"), sep = "_")# %>% 
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
                         y_lab = "Number of identified peptides",
                         subtitle_txt = "")
  
  barplt_df_ecoli <- filtered_abundances_ecoli %>% 
    pivot_longer(cols = starts_with("E3"), 
                 values_to = "intensity",
                 names_to = "sample_ids",
                 values_drop_na = T) %>%
    separate(sample_ids, into = c("Exp_id","Sample_id", "Rep_id"), sep = "_") %>%
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
                         y_lab = "Number of identified peptides",
                         subtitle_txt = "")

  
  ## DENSITY PLOT OF BEFORE IMPUTATION 
  
  sample_size <- 5
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
  

  quant_peptides_ECOLI_density_plot <- filtered_abundances_ecoli %>%
    select(!starts_with("E")) %>%
    bind_cols(abundances_ecoli_rowMeans) %>%
    tibble() %>%
    rename_with(~ paste0("mean_abun",1:5), matches("^row")) %>%
    pivot_longer(cols = starts_with("mean"), 
                 values_to = "intensity",
                 names_to = "sample_ids",
                 values_drop_na = T) %>%
    mutate(sample_id_seq = paste(sequence, sample_ids, sep = "_")) %>%
    filter(duplicated(sample_id_seq)==FALSE) %>%
    mutate(species="ECOLI") 
  
  ecoli_density_plot<- quant_peptides_ECOLI_density_plot %>%
    select(contains(c("sample_ids","intensity","species"))) 
  
  
  colnames(abundances_rowMeans) <- paste0("mean_abun",1:sample_size)
  
  #### MEAN ABUNDANCE RATIO WITH  DENSITY PLOT ####
      ### BEFORE IMPUTATION ###
  quant_phospho_density_plot <- filtered_abundances %>%
    select(!starts_with("E")) %>%
    bind_cols(abundances_rowMeans) %>% 
    #rename_with(~ paste0("mean_abun",1:5), matches("^row")) %>%
    tibble() %>% #mutate(pep_with_pos = sequence) %>% ###  At this stage, no need for phospho-position
    mutate(species="MOUSE") %>%                         # We only count number of sequence here.   
    pivot_longer(cols = starts_with("mean"), 
                 names_to = "sample_ids",
                 values_to = "intensity",
                 values_drop_na = T)
  
  density_df <-quant_phospho_density_plot %>%
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
             fill_lab = "species",
             subtitle_txt = "")

  # Calculate 1 percent quantile of each sample
  
  impute_values <- apply(abundances_for_impute, 2 , quantile , probs = 0.01 , na.rm = TRUE )
  abundances_for_impute_mouse <- filtered_abundances %>% select(exp_design)
  # Impute missing values
  for (j in 1:length(impute_values)){
    # Number NA
    num_NA <- length(abundances_for_impute_mouse[,j][is.na(abundances_for_impute_mouse[,j])])
    
    abundances_for_impute_mouse[,j][is.na(abundances_for_impute_mouse[,j])] <- impute_values[j]
    
    # After imputation number of imputed values
    num_imp <-length(abundances_for_impute_mouse[,j][(abundances_for_impute_mouse[,j]==impute_values[j])])
    
    # This is verification of imputation is done successfully
    # Because we expect to see that number of imputed values should be the same amount as number of NA
    print(setequal(num_NA,num_imp))
    #print(num_NA)
    #print(num_imp)
  }
  
  
  
  # Take log10 
  #log_10_abundances_only_for_impute <- log10(abundances_only_for_impute)
  #colnames(log_10_abundances_only_for_impute) <- paste0(colnames(abundances_only_for_impute),"_log10")
  
  ## No need is right now. I will cont with "log_10_abundances_only_for_impute" for rowMeans and FC
  #log_10_filtered_abundances <- cbind(filtered_abundances,log_10_abundances_only_for_impute)
  
  # Take mean of triplicates of each sample 
  # Ask sample_size additional parameter
  
  filtered_abundances_rowMeans <- NULL
  filtered_abundances_log10 <- NULL
  
  for (k in 1:sample_size){
    # If separate version of row means is not needed, it can be commented later.
    # Separate row Means can be collected in temp object to merge in "log_10_filtered_abundances_rowMeans"
    assign(paste0("aft_imp_abundances_A",k),as.data.frame(rowMeans(abundances_for_impute_mouse  %>% select(contains(paste0("A",k))))))
    filtered_abundances_rowMeans<- bind_cols(filtered_abundances_rowMeans,get(paste0("aft_imp_abundances_A",k)))
    
    assign(paste0("log10_abundances_A",k),as.data.frame(log10(abundances_for_impute_mouse  %>% select(contains(paste0("A",k)))))) 
    filtered_abundances_log10 <- bind_cols(filtered_abundances_log10, get(paste0("log10_abundances_A",k)))
    
  }
  
  colnames(filtered_abundances_rowMeans) <- paste0("mean_abundances_aft_imp_A",1:sample_size)
  
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
  
  final_imputed_data <- cbind(filtered_abundances,filtered_abundances_rowMeans,filtered_abundances_log10)
  
  #write.table(final_imputed_data, file = "final_imputed_normalized_data_PAL _T_cell_Exp3_( 5 conc 3reps)_NoFAIMS_DDA_with_cont_230206_2023-02-07_0947.txt",sep = "\t",row.names = F)
  
  df_mean_ab_after_impt <- final_imputed_data %>% 
    select(starts_with("mean_")| contains("common_col_for_merging")) %>%
    tibble() %>% 
    pivot_longer(cols = starts_with("mean"),
                 names_to = "Mean_abundance",
                 values_to = "values")
    
  df_FC_ratio_after_impt <- final_imputed_data %>% 
    select(starts_with("exp_")| contains("pep_with_pos")) %>%
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
                   fill_lab = "Sample Names",
                   subtitle_txt = "")
  
  p5 <- gg_raincloud(data_set = df_mean_ab_after_impt,
                     x_df = df_mean_ab_after_impt$Mean_abundance,
                     y_df = df_mean_ab_after_impt$values,
                     fill_df = df_mean_ab_after_impt$Mean_abundance,
                     header = "Distribution of mean abundance of every sample after imputation",
                     x_lab = "Sample Names",
                     y_lab = " Density of log10(Mean Abundance)",
                     fill_lab = "Sample Names",
                     caption_lab = "",
                     subtitle_txt = "")
  
  
  
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
                            fill_lab = "Sample Names",
                            subtitle_txt = "")
  
  
  ### HALF-BOX-PLOT & HALF-SCATTER-PLOT: Experimental Quantity Ratio of Synthetic Peptides  
  library(gghalves)
  
  p7 <- gg_half_boxplt_exp_ratio(data_set = df_FC_ratio_after_impt, 
                                 x_df = df_FC_ratio_after_impt$exp_FC,
                                 y_df = df_FC_ratio_after_impt$values,
                                 fill_df = df_FC_ratio_after_impt$exp_FC,
                                 header="Experimental Quantity Ratio of T-cell Phospho Peptides",
                                 x_lab="Sample Names",
                                 y_lab="Abundance Ratios",
                                 fill_lab = "Sample Names",
                                 subtitle_txt = "")
  
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
                            trim=TRUE,
                            subtitle_txt = "")
  
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
  ## TODO: ADD LIMMA
  stat_analysis <- final_imputed_data %>%
    select(pep_with_pos,spectrum_title,accession, starts_with("log10_"))
  
  rownames(stat_analysis) <- paste0(stat_analysis$pep_with_pos,"@",stat_analysis$spectrum_title,"@",stat_analysis$accession,"@",(1:nrow(stat_analysis)))
  

  all_pvalues <- NULL
  
  for (i in 2:sample_size){
    p_values_tmp <- NULL
    
    for(j in 1:dim(stat_analysis)[1]){
      
      if(test_type=="t.test"){
        p_values_tmp[j] <- ttest_func(select(stat_analysis,contains("A1"))[j,],     ### FOR DIFFERENT KIND OF EXP SETUP, 
                                      select(stat_analysis,contains(paste0("A",i)))[j,])  ## It should be defined as an input.
      
      }else if(test_type=="wilcoxon"){
        p_values_tmp[j] <- wilcox_func(select(stat_analysis,contains("A1"))[j,],     ### FOR DIFFERENT KIND OF EXP SETUP, 
                                      select(stat_analysis,contains(paste0("A",i)))[j,])  ## It should be defined as an input.
      }
      
      
    }
    p_values_tmp <- as.data.frame(p_values_tmp)
    colnames(p_values_tmp) <- paste0("A1_","vs_","A",i)
    all_pvalues <- bind_cols(all_pvalues,p_values_tmp)
    rm(p_values_tmp)
  }
  
  
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
  num_reps <- 3
  sample_size <- length(exp_design) / num_reps
  sample_names <- paste0("A",1:sample_size)
  comparisons <- NULL
  for (i in 1:sample_size){
    tmp <- paste0(sample_names[1], "_vs_",sample_names[i])
    comparisons[i] <- tmp
    rm(tmp)
  }
  
  if(test_type=="limma"){
    library(limma)
    
    design_matrix <- model.matrix(~ 0 + factor(rep(sample_names, each = 3)))
    colnames(design_matrix) <-comparisons
    data_limma <- limma::eBayes(limma::lmFit(stat_analysis[,3:18], design_matrix))
    
    ###  adjust.method argument in the topTable function does not directly control 
    # the target false discovery rate (FDR) threshold.
    
    all_pvalues <- topTable(data_limma,number = Inf, adjust.method = "BH",p.value = 1)
    # For A1 vs. all comparison
    
    for (j in 2:length(comparisons)){
      assign(paste0("results_",comparisons[j]),
             topTable(data_limma, coef = comparisons[j],
                      number = Inf, adjust.method = "BH",p.value = 1))
      
    }
    # Filter the results based on the FDR threshold
    # desired_fdr_threshold <- 0.05
    # significant_results_A1_vs_A2 <- subset(results_A1_vs_A2, adj.P.Val <= desired_fdr_threshold)
  
  }
  
  necessary_cols <- as.data.frame(rownames(all_pvalues))
  
  stat_df <- necessary_cols %>%
    rename(mult_col = 1) %>%
    separate(mult_col, c("pep_with_pos", "spectrum_title","accession","id"),sep = "@") %>%
    select(-accession) %>%
    bind_cols(all_pvalues) %>% tibble() ## tibble() will delete rownames! 
                                        ## Be careful that merging correct rownames to pvalues
  
  
  
  
  

}

