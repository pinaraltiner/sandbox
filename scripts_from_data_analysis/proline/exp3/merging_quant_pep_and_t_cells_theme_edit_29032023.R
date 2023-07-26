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

num_reps <- 3
sample_size <- length(exp_design) / num_reps
sample_names <- paste0("A",1:sample_size)
comparisons <- NULL
for (i in 1:sample_size){
  tmp <- paste0(sample_names[1], "_vs_",sample_names[i])
  comparisons[i] <- tmp
  rm(tmp)
}

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
  filtered_abundances_log10_rowMeans <- NULL
  
  for (k in 1:sample_size){
    # If separate version of row means is not needed, it can be commented later.
    # Separate row Means can be collected in temp object to merge in "log_10_filtered_abundances_rowMeans"
    assign(paste0("aft_imp_abundances_A",k),as.data.frame(rowMeans(abundances_for_impute_mouse  %>% select(contains(paste0("A",k))))))
    filtered_abundances_rowMeans<- bind_cols(filtered_abundances_rowMeans,get(paste0("aft_imp_abundances_A",k)))
    
    assign(paste0("log10_abundances_A",k),as.data.frame(log10(abundances_for_impute_mouse  %>% select(contains(paste0("A",k)))))) 
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
  
  final_imputed_data <- cbind(filtered_abundances,filtered_abundances_rowMeans,filtered_abundances_log10,filtered_abundances_log10_rowMeans)
  
  #write.table(final_imputed_data, file = "final_imputed_normalized_data_PAL _T_cell_Exp3_( 5 conc 3reps)_NoFAIMS_DDA_with_cont_230206_2023-02-07_0947.txt",sep = "\t",row.names = F)
  
  df_mean_ab_after_impt <- final_imputed_data %>% 
    select(contains("aft_imp") | contains("pep_with_pos")) %>%
    tibble() %>% 
    pivot_longer(cols = contains("aft_imp"),
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
    select(pep_with_pos,spectrum_title,accession, starts_with("log10_") | starts_with("mean_log10_") | starts_with("exp_FC"))
  
  rownames(stat_analysis) <- paste0(stat_analysis$pep_with_pos,"@",stat_analysis$spectrum_title,"@",stat_analysis$accession,"@",(1:nrow(stat_analysis)))
  if(test_type== "t.test" | test_type== "wilcoxon"){
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
      colnames(p_values_tmp) <- paste0("pvalues_A1","/","A",i)
      all_pvalues <- bind_cols(all_pvalues,p_values_tmp)
      rm(p_values_tmp)
      
    }
    all_pvalues_common_col <- stat_analysis %>% 
      select(pep_with_pos, spectrum_title, accession) %>%
      bind_cols(all_pvalues) %>% 
      pivot_longer(cols = starts_with("pvalues_"), values_to = "pvalues", names_to ="p_ratios") %>%
      separate(p_ratios, into = c("tmp","ratio"),sep = "_") %>%
      select(!tmp) %>%
      mutate(common_col = paste(pep_with_pos,spectrum_title,accession,ratio,1:((sample_size-1)*nrow(stat_analysis)),sep="@"))
    
    ### MERGING I: All used columns are merged and used to combine pvalues and ratios
          ## Before changing the shape of data 
          ## This eliminates any duplication and mismatching btw dataframes
          ## We are 100% sured that all pvalues are associated with its ratio.
      ## Shape of volcano plot looks quite weird esspecially in the A1vsA5.
    merge_stat_df <- stat_analysis %>%
      select(pep_with_pos, spectrum_title, accession, starts_with("exp_FC")) %>%
      bind_cols(all_pvalues) %>%
      pivot_longer(cols = starts_with("exp_FC"), values_to = "fold_change_values", names_to ="fold_change_ratios") %>%
      separate(fold_change_ratios, into = c("tmp","tmp1","ratio"),sep = "_") %>%
      select(!c(tmp,tmp1)) %>%
      mutate(common_col = paste(pep_with_pos,spectrum_title,accession,ratio,1:((sample_size-1)*nrow(stat_analysis)),sep="@")) %>%
      left_join(all_pvalues_common_col, by="common_col")
    
    ### MERGING II: Previous version of combining pvalues and ratios
        ## No errors were appeared when I double pivot_longer()
        ## But there is no common column between those values so that
        ## We cannot know which pvalues correspond to which ratio (MAYBE WE DON'T NEED IT)
    ## Shape of volcano plot looks better when I used this way.
      merge_stat_df1 <- stat_analysis %>%
        select(pep_with_pos, spectrum_title, accession, starts_with("exp_FC")) %>%
        bind_cols(all_pvalues) %>%
        pivot_longer(cols = starts_with("exp_FC"), values_to = "fold_change_values", names_to ="fold_change_ratios") %>%
        pivot_longer(cols = starts_with("pvalues_"), values_to = "pvalues", names_to ="p_ratios")
      
      
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
        assign(paste0("df_A1vsA",i),stat_analysis %>% select(1:3 | contains("log10_E3_A1") | contains(paste0("log10_E3_A",i))))
        # First, linear model was built
        assign(paste0("fit",i) ,lmFit(get(paste0("df_A1vsA",i))[,3:9], design_matrix))
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
        rename_with(.col=9, ~ "ratios") %>%
        separate(mult_col, into = c("pep_with_pos", "spectrum_title","accession","id"),sep = "@")
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
    
    actual_ratio <- c(2,10,20,100)
    actual_ratio <- data.frame(Ratio_col=unique(merge_stat_df_pivot$fold_change_ratios),-1*(log2(c(2,10,20,100))))
    
    
    ggplot(merge_stat_df_final ,aes(x = merge_stat_df_final$logFC, y = -log10(merge_stat_df_final$P.Value), color=merge_stat_df_final$ratios)) +
      geom_point( size = 2.5) +
      facet_wrap(~ratios) +
    #geom_line(aes(color = new_col_coloring), size = 1) +  # Add color aesthetic to geom_line()
    # scale_color_manual(values = c("ISO-REF" = "#000000", "unexpected" = "#999999",
    #                               "Others_multi_A1-A2_Ratio" = "#CC79A7",
    #                               "Others_mono_A1-A2_Ratio" = "#CC79A7",
    #                               "Others_multi_A1-A3_Ratio" = "#E69F00",
    #                               "Others_mono_A1-A3_Ratio" = "#E69F00",
    #                               "Others_multi_A1-A4_Ratio" = "#56B4E9",
    #                               "Others_mono_A1-A4_Ratio" = "#56B4E9",
    #                               "Others_multi_A1-A5_Ratio" = "#009E73",
    #                               "Others_mono_A1-A5_Ratio" = "#009E73"),
    #                    labels = c('Non-variant', 'Variant non-isomeric A1 vs A2',
    #                               'Variant non-isomeric A1 vs A3',
    #                               'Variant non-isomeric A1 vs A4',
    #                               'Variant non-isomeric A1 vs A5',
    #                               'Variant isomeric A1 vs A2',
    #                               'Variant isomeric A1 vs A3',
    #                               'Variant isomeric A1 vs A4',
    #                               'Variant isomeric A1 vs A5',
    #                               'Unexpected')) +
    #scale_shape_manual(values = c(16, 15, 12, 17),
                       #labels = comparisons) +
      scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 0.8)) +
      scale_x_continuous(limits = c(-3,3),breaks = seq(-3, 3, by = 0.8)) +
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
      #expand_limits(x = 0, y = 0) +
      #geom_vline(data = actual_ratio, aes(xintercept = actual_ratio$X.1....log2.c.2..10..20..100..., size = 1, show.legend = FALSE)) + #color=c("#CC79A7","#E69F00","#56B4E9","#009E73")
      #geom_hline(data = log10_p_thresholds, aes(yintercept = log10_p_thresholds$X.log10.p_thresholds.),color=c("#CC79A7","#E69F00","#56B4E9","#009E73"), size = 1, linetype = 2, show.legend = FALSE)+ 
      labs(title = "Experiment 3 - DDA TIMS-TOF processed by Proline", color = "Classes", shape="Type") 
    
    
  }

  
  
  
  



