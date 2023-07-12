#library(PhosR)
library(stringr)
library(dplyr)
library(data.table)
library(openxlsx)
library(tidyr)
###############################################
source("D:/dev/Desktop_copy/PHD/data_analysis/scripts/ggplot_functions.R")

final_proline_pep_quant_analysis <- function(file_path,
                                             file_name,
                                             sheet_name,
                                             theo_file_path,
                                             theo_file_name,
                                             sheet_theo_name){
  quant_peptides <- read.xlsx(paste0(file_path,file_name), sheet = sheet_name)
  selected_cols <- colnames(quant_peptides[c(1:4,15)])
  
  syn_phospho_pep <- quant_peptides %>% 
    filter(grepl("Phospho", ptm_protein_positions)) %>% 
    filter(grepl("_HUMAN", accession)) %>%
    select(contains(selected_cols) | starts_with("abundance_"))
  
  
  pep_list_w_theo_quant <- read.xlsx(paste0(theo_file_path, theo_file_name), sheet = sheet_theo_name)
  pep_list_w_theo_quant <- pep_list_w_theo_quant[,-1]
  
  
  # Extraction of phospho positions from quant peptides object
  phospho_ptm_pos <- lapply(syn_phospho_pep$ptm_protein_positions, function(each_ptm_protein_positions) {
    
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
  common_col_exp_quant <- as.data.frame(paste(syn_phospho_pep$sequence, phospho_ptm_pos_df, sep = "_"))
  colnames(common_col_exp_quant) <- "common_col_for_merging"
  
  # Bind it to the quant data
  syn_phospho_pep_new <- cbind(common_col_exp_quant,syn_phospho_pep)
  
  # Creation of common column merging peptide sequence and phospho positions -> theoretical data
  common_col_theo_quant <- as.data.frame(paste(pep_list_w_theo_quant$Phospopeptide.sequence,
                                               pep_list_w_theo_quant$modified.position.in.peptide, sep = "_"))
  colnames(common_col_theo_quant) <- "common_col_for_merging"
  pep_list_w_theo_quant_new <- cbind(common_col_theo_quant,pep_list_w_theo_quant)
  
  # Combine data frame (we will continue with this for further step)
  df_merge <- syn_phospho_pep_new %>% 
    left_join(pep_list_w_theo_quant_new,by="common_col_for_merging")
  
  ## ALWAYS KEEP IT COMMENTED TO PREVENT OVERWRITE
  #write.xlsx(df_merge,file = "D:/dev/Desktop_copy/PHD/wet_lab_experiments/experiment2_quantification_peptide_level/merging_exp_syn_phospho_pep_and_pep_list_theo_quant_exp2_before_stast_analysis_08012023.xlsx")
  
  ####################
  
  # Sorting columns A1 to A5 and R1 to R3
  sorted_colnames <- c("psm_count","abundance") # Abundance contains both w/wo raw_abundance
  for (i in 1:length(sorted_colnames)){
    assign(paste0("tmp",i), sort(colnames(select(df_merge,
                                                 contains(sorted_colnames[i])))))
  }
  
  sort_col <- c(tmp1,tmp2)
  df_merge_sorted_col <- df_merge[, sort_col]
  
  df_stat_analysis <- df_merge %>% select(-starts_with("abundance")) %>% mutate(df_merge_sorted_col)
  
  ### If row count is needed, it can be extracted from here:
  #abundances_with_row_count <- cbind(abundances,rowSums(!is.na(abundances)))
  
  # Nothing is changed -> If all replicates of one sample are NA, we will exclude them for the imputation
  filtered_abundances<-df_stat_analysis[rowSums(!is.na(select(df_stat_analysis,starts_with("abundance_"))))>0,]
  
  abundances_only_for_impute <- df_stat_analysis %>% 
    select(starts_with("abundance_"))
  
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
  sample_size <- 5
  filtered_abundances_rowMeans <- NULL
  for (k in 1:sample_size){
    # If separate version of row means is not needed, it can be commented later.
    # Separate row Means can be collected in temp object to merge in "log_10_filtered_abundances_rowMeans"
    assign(paste0("abundances_A",k),as.data.frame(rowMeans(abundances_only_for_impute  %>% select(contains(paste0("A",k))) %>%
                                                             select(starts_with("abundance_")))))
    filtered_abundances_rowMeans<- bind_cols(filtered_abundances_rowMeans,get(paste0("abundances_A",k)))
    
  }
  
  colnames(filtered_abundances_rowMeans) <- paste0("mean_abundances_A",1:sample_size)
  
  # Calculate Fold Change by keeping A1 constant (mean(S1)/mean(S2), etc.)
  cols <- ncol(filtered_abundances_rowMeans)
  for(An in 2:cols){
    filtered_abundances_rowMeans[,paste0("exp_FC_A1/A",An)] <- filtered_abundances_rowMeans[,1]/filtered_abundances_rowMeans[,An]
    
  }
  # To calculate all binary combination in the data frame
  #mat <- do.call(cbind, lapply(cols, function(xj) 
  #  sapply(cols, function(xi) (filtered_abundances_rowMeans[, xj]/(filtered_abundances_rowMeans[, xj])))))
  #colnames(mat) <-  outer(names(filtered_abundances_rowMeans), names(filtered_abundances_rowMeans), paste0)
  
  final_imputed_normalized_data <- cbind(df_stat_analysis,filtered_abundances_rowMeans)
  
  unexpectedly_id_peps <- final_imputed_normalized_data %>% 
    filter_at(vars(Pool), all_vars(is.na(.)))
  
  col_sel_stat_analysis <- c("common_col_for_merging",
                             "pool",
                             "A1-A2_Ratio",
                             "A1-A3_Ratio",
                             "A1-A3_Ratio",
                             "A1-A4_Ratio",
                             "A1-A5_Ratio",
                             "exp_FC_A1/A2",
                             "exp_FC_A1/A3",
                             "exp_FC_A1/A4",
                             "exp_FC_A1/A5")
  
  #### Distribution of experimental quantitative ratios with two different pools
  
  library(ggplot2)
  
  ## This can be used as an object name
  df_for_figure <- final_imputed_normalized_data %>% 
    select(contains(col_sel_stat_analysis)) %>%
    pivot_longer(cols = col_sel_stat_analysis[3:11], names_to = "Theo_Exp", values_to = "Ratios",values_drop_na = T) %>%
    #select(!contains(c("common_col_for_merging"))) %>%
    filter(grepl("exp_FC",Theo_Exp)) %>%
    filter_at(vars(Pool), all_vars(!is.na(.)))
  
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
  
  figure_with_mean_abundance <- final_imputed_normalized_data %>%
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

