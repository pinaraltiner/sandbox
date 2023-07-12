#library(PhosR)
library(stringr)
library(dplyr)
library(data.table)
library(openxlsx)
library(ggplot2)
###############################################

quant_peptides <- read.xlsx("D:/dev/Desktop_copy/PHD/wet_lab_experiments/CD4_T_cells_experiment/quant_pep_06022023/PAL _T_cell_Exp3_( 5 conc 3reps)_NoFAIMS_DDA_with_cont_230206_2023-02-07_0947.xlsx", sheet = "Quantified peptide ions")

quant_peptides_with_all <- quant_peptides %>% 
  filter(grepl("_MOUSE",accession)) %>% 
  filter(grepl("Phospho",modifications)) 

quant_peptides_s1 <- na.omit(quant_peptides_with_all[,c(1:6,15,17,95:97)])
quant_peptides_s2 <- na.omit(quant_peptides_with_all[,c(1:6,15,17,98:100)])
quant_peptides_s3 <- na.omit(quant_peptides_with_all[,c(1:6,15,17,101:103)])
quant_peptides_s4 <- na.omit(quant_peptides_with_all[,c(1:6,15,17,104:106)])
quant_peptides_s5 <- na.omit(quant_peptides_with_all[,c(1:6,15,17,107:109)])

length_id_phospho_pep_all_samples <- data.frame(c(dim(quant_peptides_s1)[1],
                                       dim(quant_peptides_s2)[1],
                                       dim(quant_peptides_s3)[1],
                                       dim(quant_peptides_s4)[1],
                                       dim(quant_peptides_s5)[1]),
                                       c("M1 170ng/µL",
                                                     "M2 85ng/µL",
                                                     "M3 17ng/µL",
                                                     "M4 8.5ng/µL",
                                                     "M5 1.7ng/µL"),
                                       c(1:5))

ggplot(length_id_phospho_pep_all_samples,
            aes(x=length_id_phospho_pep_all_samples$c..M1.170ng.µL....M2.85ng.µL....M3.17ng.µL....M4.8.5ng.µL...,
                y=length_id_phospho_pep_all_samples$c.dim.quant_peptides_s1..1...dim.quant_peptides_s2..1...dim.quant_peptides_s3..1...,
                fill=length_id_phospho_pep_all_samples$c..M1.170ng.µL....M2.85ng.µL....M3.17ng.µL....M4.8.5ng.µL...)) + 
  geom_bar(stat = "identity") + 
  theme_light() + scale_y_continuous(breaks = seq(0, 4200, by = 80)) +
  theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
        axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
        plot.title = element_text(size=20),
        legend.title=element_text(size=15),
        axis.text=element_text(size=15,angle = 0),
        axis.title=element_text(size=15)
  ) +   ggtitle("Number of identified T cell enriched phospho-peptides") +
  scale_x_discrete(labels=c("M1 170ng/µL",
                            "M2 85ng/µL",
                            "M3 17ng/µL",
                            "M4 8.5ng/µL",
                            "M5 1.7ng/µL")) +
  labs(x="Sample Names",y="Number of phospho-peptides", fill="Sample name")+
  scale_fill_brewer(palette="Set1")



#pep_list_w_theo_quant <- read.xlsx("D:/dev/Desktop_copy/PHD/wet_lab_experiments/experiment2_quantification_peptide_level/Synthetic peptides list_theo_conc_added_080122.xlsx", sheet = "ISO-ref and OTHER with FC")
#pep_list_w_theo_quant <- read.xlsx("D:/dev/Desktop_copy/PHD/wet_lab_experiments/experiment2_quantification_peptide_level/Synthetic peptides list -060122.xlsx", sheet = "ISO-ref and OTHER with FC")

#pep_list_w_theo_quant <- pep_list_w_theo_quant[,-1]

# Extraction of phospho positions from quant peptides object
phospho_ptm_pos <- lapply(quant_peptides$ptm_protein_positions, function(each_ptm_protein_positions) {
  
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
common_col_exp_quant <- as.data.frame(paste(quant_peptides$sequence, phospho_ptm_pos_df, sep = "_"))
colnames(common_col_exp_quant) <- "common_col_for_merging"

# Keep actual colnames to add after merging (because it contains duplicate col_names which prevents merging btw two df)

quant_peptides_new <- cbind(common_col_exp_quant,quant_peptides)

####################


### If row count is needed, it can be extracted from here:
#abundances_with_row_count <- cbind(abundances,rowSums(!is.na(abundances)))

# Nothing is changed
filtered_abundances<-quant_peptides_new[rowSums(!is.na(select(quant_peptides_new,starts_with("abundance_"))))>0,]

abundances_only_for_impute <- quant_peptides_new %>% select(starts_with("abundance_"))

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
  print(num_NA)
  print(num_imp)
  print(impute_values[j])
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

final_imputed_normalized_data <- cbind(filtered_abundances,filtered_abundances_rowMeans)

final_imputed_normalized_data_tcell <- final_imputed_normalized_data %>% 
  filter(grepl("_MOUSE",accession)) %>% 
  filter(grepl("Phospho",modifications))
write.table(final_imputed_normalized_data, file = "final_imputed_normalized_data_PAL _T_cell_Exp3_( 5 conc 3reps)_NoFAIMS_DDA_with_cont_230206_2023-02-07_0947.txt",sep = "\t",row.names = F)

melt_filtered_abundances_rowMeans_FC <- melt(final_imputed_normalized_data_tcell[,116:119])

melt_filtered_abundances_rowMeans <- melt(filtered_abundances_rowMeans[,1:5])
ggplot(data = melt_filtered_abundances_rowMeans, aes(log2(value),color=variable)) + geom_density()
  geom_violin(colour = "black",size = 0.8) + #bw(0.6)  + #geom_point() +
  theme_bw() +
  theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
        axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
        plot.title = element_text(size=20),
        legend.title=element_text(size=15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15)
  ) +   ggtitle("Experimental Abundances of Phosphopeptides from T cells (log2)") +
  scale_x_discrete(labels=c("A1","A2","A3","A4","A5")) + #scale_y_continuous(breaks = seq(10,34,2)) +
  labs(x="Sample Names",y="Abundance Means", color="Pool Names") +
  scale_fill_brewer(palette="Set1")


ggplot(melt_filtered_abundances_rowMeans_FC,aes(x =variable , y =log2(value),fill = variable)  ) +
  geom_violin(colour = "black",size = 0.8) + #bw(0.6)  + #geom_point() +
  theme_bw() +
  theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
        axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
        plot.title = element_text(size=20),
        legend.title=element_text(size=15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15)
  ) +   ggtitle("Experimental Quantity Ratio of Phosphopeptides from T cells (log2)") +
  scale_x_discrete(labels=c("A1/A2 \n (1.00)","A1/A3 \n (3.21)","A1/A4 \n (4.32)","A1/A5 \n(6.64)")) + scale_y_continuous(breaks = seq(-10,10,2)) +
  labs(x="Sample Names",y="Abundance Ratios", color="Pool Names") +
  scale_fill_brewer(palette="Set1")




############################ # Do t-test # ############################

# The code below does t-test for each row. Because of that, multiple test correction (like BH, Bonferoni)
p_values_for_all_ratio <- NULL
for(k in 2:sample_size){
  
  assign(paste0("t_test_res_A1_to_A",k) , lapply(na.omit(correct_identifed_peps$mean_abundances_A1),
                                                 na.omit(correct_identifed_peps[,paste0("mean_abundances_A",k)], 
                                                         function(x,y) t.test(x,y,alternative = "two.sided", var.equal = TRUE))))
  assign(paste0("p_values_A1_vs_A",k), cbind(p_values_for_all_ratio, get(paste0("t_test_res_A1_to_A",k))$p.value))
  
}

#adjusted_p_values <- as.data.frame(p.adjust(p_values_A1_vs_A5, method = "BH", n = length(p_values_A1_vs_A5)))
adjusted_p_values <- lapply(p_values_for_all_ratio, function(x) p.adjust (x, method = "BH", n = length(x)))

#plot(correct_identifed_peps$mean_abundances_log10_A1,correct_identifed_peps$mean_abundances_log10_A5, pch = 16, col = "blue")
#abline(h = mean(na.omit(correct_identifed_peps$mean_abundances_log10_A1)) - mean(na.omit(correct_identifed_peps$mean_abundances_log10_A5)), col = "red")

boxplot(correct_identifed_peps$mean_abundances_log10_A1,correct_identifed_peps$mean_abundances_log10_A5)
ggplot(correct_identifed_peps, aes(x = mean_abundances_log10_A1, y = mean_abundances_log10_A5)) +
  geom_point(aes(color = "red", size = 5))# +
  #scale_color_discrete(name = "") +
  #scale_size_discrete(name = "Size")
# 
# Take -log10() of results
ggplot(log_10_filtered_abundances_rowMeans_A1_A5,aes(x=log_10_filtered_abundances_rowMeans_A1_A5$mean_abundances_log10_A1,
                                                     y =log_10_filtered_abundances_rowMeans_A1_A5$mean_abundances_log10_A5,
                                                     )) +
  geom_point()+
  #facet_wrap(vars(df_for_figure$common_col_for_merging))  + 
  #geom_smooth(method = "lm", colour = "green", fill = "green") +
  theme_light()

library(gginference)
ggttest(t.test(na.omit(correct_identifed_peps$mean_abundances_log10_A1),na.omit(correct_identifed_peps$mean_abundances_log10_A5), alternative = "two.sided", var.equal = TRUE))



df_for_figure_exp <- correct_identifed_peps[,c(11,74:77)]
df_for_figure <- correct_identifed_peps[,c(11:15,74:77)]
melt_df_for_figure_exp <- melt(df_for_figure_exp)

df_for_figure_theo <- correct_identifed_peps[,c(11:15)]
melt_df_for_figure_theo <- melt(df_for_figure_theo)

p1 <- ggplot(melt_df_for_figure_theo,aes(x =melt_df_for_figure_theo$variable , y =log2(melt_df_for_figure_theo$value),fill = melt_df_for_figure_theo$Pool)  ) +
  geom_boxplot() +
  theme_light() + #scale_y_continuous(limits = c(-60, 60),breaks = seq(-60, 60, by = 20)) +
  theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
        axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
        plot.title = element_text(size=20),
        legend.title=element_text(size=15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15)
  ) +   stat_boxplot(geom = "errorbar") + ggtitle("Theoretical Quantity Ratio of Synthetic Peptides") +
  scale_x_discrete(labels=c("A1/A2","A1/A2","A1/A4","A1/A5")) +
  labs(x="Sample Names",y="Abundance Ratios", color="Pool Names")+
  scale_fill_brewer(palette="Set1")

p2 <-ggplot(melt_df_for_figure_exp,aes(x =melt_df_for_figure_exp$variable , y =log2(melt_df_for_figure_exp$value),fill = melt_df_for_figure_exp$Pool)  ) +
  geom_boxplot() +
  theme_light() +
  theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
        axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
        plot.title = element_text(size=20),
        legend.title=element_text(size=15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15)
  ) +   stat_boxplot(geom = "errorbar") + ggtitle("Experimental  Quantity Ratio of Synthetic Peptides") +
  scale_x_discrete(labels=c("A1/A2","A1/A2","A1/A4","A1/A5")) +
  labs(x="Sample Names",y="Abundance Ratios", color="Pool Names") +
  scale_fill_brewer(palette="Set1")

p3 <-ggplot(melt_df_for_figure_exp,aes(x =melt_df_for_figure_exp$variable , y =melt_df_for_figure_exp$value,color = melt_df_for_figure_exp$Pool)  ) +
  geom_point() +
  theme_light() +
  theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
        axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
        plot.title = element_text(size=20),
        legend.title=element_text(size=15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15)
  ) +   ggtitle("Experimental Quantity Ratio of Synthetic Peptides") +
  scale_x_discrete(labels=c("A1/A2","A1/A3","A1/A4","A1/A5")) +
  labs(x="Sample Names",y="Abundance Ratios", color="Pool Names") +
  scale_color_brewer(palette="Set1")

p4 <-ggplot(melt_df_for_figure_theo,aes(x =melt_df_for_figure_theo$variable , y =melt_df_for_figure_theo$value,color = melt_df_for_figure_theo$Pool)  ) +
  geom_point() +
  theme_light() +
  theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
        axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
        plot.title = element_text(size=20),
        legend.title=element_text(size=15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15),
  ) +   ggtitle("Theoretical Quantity Ratio of Synthetic Peptides") +
  scale_x_discrete(labels=c("A1/A2","A1/A3","A1/A4","A1/A5")) +
  labs(x="Sample Names",y="Abundance Ratios", color="Pool Names") +
  scale_color_brewer(palette="Set1")

p1
library(ggpubr)
ggarrange(p1, p2, common.legend = TRUE, legend="right")
ggarrange(p4, p3, common.legend = TRUE, legend="right")
