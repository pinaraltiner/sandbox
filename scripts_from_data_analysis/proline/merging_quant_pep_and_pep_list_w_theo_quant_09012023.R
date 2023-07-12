library(PhosR)
library(stringr)
library(dplyr)
library(data.table)
library(openxlsx)

###############################################

quant_peptides <- read.xlsx("D:/dev/Desktop_copy/PHD/wet_lab_experiments/experiment2_quantification_peptide_level/PAL - 2023 01 04 - Phosphopeptides exp2 ( 5 conc 3 reps) DDA_2023-01-06_1032.xlsx", sheet = "Quantified peptide ions")

pep_list_w_theo_quant <- read.xlsx("D:/dev/Desktop_copy/PHD/wet_lab_experiments/experiment2_quantification_peptide_level/Synthetic peptides list_theo_conc_added_080122.xlsx", sheet = "ISO-ref and OTHER with FC")
#pep_list_w_theo_quant <- read.xlsx("D:/dev/Desktop_copy/PHD/wet_lab_experiments/experiment2_quantification_peptide_level/Synthetic peptides list -060122.xlsx", sheet = "ISO-ref and OTHER with FC")

pep_list_w_theo_quant <- pep_list_w_theo_quant[,-1]

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
actual_colname_quant_peptides <- colnames(quant_peptides)
colnames(quant_peptides)<- c(1:length(colnames(quant_peptides)))

quant_peptides_new <- cbind(common_col_exp_quant,quant_peptides)

# Creation of common column merging peptide sequence and phospho positions -> theoretical data
common_col_theo_quant <- as.data.frame(paste(pep_list_w_theo_quant$Phospopeptide.sequence, pep_list_w_theo_quant$modified.position.in.peptide, sep = "_"))
colnames(common_col_theo_quant) <- "common_col_for_merging"
pep_list_w_theo_quant_new <- cbind(common_col_theo_quant,pep_list_w_theo_quant)

# Combine data frame (we will continue with this for further step)
df_merge <- quant_peptides_new %>% left_join(pep_list_w_theo_quant_new,by="common_col_for_merging")

# Paste the actual col names of quant_peptides (exp data) after merging two of them

setnames(df_merge, old = colnames(df_merge[2:length(colnames(quant_peptides_new))]), new =actual_colname_quant_peptides)

write.xlsx(df_merge,file = "D:/dev/Desktop_copy/PHD/wet_lab_experiments/experiment2_quantification_peptide_level/merging_exp_quant_peptides_and_pep_list_theo_quant_exp2_before_stast_analysis_08012023.xlsx")
####################

# Sorting columns A1 to A5 and R1 to R3
sorted_colnames <- c("psm_count","abundance") # Abundance contains both w/wo raw_abundance
for (i in 1:length(sorted_colnames)){
  assign(paste0("tmp",i), sort(colnames(select(df_merge,contains(sorted_colnames[i])))))
}

sort_col <- c(tmp1,tmp2)
df_merge_sorted_col <- df_merge[, sort_col]

# This object will be used for statistical analysis
df_merge_stat_analysis <- cbind(df_merge$sequence, df_merge$modifications, df_merge$common_col_for_merging, df_merge$accession, df_merge_sorted_col)
colnames(df_merge_stat_analysis)[1:4] <- c("sequence","modifications","common_col_for_merging", "accession")

### If row count is needed, it can be extracted from here:
#abundances_with_row_count <- cbind(abundances,rowSums(!is.na(abundances)))

# Nothing is changed
filtered_abundances<-df_merge_stat_analysis[rowSums(!is.na(select(df_merge_stat_analysis,starts_with("abundance_"))))>0,]

abundances_only_for_impute <- df_merge_stat_analysis %>% select(starts_with("abundance_"))

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

final_imputed_normalized_data <- cbind(filtered_abundances,filtered_abundances_rowMeans)


# Some peptides has more than PSM in that case theo values repeated while only exp values are different.
# Because of that it has 268 rows
correct_identifed_peps <- pep_list_w_theo_quant_new %>% left_join(final_imputed_normalized_data,by="common_col_for_merging")

df_for_figure <- correct_identifed_peps[,c(1,11:20,69:77)]
library(reshape2)
melt_df_for_figure <- melt(df_for_figure[,c(2,12:16)])

library(ggplot2)
ggplot(df_for_figure,aes(x =df_for_figure$mean_abundances_log10_A1 , y = df_for_figure$mean_abundances_log10_A2, color=df_for_figure$Pool)) +
  geom_point()+
  #facet_wrap(vars(df_for_figure$common_col_for_merging))  + 
  #geom_smooth(method = "lm", colour = "green", fill = "green") +
  theme_light() #+
  geom_errorbar(aes(ymin=(df_for_figure$`A1_others_conc_fmol/ul`-df_for_figure$mean_abundances_log10_A1), 
                    ymax=(df_for_figure$`A2_others_conc_fmol/ul`-df_for_figure$mean_abundances_log10_A2)), width=.2,
                position=position_dodge(0.05))


ggplot(melt_df_for_figure , aes(x=value,fill=Pool)) +
  geom_density(alpha=0.8)

# Do t-test
# The code below does t-test for each row. Because of that, multiple test correction (like BH, Bonferoni)

t_test_res <- t.test(t.test(na.omit(correct_identifed_peps$mean_abundances_log10_A1),na.omit(correct_identifed_peps$mean_abundances_log10_A5), alternative = "two.sided", var.equal = TRUE))
p_values_A1_vs_A5 <- t_test_res$p.value

adjusted_p_values <- as.data.frame(p.adjust(p_values_A1_vs_A5, method = "BH", n = length(p_values_A1_vs_A5)))


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


library(ggpubr)
ggarrange(p1, p2, common.legend = TRUE, legend="right")
ggarrange(p4, p3, common.legend = TRUE, legend="right")


iso_exp <- df_for_figure_exp %>% filter(grepl("ISO-REF",Pool))
melt_iso_theo <- melt(iso_exp)
ggplot(melt_iso_theo,aes(x =melt_iso_theo$variable , y =log2(melt_iso_theo$value),color = melt_iso_theo$Pool)  ) +
  geom_violin()  + #geom_point() +
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


other_exp <- df_for_figure_exp %>% filter(grepl("Others",Pool))
melt_other_exp <- melt(other_exp)
ggplot(melt_other_exp,aes(x =melt_other_exp$variable , y =log2(melt_other_exp$value),color = melt_other_exp$Pool)  ) +
  geom_violin() + #geom_point() +
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


