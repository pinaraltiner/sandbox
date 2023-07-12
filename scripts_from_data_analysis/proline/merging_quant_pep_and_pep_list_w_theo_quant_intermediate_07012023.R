library(PhosR)
library(stringr)
library(dplyr)
data("phospho_L6_ratio")
ppe <- PhosphoExperiment(assays = list(Quantification = as.matrix(phospho.L6.ratio)))


phospho_hepatocyte_raw <- read.delim("D:/dev/R/phosR/PhosR_STAR_Protocols-main/PhosR_STAR_Protocols-main/Data/PXD001792_raw_hepatocyte.txt", header=TRUE)

###############################################
library(openxlsx)
quant_peptides <- read.xlsx("D:/dev/Desktop_copy/PHD/wet_lab_experiments/experiment2_quantification_peptide_level/PAL - 2023 01 04 - Phosphopeptides exp2 ( 5 conc 3 reps) DDA_2023-01-06_1032.xlsx", sheet = "Quantified peptide ions")

pep_list_w_theo_quant <- read.xlsx("D:/dev/Desktop_copy/PHD/wet_lab_experiments/experiment2_quantification_peptide_level/Synthetic peptides list -060122.xlsx", sheet = "ISO-ref and OTHER with FC")
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
library(data.table)
setnames(df_merge, old = colnames(df_merge[2:length(colnames(quant_peptides_new))]), new =actual_colname_quant_peptides)

write.xlsx(df_merge,file = "D:/dev/Desktop_copy/PHD/wet_lab_experiments/experiment2_quantification_peptide_level/merging_exp_quant_peptides_and_pep_list_theo_quant_exp2_before_stast_analysis.xlsx")
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


abundances <- df_merge_stat_analysis %>% select(starts_with("abundance_"))

# If row count is needed, it can be extracted from here:
#abundances_with_row_count <- cbind(abundances,rowSums(!is.na(abundances)))

# Nothing is changed
abundances_filtered<-abundances[rowSums(!is.na(abundances))>0,]

# Calculate 1 percent quantile of each sample
impute_values <- apply(abundances, 2 , quantile , probs = 0.01 , na.rm = TRUE )

# Impute missing values
for (j in 1:length(impute_values)){
  # Number NA 
  num_NA <- length(abundances_filtered[,j][is.na(abundances_filtered[,j])])
  
  abundances_filtered[,j][is.na(abundances_filtered[,j])] <- impute_values[j]
  
  # After imputation number of imputed values
  num_imp <-length(abundances_filtered[,j][(abundances_filtered[,j]==impute_values[j])])
  
  # This is verification of imputation is done successfully
  # Because we expect to see that number of imputed values should be the same amount as number of NA
  print(setequal(num_NA,num_imp))
  #print(num_NA)
  #print(num_imp)
}


# take log10 
log_10_filtered_abundances <- log10(abundances_filtered)

# Take mean of triplicates of each sample 

# Calculate Fold Change (mean(S1)/mean(S2), etc.)

# Do t-test

# Take -log10() of results
