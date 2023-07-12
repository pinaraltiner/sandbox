library(stringr)
library(dplyr)
library(data.table)
#library(openxlsx)
library(ggplot2)
library(reshape2)
library(readr)
library(tibble)
library(purrr)
library(readr)
library(tidyr)
###############################################
## MAXQUANT DATA ANALYSIS

### Source code was taken from here: https://rdrr.io/github/singjc/mstools/src/R/getModificationPosition.R
  ## The code was modified based on what I want and based on software input tyoe (DIANN-and MaxQuant)
source("D:/dev/Desktop_copy/PHD/data_analysis/scripts/getModificationPositionMQ_func.R")

setwd("D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/MQ_data_analysis/exp2_wo_FAIMS/exp2_wo_FAIMSwith_MBR/")

exp2_wo_faims <- read_tsv("evidence.txt")


theo_file_path <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/"
theo_file_name <- "Synthetic peptides list_theo_conc_added_pool_id_iso_count.xlsx"
sheet_theo_name <- "ISO-ref and OTHER with FC"

pep_list_w_theo_quant <- read.xlsx(paste0(theo_file_path, theo_file_name), sheet = sheet_theo_name)
pep_list_w_theo_quant <- pep_list_w_theo_quant[,-1]

common_col_theo_quant <- as.data.frame(paste(pep_list_w_theo_quant$Phosphopeptide.sequence,
                                             pep_list_w_theo_quant$modified.position.in.peptide, sep = "_"))
colnames(common_col_theo_quant) <- "common_col_for_merging"
pep_list_w_theo_quant_new <- cbind(common_col_theo_quant,pep_list_w_theo_quant)

## NECESSARY FOR COUNTING NUM OF IDENTIFIED ECOLIPEPTIDES
id_ecoli_pep <- exp2_wo_faims %>% 
  filter(grepl(83333,`Taxonomy IDs`)) %>% 
  group_by(Sequence) %>% 
  summarise(num_id_ecoli_pep = n(), mods=unique(Modifications))


## NECESSARY FOR COUNTING NUM OF IDENTIFIED PHOSPHO PEPTIDES
id_syn_phospho_pep <- exp2_wo_faims %>% 
  filter(grepl(9606,`Taxonomy IDs`)) %>% 
  filter(grepl("Phospho",Modifications)) 



df2 <- apply(id_syn_phospho_pep[,"Modified sequence"],1,getModificationPosition_MQ)

results1 <- map_dfr(df2, ~ enframe(.x)) %>%
  filter(grepl("modification_",name)| grepl("pep_seq", name)) %>%
  mutate(value = map_chr(value, str_c, collapse="&")) %>%
  mutate(mods=case_when(grepl("Phospho (STY)",fixed = T,name) ~ "phospho",
                        grepl("Oxidation (M)",fixed = T,name) ~ "Oxidation",
                        grepl("(Acetyl (Protein N-term))",fixed = T,name) ~ "N-term_Acetyl",
                        TRUE ~ ""))

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
  left_join(filter(results_with_index, name == "pep_seq"), by = "id") %>% 
  ungroup() %>%
  ## SELECTING USEFUL COLUMNS
  select(c(name.x,value.x,mods.x,value.y))


result_with_common_col <- reshaped_results %>%
  filter(grepl("modification_(Phospho (STY))",fixed = T,name.x)) %>%
  mutate(common_col_for_merging = paste(value.y,value.x, sep = "_")) %>%
  select(!value.y) %>%
  rename("mods_id" = "name.x",
         "phospho_positions" ="value.x",
         "all_mods_with_mod_type" = "mods.x")


#getModificationPosition_MQ(mod_seq = id_syn_pep$`Modified sequence`,F)

final_results_with_common_col <- mutate(result_with_common_col,id_syn_phospho_pep)


# Combine data frame (we will continue with this for further step)
df_merge <- final_results_with_common_col %>% 
  left_join(pep_list_w_theo_quant_new,by="common_col_for_merging")

sep_df_merge <- df_merge %>% separate(Experiment, 
                                      c('Exp_id', 'Sample_id', 'Rep_id'),sep = "-",remove = FALSE)

quantile_values <- sep_df_merge %>%
  group_by(Sample_id) %>% 
  summarise(first_quantile = quantile(Intensity,probs=0.01, na.rm=TRUE))


imputed_abun_df_merge <- sep_df_merge %>% 
  mutate(Intensity = case_when(
    Sample_id == "A1" ~ if_else(is.na(Intensity), as.numeric(quantile_values[1,2]), Intensity),
    Sample_id == "A2" ~ if_else(is.na(Intensity), as.numeric(quantile_values[2,2]), Intensity),
    Sample_id == "A3" ~ if_else(is.na(Intensity), as.numeric(quantile_values[3,2]), Intensity),
    Sample_id == "A4" ~ if_else(is.na(Intensity), as.numeric(quantile_values[4,2]), Intensity),
    Sample_id == "A5" ~ if_else(is.na(Intensity), as.numeric(quantile_values[5,2]), Intensity)
    ))

write.table(imputed_abun_df_merge, file = "MQ_exp2_woFAIMS_evidences_after_imputation.txt",sep = "\t")
#filtered_abundances<-final_results_with_common_col[rowSums(!is.na(select(quant_peptides_with_all,starts_with("abundance_"))))>0,]

test <- imputed_abun_df_merge %>% 
  group_by(common_col_for_merging, Charge) %>%
  slice(which.max(Intensity)) %>%
  ungroup()
  #summarise(max_row_sum = max(row_sum, na.rm = TRUE)) %>%
  

toy_data <- read_tsv(file = "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/comparision_DDA_noFAIMS_mq_proline_pd/toy_data_MQ.txt")


final_selected_abun_df_merge <- imputed_abun_df_merge  %>%  
  group_by(common_col_for_merging, Experiment) %>% 
  summarise(Intensity = max(Intensity)) %>%
  ungroup()

########################


final_data <-final_selected_abun_df_merge  %>% 
  left_join(imputed_abun_df_merge,by=c("common_col_for_merging","Experiment","Intensity"),keep=FALSE)


### TODO: COMPLETE FOR THE EXTRACTION OF MEAN ABUNDANCES WHICH WILL BE NEEDED for VOLCANO PLOT.
final_data %>% 
  filter(grepl("A1",Experiment)) %>% 
  group_by(common_col_for_merging,Sample_id) %>%
  summarise(mean_intensity_A1 = mean(Intensity))



final_data_log10_abun <- final_selected_abun_df_merge %>%
  mutate(log10_int = log10(Intensity))




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




A1_final_data_log10_abun <- final_data_log10_abun %>% 
  filter(grepl("A1",Experiment))

A2_final_data_log10_abun <- final_data_log10_abun %>% 
  filter(grepl("A",Experiment))

data_wide <- data %>% spread(final_data_log10_abun log10_int)



for (i in 1:dim(common_col_theo_quant)[1]){
  
}

# iterate through each column of the data set and calculate t-tests

p_values_12 <- NULL
p_values_13 <- NULL
p_values_14 <- NULL
p_values_15 <- NULL


for(j in 1:dim(stat_analysis)[1]){
  
  p_values_12[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,5:7])
  p_values_13[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,8:10])
  p_values_14[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,11:13])
  p_values_15[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,14:16])
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
p_thresholds










intensity_match <- final_results_with_common_col %>%
  select(common_col_merging,Intensity,Experiment) %>% 
  tibble() %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = "Experiment",
              values_from = "Intensity")



  
  
  
  
  
  