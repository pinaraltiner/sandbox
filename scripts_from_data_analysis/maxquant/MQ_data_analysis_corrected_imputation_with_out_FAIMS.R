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
library(openxlsx)
###############################################
## MAXQUANT DATA ANALYSIS

### Source code was taken from here: https://rdrr.io/github/singjc/mstools/src/R/getModificationPosition.R
  ## The code was modified based on what I want and based on software input tyoe (DIANN-and MaxQuant)
source("D:/dev/Desktop_copy/PHD/data_analysis/scripts/getModificationPositionMQ_func.R")
source("D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/comparision_DDA_noFAIMS_mq_proline_pd/roc_curve_generation_proline.R")

#setwd("D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/MQ_data_analysis/exp2_wo_FAIMS/exp2_wo_FAIMSwith_MBR/")
setwd("D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/MQ_data_analysis/with_FAIMS_MBR")
#exp2_wo_faims <- read_tsv("evidence.txt")
exp2_with_faims <- read_tsv("evidence.txt")

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
id_syn_phospho_pep <- exp2_with_faims %>% 
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

final_results_max_abun <- final_results_with_common_col %>% 
  group_by(common_col_for_merging,Experiment) %>%
  slice(which.max(Intensity)) %>%
  ungroup()

df_wide <- final_results_max_abun %>% 
  select(common_col_for_merging,Experiment, Intensity) %>%
  pivot_wider(names_from = Experiment,values_from = Intensity)


abundances_for_impute <- df_wide %>% 
  select(order(colnames(df_wide))) %>% 
  select(-common_col_for_merging)
  
experiment_name <- colnames(abundances_for_impute)


#######################################
triplicate_indx <- c(1,3,4,6,7,9,10,12,13,15)
impute_values_2NA <- NULL
impute_values_1NA <- NULL
tmp_2NA <- NULL
tmp_1NA <- NULL
#imputed_abundances <- NULL
for (i in 1:5){
  tmp_1NA <- quantile (abundances_for_impute[,experiment_name[triplicate_indx[1]:triplicate_indx[2]]],
                       probs = 0.01 , na.rm = TRUE )
  tmp_2NA <- quantile (abundances_for_impute[,experiment_name[triplicate_indx[1]:triplicate_indx[2]]],
                       probs = 0.001 , na.rm = TRUE )
  impute_values_2NA[i] <- as.numeric(tmp_2NA)
  impute_values_1NA[i] <- as.numeric(tmp_1NA)
  
  triplicate_indx <- triplicate_indx[-c(1:2)]
  
  # tmp1 <- abundances_for_impute %>% select(contains(paste0("A",i))) %>%
  #   mutate(across(contains(paste0("A",i)), ~ifelse(is.na(.), impute_values[i], .)))
  # imputed_abundances <- bind_cols(imputed_abundances, tmp1)
}

df_long <- df_wide %>%
  pivot_longer(cols = -common_col_for_merging, names_to = "column", values_to = "value") %>%
  separate(column, into = c("col_group", "col_index", "row_index"), sep = "-")


# Count the number of NA values in each row for each column group
na_counts <- df_long %>%
  group_by(common_col_for_merging, col_index) %>%
  summarise(na_count = sum(is.na(value))) 
#pivot_wider(names_from = col_index, values_from = na_count)

# Replace the NA values based on the count
df_result <- df_long %>%
  left_join(na_counts, by = c("common_col_for_merging", "col_index")) %>%
  mutate(value = case_when(
    na_count == 3 ~ 0,
    
    na_count == 2 & col_index == "A1" ~ impute_values_2NA[1],
    na_count == 1 & col_index == "A1" ~ impute_values_1NA[1],
    
    na_count == 2 & col_index == "A2" ~ impute_values_2NA[2],
    na_count == 1 & col_index == "A2" ~ impute_values_1NA[2],
    
    na_count == 2 & col_index == "A3" ~ impute_values_2NA[3],
    na_count == 1 & col_index == "A3" ~ impute_values_1NA[3],
    
    na_count == 2 & col_index == "A4" ~ impute_values_2NA[4],
    na_count == 1 & col_index == "A4" ~ impute_values_1NA[4],
    
    na_count == 2 & col_index == "A5" ~ impute_values_2NA[5],
    na_count == 1 & col_index == "A5" ~ impute_values_1NA[5],
    TRUE ~ value
  )) %>%
  select(-na_count) %>% unite(new_col,col_group,col_index, row_index,sep = "_") %>%
  pivot_wider(names_from = new_col, values_from = value) %>%
  #select(-col_group) %>% 
  filter(rowSums(across(where(is.numeric)))!=0) %>%
  select(order(colnames(df_wide)))

## TEST DATA for imputation part: 
#toy_data <- read_tsv(file = "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/comparision_DDA_noFAIMS_mq_proline_pd/toy_data_MQ.txt")


#######################################
#write.table(imputed_abun_df_merge, file = "MQ_exp2_woFAIMS_evidences_after_imputation.txt",sep = "\t")




# Combine data frame (we will continue with this for further step)
df_merge <- df_result %>% 
  left_join(pep_list_w_theo_quant_new,by="common_col_for_merging")

#df_merge <- final_results_max_abun  %>% inner_join(df_result, 
#                                    by="common_col_for_merging") %>% 
#  left_join(pep_list_w_theo_quant_new,by="common_col_for_merging")

# sep_df_merge <- df_merge %>% separate(Experiment, 
#                                       c('Exp_id', 'Sample_id', 'Rep_id'),sep = "-",remove = FALSE)
# final_data <-final_selected_abun_df_merge  %>% 
#   left_join(imputed_abun_df_merge,by=c("common_col_for_merging","Experiment","Intensity"),keep=FALSE)


########################
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
  rename_with(~ paste0("log10_", .x), everything()) %>%
  mutate_all(function(x) ifelse(is.na(x), 0, x))  %>%
  mutate_all(function(x) ifelse(is.infinite(x), 0, x))



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

#test <- final_imputed_data

unexpectedly_id_peps <- final_imputed_data %>% 
  filter_at(vars(Pool), all_vars(is.na(.)))


output_path <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/comparision_DDA_noFAIMS_mq_proline_pd/"
#write.xlsx(unexpectedly_id_peps,file = paste0(output_path,"Corrected_imputation_proline_unexpectedly_identified_phosphopeptides_exp2_noFAIMS_DDA.xlsx"))

###########################################

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
  
  if (round((sum(stat_analysis[j,2:4])/3),3) !=  round(stat_analysis[j,2:4][1],3) | round((sum(stat_analysis[j,5:7])/3),3) !=  round(stat_analysis[j,5:7][1],3)){
    
    p_values_12[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,5:7])
    
  }else{
    p_values_12[j] <- NA
    
  }
  if (round((sum(stat_analysis[j,2:4])/3),3) !=  round(stat_analysis[j,2:4][1],3) | round((sum(stat_analysis[j,8:10])/3),3) !=  round(stat_analysis[j,8:10][1],3)){

    p_values_13[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,8:10])

  }else{
    p_values_13[j] <- NA
  }
  if (round((sum(stat_analysis[j,2:4])/3),3) !=  round(stat_analysis[j,2:4][1],3) | round((sum(stat_analysis[j,11:13])/3),3) !=  round(stat_analysis[j,11:13][1],3)){


    p_values_14[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,11:13])

  }else{
    p_values_14[j] <- NA

  }
  if (round((sum(stat_analysis[j,2:4])/3),3) !=  round(stat_analysis[j,2:4][1],3) | round((sum(stat_analysis[j,14:16])/3),3) !=  round(stat_analysis[j,14:16][1],3)){

    p_values_15[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,14:16])
  }else{
    p_values_15[j] <- NA

  }
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
  p_values_12[j,"test12"] <- (p_values_12[j,"rank12"]/test_num)*0.1
  p_values_12[j,"test_bool"] <- p_values_12[j,"test12"] > p_values_12[j,"p_values_12"]
  
  
  p_values_13 <- p_values_13[order(p_values_13$p_values_13),]
  colnames(p_values_13)[1] <- "common_col_for_merging"
  p_values_13["rank13"] <- 1:test_num
  p_values_13[j,"test13"] <- (p_values_13[j,"rank13"]/test_num)*0.1
  p_values_13[j,"test_bool"] <- p_values_13[j,"test13"] > p_values_13[j,"p_values_13"]
  
  
  p_values_14 <- p_values_14[order(p_values_14$p_values_14),]
  colnames(p_values_14)[1] <- "common_col_for_merging"
  p_values_14["rank14"] <- 1:test_num
  p_values_14[j,"test14"] <- (p_values_14[j,"rank14"]/test_num)*0.1
  p_values_14[j,"test_bool"] <- p_values_14[j,"test14"] > p_values_14[j,"p_values_14"]
  
  
  p_values_15 <- p_values_15[order(p_values_15$p_values_15),]
  colnames(p_values_15)[1] <- "common_col_for_merging"
  p_values_15["rank15"] <- 1:test_num
  p_values_15[j,"test15"] <- (p_values_15[j,"rank15"]/test_num)*0.1
  p_values_15[j,"test_bool"] <- p_values_15[j,"test15"] > p_values_15[j,"p_values_15"]
  
}


#######
## Without FAIMS
#p_thresholds<- c(0.0194189692,1.492922e-02,2.983289e-02,2.558849e-02)

## with FAIMS
p_thresholds <- c(0.0005494505,0.002197802,0.0010989011,7.825033e-04)
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
           mutate(ls_correctness = ifelse(test_bool == FALSE, 0, 1)))
  
  assign(paste0("roc_data",i),compute_roc_curve(get(paste0("complete_p_values_",i)),flag = "Others",expected = 145 ))
  #assign(paste0("roc_data_unexpected",i),compute_roc_curve(get(paste0("complete_p_values_",i)),flag = "unexpected",expected = 100))
  #assign(paste0("roc_data_isoref",i),compute_roc_curve(get(paste0("complete_p_values_",i)),flag = "ISO-REF",expected = 37 ))
}


## ROC CURVE GENERATION ### 
roc_data121 <- cbind(roc_data12, "A1 vs A2","Others", "MaxQuant")
roc_data131 <- cbind(roc_data13, "A1 vs A3","Others", "MaxQuant")
roc_data141 <- cbind(roc_data14, "A1 vs A4","Others", "MaxQuant")
roc_data151 <- cbind(roc_data15, "A1 vs A5","Others", "MaxQuant")

colnames(roc_data121)[c(4:6)] <- c("Comparison","Pool","Software")
colnames(roc_data131)[c(4:6)] <- c("Comparison","Pool","Software")
colnames(roc_data141)[c(4:6)] <- c("Comparison","Pool","Software")
colnames(roc_data151)[c(4:6)] <- c("Comparison","Pool","Software")



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



roc_final <- rbind(roc_data121,roc_data131,roc_data141,roc_data151)
write.table(roc_final,file="roc_curve_data_exp2_noFAIMS_MQ.txt", sep = "\t", row.names = F)
##TODO: make it more professional
ggplot(roc_final,aes(x=fdp,y=tpr,color=Comparison,linetype=Pool)) + geom_line(size=1) +
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
  #facet_wrap(~Ratio_col,scales = "free_x") +
  labs(title = "Experiment 2 - DDA no FAIMS processed by MaxQuant")

##TODO: make it more professional
ggplot(roc_data12,aes(x=fdp,y=tpr)) + geom_line()


#######
merging_all_pvalues <- final_imputed_data %>% 
  select(col_sel_stat_analysis) %>%
  left_join(p_values_12,by="common_col_for_merging") %>%
  left_join(p_values_13,by="common_col_for_merging") %>%
  left_join(p_values_14,by="common_col_for_merging") %>%
  left_join(p_values_15,by="common_col_for_merging")

merging_all_pvalues <- merging_all_pvalues %>%
  mutate(Pool=replace_na(Pool,"unexpected"))

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
  unite(Pool_new, Pool, isomericity,sep = "_") 
  #mutate(exp_FC_value = ifelse(is.infinite(exp_FC_value), 0, exp_FC_value))

ggplot(na.omit(volcano_final1),aes(x=log2(exp_FC_value),y=-log10(pvalues_value),color=Ratio_col)) +
  geom_point(aes(shape=Pool_new),size=2.5) + ###### CHANGE NA to UNEXPECTED to MAKE THEM VISIBLE
  #geom_label( 
  #  data=volcano_final %>% filter(grepl("Others",Pool) & grepl("mono", isomericity)), # Filter data first
  #  aes(label=isomericity)
  # ) +
  #scale_x_continuous(breaks=seq() + 
  scale_y_continuous(breaks=seq(0, max(-log10(na.omit(volcano_final$pvalues_value))), length.out = 21)) +
  #scale_x_continuous(breaks=seq(0, max(log2(na.omit(volcano_final$exp_FC_value)) +10), length.out = 21)) +
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
  #facet_wrap(~Ratio_col,scales = "free_x") +
  geom_vline(data=actual_ratio, aes(xintercept =actual_ratio$log2.c.2..10..20..100..,color=Ratio_col), size=1) + 
  geom_hline(data=log10_p_thresholds, aes(yintercept=log10_p_thresholds$X.log10.p_thresholds.,color=ratio_col),size=1,linetype=2)+
  labs(title = "Experiment 2 - DDA no FAIMS processed by MaxQuant")


######



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
  
  if ((sum(stat_analysis[j,2:4])/3) !=  stat_analysis[j,2:4][1] | (sum(stat_analysis[j,5:7])/3) !=  stat_analysis[j,5:7][1]){
    
    p_values_12[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,5:7])
    #p_values_13[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,8:10])
    #p_values_14[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,11:13])
    #p_values_15[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,14:16])
  }else{
    
  }
  if ((sum(stat_analysis[j,2:4])/3) !=  stat_analysis[j,2:4][1] | (sum(stat_analysis[j,8:10])/3) !=  stat_analysis[j,8:10][1]){
    
    p_values_13[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,8:10])
    
  }else{
    
  } 
  if ((sum(stat_analysis[j,2:4])/3) !=  stat_analysis[j,2:4][1] | (sum(stat_analysis[j,11:13])/3) !=  stat_analysis[j,11:13][1]){
    
    
    p_values_14[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,11:13])
    
  }else{
    
  }
  if ((sum(stat_analysis[j,2:4])/3) !=  stat_analysis[j,2:4][1] | (sum(stat_analysis[j,14:16])/3) !=  stat_analysis[j,14:16][1]){
    
    p_values_15[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,14:16])
  }else{
    
  } 
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
  p_values_12[j,"test12"] <- (p_values_12[j,"rank12"]/test_num)*0.1
  p_values_12[j,"test_bool"] <- p_values_12[j,"test12"] > p_values_12[j,"p_values_12"]
  
  
  p_values_13 <- p_values_13[order(p_values_13$p_values_13),]
  colnames(p_values_13)[1] <- "common_col_for_merging"
  p_values_13["rank13"] <- 1:test_num
  p_values_13[j,"test13"] <- (p_values_13[j,"rank13"]/test_num)*0.1
  p_values_13[j,"test_bool"] <- p_values_13[j,"test13"] > p_values_13[j,"p_values_13"]
  
  
  p_values_14 <- p_values_14[order(p_values_14$p_values_14),]
  colnames(p_values_14)[1] <- "common_col_for_merging"
  p_values_14["rank14"] <- 1:test_num
  p_values_14[j,"test14"] <- (p_values_14[j,"rank14"]/test_num)*0.1
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



  
  
  
  
  
  