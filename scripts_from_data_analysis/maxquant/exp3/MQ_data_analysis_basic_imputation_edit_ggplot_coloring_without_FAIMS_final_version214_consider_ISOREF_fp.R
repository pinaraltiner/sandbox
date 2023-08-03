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
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/get_modification_func/getModificationPositionMQ_func_edit.R")
#source("D:/dev/Desktop_copy/PHD/data_analysis/scripts/getModificationPositionMQ_func_change_condition_current_mod_sequence.R")
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/roc_curve/roc_curve_generation_proline_edit.R")
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/ggplot/ggplot_functions.R")

file_path <- "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_3/MQ_214_optimize_param_noFAIMS/"
#"D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/MQ_214_optimize_param_noFAIMS/"
#file_path <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/MQ_data_analysis/with_FAIMS_MBR/min_mod_score_0/"
#file_path <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/MQ_data_analysis/exp2_wo_FAIMS/rerun_using_1634/"
output_path <- "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_3/MQ_214_optimize_param_noFAIMS/figures/"
setwd(file_path)
exp3_wo_faims <- read_tsv("evidence.txt")

# experiment_name <- c("E2-M1-R1",
#                      "E2-M1-R2",
#                      "E2-M1-R3",
#                      "E2-M2-R1",
#                      "E2-M2-R2",
#                      "E2-M2-R3",
#                      "E2-M3-R1",
#                      "E2-M3-R2",
#                      "E2-M3-R3",
#                      "E2-M4-R1",
#                      "E2-M4-R2",
#                      "E2-M4-R3",
#                      "E2-M5-R1",
#                      "E2-M5-R2",
#                      "E2-M5-R3")

exp_design <- c("E3-A1-R1","E3-A1-R2",
                "E3-A1-R3",
                "E3-A2-R1",
                "E3-A2-R2",
                "E3-A2-R3",
                "E3-A3-R1",
                "E3-A3-R2",
                "E3-A3-R3",
                "E3-A4-R1",
                "E3-A4-R2",
                "E3-A4-R3",
                "E3-A5-R1",
                "E3-A5-R2",
                "E3-A5-R3")



num_reps <- 3
sample_size <- length(exp_design) / num_reps
sample_names <- paste0("A",1:sample_size)
comparisons <- NULL
for (i in 1:sample_size){
  tmp <- paste0(sample_names[1], "_vs_",sample_names[i])
  comparisons[i] <- tmp
  rm(tmp)
}

acquisiton_type <- c("DDA Exploris with FAIMS","DDA Exploris no FAIMS","DDA TIMS-TOF")
exp_id <- 1:3
software_name <- c("Proline", "MaxQuant", "PD")

### IMPUTATION

# Calculate 1 percent quantile of each sample

imputed_values <- exp3_wo_faims %>% 
  group_by(Experiment) %>% 
  summarise(first_quantile=quantile(Intensity,probs=0.01,na.rm=TRUE))
imputed_values_vec <- as.vector(imputed_values$first_quantile)
# Impute missing values


  
# abundances_for_impute <- exp3_wo_faims %>% select(Experiment,Intensity) %>%
#   group_by(Experiment) %>%
#   mutate(imputed_intensity=case_when(grepl("A1-R1",Experiment) ~imputed_values_vec[1],
#                                        grepl("A1-R2",Experiment) ~imputed_values_vec[2],
#                                        grepl("A1-R3",Experiment) ~imputed_values_vec[3],
#                                        grepl("A2-R1",Experiment) ~imputed_values_vec[4],
#                                        grepl("A2-R2",Experiment) ~imputed_values_vec[5],
#                                        grepl("A2-R3",Experiment) ~imputed_values_vec[6],
#                                        grepl("A3-R1",Experiment) ~imputed_values_vec[7],
#                                        grepl("A3-R2",Experiment) ~imputed_values_vec[8],
#                                        grepl("A3-R3",Experiment) ~imputed_values_vec[9],
#                                        grepl("A4-R1",Experiment) ~imputed_values_vec[10],
#                                        grepl("A4-R2",Experiment) ~imputed_values_vec[11],
#                                        grepl("A4-R3",Experiment) ~imputed_values_vec[12],
#                                        grepl("A5-R1",Experiment) ~imputed_values_vec[13],
#                                        grepl("A5-R2",Experiment) ~imputed_values_vec[14],
#                                        grepl("A5-R3",Experiment) ~imputed_values_vec[15],
#                                        TRUE ~ Intensity))
  

###### CHANGE THE NAME OF THE OBJECT LATER ### THERE WAS NO IMPUTATION HERE ######
  phospho_mouse <- exp3_wo_faims %>%
  #select(Experiment,Intensity) %>% 
  #mutate(abundances_for_impute) %>% 
  filter(grepl(10090,`Taxonomy IDs`))%>%
  #filter(grepl("_MOUSE",Proteins)) %>% 
  filter(grepl("Phospho",Modifications)) %>%
  mutate(values = paste(id,`Protein group IDs`,sep="@"))

#id_syn_phospho_pep <- phospho_mouse

#write.table(id_syn_phospho_pep,file="exp2_no_faims_MQ_evidences_filtered_human_phospho.tsv",sep="\t",row.names = FALSE)

# ## NECESSARY FOR COUNTING NUM OF IDENTIFIED ECOLIPEPTIDES
# id_ecoli_pep <- exp2_wo_faims %>% 
#   filter(grepl(83333,`Taxonomy IDs`)) %>% 
#   group_by(Sequence) %>% 
#   summarise(num_id_ecoli_pep = n(), mods=unique(Modifications))

###### WARNING: ALWAYS SELECT MODIFIED PEPTIDES OTHERWISE THE FINAL RESULT OF THE 
    ## FUNCTION WILL BE THE SAME NUMBER OF LINES.

df2 <- apply(as.data.frame(phospho_mouse[,"Modified sequence"]),1,getModificationPosition_MQ)

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

final_results_with_common_col <- mutate(result_with_common_col,phospho_mouse)

### ELIMINATION OF MULTI-CHARGE PEPTIDES (ask for whether necessary for Exp3 or not) ###
### IT IS DEFINETELY NECESSARY FOR PIVOT_WIDER()

barplt_df <- final_results_with_common_col %>%
  select(Sequence,Experiment,Intensity,common_col_for_merging,`Taxonomy IDs`, Modifications,id, `Protein group IDs`) %>%
  separate(Experiment, into = c("Exp_id","Sample_id","Rep_id"),sep = "-",remove = F) %>%
  mutate(sample_rep_id_seq = paste(common_col_for_merging, Sample_id,Rep_id, sep = "_"))%>%
  group_by(sample_rep_id_seq,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
  slice(which.max(Intensity)) %>%
  ungroup()


barplt_df_ecoli <- exp3_wo_faims %>% 
  filter(grepl(83333,`Taxonomy IDs`)) %>%
  select(Sequence,Experiment,Intensity,`Taxonomy IDs`, Modifications,id, `Protein group IDs`) %>%
  separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F) %>%
  mutate(sample_rep_id_seq = paste(Sequence, Sample_id,Rep_id, sep = "_")) %>%
  group_by(sample_rep_id_seq,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
  slice(which.max(Intensity)) %>%
  ungroup()




p1 <- gg_barplt_id_pep_count(data_set = barplt_df,
                             x_df = barplt_df$Sample_id,
                             fill_df = barplt_df$Rep_id,
                             ymax = 5000,
                             header = "Total number of quantified phospho-site across each sample",
                             caption_lab = "NA values are removed.",
                             x_lab = "Sample id",
                             fill_lab =  "Sample id",
                             y_lab = "Number of identified peptides",
                             subtitle_txt = "")




  

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

## When Sequence was used for grouping, 2169 more rows were deleted. 
## That's why I kept (common_col_for_merging) peptides with phospho positions.
# tmp1  <-final_results_with_common_col %>% 
#   group_by(common_col_for_merging,Experiment) %>%
#   slice(which.max(Intensity)) %>%
#   ungroup()

#### df_reversible is a reversible object to retrieve any features from the very beginning,
  ## It can be done by using this object thanks to unique_row_name column (in the barplot_df)
#### ONLY THING THAT SHOULD BE DONE IS TO USE PIVOT_LONGER() WITH "df_reversible" THEN USE left_join()
  ## unique_row_name is the combination of id and protein group id column from the evidence.txt
################################################################################
df_reversible <- barplt_df %>% 
  select(common_col_for_merging,Experiment, Intensity,id, `Protein group IDs`) %>%
  mutate(unique_row_name = paste(id,`Protein group IDs`,sep="@")) %>%
  select(-id,-`Protein group IDs`) %>%
  pivot_wider(names_from = Experiment,values_from = c(unique_row_name,Intensity)) %>%
  mutate(species = "MOUSE")

df_wide <- df_reversible %>% 
    select(common_col_for_merging, species,contains("Intensity"))
    
df_reversible_ecoli <- barplt_df_ecoli %>% 
  select(Sequence,Experiment, Intensity,id, `Protein group IDs`) %>%
  mutate(unique_row_name = paste(id,`Protein group IDs`,sep="@")) %>%
  select(-id,-`Protein group IDs`) %>%
  pivot_wider(names_from = Experiment,values_from = c(unique_row_name,Intensity)) %>%
  rename_with(~ paste0("common_col_for_merging"), matches("^Seq")) %>%
  mutate(species = "ECOLI")

df_wide_ecoli <- barplt_df_ecoli %>%
  select(Sequence,Experiment,Intensity) %>%
  pivot_wider(names_from = Experiment,values_from = Intensity) %>%
  rename_with(~ paste0("common_col_for_merging"), matches("^Seq")) %>%
  mutate(species = "ECOLI") %>%
  relocate(all_of(exp_design), .after = where(is.character)) %>%
  rename_with(~ paste0("Intensity_", .x), matches(exp_design))

#### EXAMPLE OF HOW TO MERGE DATA AFTER IMPUTATION ####

## final_imputated_data contains all abundance variety that we need
   ## abundance, log10(abundance), rowMeans(abundance),
## However, merging should be done one by one.
## This is the example of extraction of log10 abundances.

test <- final_imputed_data %>% filter(grepl("MOUSE",species)) %>%
  select(1:2,contains("log10")) %>% 
  pivot_longer(cols = starts_with("log10"),names_to = "Experiment",values_to = "log10_Intensity") 

tmp <- df_reversible %>% select(1:2,contains("unique")) %>%
  pivot_longer(cols = starts_with("unique"),names_to = "unique_ids",values_to = "values") %>%
  bind_cols(test) %>%
  left_join(phospho_mouse,by="values")

####################################################

filtered_abundances<-df_wide[rowSums(!is.na(select(df_wide,contains(exp_design))))>0,]
filtered_abundances_ecoli<-df_wide_ecoli[rowSums(!is.na(select(df_wide_ecoli,contains(exp_design))))>0,]

abundances_rowMeans <- NULL
abundances_ecoli_rowMeans <- NULL

for (k in 1:sample_size){
  # If separate version of row means is not needed, it can be commented later.
  # Separate row Means can be collected in temp object to merge in "log_10_filtered_abundances_rowMeans"
  assign(paste0("abundances_A",k),as.data.frame(rowMeans(filtered_abundances  %>%
                                                           select(contains(paste0("-A",k)))))) #%>%
  #select(starts_with("abundance_")))))
  abundances_rowMeans<- bind_cols(abundances_rowMeans,get(paste0("abundances_A",k)))
  
  
  assign(paste0("ecoli_abundances_A",k),as.data.frame(rowMeans(filtered_abundances_ecoli  %>%
                                                                 select(contains(paste0("-A",k))))))# %>%
  #select(starts_with("abundance_")))))  
  
  abundances_ecoli_rowMeans<- bind_cols(abundances_ecoli_rowMeans,get(paste0("ecoli_abundances_A",k)))
  
}

quant_peptides_ECOLI_density_plot <- filtered_abundances_ecoli %>%
  select(!starts_with("Intensity")) %>%
  bind_cols(abundances_ecoli_rowMeans) %>%
  tibble() %>%
  rename_with(~ paste0("mean_abun",1:5), matches("^row")) %>%
  pivot_longer(cols = starts_with("mean"), 
               values_to = "intensity",
               names_to = "sample_ids",
               values_drop_na = T)

ecoli_density_plot<- quant_peptides_ECOLI_density_plot %>%
  select(contains(c("sample_ids","intensity","species"))) 


colnames(abundances_rowMeans) <- paste0("mean_abun",1:sample_size)
#### MEAN ABUNDANCE RATIO WITH  DENSITY PLOT ####
### BEFORE IMPUTATION ###
quant_phospho_density_plot <- filtered_abundances %>%
  select(!starts_with("Intensity")) %>%
  bind_cols(abundances_rowMeans) %>% 
  #rename_with(~ paste0("mean_abun",1:5), matches("^row")) %>%
  tibble() %>% #mutate(pep_with_pos = sequence) %>% ###  At this stage, no need for phospho-position
               # We only count number of sequence here.   
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

### Counting NAs across either rows or columns can be achieved by using the apply() function.
 ## This function takes three arguments:
  # X is the input matrix, 
  # MARGIN is an integer, => ## MARGIN = 1 means to apply the function across rows and MARGIN = 2 across columns.
  # FUN is the function to apply to each row or column. 
library(kableExtra)
na_phospho_mouse <- apply(X = is.na(filtered_abundances), MARGIN = 2, FUN = sum)
na_ecoli <- apply(X = is.na(filtered_abundances_ecoli), MARGIN = 2, FUN = sum)

na_table <- bind_rows(na_phospho_mouse,na_ecoli)
na_table$species <- c("MOUSE","ECOLI")
na_table$total <- c(dim(filtered_abundances)[1],dim(filtered_abundances_ecoli)[1])

p4 <- na_table %>% select(-common_col_for_merging) %>%
  kbl(caption = "Number of NA values across all samples") %>%
  kable_material(c("striped", "hover")) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left", font_size = 12) %>%
  kable_minimal(full_width = F) %>%
  footnote(general = "This table was created after elimination of multiple charges by selecting either phospho-sites of mouse or background sequences \n that has the highest abundace.",
           # number = c("Footnote 1; ", "Footnote 2; "),
           # alphabet = c("Footnote A; ", "Footnote B; "),
           # symbol = c("Footnote Symbol 1; ", "Footnote Symbol 2")
           footnote_as_chunk = T, title_format = c("italic", "underline")
  )

for (j in 1:length(imputed_values_vec)){
  # Number NA
  #num_NA <- length(filtered_abundances[,j+1][is.na(filtered_abundances[,j+1])])
  
  filtered_abundances[,j+2][is.na(filtered_abundances[,j+2])] <- imputed_values_vec[j]
  filtered_abundances_ecoli[,j+2][is.na(filtered_abundances_ecoli[,j+2])] <- imputed_values_vec[j]
  
  # After imputation number of imputed values
  #num_imp <-length(filtered_abundances[,j+1][(filtered_abundances[,j+1]==imputed_values_vec[j])])
  
  # This is verification of imputation is done successfully
  # Because we expect to see that number of imputed values should be the same amount as number of NA
  #print(setequal(num_NA,num_imp))
  #print(num_NA)
  #print(num_imp)
}


#  ## CONDITIONAL IMPUTATION ##
# 
# #######################################
# triplicate_indx <- c(1,3,4,6,7,9,10,12,13,15)
# impute_values_2NA <- NULL
# impute_values_1NA <- NULL
# tmp_2NA <- NULL
# tmp_1NA <- NULL
# #imputed_abundances <- NULL
# for (i in 1:5){
#   tmp_1NA <- quantile (abundances_for_impute[,experiment_name[triplicate_indx[1]:triplicate_indx[2]]],
#                        probs = 0.01 , na.rm = TRUE )
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
# df_long <- df_wide %>%
#   pivot_longer(cols = -common_col_for_merging, names_to = "column", values_to = "value") %>%
#   separate(column, into = c("col_group", "col_index", "row_index"), sep = "-")
# 
# 
# # Count the number of NA values in each row for each column group
# na_counts <- df_long %>%
#   group_by(common_col_for_merging, col_index) %>%
#   summarise(na_count = sum(is.na(value))) 
# #pivot_wider(names_from = col_index, values_from = na_count)
# 
# # Replace the NA values based on the count
# df_result <- df_long %>%
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
#   filter(rowSums(across(where(is.numeric)))!=0) %>%
#   select(order(colnames(df_wide)))

## TEST DATA for imputation part: 
#toy_data <- read_tsv(file = "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/comparision_DDA_noFAIMS_mq_proline_pd/toy_data_MQ.txt")
#######################################



########################

filtered_abundances_rowMeans <- NULL
filtered_abundances_log10 <- NULL
filtered_abundances_log10_rowMeans <- NULL

filtered_all <- bind_rows(filtered_abundances,filtered_abundances_ecoli)

for (k in 1:sample_size){
  # If separate version of row means is not needed, it can be commented later.
  # Separate row Means can be collected in temp object to merge in "log_10_filtered_abundances_rowMeans"
  
  #### COMMENTED CODES ARE CORRESPOND TO IMPUTATION FOR ONLY MOUSE ####
  
  # assign(paste0("aft_imp_abundances_A",k),as.data.frame(rowMeans(abundances_for_impute_mouse  %>% select(contains(paste0("A",k))))))
  # filtered_abundances_rowMeans<- bind_cols(filtered_abundances_rowMeans,get(paste0("aft_imp_abundances_A",k)))
  # 
  # assign(paste0("log10_abundances_A",k),as.data.frame(log10(abundances_for_impute_mouse  %>% select(contains(paste0("A",k)))))) 
  # filtered_abundances_log10 <- bind_cols(filtered_abundances_log10, get(paste0("log10_abundances_A",k)))
  # 
  # assign(paste0("rowMean_log10_abundances_A",k),as.data.frame(rowMeans(get(paste0("log10_abundances_A",k)) %>% select(contains(paste0("A",k))))))
  # filtered_abundances_log10_rowMeans<- bind_cols(filtered_abundances_log10_rowMeans,get(paste0("rowMean_log10_abundances_A",k)))
  
  assign(paste0("aft_imp_abundances_A",k),as.data.frame(rowMeans(filtered_all  %>% select(contains(paste0("A",k))))))
  filtered_abundances_rowMeans<- bind_cols(filtered_abundances_rowMeans,get(paste0("aft_imp_abundances_A",k)))
  
  assign(paste0("log10_abundances_A",k),as.data.frame(log10(filtered_all  %>% select(contains(paste0("A",k))))))
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

final_imputed_data <- cbind(filtered_all, filtered_abundances_rowMeans,filtered_abundances_log10) #filtered_abundances_log10_rowMeans


df_mean_ab_after_impt <- final_imputed_data %>% 
  select(contains("aft_imp"), species) %>%
  tibble() %>% 
  pivot_longer(cols = contains("aft_imp"),
               names_to = "Mean_abundance",
               values_to = "values")

df_FC_ratio_after_impt <- final_imputed_data %>% 
  select(starts_with("exp_"),species) %>%
  tibble() %>% 
  pivot_longer(cols = starts_with("exp_"),
               names_to = "exp_FC",
               values_to = "values")

p5 <- gg_density(data_set = df_mean_ab_after_impt, 
                 x_df = df_mean_ab_after_impt$values,
                 fill_df = df_mean_ab_after_impt$Mean_abundance,
                 color_df = df_mean_ab_after_impt$species,
                 header="Distribution of mean abundance of every sample after imputation",
                 facet_df = "Mean_abundance",
                 x_lab = "log10(values)",
                 color_lab= "",
                 fill_lab = "Sample Names",
                 subtitle_txt = "")

# px <- gg_raincloud(data_set = df_mean_ab_after_impt,
#                    x_df = df_mean_ab_after_impt$Mean_abundance,
#                    y_df = df_mean_ab_after_impt$values,
#                    fill_df = df_mean_ab_after_impt$Mean_abundance,
#                    header = "Distribution of mean abundance of every sample after imputation",
#                    x_lab = "Sample Names",
#                    y_lab = " Density of log10(Mean Abundance)",
#                    fill_lab = "Sample Names",
#                    caption_lab = "",
#                    subtitle_txt = "")



p6 <- gg_density(data_set = df_FC_ratio_after_impt,
                 x_df = df_FC_ratio_after_impt$values,
                 fill_df = df_FC_ratio_after_impt$exp_FC,
                 color_df = df_FC_ratio_after_impt$species,
                 header="Distribution of Fold change Ratio of every sample after imputation",
                 facet_df = "exp_FC",
                 x_lab = "log10(values)",
                 color_lab= "",
                 fill_lab = "Sample Names",
                 subtitle_txt = "")


### BOX-PLOT: Experimental Quantity Ratio of Phospho Peptides  

p7 <- gg_boxplt_exp_ratio(data_set = df_FC_ratio_after_impt, 
                          x_df = df_FC_ratio_after_impt$exp_FC,
                          y_df = df_FC_ratio_after_impt$values,
                          fill_df = df_FC_ratio_after_impt$species,
                          header="Experimental Quantity Ratio of T-cell Phospho Peptides",
                          x_lab="Sample Names",
                          y_lab="Abundance Ratios",
                          fill_lab = "Sample Names",
                          subtitle_txt = "")



### HALF-BOX-PLOT & HALF-SCATTER-PLOT: Experimental Quantity Ratio of Synthetic Peptides  
library(gghalves)

p8 <- gg_half_boxplt_exp_ratio(data_set = df_FC_ratio_after_impt, 
                               x_df = df_FC_ratio_after_impt$exp_FC,
                               y_df = df_FC_ratio_after_impt$values,
                               fill_df = df_FC_ratio_after_impt$species,
                               header="Experimental Quantity Ratio of T-cell Phospho Peptides with Background",
                               x_lab="Sample Names",
                               y_lab="Abundance Ratios",
                               fill_lab = "Sample Names",
                               subtitle_txt = "")

### VIOLIN-PLOT: Experimental Quantity Ratio of Synthetic Peptides   

### TODO: fix y scaling without trimming 
p9 <- gg_violin_exp_ratio(data_set = df_FC_ratio_after_impt, 
                          x_df = df_FC_ratio_after_impt$exp_FC,
                          y_df = df_FC_ratio_after_impt$values,
                          fill_df = df_FC_ratio_after_impt$species,
                          header="Experimental Quantity Ratio of T-cell Phospho Peptides with Background",
                          x_lab="Sample Names",
                          y_lab="Abundance Ratios",
                          fill_lab = "Sample Names",
                          trim=TRUE,
                          subtitle_txt = "")

# test <- final_imputed_data %>% 
#   select(common_col_for_merging, contains(exp_design) | starts_with("log10") | starts_with("exp_FC")) %>%
#   pivot_longer(cols = starts_with("E3"), values_to = "abundances", names_to = "Sample_names1" ) %>%
#   pivot_longer(cols = starts_with("log10"), values_to = "mean_log10_abundances", names_to = "Sample_names2" ) %>%
#   pivot_longer(cols = starts_with("exp_FC"), values_to = "row_mean_abundances", names_to = "Sample_names3")
# 
# test <- final_imputed_data %>% 
#   select(common_col_for_merging, starts_with("log10")) %>%
#   pivot_longer(cols = starts_with("log10"), values_to = "mean_log10_abundances", names_to = "Sample_names2" )
# 
# 
# test1 <- final_imputed_data %>% 
#   select(common_col_for_merging, contains(exp_design)) %>%
#   pivot_longer(cols = starts_with("E3"), values_to = "abundances", names_to = "Sample_names1" )
# 
# 
# ### COMBINE IMPUTED INTENSITIES WITH THE REST OF THE DF
# final_results_max_abun <- final_results_with_common_col %>% 
#      group_by(common_col_for_merging,Experiment) %>%
#      slice(which.max(Intensity)) %>%
#      ungroup()
# 
# df_merge <- df_merge %>% 
#   pivot_longer(cols = -common_col_for_merging, 
#                names_to = "Experiment", 
#                values_to = "Imputed_intesity") %>%
#   left_join(final_results_max_abun, 
#             by=c("Experiment","common_col_for_merging"))
# To calculate all binary combination in the data frame
#mat <- do.call(cbind, lapply(cols, function(xj) 
#  sapply(cols, function(xi) (filtered_abundances_rowMeans[, xj]/(filtered_abundances_rowMeans[, xj])))))
#colnames(mat) <-  outer(names(filtered_abundances_rowMeans), names(filtered_abundances_rowMeans), paste0)

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


rownames(final_imputed_data) <- paste0(final_imputed_data$common_col_for_merging,"@",final_imputed_data$species,"@",(1:nrow(final_imputed_data)))
if(test_type== "t.test" | test_type== "wilcoxon"){
  all_pvalues <- NULL
  
  for (i in 2:sample_size){
    p_values_tmp <- NULL
    
    for(j in 1:dim(final_imputed_data)[1]){
      
      if(test_type=="t.test"){
        p_values_tmp[j] <- ttest_func(select(final_imputed_data,contains("A1-") & contains("log10_"))[j,],     ### FOR DIFFERENT KIND OF EXP SETUP, 
                                      select(final_imputed_data,contains(paste0("A",i,"-")) & contains("log10_"))[j,] ) ## It should be defined as an input.
        
      }else if(test_type=="wilcoxon"){
        p_values_tmp[j] <- wilcox.test(select(final_imputed_data,contains("A1-") & contains("log10_"))[j,],     ### FOR DIFFERENT KIND OF EXP SETUP, 
                                       select(final_imputed_data,contains(paste0("A",i,"-")) & contains("log10_"))[j,])     ## It should be defined as an input.
      }
      
      
    }
    p_values_tmp <- as.data.frame(p_values_tmp)
    colnames(p_values_tmp) <- paste0("pvalues_A1","/","A",i)
    all_pvalues <- bind_cols(all_pvalues,p_values_tmp)
    rm(p_values_tmp)
    
  }
  all_pvalues_common_col <- all_pvalues %>% 
    bind_cols(final_imputed_data$common_col_for_merging,final_imputed_data$species) %>%
    pivot_longer(cols = starts_with("pvalues_"), values_to = "pvalues", names_to ="p_ratios") %>%
    separate(p_ratios, into = c("tmp","ratio"),sep = "_") %>%
    select(!tmp) %>%
    rename_with(.col =1 , ~"common_col_for_merging") %>%
    rename_with(.col=2, ~ "species") %>%
    mutate(common_col = paste(common_col_for_merging,species,1:((sample_size-1)*nrow(final_imputed_data)),sep="@"))
  
  
  ## THE BEST WAY TO DO is this:
  merge_stat_df <- final_imputed_data %>%
    select(common_col_for_merging, species, starts_with("exp_FC")) %>%
    pivot_longer(cols = starts_with("exp_FC"), values_to = "fold_change_values", names_to ="fold_change_ratios") %>%
    separate(fold_change_ratios, into = c("tmp","tmp1","ratio"),sep = "_") %>%
    select(!c(tmp,tmp1)) %>%
    mutate(common_col = paste(common_col_for_merging,species,1:((sample_size-1)*nrow(final_imputed_data)),sep="@")) %>%
    bind_cols(all_pvalues_common_col$ratio,all_pvalues_common_col$pvalues) %>%
    rename_with(.col =6 , ~"ratio1") %>%
    rename_with(.col=7, ~ "pvalues")
  
  
  ### MERGING I: All used columns are merged and used to combine pvalues and ratios
  ## Before changing the shape of data 
  ## This eliminates any duplication and mismatching btw dataframes
  ## We are 100% sured that all pvalues are associated with its ratio.
  ## Shape of volcano plot looks quite weird esspecially in the A1vsA5.
  
  ## THIS IS CORRECT BUT UNNECESSARILY LONG WAY
  # merge_stat_df <- final_imputed_data %>%
  #   select(common_col_for_merging, species, starts_with("exp_FC")) %>%
  #   bind_cols(all_pvalues) %>%
  #   pivot_longer(cols = starts_with("exp_FC"), values_to = "fold_change_values", names_to ="fold_change_ratios") %>%
  #   separate(fold_change_ratios, into = c("tmp","tmp1","ratio"),sep = "_") %>%
  #   select(!c(tmp,tmp1)) %>%
  #   mutate(common_col = paste(common_col_for_merging,species,1:((sample_size-1)*nrow(final_imputed_data)),sep="@")) %>%
  #   left_join(all_pvalues_common_col, by="common_col")
  # 
  ### MERGING II: Previous version of combining pvalues and ratios
  ## No errors were appeared when I double pivot_longer()
  ## But there is no common column between those values so that
  ## We cannot know which pvalues correspond to which ratio (MAYBE WE DON'T NEED IT)
  ## Shape of volcano plot looks better when I used this way.
  ## IT IS NOT CORRECT
  # merge_stat_df1 <- final_imputed_data %>%
  #   select(common_col_for_merging, species, starts_with("exp_FC")) %>%
  #   bind_cols(all_pvalues) %>%
  #   pivot_longer(cols = starts_with("exp_FC"), values_to = "fold_change_values", names_to ="fold_change_ratios") %>%
  #   pivot_longer(cols = starts_with("pvalues_"), values_to = "pvalues", names_to ="p_ratios")
  # 
  
  
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
    stat_analysis <- final_imputed_data %>%
      select(common_col_for_merging,species, starts_with("log10_") | starts_with("mean_log10_") | starts_with("exp_FC")) #spectrum_title
    
    
  design_matrix <- model.matrix(~factor(c(rep(2,num_reps),rep(1,num_reps))))
  merge_stat_df <-NULL
  for ( i in 2:sample_size){
    # Change only the colname iteratively makes fit to every comparison
    colnames(design_matrix) <- c("Intercept", paste0("A1-A",i))
    #print(colnames(design_matrix))
    # Col selection for each comparison
    assign(paste0("df_A1vsA",i),stat_analysis %>% select(1:2 | contains("log10_Intensity_E3-A1") | contains(paste0("log10_Intensity_E3-A",i))))
    # First, linear model was built
    assign(paste0("fit",i) ,lmFit(get(paste0("df_A1vsA",i))[,3:8], design_matrix))
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
    rename_with(.col=8, ~ "ratios") %>%
    separate(mult_col, into = c("sequence","accession","id"),sep = "@")
    #separate(accession, into = c("uniprot_id","species"),sep = "_")
  
  
  
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


p10 <- ggplot(merge_stat_df ,aes(x =log2(fold_change_values), y = -log10(pvalues), color=species)) +
  geom_point(size = 1) + #, aes(shape=merge_stat_df_final$species)
  facet_wrap(~ratio) +
  #scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 0.8)) +
  #scale_x_continuous(limits = c(-3,3),breaks = seq(-3, 3, by = 0.8)) +
  scale_color_brewer(palette = "Set1") +
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
  labs(title = paste0("Experiment - ", exp_id[3], acquisiton_type[2], " data processed by ", software_name[2]), subtitle = "T-test was used")
#expand_limits(x = 0, y = 0) +
#geom_vline(data = actual_ratio, aes(xintercept = actual_ratio$X.1....log2.c.2..10..20..100..., size = 1, show.legend = FALSE)) + #color=c("#CC79A7","#E69F00","#56B4E9","#009E73")
#geom_hline(data = log10_p_thresholds, aes(yintercept = log10_p_thresholds$X.log10.p_thresholds.),color=c("#CC79A7","#E69F00","#56B4E9","#009E73"), size = 1, linetype = 2, show.legend = FALSE)+ 
#

# gg_volcano(data_set = merge_stat_df,
#                 x_df = log2(merge_stat_df$fold_change_values),
#                 y_df = -log10(merge_stat_df$pvalues),
#                 color_df = merge_stat_df_final$species,
#                 facet_df = "ratio.x",
#                 header ="", #paste0("Experiment ", exp_id[3], acquisiton_type[2], " data processed by ", software_name[1]),
#                 color_lab = "Classes",
#                 subtitle_txt = "with adjusted p values")

p11 <- gg_volcano(data_set = merge_stat_df_final,
                  x_df = merge_stat_df_final$logFC,
                  y_df = -log10(merge_stat_df_final$P.Value),
                  color_df = merge_stat_df_final$accession,
                  facet_df = "ratios",
                  header = paste0("Experiment - ", exp_id[3], acquisiton_type[2], " data processed by ", software_name[1]), ## Specifications of the header will be asked as an input.
                  color_lab = "Classes",
                  subtitle_txt = "limma was used")    

# ggplot(data=merge_stat_df1,aes(x=log2(merge_stat_df1$fold_change_values),y=-log10(merge_stat_df1$pvalues),color=species)) +
#   geom_point() + facet_wrap(~fold_change_ratios)
# 





sapply(10:11,function(x) ggsave(filename = paste0("p",x,".tiff"),
                                width = 50, height = 40, 
                                path = paste0(file_path,"/outputs_with_new_script/"),
                                units = "cm",
                                get(paste0("p",x)),
                                device = "tiff", #".svg"
))




