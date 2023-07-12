library(stringr)
library(dplyr)
library(data.table)
#library(openxlsx)
library(ggplot2)
library(reshape2)
library(readr)
library(tibble)
#library(purrr)
library(readr)
library(tidyr)
###############################################
## MAXQUANT DATA ANALYSIS

source("D:/dev/Desktop_copy/PHD/data_analysis/scripts/getModificationPositionMQ_func.R")

setwd("D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/MQ_data_analysis/exp2_wo_FAIMS/txt")

exp2_wo_faims <- read_tsv("evidence.txt")

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
  mutate(common_col_merging = paste(value.y,value.x, sep = "_")) %>%
  select(!value.y) %>%
  rename("mods_id" = "name.x",
         "phospho_positions" ="value.x",
         "all_mods_with_mod_type" = "mods.x")


#getModificationPosition_MQ(mod_seq = id_syn_pep$`Modified sequence`,F)

final_results_with_common_col <- mutate(result_with_common_col,id_syn_phospho_pep)

#filtered_abundances<-final_results_with_common_col[rowSums(!is.na(select(quant_peptides_with_all,starts_with("abundance_"))))>0,]

intensity_match <- final_results_with_common_col %>%
  select(common_col_merging,Intensity,Experiment) %>% 
  tibble() %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = "Experiment",
              values_from = "Intensity")





pep_list_with_pool_id <- read.xlsx(paste0("D:/dev/Desktop_copy/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/",
                                          "pep_list_ordered_in_pool_id_with_pos.xlsx"), sheet = "final_pep_list_with_pos")

pep_list_with_pool_id[,"common_col_merging"] <- paste(pep_list_with_pool_id$sequence,
                                                      pep_list_with_pool_id$modified.position.in.peptide,
                                                      sep="_")

  
  
  
  
  
  
  
  