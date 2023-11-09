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
library(purrr)
###############################################
## MAXQUANT DATA ANALYSIS

### Source code was taken from here: https://rdrr.io/github/singjc/mstools/src/R/getModificationPosition.R
## The code was modified based on what I want and based on software input tyoe (DIANN-and MaxQuant)
#source("D:/dev/Desktop_copy/PHD/data_analysis/scripts/getModificationPositionMQ_func.R")
#source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/get_modification_func/getModificationPositionMQ_func_edit.R")
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/get_modification_func/getModificationPosition_general_change_condition_current_mod_sequence_MQ_Spectronaut.R")
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/ggplot/ggplot_functions.R")

file_path = "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_1/MaxQuant/MQ_v2.1.4_default_params/"
selected_species = "HUMAN"
background_species ="ECOLI"
all_dirs <- list.files(file_path)
acquisiton_type <- c("DDA_with_FAIMS","DDA_with_FAIMS", "DDA_no_FAIMS","DDA_no_FAIMS")
software_name <-"MaxQuant"
exp_id<-1
#all_files <- list.files(paste0(file_path,all_dirs[i],"/"),pattern = ".xlsx")

#comb_result <- NULL
#comb_ecoli <- NULL
for(i in 1:length(all_dirs)){
  
  assign(paste0("tmp"), read.table(file=paste0(file_path,all_dirs[i],"/evidence.txt"),sep = "\t",header = T))
  
  phospho_tmp <- tmp %>% 
    filter(grepl("Phospho",Modified.sequence) & grepl(selected_species,Proteins)) %>%
    drop_na(Modified.sequence) %>%
    mutate(all_dirs[i]) %>%
    mutate(acq_type=acquisiton_type[i]) %>%
    mutate(soft_name=software_name)
    #mutate(phospho_pos = proline_phospho_pos_extraction(Modifications)) %>%
    #mutate(pep_with_pos = paste(Sequence,phospho_pos,sep = "_")) %>%
    #distinct(pep_with_pos,.keep_all = TRUE)
  
  #assign(paste0(all_dirs[i]), phospho_tmp)
  #assign(paste0("comb_result"),bind_rows(comb_result, phospho_tmp))
  
  if(background_species == "Escherichia coli" |background_species == "ECOLI"){
    
    ecoli_tmp <- tmp %>% 
      filter(!grepl(selected_species,Proteins) & !grepl("CON_", Proteins)) %>%
      mutate(Experiment = ifelse(Raw.file == "OXPAL230421_08_-45","E1-M2-coli-inj1",Experiment)
      ) %>%
      mutate(all_dirs[i]) %>%
      mutate(acq_type=acquisiton_type[i]) %>%
      mutate(soft_name=software_name) %>%
      group_by(Sequence,Experiment) %>%
      slice(which.max(Intensity)) %>%
      ungroup() %>%
      #distinct(Sequence,.keep_all = TRUE) %>%
      separate(Experiment,into =c("Exp_id","Sample_id","tmp","inj_id"),sep = "-",remove = F) %>%
      mutate(new_col=paste(Exp_id,Sample_id,sep = "_"))
      
    ecoli_seq_dist <- quant_peptides %>% 
      select(Sequence,Modifications,Proteins) %>%
      filter(!grepl(selected_spcies, Proteins) & !grepl("CON__", Proteins)) %>%
      
      #filter(grepl(background_species,Proteins)) %>%
      distinct(Sequence,.keep_all = T) %>%
      mutate(species=background_species) 
    
    
    
    plot6 <- gg_barplt_id_pep_count(data_set = ecoli_tmp,
                                    x_df = ecoli_tmp$Raw.file,
                                    fill_df = ecoli_tmp$new_col,
                                    ymax = 20000,
                                    header = paste("Total number of identified synthetic phospho-sequences across each pool",sep=" "),
                                    caption_lab = "Wrong Sequences were removed for the futher analysis.",
                                    x_lab = "Sample id",
                                    fill_lab =  "Sample id",
                                    y_lab = "Number of identified Sequence",
                                    subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) + 
      scale_fill_brewer(palette = "Dark2") + theme(axis.text.x = element_text(angle = 90))
    
    write.table(comb_ecoli,file=paste(file_path, software_name,
                                      "experiment",
                                      exp_id,acquisiton_type,
                                      "merge_identified_ecoli_sequences.tsv",sep = "_"),
                sep = "\t",col.names = T,row.names = F)
    
    
    #assign(paste0("comb_ecoli"), bind_rows(comb_ecoli,ecoli_tmp))
  }else{
    
  }
  
  df2 <- apply(X = as.data.frame(phospho_tmp[,"Modified.sequence"]),1,function(x){getModificationPosition_general(mod_seq = x,software_name = "MQ_214")})
  
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
  
  final_results_with_common_col <- mutate(result_with_common_col,phospho_tmp)
  
  
  rm(tmp,phospho_tmp,ecoli_tmp)
  
}
exp1_wo_faims_wo_Ecoli <- read_tsv("evidence.txt")
#exp1_wo_faims_with_Ecoli <- read_tsv("evidence.txt")

# ## NECESSARY FOR COUNTING NUM OF IDENTIFIED ECOLIPEPTIDES
# id_ecoli_pep <- exp2_wo_faims %>% 
#   filter(grepl(83333,`Taxonomy IDs`)) %>% 
#   group_by(Sequence) %>% 
#   summarise(num_id_ecoli_pep = n(), mods=unique(Modifications))


## NECESSARY FOR COUNTING NUM OF IDENTIFIED PHOSPHO PEPTIDES
id_syn_phospho_pep <- exp1_wo_faims_wo_Ecoli %>% 
  filter(grepl(9606,`Taxonomy IDs`)) %>% 
  filter(grepl("Phospho",Modifications)) 



#filtered_abundances<-final_results_with_common_col[rowSums(!is.na(select(quant_peptides_with_all,starts_with("abundance_"))))>0,]



sheet_theo_name <- "final_pep_list_with_pos"
theo_path = "D:/dev/Desktop_copy/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/"
theo_file = "pep_list_ordered_in_pool_id_with_pos.xlsx"
pep_list_w_theo <- read.xlsx(paste0(theo_path, theo_file), sheet = sheet_theo_name)
#pep_list_w_theo <- pep_list_w_theo[,-1]


common_col_theo <- as.data.frame(paste(pep_list_w_theo$Sequence,
                                             pep_list_w_theo$modified.position.in.peptide, sep = "_"))
colnames(common_col_theo) <- "common_col_merging"
pep_list_w_theo_new <- cbind(common_col_theo,pep_list_w_theo)

#####
# Combine data frame (we will continue with this for further step)
df_merge_phosphosite <- final_results_with_common_col %>% 
  left_join(pep_list_w_theo_new,by="common_col_merging")
#####

mapping_path <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_1/MaxQuant/"
mapping_file <- "exp1_wo_FAIMS_wo_Ecoli_mapping.txt"
## MAPPING_PATH
raw_files_order <- as.data.frame(read_tsv(paste0(mapping_path,mapping_file)))


### EXPERIMENTAL DATA SORTING BASED ON POOL ID
for (i in 1:dim(raw_files_order)[1]){
  
  query <- paste0("_", strsplit(raw_files_order[i,"raw_file"], "_")[[1]][2])
  assign(paste0(raw_files_order[i,"pool_id"],query), final_results_with_common_col %>% filter(grepl(query,`Raw file`)))
  
  rm(query)
  
}
pool_size <- 8
## TO MERGE RAWS FILES ARE ASSOCIATED TO THE SAME POOL_ID
for (j in 1:pool_size){
  
  object_list <- mget(ls(pattern=paste0("pool",j)))
  assign(paste0("pool",j,"_all"), do.call(rbind, object_list))
  
  rm(object_list)
  
  # new_col_name <- "common_col_merging"
  # pool_j_all <- get(paste0("pool", j, "_all"))
  # pool_j_all[[new_col_name]] <- paste(pool_j_all$Stripped.Sequence, pool_j_all$value, sep = "_")
  #assign(paste0("pool", j, "_all"), pool_j_all)
  
  ## WITHOUT INSIDE FOR LOOP
  #pool1_all[,"common_col_merging"] <- paste(pool1_all$Stripped.Sequence,
  #pool1_all$value, sep="_")
  
  tmp <- pep_list_w_theo_new %>%
    filter(grepl(paste0("pool",j), pool_id)) %>%
    ## DON'T WORRY ABOUT WARNING MESSAGE,
    ## DUPLICATES WILL BE ELIMINATED IN THE NEXT STEP
    full_join(get(paste0("pool",j,"_all")),by="Sequence") ## TODO:TRY left_join() without filtering duplicates
  
  # ## REMOVE MULTIPLE PHOSPHOSITE THAT ARE RELATED TO THE SAME PEPTIDE
  assign(paste0("pool",j,"merged"),tmp %>%
           filter(duplicated(Sequence) == FALSE))
  
  ## UNEXPECTEDLY IDENTIFIED - FALSE POSITIVES
  assign(paste0("unexp_id_pool",j),get(paste0("pool",j,"merged")) %>% 
           filter(is.na(`Raw file`) == FALSE & is.na(Well.position) == TRUE) %>%
           mutate(type=paste0("unexpected_pool",j)))
  
  ## CORRECTLY IDENTIFIED - TRUE POSTIVIES
  assign(paste0("correct_id_pool",j),get(paste0("pool",j,"merged")) %>% 
           filter(is.na(`Raw file`) == FALSE & is.na(Well.position) == FALSE)%>%
           mutate(type=paste0("correct_pool",j)))
  
  ### MISSED IDENTIFIED - FALSE NEGATIVES
  assign(paste0("missed_id_pool",j),get(paste0("pool",j,"merged")) %>% 
           filter(is.na(`Raw file`) == TRUE & is.na(Well.position) == FALSE)%>%
           mutate(type=paste0("missed_pool",j)))
  
}

type_of_id <- c("correct_id_pool","missed_id_pool","unexp_id_pool")  

for (k in 1:length(type_of_id)){
  
  assign(paste0("final_",type_of_id[k],"list"),mget(ls(pattern=paste0(type_of_id[k]))))
  assign(paste0("ffinal_",type_of_id[k]), do.call(rbind, get(paste0("final_",type_of_id[k],"list"))))
  #assign(paste0("ffinal_",type_of_id[k]),cbind(get(paste0("ffinal_",type_of_id[k])),paste0(type_of_id[k])))
  
}


final_all_id_list <- mget(ls(pattern=paste0("ffinal")))
final_all_id <- do.call(rbind, final_all_id_list)


output_path <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_1/MaxQuant/exp1_wo_FAIMS_wo_Ecoli_OXPAL230121/modified_version_of_evidence/"
output_file <- "modified_processable_version_evidence.txt"
write_tsv(final_all_id,paste0(output_path,output_file))

library(viridis)
library(ggplot2)
final_all_id  %>% 
  group_by(type) %>% 
  summarise(n()) %>%
  separate(type,c('type_name','pool_id')) %>%
  ggplot( aes(fill=type_name,x=pool_id,y=`n()`)) + geom_bar(stat="identity",position = "dodge")+
  geom_text(aes(label=`n()`),colour = "gray", size = 4,position = position_dodge(0.9),vjust=1.6) +
  theme_minimal() +
  theme(legend.text = element_text(size=15), 
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size=20),
        legend.title=element_text(size=15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15)) + 
  scale_fill_viridis(discrete = TRUE,option = "viridis")







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








