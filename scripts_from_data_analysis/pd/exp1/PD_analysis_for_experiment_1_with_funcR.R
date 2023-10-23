#library(PhosR)
library(stringr)
library(dplyr)
library(data.table)
library(openxlsx)
library(tidyr)
library(ggplot2)
###############################################
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/ggplot/ggplot_functions.R")
#source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/roc_curve/new_roc_curve_generation_with_custom_threshold.R")

file_path <- "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_1/PD/target_decoy_no_FAIMS_woEcoli/" #target_decoy_withFAIMS_woEcoli
#file_name <-"PAL _Phosphopeptides exp1_wo_Ecoli_noFAIMS_DDA_230208_2023-05-23_1650.xlsx"
theo_file_path= "D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/"
theo_file_name="Synthetic peptides list_theo_conc_corrected_isomericity.xlsx"
sheet_theo_name="ISO-ref and OTHER with FC"

exp_id=1
acquisiton_type="target_decoy_withFAIMS_woEcoli" #"target_decoy_withFAIMS_woEcoli"
software_name="Proteome Discoverer"

final_pd_pep_quant_analysis_exp1 <- function(file_path,
                                                 file_name,
                                                 
                                                 theo_file_path,
                                                 theo_file_name,
                                                 sheet_theo_name,
                                                
                                                 #selected_spcies,
                                                 exp_id,
                                                 #exp_design,
                                                 
                                                 acquisiton_type,
                                                 software_name,
                                                 
                                                 
                                                 
                                                 subtitle){
  
  proline_phospho_pos_extraction <- function(df){
    # Extraction of phospho positions from quant peptides object
    phospho_ptm_pos <- lapply(df, function(each_ptm_protein_positions) {
      
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
    
    return(phospho_ptm_pos_df)
    
  }
  
  all_dirs <- list.files(file_path,pattern = '.txt')
  
  #all_files <- list.files(paste0(file_path,all_dirs[i],"/"),pattern = ".xlsx")
  comb_result <- NULL
  for(i in 1:length(all_dirs)){
    
    #for(j in 1:length(all_files)){
    
      assign(paste0("tmp"), read.table(file=paste0(file_path,all_dirs[i]),sep = "\t",header = T))
      phospho_tmp <- tmp %>% 
        filter(grepl("Phospho",Modifications)) %>%
        mutate(all_dirs[i]) %>%
        mutate(phospho_pos = proline_phospho_pos_extraction(Modifications)) %>%
        mutate(pep_with_pos = paste(Sequence,phospho_pos,sep = "_")) %>%
        distinct(pep_with_pos,.keep_all = TRUE)
      
      assign(paste0(all_dirs[i]), phospho_tmp)
      assign(paste0("comb_result"),bind_rows(comb_result, phospho_tmp))
      rm(tmp,phospho_tmp)
  #}
  
    }
  


  pep_list_w_theo <- read.xlsx(paste0(theo_file_path, theo_file_name), sheet = sheet_theo_name)
  
  pep_list_w_theo_quant <- pep_list_w_theo[,-1]
  common_col_theo_quant <- as.data.frame(paste(pep_list_w_theo_quant$Phosphopeptide.sequence,
                                               pep_list_w_theo_quant$modified.position.in.peptide, sep = "_"))
  colnames(common_col_theo_quant) <- "pep_with_pos"
  pep_list_w_theo_quant_new <- cbind(common_col_theo_quant,pep_list_w_theo_quant)
  
  pep_list_w_theo_unique <- pep_list_w_theo %>% 
    distinct(Phospopeptide.sequence,.keep_all = TRUE) %>%
    rename(Sequence = Phospopeptide.sequence) %>% 
    select(Sequence,Pool) %>%
    rename(Pool_for_seq_merge=Pool)

  comb_result2 <- comb_result %>% select(Sequence, pep_with_pos,
                                         ptmRS.Best.Site.Probabilities,
                                         Spectrum.File,`all_dirs[i]`, First.Scan,
                                         Confidence,Intensity,
                                         Protein.Accessions) %>%
    left_join(pep_list_w_theo_quant_new,by="pep_with_pos") %>%
    drop_na(pool_id)
  
  write.table(comb_result2,file=paste(file_path, software_name,
                                      "experiment",
                                      exp_id,acquisiton_type,
                                      "merge_theo_list_with_identified_phospho_sites.tsv",sep = "_"),
              sep = "\t",col.names = T,row.names = F)
  # comb_result2$new_cat <- factor(comb_result2$Spectrum.File, levels = c("OXPAL221027_01.raw","OXPAL221030_02.raw",
  #                                                                       "OXPAL221030_06.raw","OXPAL221030_10.raw",
  #                                                                       "OXPAL221027_11.raw","OXPAL221030_14.raw",
  #                                                                       "OXPAL221030_18.raw","OXPAL221030_22.raw",
  #                                                                       "OXPAL221030_26.raw"))
  # 
  
  comb_result_seq <- comb_result %>% distinct(Sequence,.keep_all = TRUE) %>%
    select(Sequence, pep_with_pos,ptmRS.Best.Site.Probabilities,
                                            Spectrum.File,`all_dirs[i]`, First.Scan,
                                            Confidence,Marked.as,Intensity,
                                            Protein.Accessions) %>%
    
    left_join(pep_list_w_theo_unique,by="Sequence")
  
  ggplot(comb_result2, aes(x=Spectrum.File,fill=pool_id)) +
    geom_bar() + theme_bw() +
    geom_text(color="black",aes(label=after_stat(count)),
              show.legend = F,stat = "count",size=8,
              position = position_stack(vjust = 0.5)) + 
    scale_fill_brewer(palette = "Dark2") +
    labs(title=paste("Number of Identified phospho-sites for each pool \n","Experiment", exp_id, acquisiton_type))
  
  
  all_seq <- quant_peptides_cor_abun %>% 
    filter(grepl("Phospho",modifications) & grepl(selected_spcies,accession)) %>%
    #distinct(sequence, .keep_all = T) %>%
    separate(accession, into = c("protein","species"),remove = F,sep="_") %>%
    full_join(pep_list_w_theo_unique,by="sequence") %>%
    mutate_at("Pool_for_seq_merge", ~replace_na(.,"Unexpected")) %>%
    mutate(Pool_for_seq_merge= ifelse(is.na(species),"missing",Pool_for_seq_merge)) %>%
    filter(!grepl("Unexpected",Pool_for_seq_merge)) %>%
    bind_rows(ecoli_seq) %>% 
    mutate(Pool_for_seq_merge= ifelse(is.na(Pool_for_seq_merge),background_species,Pool_for_seq_merge))
  #mutate(acq_type=acquisiton_type) %>%
  #mutate(soft_name=software_name)
  
  all_seq_syn <- all_seq %>%
    select(sequence, modifications, Pool_for_seq_merge) %>%
    distinct(sequence, .keep_all = T) %>%
    bind_rows(ecoli_seq_dist) %>%
    mutate(Pool_for_seq_merge= ifelse(is.na(Pool_for_seq_merge),background_species,Pool_for_seq_merge)) %>%
    mutate(acq_type=acquisiton_type) %>%
    mutate(soft_name=software_name)
  
  p13 <- gg_barplt_id_pep_count(data_set = all_seq_syn,
                                x_df = all_seq_syn$Pool_for_seq_merge,
                                fill_df = all_seq_syn$Pool_for_seq_merge,
                                ymax = 20000,
                                header = paste("Total number of identified phosphorylated", selected_spcies,"and", background_species,"across each sample",sep=" "),
                                caption_lab = "NA values are removed.",
                                x_lab = "Sample id",
                                fill_lab =  "Sample id",
                                y_lab = "Number of identified peptides",
                                subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
  write.table(all_seq_syn, file=paste0(file_path,"Experiment2",software_name,"_number_of_unique_sequence_for_each_species.txt"),sep = "\t",col.names = T,row.names = F)
  #################################################
  
  
  
}  
  
  
# Creation of common column merging peptide sequence and phospho positions -> experimental data
common_col_exp <- as.data.frame(paste(id_syn_phospho_pep$sequence, phospho_ptm_pos_df, sep = "_"))
colnames(common_col_exp) <- "common_col_for_merging"

# Bind it to the quant data
syn_phospho_pep_new <- cbind(common_col_exp,id_syn_phospho_pep)


#### THEORETICAL SEQUENCE DATA ####
sheet_theo_name <- "final_pep_list_with_pos"
theo_path = "D:/dev/Desktop_copy/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/"
theo_file = "pep_list_ordered_in_pool_id_with_pos.xlsx"
pep_list_w_theo <- read.xlsx(paste0(theo_path, theo_file), sheet = sheet_theo_name)
#pep_list_w_theo <- pep_list_w_theo[,-1]


common_col_theo <- as.data.frame(paste(pep_list_w_theo$Sequence,
                                       pep_list_w_theo$modified.position.in.peptide, sep = "_"))
colnames(common_col_theo) <- "common_col_merging"
pep_list_w_theo_new <- cbind(common_col_theo,pep_list_w_theo)
colnames(pep_list_w_theo_new)[3] <- "sequence"

#####
#
final_results_with_common_col <- mutate(result_with_common_col,id_syn_phospho_pep)
#
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
  assign(paste0(raw_files_order[i,"pool_id"],query), syn_phospho_pep_new %>% filter(grepl(query,spectrum_title)))
  
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
    full_join(get(paste0("pool",j,"_all")),by="sequence") ## TODO:TRY left_join() without filtering duplicates
  
  # ## REMOVE MULTIPLE PHOSPHOSITE THAT ARE RELATED TO THE SAME PEPTIDE
  assign(paste0("pool",j,"merged"),tmp %>%
           filter(duplicated(sequence) == FALSE))
  
  ## UNEXPECTEDLY IDENTIFIED - FALSE POSITIVES
  assign(paste0("unexp_id_pool",j),get(paste0("pool",j,"merged")) %>% 
           filter(is.na(spectrum_title) == FALSE & is.na(Well.position) == TRUE) %>%
           mutate(type=paste0("unexpected_pool",j)))
  
  ## CORRECTLY IDENTIFIED - TRUE POSTIVIES
  assign(paste0("correct_id_pool",j),get(paste0("pool",j,"merged")) %>% 
           filter(is.na(spectrum_title) == FALSE & is.na(Well.position) == FALSE)%>%
           mutate(type=paste0("correct_pool",j)))
  
  ### MISSED IDENTIFIED - FALSE NEGATIVES
  assign(paste0("missed_id_pool",j),get(paste0("pool",j,"merged")) %>% 
           filter(is.na(spectrum_title) == TRUE & is.na(Well.position) == FALSE)%>%
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


output_path <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_1/Proline/exp1_wo_Ecoli_noFAIMS/"
output_file <- "modified_processable_version_proline_output.txt"
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
  labs(title="Number of identified peptides using Proline") + 
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








