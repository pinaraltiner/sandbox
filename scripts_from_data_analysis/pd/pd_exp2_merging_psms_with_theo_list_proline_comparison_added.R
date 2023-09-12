#library(PhosR)
library(stringr)
library(dplyr)
library(data.table)
library(openxlsx)
library(tidyr)
library(ggplot2)
library(tidyverse)


# extract_phospho_numbers <- function(input_string) {
#   phospho_part <- regmatches(input_string, gregexpr("Phospho \\[[^]]+\\]", input_string))
#   
#   if (length(phospho_part) == 0) {
#     return(list(NULL, NULL))
#   }
#   
#   phospho_text <- phospho_part[[1]]
#   letter_numbers <- gregexpr("[A-Z](\\d+)", phospho_text)
#   extracted_letters <- regmatches(phospho_text, letter_numbers)[[1]]
#   extracted_values <- gregexpr("\\((\\d+(?:\\.\\d+)?)\\)", phospho_text)
#   extracted_values <- regmatches(phospho_text, extracted_values)[[1]]
#   
#   extracted_values <- gsub("\\(|\\)", "", extracted_values)  # Remove parentheses
#   extracted_letters <- gsub("[A-Z]", "", extracted_letters)  # Remove letters
#   
#   combined_results <- list(paste(extracted_letters, collapse = "&"), paste(extracted_values, collapse = "&"))
#   return(combined_results)
# }
######################
file_path <- "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/PD_data_analysis/target_decoy_no_FAIMS/"
file_name <- "Multiconsensus_Exp2_TargetDecoy_woFAIMS_230612_FerriesPhosphoMarkers_PSMs.txt"
selected_spcies="HUMAN"
background_species= "ECOLI"
theo_file_path= "D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/"
theo_file_name="Synthetic peptides list_theo_conc_corrected_pool_id_iso_count_final.xlsx"
sheet_theo_name = "ISO-ref and OTHER with FC"

setwd(file_path)
quant_peptides <- read.table(file_name, sep = "\t", header = T)

pep_list_w_theo_quant <- read.xlsx(paste0(theo_file_path, theo_file_name), sheet = sheet_theo_name)
pep_list_w_theo_quant <- pep_list_w_theo_quant[,-1]

common_col_theo_quant <- as.data.frame(paste(pep_list_w_theo_quant$Phosphopeptide.sequence,
                                             pep_list_w_theo_quant$modified.position.in.peptide, sep = "_"))
colnames(common_col_theo_quant) <- "pep_with_pos"
pep_list_w_theo_quant_new <- cbind(common_col_theo_quant,pep_list_w_theo_quant)


quant_phospho_peptides <-quant_peptides %>%
  filter(grepl("Homo sapiens",Master.Protein.Descriptions) & 
           grepl("Phospho",Modifications)) %>%#& 
           #!grepl("positions not distinguishable", Modification.Pattern)) %>%
  select(Sequence,
         Modifications,
         Original.Precursor.Charge,
         #Number.of.PSMs,
         Delta.Score,
         Spectrum.File,
         First.Scan,
         Last.Scan,
         Ions.Score,
         Quan.Info,
         ptmRS.Binomial.Peptide.Score,
         ptmRS.Isoform.Confidence.Probability,
         ptmRS.Best.Site.Probabilities,
         ptmRS.Phospho.Site.Probabilities,
         RT.in.min,
         PSM.Ambiguity,
         Master.Protein.Descriptions,
         Protein.Accessions,
         Marked.as,
         Spectrum.File, starts_with("Abundances.Normalized")) %>%
  rowwise() %>%
  separate(ptmRS.Best.Site.Probabilities, into = c("first","second","third",
                                                   'fourth','fifth',"sixth",
                                                   'seventh'),sep = ";",remove = F) %>%
  separate(first,into = c("position","tmp"),sep = "\\(",remove = F) %>%
  separate(position, into =c("aa","pos"),sep = c(1,2)) %>%
  mutate(pep_with_pos=paste(Sequence,pos,sep="_"))  %>%
  full_join(pep_list_w_theo_quant_new,by="pep_with_pos") %>%
  separate(Spectrum.File, into = c("file_name","tmp"),sep = "\\.") %>%
  select(!tmp) %>%
  mutate(new_col_merging = paste(First.Scan,file_name,sep = "@")) #Sequence


#write.table(quant_phospho_peptides,file = paste0(file_path,"exp2_psm_merge_theo_list_pd_dda_no_faims_target_decoy.txt"),sep = "\t",col.names = T,row.names = F)


file_path <- "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/Proline_investigating_missing_values_compared_to_PD/"

file_list <- list.files(file_path)
file_list <- file_list[-c(1,2)]
final_df <- NULL

for (i in 1:length(file_list)){
  
  file_name_list <- paste0(file_list[i], "/", list.files(paste0(file_path,file_list[i])))
  
  for (j in 1:length(file_name_list)){
    tmp <- read.xlsx(paste0(file_path,file_name_list[j]),sheet = "Peptide Match")
    #tmp$file_name <- file_list[j]
    assign(paste0("final_df"), rbind(final_df, tmp))
           
  }
  
}

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

#write.table(final_df,file = paste0(file_path,"exp2_ms_queries_missing_Proline_dda_no_faims_target_decoy.txt"),sep = "\t",col.names = T,row.names = F)

ptms <- proline_phospho_pos_extraction(df = final_df$PTMs)
test <- cbind(final_df, paste(final_df$Peptide, ptms,sep = "_"))
colnames(test)[16] <- "pep_with_pos"
colnames(test)[15] <- "spectrum_title"

file_path <- "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/exp2_re_injection/"
file_name <- "PAL _Phosphopeptides exp2 ( 5 conc 3reps) DDA_230117 with Design_2023-06-07_1003.xlsx"


proline_actual_df <- read.xlsx(paste0(file_path,file_name),sheet = "Best PSM from protein sets")

merge_msquery <- test %>% left_join(proline_actual_df,by="spectrum_title") 


test1 <- merge_msquery %>% filter(is.na(sequence)) %>%
  separate(spectrum_title,
                           into = c("first_cycle","last_cycle","first_scan",
                                    "last_scan","first_time","last_time",
                                    "raw_file"),sep = ";") %>%
  separate(first_scan, into = c("tmp","first_scan1"),sep = ":") %>%
  separate(raw_file, into = c("tmp","file_name"),sep = ":") %>%
  mutate(new_col_merging =  paste(first_scan1,file_name,sep = "@")) %>% #Peptide
  left_join(quant_phospho_peptides,by="new_col_merging")

write.table(test1,file = paste0(file_path,"NEW_exp2_combining_missing_Proline_to_pd_dda_no_faims_target_decoy.txt"),sep = "\t",col.names = T,row.names = F)



write.table(merge_msquery,file = paste0(file_path,"exp2_combining_Proline_MSQuery_main_output.txt"),sep = "\t",col.names = T,row.names = F)















