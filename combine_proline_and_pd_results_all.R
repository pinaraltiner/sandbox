#
library(dplyr)
library(tidyr)
library(stringr)
############### FUNCTIONS TO EXTRACT POSITIONS FOR PROLINE AND PD ###############

extracted_characters <- function(df){
  pattern <- "(\\d+)\\[H O\\(3\\) P\\]"
  
  matches <- NULL
  
  for (i in 1:nrow(df)) {
    
    tmp <- str_match_all(df[i,], pattern)
    tmp <-as.data.frame(tmp)[,2]
    
    if(length(tmp) > 1){
      
      matches[i] <- paste(tmp, collapse = "&")
      
    }else{
      matches[i] <- tmp
    }
    
    rm(tmp)
  }
  return(matches)
}

pd_psm_pos_extract <- function(string){
  components <- unlist(strsplit(string, ";"))
  phospho_numbers <- c()
  
  # Loop through the components to find and extract numbers before "(Phospho)"
  for (comp in components) {
    if (grepl("\\(Phospho\\)", comp)) {
      number <- as.numeric(sub(".*?(\\d+)\\(Phospho\\).*", "\\1", comp))
      phospho_numbers <- c(phospho_numbers, number)
    }
  }
  
  # Combine the numbers with "&" if there are multiple
  if (length(phospho_numbers) > 0) {
    combined_numbers <- paste(phospho_numbers, collapse = "&")
    #print(combined_numbers)
    return(combined_numbers)
  } else {
    print("No '(Phospho)' components found.")
    return(NA)
  }
}
########################################################################### 

###############  PD DATA ###############   
pd_psms_path <- "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/PD_data_analysis/target_decoy_no_FAIMS/"
pd_psms <- read.table(paste0(pd_psms_path,"Multiconsensus_Exp2_TargetDecoy_woFAIMS_230612_FerriesPhosphoMarkers_PSMs.txt"), sep = "\t",header = T)

pd_psms_human <-pd_psms %>% 
  filter(grepl("Homo sapiens", Protein.Descriptions) & grepl("Phospho",Modifications)) %>%
  select(Sequence,
         Confidence,
         Master.Protein.Accessions,
         Modifications,
         Charge,
         First.Scan,
         Spectrum.File,
         Delta.Score,
         Ions.Score,
         contains("ptmRS"),
         Intensity) %>%
  mutate(extracted_number=sapply(Modifications,pd_psm_pos_extract)) %>%
  mutate(pep_with_pos = paste(Sequence,extracted_number,sep = "_")) %>%
  mutate(common_col=paste(Spectrum.File,First.Scan,sep = "@"))

############### PROLINE DATA ############### 
file_path <- "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/"
proline_data <- read.csv(file = paste0(file_path,"Proline_MS_Query_Exp2_noFAIMS_DDA.txt"),sep = ",")

phospho_data <- proline_data %>% 
  filter(grepl("P",ptm_string))

#test <- phospho_data[1:100,"ptm_string"]

result <- extracted_characters(df=as.data.frame(phospho_data$ptm_string))

result <- as.data.frame(result)

final_df <- phospho_data %>%
  bind_cols(result) %>% 
  mutate(pep_with_pos = paste(sequence,result,sep = "_")) %>%
  mutate(common_col=paste(raw_file_identifier,first_scan,sep = "@"))















