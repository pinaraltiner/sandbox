#library(PhosR)
library(stringr)
library(dplyr)
library(data.table)
library(openxlsx)
library(tidyr)
library(ggplot2)
###############################################
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/ggplot/ggplot_functions.R")

file_path="D:/dev/Pinar/PHD/wet_lab_experiments/DIA_data_analysis/experiment_2/Spectronaut/trypsinP/"
acquisiton_types <- c("DIA Exploris no FAIMS") #DIA TIMS-TOF will be added later.
experiment_name <- c("E2-A1-R1","E2-A1-R2",
                     "E2-A1-R3",
                     "E2-A2-R1",
                     "E2-A2-R2",
                     "E2-A2-R3",
                     "E2-A3-R1",
                     "E2-A3-R2",
                     "E2-A3-R3",
                     "E2-A4-R1",
                     "E2-A4-R2",
                     "E2-A4-R3",
                     "E2-A5-R1",
                     "E2-A5-R2",
                     "E2-A5-R3")

exp_design <- experiment_name
final_spectronaut_pep_quant_analysis_syn <- function(file_path,
                                                 file_name,
                                                 sheet_name,
                                                 theo_file_path,
                                                 theo_file_name,
                                                 sheet_theo_name,
                                                 background_species,
                                                 selected_spcies,
                                                 exp_id,
                                                 exp_design,
                                                 fdr_threshold,
                                                 acquisiton_type,
                                                 software_name,
                                                 test_type,
                                                 num_reps,
                                                 actual_ratio,
                                                 subtitle){
    sample_size <- length(exp_design) / num_reps
    sample_names <- paste0("A",1:sample_size)
    comparisons <- NULL
    for (i in 1:sample_size){
        tmp <- paste0(sample_names[1], "/",sample_names[i])
        comparisons[i] <- tmp
        rm(tmp)
    }
    comparisons <- comparisons[-1]
    
    quant_peptides <- read_tsv(paste0(file_path,file_name))
    
    pep_list_w_theo_quant <- read.xlsx(paste0(theo_file_path, theo_file_name), sheet = sheet_theo_name)
    pep_list_w_theo_quant <- pep_list_w_theo_quant[,-1]

    
    
    
    
    quant_phospho <- quant_peptides %>% filter(grepl("HUMAN", PG.ProteinLabel)) %>%
        filter(grepl("Phospho",EG.PrecursorId))
    # Function to find max value and corresponding position
    find_max_value_and_pos <- function(ptm_prob, ptm_pos) {
        prob_values <- as.numeric(unlist(strsplit(ptm_prob, ";")))
        pos_values <- as.numeric(unlist(strsplit(ptm_pos, ";")))
        max_prob <- max(prob_values)
        max_pos <- pos_values[which(prob_values == max_prob)]
        return(paste(max_prob, max_pos, sep = "_"))
    }
    
    # Apply the function to each row
    max_prob_and_pos <- apply(quant_phospho, 1, function(row) {
        find_max_value_and_pos(row["EG.PTMProbabilities [Phospho (STY)]"], row["EG.PTMPositions [Phospho (STY)]" ])
    })
    
    
    
    
    
    find_max_value_and_pos <- function(ptm_prob, ptm_pos) {
        prob_values <- as.numeric(unlist(strsplit(ptm_prob, ";")))
        pos_values <- as.numeric(unlist(strsplit(ptm_pos, ";")))
        max_prob <- max(prob_values)
        max_pos <- pos_values[which(prob_values == max_prob)]
        
        result <- list()
        
        if (length(max_pos) == 1) {
            result$position <- paste(max_prob, max_pos, sep = "_")
            result$undistinguishable <- NA
        } else {
            for (i in 1:length(max_pos)){
                result<-list.append(paste(max_prob, max_pos, sep = "_"))
            }
            
            result$undistinguishable <-"undistinguishable"
        }
        
        return(result)
    }
    
    
    
    test_result <-apply(quant_phospho, 1, function(row) {
        find_max_value_and_pos(row["EG.PTMProbabilities [Phospho (STY)]"], row["EG.PTMPositions [Phospho (STY)]" ])
    })
    
    test_result2 <- lapply(lapply(test_result, unlist), function(x) {
        names(x)[names(x) == "parent_categories"] <- "parent_categories1"
        data.frame(t(x))
    })
    
    library(plyr)
    test_result3 <- rbind.fill(test_result2)
    
    
    test_result_df <- as.data.frame(test_result)
    test_phospho <-cbind(test_result,quant_phospho)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        