#
library(stringr)
library(dplyr)
library(data.table)
library(openxlsx)
library(tidyr)
library(ggplot2)
library(tidyverse)
###############################################
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/ggplot/ggplot_functions.R")
#source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/roc_curve/roc_curve_generation_proline_edit.R")
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/roc_curve/new_roc_curve_generation_with_custom_threshold.R")

file_path <- "D:/dev/Pinar/PHD/wet_lab_experiments/DIA_data_analysis/experiment_3/Spectronaut/"
file_name <- "20230705_125301_OXPAL230127_Exp3_woFAIMS_DIA_Report.tsv"

selected_spcies="MOUSE"
background_species= "ECOLI"
acquisiton_type="DIA no FAIMS Exploris"
subtitle = ""
fdr_threshold = 0.05
actual_ratio = c(2,10,20,100)
exp_design= c("E3-A5-R1", "E3-A4-R1",
              "E3-A3-R1",
              "E3-A2-R1",
              "E3-A1-R1",
              "E3-A5-R2",
              "E3-A4-R2",
              "E3-A3-R2",
              "E3-A2-R2",
              "E3-A1-R2",
              "E3-A5-R3",
              "E3-A4-R3",
              "E3-A3-R3",
              "E3-A2-R3",
              "E3-A1-R3")
exp_id=3
software_name="Spectronaut"
num_reps=3
test_type="limma"
loc_filter_opt =TRUE
loc_filter =0.75

final_spectronaut_pep_quant_analysis_bio <- function(file_path,
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
                                                     loc_filter_opt,
                                                     loc_filter,
                                                     acquisiton_type,
                                                     software_name,
                                                     test_type,
                                                     num_reps,
                                                     actual_ratio,
                                                     subtitle){
  
  find_max_value_and_pos <- function(ptm_count, ptm_prob, ptm_pos) {
    
    prob_values <- as.numeric(unlist(strsplit(ptm_prob, ";")))
    pos_values <- as.numeric(unlist(strsplit(ptm_pos, ";")))
    result <- list()
    dim_count <- dim(as.data.frame(str_match_all(pattern = "\\[Phospho", ptm_count)))[1]
    max_prob <- max(prob_values)
    max_pos <- pos_values[which(prob_values == max_prob)]
    
    
    if (dim_count >= 2) {
      
      index_order <- order(prob_values, decreasing = TRUE)
      
      # Sort prob_values and pos_values using the same index order
      sorted_prob_values <- prob_values[index_order]
      sorted_pos_values <- pos_values[index_order]
      
      highest_dim <- sorted_prob_values[1:dim_count]
      
      #for (i in 1:length(highest_dim)){
      result["prob"] <- max(sorted_prob_values)#paste(sorted_prob_values[1:length(highest_dim)],collapse ="&")
      result["pos"] <- paste(sorted_pos_values[1:length(highest_dim)],collapse ="&")
      result$mod <- "two_phospho"
      #result<-list.append(paste(max_prob, max_pos, sep = "_"))
      #result$position <- paste(max_prob, max_pos, sep = "_")
      #result$position1 <- paste(max_prob, max_pos, sep = "_")
      #positions = paste0("position",i)
      #probs=paste0("prob",i)
      #max_pos <- 
      #result[positions] <- paste(pos_values[which(prob_values[i+1] == highest_dim[i])],pos_values[which(prob_values[i+1] == highest_dim[i])],sep = "&")
      #append(result[[positions]], paste(pos_values[i],pos_values[i],sep = "&"))
      #result[[probs]] <- append(result[[probs]], paste(sorted_probs[i],sorted_probs[i],sep = "&"))
      #}
    }else if (length(max_pos) == 1) {
      result$prob <- max_prob
      result$pos <- max_pos
      result$mod <- "mono phospho"
    } else {
      #### IF THIS PART CREATES AN ERROR, REMOVE THE POS AND 
      #### SCORE VALUE JUST RETURNED " Non-Distinguishable"
      for (i in 1:length(max_pos)){
        #result<-list.append(paste(max_prob, max_pos, sep = "_"))
        #result$position <- paste(max_prob, max_pos, sep = "_")
        #result$position1 <- paste(max_prob, max_pos, sep = "_")
        probs = paste0("prob",i)
        poses = paste0("pos",i)
        result[[probs]] <- append(result[[probs]], paste(max_pos[i]))
        result[[poses]] <- append(result[[poses]], paste(max_prob[i]))
        
        result$mod <- "non-distinguishable"
      }
    } 
    return(result)
  }
  
  sample_size <- length(exp_design) / num_reps
  sample_names <- paste0("A",1:sample_size)
  comparisons <- NULL
  for (i in 1:sample_size){
    tmp <- paste0(sample_names[1], "/",sample_names[i])
    comparisons[i] <- tmp
    rm(tmp)
  }
  comparisons <- comparisons[-1]
  
  ## EXPERIMENTAL DATA
  quant_peptides <- read_tsv(paste0(file_path,file_name))
  
  ##IMPUTED DATA
  ### IF GROUP_BY() does not work try to detach the plyr library: detach("package:plyr", unload = TRUE)
  # https://stackoverflow.com/questions/26923862/why-are-my-dplyr-group-by-summarize-not-working-properly-name-collision-with
  
  quant_peptides_with_cond <- quant_peptides %>% 
    mutate(Experiment=paste0(R.Condition,"-R",R.Replicate)) %>%
    rename("Intensity"= "EG.TotalQuantity (Settings)")
  
  imputed_values <- quant_peptides_with_cond  %>%
    group_by(Experiment) %>% 
    summarise(first_quantile=quantile(Intensity,probs=0.01,na.rm=TRUE))
  
  imputed_values_vec <- as.vector(imputed_values$first_quantile)
  
  ecoli_seq <- quant_peptides %>% 
    select(PEP.GroupingKey, EG.PrecursorId, PG.ProteinLabel) %>%
    filter(grepl(background_species,PG.ProteinLabel)) %>%
    mutate(species=background_species) 

  
  quant_phospho <- quant_peptides_with_cond %>% filter(grepl(selected_spcies, PG.ProteinLabel)) %>%
    filter(grepl("Phospho",EG.PrecursorId))
  
  # Apply the function to each row
  max_prob_and_pos <- apply(quant_phospho, 1, function(row) {
    find_max_value_and_pos(ptm_prob = row["EG.PTMProbabilities [Phospho (STY)]"],
                           ptm_pos = row["EG.PTMPositions [Phospho (STY)]" ],
                           ptm_count = row["EG.PrecursorId"])
  })
  
  # Determine the maximum number of columns
  max_length <- max(sapply(max_prob_and_pos, function(x) length(unlist(x))))
  
  # Create a data frame with the determined number of columns
  df <- data.frame(matrix(NA, ncol = max_length))
  
  # Fill the data frame with values from the list
  ##TODO: convert it to apply func instead of for loop
  for (i in 1:length(max_prob_and_pos)) {
    element <- max_prob_and_pos[[i]]
    df[i, 1:length(element)] <- unlist(element)
  }
  
  # Rename the columns as needed
  colnames(df) <- c("ptm_score","ptm_position","ptm_type")
  
  quant_phospho_combined <- quant_phospho %>% 
    bind_cols(df) %>% 
    #filter(!grepl("undistinguishable",undistinguishable)) %>%
    #separate(position1, into = c("phospho_score","phospho_pos"),sep = "_") %>%
    mutate(pep_with_pos=paste(PEP.GroupingKey,ptm_position,sep = "_"))
  
  barplt_df <- quant_phospho_combined %>%
    select(PEP.GroupingKey,Experiment,Intensity,PG.ProteinLabel,ptm_score) %>%
    separate(Experiment, into = c("Exp_id","Sample_id","Rep_id"),sep = "-",remove = F) %>%
    #mutate(sample_rep_id_seq = paste(pep_with_pos, Sample_id,Rep_id, sep = "_")) %>%
    group_by(PEP.GroupingKey,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    slice(which.max(Intensity)) %>%
    ungroup()
  
  if(loc_filter_opt == TRUE){
    site_prob <- barplt_df %>% 
      select(ptm_score, PEP.GroupingKey,Experiment) %>% 
      pivot_wider(names_from = "Experiment",
                  values_from = "ptm_score") %>% 
      rowwise() %>%
      mutate(max_value = max(c_across(all_of(exp_design)), na.rm = TRUE)) %>%
      select(!exp_design)
    
    barplt_df_wide <- barplt_df %>%  ## If you select "charge" column, it will bring multiple rows for one seq
      select(PEP.GroupingKey,Experiment,Intensity,PG.ProteinLabel) %>%
      pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
      left_join(site_prob, by="PEP.GroupingKey") %>%
      mutate(species=selected_spcies) %>%
      filter(max_value >= loc_filter) #%>%
    #select(!PEP.GroupingKey) #%>%
    #relocate(pep_with_pos, .after = Sequence) 
    #relocate(exp_design,.after = "pep_with_pos")
    
  }else{
    barplt_df_wide <- barplt_df %>%  ## If you select "charge" column, it will bring multiple rows for one seq
      select(PEP.GroupingKey,Experiment,Intensity,pep_with_pos,PG.ProteinLabel) %>%
      pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
      mutate(species=selected_spcies) %>% 
      #relocate(pep_with_pos, .after = Sequence)
      relocate(exp_design,.after = "pep_with_pos")
  }
  
  p1 <- gg_barplt_id_pep_count(data_set = barplt_df,
                               x_df = barplt_df$Sample_id,
                               fill_df = barplt_df$Rep_id,
                               ymax = 20000,
                               header = "Total number of quantified phospho-site across each sample",
                               caption_lab = "NA values are removed.",
                               x_lab = "Sample id",
                               fill_lab =  "Sample id",
                               y_lab = "Number of identified peptides",
                               subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
  barplt_df_ecoli <- quant_peptides_with_cond %>%
    filter(grepl(background_species, PG.ProteinLabel)) %>%
    select(PEP.GroupingKey,Experiment,Intensity,PG.ProteinLabel) %>%
    separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F) %>%
    mutate(sample_rep_id_seq = paste(PEP.GroupingKey, Sample_id,Rep_id, sep = "_")) %>%
    group_by(sample_rep_id_seq,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    slice(which.max(Intensity)) %>%
    ungroup()
  
  p2 <- gg_barplt_id_pep_count(data_set = barplt_df_ecoli,
                               x_df = barplt_df_ecoli$Sample_id,
                               fill_df = barplt_df_ecoli$Rep_id,
                               ymax = 20500,
                               header = "Total number of quantified Ecoli across each sample",
                               caption_lab = "NA values are removed.",
                               x_lab = "Sample id",
                               fill_lab =  "Sample id",
                               y_lab = "Number of identified peptides",
                               subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
  
 
  
  barplt_df_ecoli_wide <- barplt_df_ecoli %>%
    select(PEP.GroupingKey,Experiment,Intensity,PG.ProteinLabel)%>%
    pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
    mutate(species=background_species) %>% relocate(exp_design,.after = "PG.ProteinLabel")
  
  # Nothing is changed
  filtered_abundances<-barplt_df_wide[rowSums(!is.na(select(barplt_df_wide,starts_with(exp_design))))>0,]
  filtered_abundances_ecoli <-barplt_df_ecoli_wide[rowSums(!is.na(select(barplt_df_ecoli_wide,starts_with(exp_design))))>0,]
  
  
 
  
  
  
  
  
  barplt_df <- quant_phospho %>%
    select(PEP.GroupingKey,Experiment,Intensity,PG.ProteinLabel,ptm_score) %>%
    separate(Experiment, into = c("Exp_id","Sample_id","Rep_id"),sep = "-",remove = F) %>%
    #mutate(sample_rep_id_seq = paste(pep_with_pos, Sample_id,Rep_id, sep = "_")) %>%
    group_by(pep_with_pos,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    slice(which.max(Intensity)) %>%
    ungroup()
  
  
  
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  