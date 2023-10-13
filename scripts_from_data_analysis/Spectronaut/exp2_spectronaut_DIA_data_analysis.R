#library(PhosR)
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
                                                 loc_filter_opt,
                                                 loc_filter,
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
    #### FUNCTION FOR MERGING POSITION AND SCORING INFO
    #### FUNCTION FOR MERGING POSITION AND SCORING INFO
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
    
    ## THEORETICAL PEPTIDE LIST
    pep_list_w_theo <- read.xlsx(paste0(theo_file_path, theo_file_name), sheet = sheet_theo_name)
    pep_list_w_theo_quant <- pep_list_w_theo[,-1]
    
    common_col_theo_quant <- as.data.frame(paste(pep_list_w_theo_quant$Phosphopeptide.sequence,
                                                 pep_list_w_theo_quant$modified.position.in.peptide, sep = "_"))
    colnames(common_col_theo_quant) <- "pep_with_pos"
    pep_list_w_theo_quant_new <- cbind(common_col_theo_quant,pep_list_w_theo_quant)
    
    #################################################
    pep_list_w_theo_unique <- pep_list_w_theo %>% select(Phospopeptide.sequence, Pool) %>%
      distinct(Phospopeptide.sequence,.keep_all = TRUE) %>%
      rename(PEP.GroupingKey = Phospopeptide.sequence) %>%
      rename(Pool_for_seq_merge=Pool)
    
    ecoli_seq <- quant_peptides %>% 
      select(PEP.GroupingKey, EG.PrecursorId, PG.ProteinLabel) %>%
      filter(grepl(background_species,PG.ProteinLabel)) %>%
      mutate(species=background_species) 
    
    ecoli_seq_dist <-  quant_peptides %>% 
      select(PEP.GroupingKey, EG.PrecursorId, PG.ProteinLabel) %>%
      filter(grepl(background_species,PG.ProteinLabel)) %>%
      distinct(PEP.GroupingKey,.keep_all = T) %>%
      mutate(species=background_species)
    
    
    all_seq <- quant_peptides_with_cond %>% 
      filter(grepl(selected_spcies,PG.ProteinLabel) & 
               grepl("Phospho",EG.PrecursorId)) %>%
      mutate(species=selected_spcies) %>%
      full_join(pep_list_w_theo_unique,by="PEP.GroupingKey") %>%
      mutate_at("Pool_for_seq_merge", ~replace_na(.,"Unexpected")) %>%
      mutate(Pool_for_seq_merge= ifelse(is.na(species),"missing",Pool_for_seq_merge)) %>%
      filter(!grepl("Unexpected",Pool_for_seq_merge)) %>%
      bind_rows(ecoli_seq) %>% 
      mutate(Pool_for_seq_merge= ifelse(is.na(Pool_for_seq_merge),background_species,Pool_for_seq_merge))
    
    
    all_seq_syn <- all_seq %>% 
      select(PEP.GroupingKey, EG.PrecursorId, PG.ProteinLabel,Pool_for_seq_merge) %>% 
      distinct(PEP.GroupingKey, .keep_all = T) %>%
      bind_rows(ecoli_seq_dist) %>%
      mutate(Pool_for_seq_merge= ifelse(is.na(Pool_for_seq_merge),background_species,Pool_for_seq_merge)) %>%
      mutate(acq_type=acquisiton_type) %>%
      mutate(soft_name=software_name)
      
    ## This is the same as p13 in DDA data analysis
    p12 <- gg_barplt_id_pep_count(data_set = all_seq_syn,
                                  x_df = all_seq_syn$Pool_for_seq_merge,
                                  fill_df = all_seq_syn$Pool_for_seq_merge,
                                  ymax = 20000,
                                  header = paste("Total number of identified phosphorylated", selected_spcies,"and", background_species,"across each sample",sep=" "),
                                  caption_lab = "NA values are removed. (p13)",
                                  x_lab = "Sample id",
                                  fill_lab =  "Sample id",
                                  y_lab = "Number of identified peptides",
                                  subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
    
    write.table(all_seq_syn, file=paste0(file_path,"Experiment2",software_name,"_number_of_unique_sequence_for_each_species.txt"),sep = "\t",col.names = T,row.names = F)
    #################################################
    ## PHOSPHO-FILTERING
    quant_phospho <- all_seq %>% filter(grepl("HUMAN", PG.ProteinLabel)) %>%
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
      select(PEP.GroupingKey,Experiment,Intensity,pep_with_pos,PG.ProteinLabel,ptm_score) %>%
      separate(Experiment, into = c("Exp_id","Sample_id","Rep_id"),sep = "-",remove = F) %>%
      #mutate(sample_rep_id_seq = paste(pep_with_pos, Sample_id,Rep_id, sep = "_")) %>%
      group_by(pep_with_pos,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
      slice(which.max(Intensity)) %>%
      ungroup()

    
    ####### ADDITIONAL PLOT TO DISPLAY MISSING and UNEXPECTED PEPTIDES ########
    df_merge_syn <- barplt_df %>%
      select(pep_with_pos,Experiment,Intensity, PG.ProteinLabel) %>% 
      pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
      full_join(pep_list_w_theo_quant_new,by="pep_with_pos") %>% 
      mutate_at("Pool", ~replace_na(.,"Unexpected")) %>%
      mutate(Pool= ifelse(is.na(PG.ProteinLabel),"missing",Pool)) %>%
      select(pep_with_pos,starts_with(exp_design),Pool) %>%
      mutate(soft_name=software_name)
    
    p11 <- gg_barplt_id_pep_count(data_set = df_merge_syn,
                                  x_df = df_merge_syn$Pool,
                                  fill_df = df_merge_syn$Pool,
                                  ymax = 20000,
                                  header = "Total number of quantified phospho-site across each sample",
                                  caption_lab = "NA values are removed.",
                                  x_lab = "Sample id",
                                  fill_lab =  "Sample id",
                                  y_lab = "Number of identified peptides",
                                  subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
    
    write.table(df_merge_syn,file = paste0(file_path,"Count_of_missing_unexpected_correct_phospho-sites_",
                                           software_name,"_Experiment",exp_id,".txt"),
                sep = "\t",row.names = F)
    #############################################################################
    
    barplt_df_ecoli <- quant_peptides_with_cond %>%
      filter(grepl(background_species, PG.ProteinLabel)) %>%
      select(PEP.GroupingKey,Experiment,Intensity,PG.ProteinLabel) %>%
      separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F) %>%
      mutate(sample_rep_id_seq = paste(PEP.GroupingKey, Sample_id,Rep_id, sep = "_")) %>%
      group_by(sample_rep_id_seq,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
      slice(which.max(Intensity)) %>%
      ungroup()
    
    
    
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
    
    if(loc_filter_opt == TRUE){
      site_prob <- barplt_df %>% 
        select(ptm_score, pep_with_pos,Experiment) %>% 
        pivot_wider(names_from = "Experiment",
                    values_from = "ptm_score") %>% 
        rowwise() %>%
        mutate(max_value = max(c_across(all_of(exp_design)), na.rm = TRUE)) %>%
        select(!exp_design)
      
      barplt_df_wide <- barplt_df %>%  ## If you select "charge" column, it will bring multiple rows for one seq
        select(PEP.GroupingKey,Experiment,Intensity,pep_with_pos,PG.ProteinLabel) %>%
        pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
        left_join(site_prob, by="pep_with_pos") %>%
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
   
    barplt_df_ecoli_wide <- barplt_df_ecoli %>%
      select(PEP.GroupingKey,Experiment,Intensity,PG.ProteinLabel)%>%
      pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
      mutate(species=background_species) %>% relocate(exp_design,.after = "PG.ProteinLabel")
    
    # Nothing is changed
    filtered_abundances<-barplt_df_wide[rowSums(!is.na(select(barplt_df_wide,starts_with(exp_design))))>0,]
    filtered_abundances_ecoli <-barplt_df_ecoli_wide[rowSums(!is.na(select(barplt_df_ecoli_wide,starts_with(exp_design))))>0,]
    
    library(kableExtra)
    na_phospho_selected <- apply(X = is.na(filtered_abundances %>% select(PEP.GroupingKey, PG.ProteinLabel,starts_with("E2"))), MARGIN = 2, FUN = sum)
    na_ecoli <- apply(X = is.na(filtered_abundances_ecoli %>% select(PEP.GroupingKey,PG.ProteinLabel,starts_with("E2"))), MARGIN = 2, FUN = sum)
    na_table <- bind_rows(na_phospho_selected,na_ecoli)
    na_table$species <- c(selected_spcies,background_species)
    na_table$total <- c(dim(filtered_abundances)[1],dim(filtered_abundances_ecoli)[1])
    
    na_table %>% select(-PEP.GroupingKey) %>%
      kbl(caption = paste("Number of NA values across all samples \n Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) %>%
      kable_material(c("striped", "hover")) %>%
      kable_styling(bootstrap_options = "striped", full_width = F, position = "left", font_size = 12) %>%
      kable_minimal(full_width = F) %>%
      footnote(general =  paste("This table was created after elimination of multiple charges by selecting either phospho-sites of", selected_spcies, "and", background_species, "sequences \n that has the highest abundace."),
               # number = c("Footnote 1; ", "Footnote 2; "),
               # alphabet = c("Footnote A; ", "Footnote B; "),
               # symbol = c("Footnote Symbol 1; ", "Footnote Symbol 2")
               footnote_as_chunk = T, title_format = c("italic", "underline")) %>%
      #as_image(width = 8) %>%
      save_kable(paste0(file_path,"outputs_with_new_script/table1.png"))
    
    filtered_abundances <- filtered_abundances %>% select(!`PEP.GroupingKey`)
    
    abundances_rowMeans <- NULL
    abundances_ecoli_rowMeans <- NULL
    for (k in 1:sample_size){
      # If separate version of row means is not needed, it can be commented later.
      # Separate row Means can be collected in temp object to merge in "log_10_filtered_abundances_rowMeans"
      assign(paste0("abundances_A",k),as.data.frame(rowMeans(filtered_abundances  %>%
                                                               select(contains(paste0("A",k)))))) #%>%
      #select(starts_with("abundance_")))))
      abundances_rowMeans<- bind_cols(abundances_rowMeans,get(paste0("abundances_A",k)))
      assign(paste0("ecoli_abundances_A",k),as.data.frame(rowMeans(filtered_abundances_ecoli  %>%
                                                                     select(contains(paste0("A",k))))))# %>%
      #select(starts_with("abundance_")))))
      abundances_ecoli_rowMeans<- bind_cols(abundances_ecoli_rowMeans,get(paste0("ecoli_abundances_A",k)))
    }
    
    #### MEAN ABUNDANCE RATIO WITH  DENSITY PLOT ####
    ### BEFORE IMPUTATION ###
    
    quant_peptides_ECOLI_density_plot <- filtered_abundances_ecoli %>%
      select(!starts_with("E")) %>%
      bind_cols(abundances_ecoli_rowMeans) %>%
      tibble() %>%
      rename_with(~ paste0("mean abundance",1:5), matches("^row")) %>%
      pivot_longer(cols = starts_with("mean"),
                   values_to = "Intensity",
                   names_to = "sample_ids",
                   values_drop_na = T) %>%
      mutate(sample_id_seq = paste(PEP.GroupingKey, sample_ids, sep = "_"))

    
    ecoli_density_plot<- quant_peptides_ECOLI_density_plot %>%
      select(contains(c("sample_ids","Intensity","species")))
    
    colnames(abundances_rowMeans) <- paste0("mean abundance",1:sample_size)
    
    quant_phospho_density_plot <- filtered_abundances %>%
      select(!starts_with("E")) %>%
      bind_cols(abundances_rowMeans) %>%
      #rename_with(~ paste0("mean_abun",1:5), matches("^row")) %>%
      tibble() %>% #mutate(pep_with_pos = sequence) %>% ###  At this stage, no need for phospho-position#
      pivot_longer(cols = starts_with("mean"),
                   names_to = "sample_ids",
                   values_to = "Intensity",
                   values_drop_na = T) %>%
      mutate(sample_id_seq = paste(pep_with_pos, sample_ids, sep = "_"))
    
    
    density_df <-quant_phospho_density_plot %>%
      select(c(sample_ids,Intensity,species)) %>%
      bind_rows(ecoli_density_plot) #%>%
    
    p3 <- gg_density(data_set = density_df,
                     x_df = density_df$Intensity,
                     fill_df = density_df$species,
                     color_df = NULL,
                     header="Distribution of mean abundance of every sample before imputation",
                     facet_df = "sample_ids",
                     x_lab = "log10(intensities)",
                     color_lab= "",
                     fill_lab = "species",
                     subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
    # Impute missing values
    for (j in 1:length(imputed_values_vec)){
      # Number NA
      #num_NA <- length(abundances_for_impute_all[,j+2][is.na(abundances_for_impute_all[,j+2])])
      filtered_abundances[,j+2][is.na(filtered_abundances[,j+2])] <- imputed_values_vec[j]
      filtered_abundances_ecoli[,j+2][is.na(filtered_abundances_ecoli[,j+2])] <- imputed_values_vec[j]
      #abundances_for_impute_all[,j+2][is.na(abundances_for_impute_all)[,j+2]] <- impute_values[j]
      # After imputation number of imputed values
      #num_imp <-length(abundances_for_impute_all[,j+2][(abundances_for_impute_all[,j+2]==impute_values[j])])
      # This is verification of imputation is done successfully
      # Because we expect to see that number of imputed values should be the same amount as number of NA
      #print(setequal(num_NA,num_imp))
      #print(num_NA)
      #print(num_imp)
    }
    
    # UNNECESSARY TO KEEP DATA BEFORE IMPUTATION
    abundances_all_aft_imputation <- filtered_abundances_ecoli %>%
      rename_with(~ paste0("pep_with_pos"), matches("^PEP")) %>%
      bind_rows(filtered_abundances)
    
    filtered_abundances_rowMeans <- NULL
    filtered_abundances_log10 <- NULL
    filtered_abundances_log10_rowMeans <- NULL
    for (k in 1:sample_size){
      assign(paste0("aft_imp_abundances_A",k),as.data.frame(rowMeans(abundances_all_aft_imputation  %>% select(contains(paste0("A",k))))))
      filtered_abundances_rowMeans<- bind_cols(filtered_abundances_rowMeans,get(paste0("aft_imp_abundances_A",k)))
      assign(paste0("log10_abundances_A",k),as.data.frame(log10(abundances_all_aft_imputation  %>% select(contains(paste0("A",k))))))
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
    
    final_imputed_data <- cbind(abundances_all_aft_imputation, filtered_abundances_rowMeans,filtered_abundances_log10,filtered_abundances_log10_rowMeans) #filtered_abundances
    final_imputed_data_syn <- final_imputed_data %>% filter(grepl(selected_spcies,species))
    final_imputed_data_ecoli <- final_imputed_data  %>% filter(!grepl(selected_spcies, species))
    
    df_merge <- final_imputed_data_syn %>%
      left_join(pep_list_w_theo_quant_new,by="pep_with_pos") %>%
      mutate_at("Pool", ~replace_na(.,"Unexpected")) %>%
      bind_rows(final_imputed_data_ecoli) %>%
      mutate_at("Pool", ~replace_na(.,background_species))
    
    df_mean_ab_after_impt <- df_merge %>% 
      select(contains("aft_imp") | contains("species"),Pool) %>%
      tibble() %>% 
      #separate(accession, into = c("uniprot_id", "species", "position"), remove = F) %>%
      pivot_longer(cols = contains("aft_imp"),
                   names_to = "Mean_abundance",
                   values_to = "values")# %>%
    #mutate_at("Pool", ~replace_na(.,background_species))
    
    
    df_FC_ratio_after_impt <- df_merge %>% 
      select(starts_with("exp_") | contains("species"),Pool) %>%
      #separate(accession, into = c("uniprot_id", "species", "position"), remove = F) %>%
      tibble() %>% 
      pivot_longer(cols = starts_with("exp_"),
                   names_to = "exp_FC",
                   values_to = "values")# %>%
    #mutate_at("Pool", ~replace_na(.,background_species))
    
    
    p4 <- gg_density(data_set = df_mean_ab_after_impt, 
                     x_df = df_mean_ab_after_impt$values,
                     fill_df = df_mean_ab_after_impt$Mean_abundance,
                     color_df = df_mean_ab_after_impt$Pool,
                     header="Distribution of mean abundance of every sample after imputation",
                     facet_df = "Mean_abundance",
                     x_lab = "log10(values)",
                     color_lab= "",
                     fill_lab = "Sample Names",
                     subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
    
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
    
    p5 <- gg_density(data_set = df_FC_ratio_after_impt,
                     x_df = df_FC_ratio_after_impt$values,
                     fill_df = df_FC_ratio_after_impt$exp_FC,
                     color_df = df_FC_ratio_after_impt$Pool,
                     header="Distribution of Fold change Ratio of every sample after imputation",
                     facet_df = "exp_FC",
                     x_lab = "log10(values)",
                     color_lab= "",
                     fill_lab = "Sample Names",
                     subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
    # 
    
    ### BOX-PLOT: Experimental Quantity Ratio of Phospho Peptides  
    
    p6 <- gg_boxplt_exp_ratio(data_set = df_FC_ratio_after_impt, 
                              x_df = df_FC_ratio_after_impt$exp_FC,
                              y_df = df_FC_ratio_after_impt$values,
                              fill_df = df_FC_ratio_after_impt$Pool,
                              header="Experimental Quantity Ratio of Phospho Peptides",
                              x_lab="Sample Names",
                              y_lab="Abundance Ratios",
                              fill_lab = "Sample Names",
                              subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
    
    
    ### HALF-BOX-PLOT & HALF-SCATTER-PLOT: Experimental Quantity Ratio of Synthetic Peptides  
    library(gghalves)
    
    p7 <- gg_half_boxplt_exp_ratio(data_set = df_FC_ratio_after_impt, 
                                   x_df = df_FC_ratio_after_impt$exp_FC,
                                   y_df = df_FC_ratio_after_impt$values,
                                   fill_df = df_FC_ratio_after_impt$Pool,
                                   header="Experimental Quantity Ratio of Phospho Peptides with Background",
                                   x_lab="Sample Names",
                                   y_lab="Abundance Ratios",
                                   fill_lab = "Sample Names",
                                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
    
    ### VIOLIN-PLOT: Experimental Quantity Ratio of Synthetic Peptides   
    
    ### TODO: fix y scaling without trimming 
    p8 <- gg_violin_exp_ratio(data_set = df_FC_ratio_after_impt, 
                              x_df = df_FC_ratio_after_impt$exp_FC,
                              y_df = df_FC_ratio_after_impt$values,
                              fill_df = df_FC_ratio_after_impt$Pool,
                              header="Experimental Quantity Ratio of Phospho Peptides with Background",
                              x_lab="Sample Names",
                              y_lab="Abundance Ratios",
                              fill_lab = "Sample Names",
                              trim=TRUE,
                              subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
    
    
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
    ## DATA FOR STATISTICAL ANALYSIS
    stat_analysis <- df_merge %>%
      select(pep_with_pos,Pool,species,isomericity,starts_with("log10_") | starts_with("mean_log10_") | starts_with("exp_FC")) %>%
      filter(!grepl(background_species, Pool))
    #filter(grepl("HUMAN",accession)) #spectrum_title
    
    rownames(stat_analysis) <- paste0(stat_analysis$pep_with_pos,"@",stat_analysis$Pool,"@",(1:nrow(stat_analysis))) #stat_analysis$spectrum_title,"@"
    library(multtest)
    if(test_type== "t.test" | test_type== "wilcoxon"){
      # Pre-allocate memory for results
      num_iterations <- sample_size - 1
      all_pvalues <- matrix(NA, nrow(stat_analysis), num_iterations)
      all_adjust_pval <- matrix(NA, nrow(stat_analysis), num_iterations)
      
      # Loop through iterations using lapply
      for (i in 2:sample_size) {
        p_values_tmp <- lapply(1:dim(stat_analysis)[1], function(j) {
          if (test_type == "t.test") {
            ttest_func(select(stat_analysis, contains("A1-") & contains("log10_"))[j,],
                       select(stat_analysis, contains(paste0("A", i, "-")) & contains("log10_"))[j,])
          } else if (test_type == "wilcoxon") {
            wilcox.test(select(stat_analysis, contains("A1-") & contains("log10_"))[j,],
                        select(stat_analysis, contains(paste0("A", i, "-")) & contains("log10_"))[j,])
          }
        })
        
        # Extract p-values from the list
        p_values_tmp <- sapply(p_values_tmp, function(x) x)
        
        # Store p-values in the pre-allocated matrix
        all_pvalues[, i - 1] <- p_values_tmp
        
        # Perform adjustment
        adjust_pval_tmp <- mt.rawp2adjp(p_values_tmp, proc = "BH", alpha = 0.05)
        qval <- data.frame(adjust_pval_tmp$adjp, adjust_pval_tmp$index)[order(adjust_pval_tmp$index), 2]
        
        # Store adjusted p-values in the pre-allocated matrix
        all_adjust_pval[, i - 1] <- qval
      }
      
      # Create data frames from matrices
      all_pvalues <- as.data.frame(all_pvalues)
      colnames(all_pvalues) <- paste0("pvalues_A1/", "A", 2:sample_size)
      rownames(all_pvalues) <- row.names(stat_analysis)
      
      all_adjust_pval <- as.data.frame(all_adjust_pval)
      colnames(all_adjust_pval) <- paste0("adjust_pval_A1/", "A", 2:sample_size)
      rownames(all_adjust_pval) <- row.names(stat_analysis)
      
      
      all_pvalues_common_col <- stat_analysis %>% 
        select(pep_with_pos, isomericity, Pool) %>% #spectrum_title
        bind_cols(all_pvalues) %>%  #all_adjust_pval
        pivot_longer(cols = starts_with("pvalues"), values_to = "pval", names_to ="pratios") %>% #adj_pval and adj_pratios
        separate(pratios, into = c("tmp","A1vs_Ai","tmp1"),sep = "_") %>%
        select(!c(tmp,tmp1)) %>%
        mutate(common_col = paste(pep_with_pos,Pool,A1vs_Ai,1:((sample_size-1)*nrow(stat_analysis)),sep="@")) #spectrum_title
      
      ## THE BEST WAY TO DO is this:
      merge_stat_df <- stat_analysis %>%
        select(pep_with_pos, Pool, isomericity,starts_with("exp_FC")) %>%
        pivot_longer(cols = starts_with("exp_FC"), values_to = "fold_change_values", names_to ="fold_change_ratios") %>%
        separate(fold_change_ratios, into = c("tmp","tmp1","ratio"),sep = "_") %>%
        select(!c(tmp,tmp1)) %>%
        #mutate(common_col = paste(pep_with_pos,Pool,1:((sample_size-1)*nrow(stat_analysis)),sep="@")) %>%
        bind_cols(all_pvalues_common_col$A1vs_Ai,all_pvalues_common_col$pval) %>%
        rename_with(.col =6 , ~"A1vs_Ai") %>%
        rename_with(.col=7, ~ "P.Value") %>%
        #filter(!grepl("ECOLI",Pool))
        #separate(accession, into = c("prot_id","species"),sep = "_")
        mutate(isomericity = ifelse(is.na(isomericity), "False Positive", isomericity)) %>%
        unite(Pool_new, Pool, isomericity,sep = "_",remove = FALSE) %>%
        unite('new_col_coloring',Pool_new,ratio,sep = "_",remove = FALSE) %>%
        mutate(new_col_coloring = if_else(grepl("ISO-REF", new_col_coloring), "ISO-REF", new_col_coloring)) %>%
        mutate(new_col_coloring = if_else(grepl("unexpected", new_col_coloring), "unexpected", new_col_coloring))
      
      
      ###############################################################################
      merge_stat_df_final <- merge_stat_df
      
      
    }else if(test_type=="limma"){
      library(limma)
      design_matrix <- model.matrix(~factor(c(rep(2,num_reps),rep(1,num_reps))))
      merge_stat_df <-NULL
      for ( i in 2:sample_size){
        # Change only the colname iteratively makes fit to every comparison
        colnames(design_matrix) <- c("Intercept", paste0("A1-A",i))
        #print(colnames(design_matrix))
        # Col selection for each comparison
        assign(paste0("df_A1vsA",i),stat_analysis %>% select(1:2 | contains("A1-") & contains("log10_") | contains(paste0("A",i,"-")) & contains("log10_")))
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
      colnames(merge_stat_df)[1] <- "common_col"
      
      merge_stat_df1 <- merge_stat_df %>% 
        separate(common_col, into = c("pep_with_pos","Pool","tmp"),sep = "@") %>%
        select(!tmp) %>%
        rename_with(.col=9, ~ "A1vs_Ai") %>%
        separate(A1vs_Ai, into = c("first","second"),sep = "-") %>%
        mutate(A1vs_Ai = paste(first,second,sep = "/")) %>%
        mutate(common_col = paste(pep_with_pos,Pool,A1vs_Ai,sep = "@")) %>%
        select(!c(first,second))
      
      merge_stat_df_final <- stat_analysis %>%
        select(pep_with_pos, Pool, isomericity,starts_with("exp_FC")) %>%
        pivot_longer(cols = starts_with("exp_FC"), values_to = "fold_change_values", names_to ="fold_change_ratios") %>%
        separate(fold_change_ratios, into = c("tmp","tmp1","ratio"),sep = "_") %>%
        select(!c(tmp,tmp1)) %>%
        #mutate(common_col = paste(common_col_for_merging,Pool,1:((sample_size-1)*nrow(stat_analysis)),sep="@")) %>%
        #bind_cols(merge_stat_df$logFC,merge_stat_df$P.Value,merge_stat_df$adj.P.Val) %>%
        mutate(common_col = paste(pep_with_pos,Pool,ratio,sep = "@")) %>%
        #rename_with(.col =8 , ~"common_col") %>%
        left_join(merge_stat_df1,by="common_col") %>%
        # rename_with(.col=7, ~ "pvalues") %>%
        #rename_with(.col=8, ~ "adj_pvalues") %>%
        
        #filter(!grepl("ECOLI",Pool))
        #separate(accession, into = c("prot_id","species"),sep = "_")
        mutate(isomericity = ifelse(is.na(isomericity), "False Positive", isomericity)) %>%
        unite(Pool_new, Pool.x, isomericity,sep = "_",remove = FALSE) %>%
        unite('new_col_coloring',Pool_new,ratio,sep = "_",remove = FALSE) %>%
        mutate(new_col_coloring = if_else(grepl("ISO-REF", new_col_coloring), "ISO-REF", new_col_coloring)) %>%
        mutate(new_col_coloring = if_else(grepl("unexpected", new_col_coloring), "unexpected", new_col_coloring)) %>%
        rename(pep_with_pos=pep_with_pos.x) %>% rename(Pool=Pool.x)
      
    }else{
      print("Statistical test could not be assessed. Check the input files!")
    }
    
    ## Generation of df -> expected abundance ratio for volcano plot
    actual_ratio_col <- merge_stat_df_final %>%
      select(A1vs_Ai) %>% distinct() %>%
      mutate(actual_ratio_val = case_when(grepl(comparisons[1],A1vs_Ai) ~actual_ratio[1],
                                          grepl(comparisons[2],A1vs_Ai) ~actual_ratio[2],
                                          grepl(comparisons[3],A1vs_Ai) ~actual_ratio[3],
                                          grepl(comparisons[4],A1vs_Ai) ~actual_ratio[4]))
    
    point_count_y_axis <- merge_stat_df_final %>%
      group_by(A1vs_Ai, new_col_coloring) %>%
      filter(P.Value < 0.05) %>% 
      count(new_col_coloring) %>% left_join(actual_ratio_col)
    
    
    ymax <- 9 + 0.5 #max(-log10(merge_stat_df_final$P.Value))
    y_decrement <- 0.50
    
    calculate_y_pos <- function(group) {
      group_length <- length(group)
      y_pos <- ymax - seq(0, by = y_decrement, length.out = group_length)
      return(y_pos)
    }
    
    # Apply the function to calculate y_pos within each group
    point_count_y_axis$y_pos <- unlist(by(point_count_y_axis$A1vs_Ai, point_count_y_axis$A1vs_Ai, calculate_y_pos))
    
    
    p9 <- ggplot(merge_stat_df_final,aes(x =log2(merge_stat_df_final$fold_change_values), y = -log10(merge_stat_df_final$P.Value))) +
      geom_point(aes(color = new_col_coloring,shape=Pool_new), size = 4) +
      #geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "red") +
      scale_fill_manual(values = c("ISO-REF" = "#000000",
                                   "Unexpected_False Positive_A1/A2"="#999999",
                                   "Unexpected_False Positive_A1/A3" ="#999999",
                                   "Unexpected_False Positive_A1/A4"="#999999",
                                   "Unexpected_False Positive_A1/A5"="#999999",
                                   "Others_multi_A1/A2" = "#CC79A7",
                                   "Others_mono_A1/A2" = "#CC79A7",
                                   "Others_multi_A1/A3" = "#E69F00",
                                   "Others_mono_A1/A3" = "#E69F00",
                                   "Others_multi_A1/A4" = "#56B4E9",
                                   "Others_mono_A1/A4" = "#56B4E9",
                                   "Others_multi_A1/A5" = "#009E73",
                                   "Others_mono_A1/A5" = "#009E73")) + 
      #geom_line(aes(color = new_col_coloring), size = 1) +  # Add color aesthetic to geom_line()
      scale_color_manual(values = c("ISO-REF" = "#000000",
                                    "Unexpected_False Positive_A1/A2"="#999999",
                                    "Unexpected_False Positive_A1/A3" ="#999999",
                                    "Unexpected_False Positive_A1/A4"="#999999",
                                    "Unexpected_False Positive_A1/A5"="#999999",
                                    "Others_multi_A1/A2" = "#CC79A7",
                                    "Others_mono_A1/A2" = "#CC79A7",
                                    "Others_multi_A1/A3" = "#E69F00",
                                    "Others_mono_A1/A3" = "#E69F00",
                                    "Others_multi_A1/A4" = "#56B4E9",
                                    "Others_mono_A1/A4" = "#56B4E9",
                                    "Others_multi_A1/A5" = "#009E73",
                                    "Others_mono_A1/A5" = "#009E73"),
                         
                         labels = c('Non-variant', 'Variant non-isomeric A1 vs A2',
                                    'Variant non-isomeric A1 vs A3',
                                    'Variant non-isomeric A1 vs A4',
                                    'Variant non-isomeric A1 vs A5',
                                    'Variant isomeric A1 vs A2',
                                    'Variant isomeric A1 vs A3',
                                    'Variant isomeric A1 vs A4',
                                    'Variant isomeric A1 vs A5',
                                    "Unexpected_False Positive A1/A2",
                                    "Unexpected_False Positive A1/A3",
                                    "Unexpected_False Positive A1/A4",
                                    "Unexpected_False Positive A1/A5")) +
      scale_shape_manual(values = c(16, 15, 12, 17),
                         labels = c('Non-variant', 'Variant non-isomeric', 'Variant isomeric', 'Unexpected')) +
      #scale_y_continuous(limits = c(0, max(-log10(merge_stat_df_final$adj.P.Val))), breaks = seq(0, max(-log10(merge_stat_df_final$adj.P.Val)), by = 0.8)) +
      #scale_x_continuous(limits = c(min(log2(merge_stat_df_final$fold_change_values)),max(log2(merge_stat_df_final$fold_change_values)))) +#facet_wrap(~ratio) +
      #scale_x_continuous(breaks = seq(from =round(min(log2(merge_stat_df_final$fold_change_values))), to=(round(max(log2(merge_stat_df_final$fold_change_values)))+2),by=1)) +
      #scale_y_continuous(breaks = seq(from =round(min(-log10(merge_stat_df_final$P.Value))), to=(round(max(-log10(merge_stat_df_final$P.Value)))+2),by=1)) +
      #scale_y_continuous(breaks = seq(0, max(-log10(volcano_final1$pvalues_value)), length.out = 21)) +
      scale_x_continuous(limits = c(-8, 10),breaks = seq(from = -8, to = 10, by = 2)) +  # Set the ticks for the x-axis
      scale_y_continuous(limits = c(0, 10),breaks = seq(from = 0, to = 10, by = 2))+  # Set the ticks for the y-axis
      
      theme_bw() +
      theme(legend.text = element_text(size = 30),
            axis.title.x = element_text(size = 30),
            axis.title.y = element_text(size = 30),
            plot.title = element_text(size = 35),
            legend.title = element_text(size = 30),
            axis.text.x = element_text(size = 30),
            axis.title = element_text(size = 30),
            axis.text.y = element_text(size = 30),
            plot.subtitle = element_text(size = 30)) +
      labs( y= "-log10(p values)", x="log2(fold change)",title = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name), subtitle = paste("Limma was used \n", subtitle)) +
      geom_vline(data = actual_ratio_col, aes(xintercept = log2(actual_ratio_val), show.legend = FALSE),color=c("#CC79A7","#E69F00","#56B4E9","#009E73"),size=2) +
      geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "red",size=2) + 
      geom_label(data = point_count_y_axis, aes(x = log2(actual_ratio_val), y = y_pos,fill=new_col_coloring, label = n),size=14, colour="white",show.legend = FALSE) 
    
    
    
    merge_stat_df_final_text <- merge_stat_df_final %>% mutate(soft_name=paste0(software_name)) %>% mutate(acq_type=paste0(acquisiton_type))
    write.table(merge_stat_df_final_text,file = paste0(file_path,"volcano_plot_",software_name,"_",acquisiton_type,".txt"),sep = 
                  "\t",col.names = T,row.names = F)
    
    df_roc <- merge_stat_df_final %>%
      select(pep_with_pos, Pool,P.Value)
    #filter(!grepl("Unexpected",Pool))
    
    ### ROC analysis custom func
    df_roc_order <- df_roc[order(df_roc$P.Value),]
    
    df_roc_func <- compute_roc_curve(df=df_roc_order, flag = "Others",expected = (4*141))
    
    p13 <- ggplot(df_roc_func, aes(y=tpr, x = fdr)) +
      geom_path(size=1.5) +
      #geom_vline(aes(xintercept=fdr)) +
      #geom_text(data=as.data.frame(result),aes(label=fdr)) +
      #scale_x_reverse() + 
      theme_bw() +
      theme(legend.text = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            plot.title = element_text(size = 25),
            legend.title = element_text(size = 20),
            axis.text.x = element_text(size = 20),
            axis.title = element_text(size = 20),
            axis.text.y = element_text(size = 20)) +
      scale_color_brewer(palette = "Dark2") +
      labs(y="True Positive Rate \n (Sensitivity)", x="False Positive Rate \n (Specificity)",
           title =  paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name), 
           subtitle = paste(subtitle,"including unexpected"), color="Pool Type")
    
    write.table(df_roc_func, file = paste0(file_path,"outputs_with_new_script/new_custom_Roc_analysis_",exp_id,"_",software_name,"_",".txt"),sep = "\t",row.names = F)
    
    #### ROC Analysis using pROC 
    
    df_roc$variant <- ifelse(df_roc$Pool == "Others", TRUE, FALSE)
    #df_roc$non_var <- ifelse(df_roc$Pool == "ISO-REF", TRUE, FALSE)
    
    library(pROC)
    # Calculate ROC curve for raw p-values
    roc_raw_variant <- roc(df_roc$variant, df_roc$P.Value)
    fpr <- as.data.frame(1 - roc_raw_variant$specificities)
    tpr_and_fpr_variant  <- cbind(roc_raw_variant$sensitivities,
                                  fpr,#roc_raw_variant$specificities,
                                  "Variant Pool")
    
    
    #roc_raw_non_var <- roc(df_roc$non_var, df_roc$P.Value)
    #tpr_and_fpr_non_var  <- cbind(roc_raw_non_var$sensitivities,
    #roc_raw_non_var$specificities,
    #"Non-variant Pool")
    
    roc_plot_df <- as.data.frame(tpr_and_fpr_variant) %>% 
      #bind_rows(as.data.frame(tpr_and_fpr_non_var)) 
      bind_cols(software_name)
    
    colnames(roc_plot_df) <- c("sensitivity", "fpr","Pool_type","Software_name")
    
    p10 <- roc_plot_df %>% #group_by(Pool_type) %>% 
      ggplot( aes(y=as.numeric(sensitivity), x = as.numeric(fpr), color=Pool_type)) +
      geom_path(size=1.5) + theme_bw() +# scale_x_reverse() +  +
      theme(legend.text = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            plot.title = element_text(size = 25),
            legend.title = element_text(size = 20),
            axis.text.x = element_text(size = 20),
            axis.title = element_text(size = 20),
            axis.text.y = element_text(size = 20)) +
      scale_color_brewer(palette = "Dark2") +
      labs(y="True Positive Rate \n (Sensitivity)", x="False Positive Rate \n (Specificity)",
           title =  paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name), 
           subtitle = paste(subtitle), color="Pool Type")
    
    write.table(roc_plot_df, file = paste0(file_path,"outputs_with_new_script/pRoc_analysis_",exp_id,"_",software_name,"_",".txt"),sep = "\t",row.names = F)
    
    sapply(1:13,function(x) ggsave(filename = paste0("p",x,".tiff"),
                                   width = 60, height = 45, 
                                   path = paste0(file_path,"/outputs_with_new_script/"),
                                   units = "cm",
                                   get(paste0("p",x)),
                                   device = "tiff", #".svg"
    ))
    
    
    }
    
    
    
    
    
    
    
    
        