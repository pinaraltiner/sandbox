#library(PhosR)
library(stringr)
library(dplyr)
library(data.table)
library(openxlsx)
library(tidyr)
library(ggplot2)
###############################################
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/ggplot/ggplot_functions.R")
#source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/roc_curve/roc_curve_generation_proline_edit.R")
# Experiment 2
 #file_path <- "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/exp2_re_injection/"
 #file_name <- "PAL _Phosphopeptides exp2 ( 5 conc 3reps) DDA_230117 with Design_2023-06-07_1003.xlsx"


#file_path <- "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/with_FAIMS/"
#file_name <- "PAL _Exp2_( 5 conc 3reps)_withFAIMS_DDA_25052023_noDesign - correct_2023-06-13_1514.xlsx"
#selected_spcies = "_HUMAN"


# # Common constant objects
# sample_size <- 5
# sheet_name <- "Best PSM from protein sets"
# #sheet_name <- "Quantified peptide ions"
# exp_id <- 3
# background_species <- "ECOLI"

final_proline_pep_quant_analysis_syn <- function(file_path,
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
    
    quant_peptides <- read.xlsx(paste0(file_path,file_name), sheet = sheet_name)
    
    pep_list_w_theo_quant <- read.xlsx(paste0(theo_file_path, theo_file_name), sheet = sheet_theo_name)
    pep_list_w_theo_quant <- pep_list_w_theo_quant[,-1]
    
    #rechecking_df <- exp2_noFAIMS_DDA_proline %>% filter(grepl("HUMAN",accession)) %>% filter(grepl("Phospho",modifications))
    #write.table(rechecking_df,file="exp2_no_faims_Proline_output_filtered_human_phospho.tsv",sep="\t",row.names = FALSE)
    
    ##### REMOVE REDUNDANCY OF ECOLI SEQUENCE
    # ecoli_pep_max_absum <- exp2_noFAIMS_DDA_proline %>% filter(grepl("ECOLI", accession)) %>%
    #   mutate(across(all_of(sort_col), ~ ifelse(is.na(.), 0, .))) %>%
    #   rowwise() %>%
    #   mutate(row_sum = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
    #   group_by(sequence) %>%
    #   #summarise(max_row_sum = max(row_sum, na.rm = TRUE)) %>%
    #   slice(which.max(row_sum))
    
    ### IMPUTATION
    abundances_for_impute <- quant_peptides %>% 
        select(starts_with("abundance_")) %>% ## spectrum_title remove it because it was not make it as rownames (has duplicates)
        rename_with(~exp_design,matches("abundance"))
    
    quant_peptides_cor_abun <- quant_peptides %>% 
        rename_with(~exp_design,matches("^abundance"))
#################################################  
    ecoli_seq <- quant_peptides_cor_abun %>% 
      select(sequence, modifications, accession) %>%
      filter(grepl(background_species,accession)) %>%
      distinct(sequence,.keep_all = T)
    
    all_seq <- quant_peptides_cor_abun %>% 
      select(sequence, modifications, accession) %>%
      filter(grepl("Phospho",modifications) & grepl(selected_spcies,accession)) %>%
      distinct(sequence, .keep_all = T) %>%
      bind_rows(ecoli_seq) %>% 
      separate(accession, into = c("protein","species"),sep="_") %>%
      mutate(acq_type=acquisiton_type) %>%
      mutate(soft_name=software_name)
    
    p13 <- gg_barplt_id_pep_count(data_set = all_seq,
                           x_df = all_seq$species,
                           fill_df = all_seq$species,
                           ymax = 20000,
                           header = paste("Total number of identified phosphorylated", selected_spcies,"and", background_species,"across each sample",sep=" "),
                           caption_lab = "NA values are removed.",
                           x_lab = "Sample id",
                           fill_lab =  "Sample id",
                           y_lab = "Number of identified peptides",
                           subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))

    write.table(all_seq, file=paste0(file_path,"Experiment2",software_name,"_number_of_unique_sequence_for_each_species.txt"),sep = "\t",col.names = T,row.names = F)
    #################################################
    quant_phospho_peptides <- quant_peptides_cor_abun %>% 
        filter(grepl(selected_spcies,accession)) %>% 
        filter(grepl("Phospho",modifications)) %>%
        select(sequence,
               modifications,
               accession,
               ptm_protein_positions,
               charge,
               spectrum_title,
               starts_with(exp_design))
    
    quant_peptides_ECOLI <- quant_peptides_cor_abun %>% 
      filter(grepl(background_species,accession)) %>% 
      select(sequence,
             modifications,
             accession,
             spectrum_title,
             starts_with(exp_design))
    
    # Extraction of phospho positions from quant peptides object
    phospho_ptm_pos <- proline_phospho_pos_extraction(quant_phospho_peptides$ptm_protein_positions)
    # Data conversion 
    #phospho_ptm_pos_df <- t(as.data.frame(phospho_ptm_pos))
    #rownames(phospho_ptm_pos_df) <- 1:length(phospho_ptm_pos_df)
    
    # Creation of common column merging peptide sequence and phospho positions -> experimental data
    common_col_exp_quant <- as.data.frame(paste(quant_phospho_peptides$sequence, phospho_ptm_pos_df, sep = "_"))
    colnames(common_col_exp_quant) <- "pep_with_pos"
    
    # Bind it to the quant data
    syn_phospho_pep_proline_new <- cbind(common_col_exp_quant,quant_phospho_peptides)
    
    # Creation of common column merging peptide sequence and phospho positions -> theoretical data
    common_col_theo_quant <- as.data.frame(paste(pep_list_w_theo_quant$Phosphopeptide.sequence,
                                                 pep_list_w_theo_quant$modified.position.in.peptide, sep = "_"))
    colnames(common_col_theo_quant) <- "pep_with_pos"
    pep_list_w_theo_quant_new <- cbind(common_col_theo_quant,pep_list_w_theo_quant)
    
    # Nothing is changed
    filtered_abundances<-syn_phospho_pep_proline_new[rowSums(!is.na(select(syn_phospho_pep_proline_new,starts_with(exp_design))))>0,]
    filtered_abundances_ecoli <-quant_peptides_ECOLI[rowSums(!is.na(select(quant_peptides_ECOLI,starts_with(exp_design))))>0,]
    
    df_id_pep <- filtered_abundances %>% 
        select(sequence,modifications,charge, pep_with_pos, starts_with(exp_design),accession) %>%
        tibble() %>%
        # rename(A1_R1= 4, # Using column index to rename the colnames
        #        A1_R2= 5,
        #        A1_R3= 6,
        #        A2_R1= 7,
        #        A2_R2= 8,
        #        A2_R3= 9,
        #        A3_R1= 10,
        #        A3_R2= 11,
        #        A3_R3= 12,
        #        A4_R1= 13,
        #        A4_R2= 14,
    #        A4_R3= 15,
    #        A5_R1= 16,
    #        A5_R2= 17,
    #        A5_R3= 18) %>%
    pivot_longer(cols = starts_with("E2"), 
                 values_to = "intensity",
                 names_to = "sample_ids",
                 values_drop_na = T) %>%
        separate(sample_ids, into = c("Exp_id","Sample_id", "Rep_id"), sep = "_",remove = F)# %>% 
    #group_by(Sample_id) %>%
    #count()
    
    ### ELIMINATION of MULTIPLE CHARGED PEPTIDES ### 
    
    #### 1st STRATEGY: without pivot_longer() func. 
    
    ### IT's an alternative way to remove duplicates from mouse dataset
    ### DOES NOT WORK BACKGROUND DF BECAUSE POSITION IS IMPORTANT FOR THE DISCRIMINATION
    
    # syn_phospho_pep_proline_unique_max_absum <- filtered_abundances %>% 
    #  #mutate(across(all_of(experiment_name), ~ ifelse(is.na(.), 0, .))) %>%
    #  rowwise() %>%
    #  mutate(row_sum = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
    #  group_by(pep_with_pos) %>%
    #  #summarise(max_row_sum = max(row_sum, na.rm = TRUE)) %>%
    #  slice(which.max(row_sum)) %>%
    #  ungroup()
    
    ### THE RESULT OF barplt_df is the same as syn_phospho_pep_proline_unique_max_absum
    
    ## 2nd STRATEGY: removing duplicates automatically
    
    # barplt_df <- df_id_pep %>%
    #   mutate(sample_rep_id_seq = paste(pep_with_pos, Sample_id,Rep_id, sep = "_"))%>%
    #   filter(duplicated(sample_rep_id_seq)==FALSE)
    
    ## 3rd STRATEGY: removing duplicates the one has lower abundances
    barplt_df <- df_id_pep %>%
        #mutate(sample_rep_id_seq = paste(common_col_for_merging, Sample_id,Rep_id, sep = "@"))%>%
        group_by(pep_with_pos,sample_ids) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
        slice(which.max(intensity)) %>%
        ungroup()
    
    ####### ADDITIONAL PLOT TO DISPLAY MISSING and UNEXPECTED PEPTIDES ########
    df_merge_syn <- barplt_df %>%
      select(pep_with_pos,sample_ids,intensity, accession) %>% 
      pivot_wider(names_from = "sample_ids",values_from = "intensity") %>%
      full_join(pep_list_w_theo_quant_new,by="pep_with_pos") %>% 
      mutate_at("Pool", ~replace_na(.,"Unexpected")) %>%
      mutate(Pool= ifelse(is.na(accession),"missing",Pool)) %>%
      select(pep_with_pos,starts_with(exp_design),Pool) %>%
      mutate(soft_name=software_name,ion_mobility=acquisiton_type)
    
    ##############################################################################
    ####### GATHERING ALL COLUMNS OF MAIN OUTPUT FROM PROLINE WITH THE CORRECT RESULTS ########
     ### This is necessary only for Proline and PD additionally to compare 
      ## the missing peptides with their scan number.
    
    merge_phospho_peptides <- quant_peptides_cor_abun %>% 
      filter(grepl(selected_spcies,accession)) %>% 
      filter(grepl("Phospho",modifications)) 
    
    merge_phospho_pos <- proline_phospho_pos_extraction(merge_phospho_peptides[,"ptm_protein_positions"])
    
    pep_with_pos_merge <- cbind(merge_phospho_peptides, paste(merge_phospho_peptides$sequence,merge_phospho_pos,sep = "_"))
    colnames(pep_with_pos_merge)[dim(pep_with_pos_merge)[2]] <- "pep_with_pos"
    
    df_merge_all_col <- pep_with_pos_merge %>% tibble() %>%
      full_join(pep_list_w_theo_quant_new,by="pep_with_pos") %>%
      mutate_at("Pool", ~replace_na(.,"Unexpected")) %>%
      mutate(Pool= ifelse(is.na(accession),"missing",Pool)) %>%
      #select(pep_with_pos,starts_with(exp_design),Pool) %>%
      mutate(soft_name=software_name,ion_mobility=acquisiton_type)
    
    ################################################################################ 
    
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
    
    write.table(df_merge_all_col,file = paste0(file_path,"Count_of_missing_unexpected_correct_phospho-sites_with_all_col_",
                                           software_name,"_Experiment",exp_id,".txt"),
                sep = "\t",row.names = F)
    #############################################################################
    
    
    ## SAME STRATEGIES ABOVE (3rd) WAS APPLIED TO BACKGROUND AS WELL
    barplt_df_ecoli <- filtered_abundances_ecoli %>% 
        select(sequence,modifications, sequence, starts_with(exp_design),accession) %>%
        pivot_longer(cols = starts_with("E2"), 
                     values_to = "intensity",
                     names_to = "sample_ids",
                     values_drop_na = T) %>%
        separate(sample_ids, into = c("Exp_id","Sample_id", "Rep_id"), sep = "_",remove = F) %>%
        mutate(sample_rep_id_seq = paste(sequence, Sample_id,Rep_id, sep = "@")) %>%
        group_by(sample_rep_id_seq,sample_ids) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
        slice(which.max(intensity)) %>%
        ungroup()
    
    
    barplt_phospho_seq <- filtered_abundances %>%   
      select(sequence,modifications, sequence, starts_with(exp_design),accession) %>%
      pivot_longer(cols = starts_with("E2"), 
                   values_to = "intensity",
                   names_to = "sample_ids",
                   values_drop_na = T) %>%
      separate(sample_ids, into = c("Exp_id","Sample_id", "Rep_id"), sep = "_",remove = F) %>%
      mutate(sample_rep_id_seq = paste(sequence, Sample_id,Rep_id, sep = "@")) %>%
      group_by(sample_rep_id_seq,sample_ids) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
      slice(which.max(intensity)) %>%
      ungroup() %>% mutate(Software_name=software_name) %>%
      mutate(Acquisition_type=acquisiton_type)
    
    write.table(barplt_phospho_seq, file = paste0(file_path,"Number_of_human_phospho_sequences_",
                                                  software_name,"_Experiment",exp_id,".txt"),
                sep = "\t",row.names = F)
    
    #colnames(abundances_for_impute) <- experiment_name
    
    # 
    # # Calculate 1 percent quantile of each sample
    # triplicate_indx <- c(1,3,4,6,7,9,10,12,13,15)
    # impute_values_2NA <- NULL
    # impute_values_1NA <- NULL
    # tmp_2NA <- NULL
    # tmp_1NA <- NULL
    # #imputed_abundances <- NULL
    # for (i in 1:5){
    #   tmp_1NA <- quantile (abundances_for_impute[,experiment_name[triplicate_indx[1]:triplicate_indx[2]]],
    #                              probs = 0.01 , na.rm = TRUE )
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
    
    p12 <- gg_barplt_id_pep_count(data_set = barplt_phospho_seq,
                                  x_df = barplt_phospho_seq$Sample_id,
                                  fill_df = barplt_phospho_seq$Rep_id,
                                  ymax = 20000,
                                  header = "Total number of quantified phospho-sequence across each sample",
                                  caption_lab = "NA values and multiple sequences are removed.",
                                  x_lab = "Sample id",
                                  fill_lab =  "Sample id",
                                  y_lab = "Number of identified peptides",
                                  subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
    
    
    ## THIS RESHAPING IS ONLY FOR ELIMINATION OF MULTIPLE PHOSPHO-SITES and ECOLI PEPTIDES
    ## ELIMINATION STEP IS NOT NECESSARY FOR ECOLI, 1st STRATEGY can be used only (this will decrease lines of code)
    
    barplt_df_wide <- barplt_df %>%  ## If you select "charge" column, it will bring multiple rows for one seq
        select(pep_with_pos,accession, sample_ids,intensity) %>%
        pivot_wider(names_from = "sample_ids",values_from = "intensity")
    
    barplt_df_ecoli_wide <- barplt_df_ecoli %>% 
        select(sequence,sample_ids, accession,intensity) %>%
        pivot_wider(names_from = "sample_ids",values_from = "intensity")
    
    library(kableExtra)
    
    na_phospho_mouse <- apply(X = is.na(filtered_abundances %>% select(sequence, accession,starts_with("E2"))), MARGIN = 2, FUN = sum)
    na_ecoli <- apply(X = is.na(filtered_abundances_ecoli %>% select(sequence,accession,starts_with("E2"))), MARGIN = 2, FUN = sum)
    
    na_table <- bind_rows(na_phospho_mouse,na_ecoli)
    na_table$species <- c(selected_spcies,background_species)
    na_table$total <- c(dim(filtered_abundances)[1],dim(filtered_abundances_ecoli)[1])
    
    na_table %>% select(-sequence) %>%
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
        save_kable(paste0(file_path,"/outputs_with_new_script/table1.png"))
    
    ## DENSITY PLOT OF BEFORE IMPUTATION 
    
    abundances_rowMeans <- NULL
    abundances_ecoli_rowMeans <- NULL
    for (k in 1:sample_size){
        # If separate version of row means is not needed, it can be commented later.
        # Separate row Means can be collected in temp object to merge in "log_10_filtered_abundances_rowMeans"
        assign(paste0("abundances_A",k),as.data.frame(rowMeans(barplt_df_wide  %>%
                                                                   select(contains(paste0("A",k)))))) #%>%
        #select(starts_with("abundance_")))))
        abundances_rowMeans<- bind_cols(abundances_rowMeans,get(paste0("abundances_A",k)))
        
        
        assign(paste0("ecoli_abundances_A",k),as.data.frame(rowMeans(barplt_df_ecoli_wide  %>%
                                                                         select(contains(paste0("A",k))))))# %>%
        #select(starts_with("abundance_")))))  
        
        abundances_ecoli_rowMeans<- bind_cols(abundances_ecoli_rowMeans,get(paste0("ecoli_abundances_A",k)))
        
    }
    
    
    quant_peptides_ECOLI_density_plot <- barplt_df_ecoli_wide %>%
        select(!starts_with("E")) %>%
        bind_cols(abundances_ecoli_rowMeans) %>%
        tibble() %>%
        rename_with(~ paste0("mean_abun",1:5), matches("^row")) %>%
        pivot_longer(cols = starts_with("mean"), 
                     values_to = "intensity",
                     names_to = "sample_ids",
                     values_drop_na = T) %>%
        mutate(sample_id_seq = paste(sequence, sample_ids, sep = "_"))
    
    ecoli_density_plot<- quant_peptides_ECOLI_density_plot %>%
        select(contains(c("sample_ids","intensity","accession"))) 
    
    
    colnames(abundances_rowMeans) <- paste0("mean_abun",1:sample_size)
    
    #### MEAN ABUNDANCE RATIO WITH  DENSITY PLOT ####
    ### BEFORE IMPUTATION ###
    quant_phospho_density_plot <- barplt_df_wide %>%
        select(!starts_with("E")) %>%
        bind_cols(abundances_rowMeans) %>% 
        #rename_with(~ paste0("mean_abun",1:5), matches("^row")) %>%
        tibble() %>% #mutate(pep_with_pos = sequence) %>% ###  At this stage, no need for phospho-position#   
        pivot_longer(cols = starts_with("mean"),
                     names_to = "sample_ids",
                     values_to = "intensity",
                     values_drop_na = T) %>%
        mutate(sample_id_seq = paste(pep_with_pos, sample_ids, sep = "_"))
    
    density_df <-quant_phospho_density_plot %>%
        select(c(sample_ids,intensity,accession)) %>%
        bind_rows(ecoli_density_plot) %>%
        separate(accession, into = c("prot_id","species","position"),sep = "_",remove = F)
    
    
    p3 <- gg_density(data_set = density_df, 
                     x_df = density_df$intensity,
                     fill_df = density_df$species,
                     color_df = NULL,
                     header="Distribution of mean abundance of every sample before imputation",
                     facet_df = "sample_ids",
                     x_lab = "log10(intensities)",
                     color_lab= "",
                     fill_lab = "species",
                     subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
    
    
    # abundances_for_before_impt <- barplt_df_wide %>%
    #     bind_rows(barplt_df_ecoli_wide) %>% select(exp_design)
    # 
    # Calculate 1 percent quantile of each sample
    impute_values <- apply(abundances_for_impute, 2 , quantile , probs = 0.01 , na.rm = TRUE )
    
    # Impute missing values
    for (j in 1:length(impute_values)){
        # Number NA
        #num_NA <- length(abundances_for_impute_all[,j+2][is.na(abundances_for_impute_all[,j+2])])
        
        barplt_df_wide[,j+2][is.na(barplt_df_wide[,j+2])] <- impute_values[j]
        barplt_df_ecoli_wide[,j+2][is.na(barplt_df_ecoli_wide[,j+2])] <- impute_values[j]
        
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
    abundances_all_aft_imputation <- barplt_df_ecoli_wide %>%
        rename_with(~ paste0("pep_with_pos"), matches("^seq")) %>%
        bind_rows(barplt_df_wide) 
    
    
    
    filtered_abundances_rowMeans <- NULL
    filtered_abundances_log10 <- NULL
    filtered_abundances_log10_rowMeans <- NULL
    
    
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
    # To calculate all binary combination in the data frame
    #mat <- do.call(cbind, lapply(cols, function(xj) 
    #  sapply(cols, function(xi) (filtered_abundances_rowMeans[, xj]/(filtered_abundances_rowMeans[, xj])))))
    #colnames(mat) <-  outer(names(filtered_abundances_rowMeans), names(filtered_abundances_rowMeans), paste0)
    
    final_imputed_data <- cbind(abundances_all_aft_imputation, filtered_abundances_rowMeans,filtered_abundances_log10,filtered_abundances_log10_rowMeans) #filtered_abundances
    
    final_imputed_data_syn <- final_imputed_data %>% filter(grepl(selected_spcies, accession))
    
    final_imputed_data_ecoli <- final_imputed_data  %>% filter(!grepl(selected_spcies, accession))
    
    df_merge <- final_imputed_data_syn %>%
        left_join(pep_list_w_theo_quant_new,by="pep_with_pos") %>% 
        mutate_at("Pool", ~replace_na(.,"Unexpected")) %>%
        bind_rows(final_imputed_data_ecoli) %>%
        mutate_at("Pool", ~replace_na(.,background_species)) 
    
    #write.table(final_imputed_data, file = "final_imputed_normalized_data_PAL _T_cell_Exp3_( 5 conc 3reps)_NoFAIMS_DDA_with_cont_230206_2023-02-07_0947.txt",sep = "\t",row.names = F)
    
    df_mean_ab_after_impt <- df_merge %>% 
        select(contains("aft_imp") | contains("accession"),Pool) %>%
        tibble() %>% 
        separate(accession, into = c("uniprot_id", "species", "position"), remove = F) %>%
        pivot_longer(cols = contains("aft_imp"),
                     names_to = "Mean_abundance",
                     values_to = "values") %>%
        mutate_at("Pool", ~replace_na(.,background_species))
    
    
    df_FC_ratio_after_impt <- df_merge %>% 
        select(starts_with("exp_") | contains("accession"),Pool) %>%
        separate(accession, into = c("uniprot_id", "species", "position"), remove = F) %>%
        tibble() %>% 
        pivot_longer(cols = starts_with("exp_"),
                     names_to = "exp_FC",
                     values_to = "values") %>%
        mutate_at("Pool", ~replace_na(.,background_species))
    
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
                              header="Experimental Quantity Ratio of T-cell Phospho Peptides",
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
                                   header="Experimental Quantity Ratio of T-cell Phospho Peptides with Background",
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
                              header="Experimental Quantity Ratio of T-cell Phospho Peptides with Background",
                              x_lab="Sample Names",
                              y_lab="Abundance Ratios",
                              fill_lab = "Sample Names",
                              trim=TRUE,
                              subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
    
    #   selected_cols <- colnames(exp2_FAIMS_DDA_proline[c(1:4,6,28)])
    # #selected_cols <- colnames(exp2_noFAIMS_DDA_proline[c(1:4,15)])
    # 
    # syn_phospho_pep_proline <- exp2_FAIMS_DDA_proline %>% 
    #   select(-starts_with("abundance_")) %>% 
    #   mutate(abundances_for_impute) %>%
    #   filter(grepl("Phospho", ptm_protein_positions)) %>% 
    #   filter(grepl("_HUMAN", accession)) %>%
    #   select(contains(selected_cols) | starts_with("abundance"))
    # 
    # 
    # # Sorting columns A1 to A5 and R1 to R3
    # sorted_colnames <- c("psm_count","abundance") # Abundance contains both w/wo raw_abundance
    # for (i in 1:length(sorted_colnames)){
    #   assign(paste0("tmp",i), sort(colnames(select(syn_phospho_pep_proline_new,
    #                                                contains(sorted_colnames[i])))))
    # }
    # 
    # sort_col <- c(tmp1,tmp2)
    # 
    # sort_syn_phospho_pep_proline_new <- syn_phospho_pep_proline_new %>%
    #   relocate(all_of(sort_col), .after = 7) %>% rename_with(~ experiment_name, all_of(sort_col)) 
    # 
    # 
    ##### REMOVE REDUNDANCY OF PHOSPHO SEQUENCES DUE TO DIFFERENT CHARGES
    # syn_phospho_pep_proline_unique_max_absum <- sort_syn_phospho_pep_proline_new %>% 
    #   #mutate(across(all_of(experiment_name), ~ ifelse(is.na(.), 0, .))) %>%
    #   rowwise() %>%
    #   mutate(row_sum = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
    #   group_by(common_col_for_merging) %>%
    #   #summarise(max_row_sum = max(row_sum, na.rm = TRUE)) %>%
    #   slice(which.max(row_sum)) %>%
    #   ungroup()
    
    ####################################################
    # #### IMPUTATION BASED ON COMPLEX CONDITIONS ####
    # 
    # df_long <- syn_phospho_pep_proline_unique_max_absum %>% 
    #   select(common_col_for_merging,all_of(experiment_name)) %>%
    #   pivot_longer(cols = -common_col_for_merging, names_to = "column", values_to = "value") %>%
    #   separate(column, into = c("col_group", "col_index", "row_index"), sep = "_")
    # #mutate(col_index = str_remove(col_index, "A"))
    # 
    # # Count the number of NA values in each row for each column group
    # na_counts <- df_long %>%
    #   group_by(common_col_for_merging, col_index) %>%
    #   summarise(na_count = sum(is.na(value))) 
    #   #pivot_wider(names_from = col_index, values_from = na_count)
    #   
    #   # Replace the NA values based on the count
    #   df_result <- df_long %>%
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
    #   filter(rowSums(across(where(is.numeric)))!=0)
    # 
    # 
    # ####################
    
    # Combine data frame (we will continue with this for further step)
    #df_merge <- df_result %>% left_join(syn_phospho_pep_proline_unique_max_absum, 
    # by="common_col_for_merging") 
    #df_merge <- syn_phospho_pep_proline_unique_max_absum %>% left_join(pep_list_w_theo_quant_new,by="common_col_for_merging")
    
    ## ALWAYS KEEP IT COMMENTED TO PREVENT OVERWRITE
    #write.xlsx(df_merge,file = "D:/dev/Pinar/PHD/wet_lab_experiments/experiment2_quantification_peptide_level/Correct_imputation_exp2proline_and_peplist_theo_before_stast_analysis.xlsx")
    
    
    # filter_at(vars(Pool), all_vars(is.na(.)))
    
    #output_path <- "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/comparision_DDA_noFAIMS_mq_proline_pd/"
    #write.xlsx(unexpectedly_id_peps,file = paste0(output_path,"Corrected_imputation_proline_unexpectedly_identified_phosphopeptides_exp2_noFAIMS_DDA.xlsx"))
    
    #test_data_one_roc <- final_imputed_data %>% select(common_col_for_merging,pool_id)
    # 
    # stat_analysis <- final_imputed_data %>%
    #   select(common_col_for_merging | starts_with("log10_E2"))
    
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
        select(pep_with_pos,Pool,accession,isomericity,starts_with("log10_") | starts_with("mean_log10_") | starts_with("exp_FC")) %>%
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
                    ttest_func(select(stat_analysis, contains("A1_") & contains("log10_"))[j,],
                               select(stat_analysis, contains(paste0("A", i, "_")) & contains("log10_"))[j,])
                } else if (test_type == "wilcoxon") {
                    wilcox.test(select(stat_analysis, contains("A1_") & contains("log10_"))[j,],
                                select(stat_analysis, contains(paste0("A", i, "_")) & contains("log10_"))[j,])
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
            separate(pratios, into = c("tmp","ratio","tmp1"),sep = "_") %>%
            select(!c(tmp,tmp1)) %>%
            mutate(common_col = paste(pep_with_pos,Pool,ratio,1:((sample_size-1)*nrow(stat_analysis)),sep="@")) #spectrum_title
        
        ## THE BEST WAY TO DO is this:
        merge_stat_df <- stat_analysis %>%
            select(pep_with_pos, Pool, isomericity,starts_with("exp_FC")) %>%
            pivot_longer(cols = starts_with("exp_FC"), values_to = "fold_change_values", names_to ="fold_change_ratios") %>%
            separate(fold_change_ratios, into = c("tmp","tmp1","ratio"),sep = "_") %>%
            select(!c(tmp,tmp1)) %>%
            #mutate(common_col = paste(pep_with_pos,Pool,1:((sample_size-1)*nrow(stat_analysis)),sep="@")) %>%
            bind_cols(all_pvalues_common_col$ratio,all_pvalues_common_col$pval) %>%
            rename_with(.col =6 , ~"ratio1") %>%
            rename_with(.col=7, ~ "pvalues") %>%
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
            assign(paste0("df_A1vsA",i),stat_analysis %>% select(1:2 | contains("A1_") & contains("log10_") | contains(paste0("A",i,"_")) & contains("log10_")))
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
    
    
    ymax <- max(-log10(merge_stat_df_final$P.Value)) + 0.5
    y_decrement <- 0.15
    
    calculate_y_pos <- function(group) {
        group_length <- length(group)
        y_pos <- ymax - seq(0, by = y_decrement, length.out = group_length)
        return(y_pos)
    }
    
    # Apply the function to calculate y_pos within each group
    point_count_y_axis$y_pos <- unlist(by(point_count_y_axis$A1vs_Ai, point_count_y_axis$A1vs_Ai, calculate_y_pos))
    
    
    p9 <- ggplot(merge_stat_df_final,aes(x =log2(merge_stat_df_final$fold_change_values), y = -log10(merge_stat_df_final$P.Value))) +
        geom_point(aes(color = new_col_coloring,shape=Pool_new), size = 2.5) +
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
        scale_x_continuous(breaks = seq(from =round(min(log2(merge_stat_df_final$fold_change_values))), to=(round(max(log2(merge_stat_df_final$fold_change_values)))+2),by=1)) +
        scale_y_continuous(breaks = seq(from =round(min(-log10(merge_stat_df_final$P.Value))), to=(round(max(-log10(merge_stat_df_final$P.Value)))+2),by=1)) +
        #scale_y_continuous(breaks = seq(0, max(-log10(volcano_final1$pvalues_value)), length.out = 21)) +
        theme_bw() +
        theme(legend.text = element_text(size = 15),
              axis.title.x = element_text(size = 15),
              axis.title.y = element_text(size = 15),
              plot.title = element_text(size = 30),
              legend.title = element_text(size = 15),
              axis.text.x = element_text(size = 15),
              axis.title = element_text(size = 15),
              axis.text.y = element_text(size = 15),
              plot.subtitle = element_text(size = 15)) +
        labs( y= "-log10(p values)", x="log2(fold change)",title = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name), subtitle = paste("Limma was used \n",subtitle)) +
        geom_vline(data = actual_ratio_col, aes(xintercept = log2(actual_ratio_val), show.legend = FALSE),color=c("#CC79A7","#E69F00","#56B4E9","#009E73"),size=1) +
        geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "red",size=1) + 
        geom_label(data = point_count_y_axis, aes(x = log2(actual_ratio_val), y = y_pos,fill=new_col_coloring, label = n),size=6, colour="white",show.legend = FALSE) 
    
    
    
    
    
    
    #### ROC Analysis
    df_roc <- merge_stat_df_final %>%
        select(pep_with_pos, Pool,P.Value)
    
    df_roc$variant <- ifelse(df_roc$Pool == "Others", TRUE, FALSE)
    df_roc$non_var <- ifelse(df_roc$Pool == "ISO-REF", TRUE, FALSE)
    
    library(pROC)
    # Calculate ROC curve for raw p-values
    roc_raw_variant <- roc(df_roc$variant, df_roc$P.Value)
    tpr_and_fpr_variant  <- cbind(roc_raw_variant$sensitivities,
                                  roc_raw_variant$specificities,
                                  "Variant Pool")
    
    
    roc_raw_non_var <- roc(df_roc$non_var, df_roc$P.Value)
    tpr_and_fpr_non_var  <- cbind(roc_raw_non_var$sensitivities,
                                  roc_raw_non_var$specificities,
                                  "Non-variant Pool")
    
    roc_plot_df <- as.data.frame(tpr_and_fpr_variant) %>% 
        bind_rows(as.data.frame(tpr_and_fpr_non_var)) %>% bind_cols(software_name)
    
    colnames(roc_plot_df) <-    c("sensitivity", "specificity","Pool_type","Software_name")
    
    
    p10 <- roc_plot_df %>% group_by(Pool_type) %>% 
        ggplot( aes(y=as.numeric(sensitivity), x = as.numeric(specificity), color=Pool_type)) +
        geom_path(size=1.5) +  scale_x_reverse() + theme_bw() +
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
    
    write.table(roc_plot_df, file = paste0(file_path,"Roc_analysis_",exp_id,"_",software_name,"_",".txt"),sep = "\t",row.names = F) #acquisiton_type ## IT WAS TOO LONG-> GIVES AN ERROR
    
    sapply(1:13,function(x) ggsave(filename = paste0("p",x,".tiff"),
                                   width = 50, height = 45, 
                                   path = paste0(file_path,"/outputs_with_new_script/"),
                                   units = "cm",
                                   get(paste0("p",x)),
                                   device = "tiff", #".svg"
    ))
    
    
    
    
    #  
    #  
    #  
    #  
    #  
    #  
    #  df_roc <- merge_stat_df_final %>%
    #      select(pep_with_pos.x, Pool.x,P.Value) %>%
    #      
    #    
    #    roc_data_all <- compute_roc_curve(complete_pvalues_all,flag = "Others",expected = 145)
    #    
    #    
    #    ## ROC CURVE GENERATION ### 
    #    roc_data121 <- cbind(roc_data12, "A1 vs A2","Others", "Proline")
    #    roc_data131 <- cbind(roc_data13, "A1 vs A3","Others", "Proline")
    #    roc_data141 <- cbind(roc_data14, "A1 vs A4","Others", "Proline")
    #    roc_data151 <- cbind(roc_data15, "A1 vs A5","Others", "Proline")
    #    
    #    colnames(roc_data121)[c(4:6)] <- c("Comparison","Pool","Software")
    #    colnames(roc_data131)[c(4:6)] <- c("Comparison","Pool","Software")
    #    colnames(roc_data141)[c(4:6)] <- c("Comparison","Pool","Software")
    #    colnames(roc_data151)[c(4:6)] <- c("Comparison","Pool","Software")
    #    
    #    roc_final <- rbind(roc_data121,roc_data131,roc_data141,roc_data151)
    #    
    #    #roc_curve_data_exp2_noFAIMS_MQ <- read.delim("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/MQ_data_analysis/exp2_wo_FAIMS/exp2_wo_FAIMSwith_MBR/roc_curve_data_exp2_noFAIMS_MQ.txt")
    #    write.table(roc_final,file = "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/exp2_re_injection/exp2_no_FAIMS_proline_roc_curve_data_consider_isoref_fp.tsv",sep = "\t",row.names = F)
    #    #roc_curve_data_pro_faims <- read.delim("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/exp2_re_injection/exp2wo_FAIMS_DDA_Proline_volcano_plot_data.tsv",sep = "\t")
    #    
    #    
    #    roc_proline_mq <- rbind(roc_curve_data_exp2_noFAIMS_MQ,roc_final)
    #    
    #    
    #    roc_data_unexpected121 <- cbind(roc_data_unexpected12, "A1 vs A2","unexcepted")
    #    roc_data_unexpected131 <- cbind(roc_data_unexpected13, "A1 vs A3","unexcepted")
    #    roc_data_unexpected141 <- cbind(roc_data_unexpected14, "A1 vs A4","unexcepted")
    #    roc_data_unexpected151 <- cbind(roc_data_unexpected15, "A1 vs A5","unexcepted")
    #    
    #    colnames(roc_data_unexpected121)[c(4,5)] <- c("Comparison","Pool")
    #    colnames(roc_data_unexpected131)[c(4,5)] <- c("Comparison","Pool")
    #    colnames(roc_data_unexpected141)[c(4,5)] <- c("Comparison","Pool")
    #    colnames(roc_data_unexpected151)[c(4,5)] <- c("Comparison","Pool")
    #    
    #    
    #    
    #    roc_final2 <- rbind(roc_data_unexpected121,roc_data_unexpected131,
    #                        roc_data_unexpected141,roc_data_unexpected151)
    #    
    #    roc_curve_combine_proline_mq_score_40 <- read.delim("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/comparision_DDA_noFAIMS_mq_proline_pd/comparison_wrong_localized_pep_MQ_Proline_with_out_FAIMS/final_figures/roc_curve_combine_proline_mq_score_40.txt")
    #    
    #    line_types <- c("solid","dashed" )
    #    
    #    ##TODO: make it more professional
    #    ggplot(roc_final,aes(x=fdp,y=tpr)) + 
    #      geom_line(size=1.2,aes(color=Comparison,linetype=Software)) +
    #      scale_color_manual(values = c("#CC79A7", "#E69F00", "#56B4E9", "#009E73"),
    #                         labels=c('A1 vs A2','A1 vs A3','A1 vs A4','A1 vs A5')) +
    #      scale_linetype_manual(values = line_types)+
    #      theme_bw() +
    #      theme(legend.text = element_text(size=15), 
    #            axis.title.x = element_text(size = 15),
    #            axis.title.y = element_text(size = 15),
    #            plot.title = element_text(size=30),
    #            legend.title=element_text(size=15),
    #            axis.text.x=element_text(size=15),
    #            axis.title=element_text(size=15),
    #            axis.text.y = element_text(size = 15)) + 
    #      expand_limits(x = 0, y = 0) +
    #      scale_y_continuous(limits = c(0,100)) +
    #      #facet_wrap(~Ratio_col,scales = "free_x") +
    #      labs(title = "Experiment 2 - DDA without FAIMS - \n Comparision of Softwares")
    #    
    # #    
    # #  # 
    # #  # # iterate through each column of the data set and calculate t-tests
    # #  # triplicate_indx <- c(1,3,4,6,7,9,10,12,13,15) +1
    # #  # 
    # # 
    # #  p_values_12 <- NULL
    # #  p_values_13 <- NULL
    # #  p_values_14 <- NULL
    # #  p_values_15 <- NULL
    # # 
    # #  ##TODO put this into another for loop to reduce redundancy of the code
    # #    for(j in 1:dim(stat_analysis)[1]){
    # # 
    # #      #if ((sum(stat_analysis[j,2:4])/3) !=  stat_analysis[j,2:4][1] | (sum(stat_analysis[j,5:7])/3) !=  stat_analysis[j,5:7][1]){
    # # 
    # #      p_values_12[j] <- ttest_func(stat_analysis[j,3:5], stat_analysis[j,6:8])
    # #      #p_values_13[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,8:10])
    # #      #p_values_14[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,11:13])
    # #      #p_values_15[j] <- ttest_func(stat_analysis[j,2:4], stat_analysis[j,14:16])
    # #      #}else{
    # # 
    # #      #}
    # #      #if ((sum(stat_analysis[j,2:4])/3) !=  stat_analysis[j,2:4][1] | (sum(stat_analysis[j,8:10])/3) !=  stat_analysis[j,8:10][1]){
    # # 
    # #        p_values_13[j] <- ttest_func(stat_analysis[j,3:5], stat_analysis[j,9:11])
    # # 
    # #      #}else{
    # # 
    # #      #}
    # #      #if ((sum(stat_analysis[j,2:4])/3) !=  stat_analysis[j,2:4][1] | (sum(stat_analysis[j,11:13])/3) !=  stat_analysis[j,11:13][1]){
    # # 
    # # 
    # #        p_values_14[j] <- ttest_func(stat_analysis[j,3:5], stat_analysis[j,12:14])
    # # 
    # #      #}else{
    # # 
    # #     # }
    # #      #if ((sum(stat_analysis[j,2:4])/3) !=  stat_analysis[j,2:4][1] | (sum(stat_analysis[j,14:16])/3) !=  stat_analysis[j,14:16][1]){
    # # 
    # #        p_values_15[j] <- ttest_func(stat_analysis[j,3:5], stat_analysis[j,15:17])
    # #     # }else{
    # # 
    # #     # }
    # #      }
    # # 
    # # 
    # #  p_values_12 <- as.data.frame(p_values_12)
    # #  p_values_13 <- as.data.frame(p_values_13)
    # #  p_values_14 <- as.data.frame(p_values_14)
    # #  p_values_15 <- as.data.frame(p_values_15)
    # # 
    # # 
    # #    for (i in 2:sample_size){
    # # 
    # #      assign(paste0("p_values_1",i),
    # #             cbind(stat_analysis$common_col_for_merging,
    # #                   as.data.frame(get(paste0("p_values_1",i)))))
    # # 
    # #    }
    # # 
    # #  
    # #  test_num <- nrow(stat_analysis)
    # #  
    # #  for (j in 1:nrow(stat_analysis)){
    # #    
    # #    p_values_12 <- p_values_12[order(p_values_12$p_values_12),]
    # #    colnames(p_values_12)[1] <- "common_col_for_merging"
    # #    p_values_12["rank12"] <- 1:test_num
    # #    p_values_12[j,"test12"] <- (p_values_12[j,"rank12"]/test_num)*0.05
    # #    p_values_12[j,"test_bool"] <- p_values_12[j,"test12"] > p_values_12[j,"p_values_12"]
    # #  
    # #    
    # #    p_values_13 <- p_values_13[order(p_values_13$p_values_13),]
    # #    colnames(p_values_13)[1] <- "common_col_for_merging"
    # #    p_values_13["rank13"] <- 1:test_num
    # #    p_values_13[j,"test13"] <- (p_values_13[j,"rank13"]/test_num)*0.05
    # #    p_values_13[j,"test_bool"] <- p_values_13[j,"test13"] > p_values_13[j,"p_values_13"]
    # #    
    # #    
    # #    p_values_14 <- p_values_14[order(p_values_14$p_values_14),]
    # #    colnames(p_values_14)[1] <- "common_col_for_merging"
    # #    p_values_14["rank14"] <- 1:test_num
    # #    p_values_14[j,"test14"] <- (p_values_14[j,"rank14"]/test_num)*0.05
    # #    p_values_14[j,"test_bool"] <- p_values_14[j,"test14"] > p_values_14[j,"p_values_14"]
    # #    
    # #    
    # #    p_values_15 <- p_values_15[order(p_values_15$p_values_15),]
    # #    colnames(p_values_15)[1] <- "common_col_for_merging"
    # #    p_values_15["rank15"] <- 1:test_num
    # #    p_values_15[j,"test15"] <- (p_values_15[j,"rank15"]/test_num)*0.05
    # #    p_values_15[j,"test_bool"] <- p_values_15[j,"test15"] > p_values_15[j,"p_values_15"]
    # #    
    # #  }
    # # ### TODO: Extract where the first FALSE was generated to use for threshold of each comparison.
    # #  #p_thresholds <- c(2.100824e-02,2.381169e-02,3.119312e-02,2.567998e-02)
    # #  p_thresholds<- c(2.045695e-02,2.359024e-02,2.751729e-02,2.703404e-02)
    # #  col_sel_stat_analysis <- c("common_col_for_merging",
    # #                             #"modifications",
    # #                             #"accession",
    # #                             "isomericity",
    # #                             "isomeric_count",
    # #                             #"pool_id",
    # #                             "Pool",
    # #                             "A1-A2_Ratio",
    # #                             "A1-A3_Ratio",
    # #                             "A1-A3_Ratio",
    # #                             "A1-A4_Ratio",
    # #                             "A1-A5_Ratio",
    # #                             "exp_FC_A1/A2",
    # #                             "exp_FC_A1/A3",
    # #                             "exp_FC_A1/A4",
    # #                             "exp_FC_A1/A5")
    # #    
    # #  
    # 
    #  # for (i in 12:15){
    #  #   assign(paste0("complete_p_values_",i), 
    #  #           final_imputed_data %>% 
    #  #             select(col_sel_stat_analysis) %>%
    #  #             left_join(get(paste0("p_values_",i)),
    #  #                       by="common_col_for_merging") %>%
    #  #             mutate(ls_correctness = ifelse(test_bool == FALSE, 0, 1)))}
    #  #   
    #  #   assign(paste0("roc_data",i),compute_roc_curve(get(paste0("complete_p_values_",i)),flag = "Others",expected = 145 ))
    #  #   #assign(paste0("roc_data_unexpected",i),compute_roc_curve(get(paste0("complete_p_values_",i)),flag = "unexpected",expected = 100))
    #  #   #assign(paste0("roc_data_isoref",i),compute_roc_curve(get(paste0("complete_p_values_",i)),flag = "ISO-REF",expected = 37 ))
    #  #   
    #  # common_colnames <-c("common_col_for_merging", "p_values","rank","test","test_bool")
    #  # 
    #  # 
    #  # 
    #  # colnames(p_values_12) <- common_colnames
    #  # colnames(p_values_13) <- common_colnames
    #  # colnames(p_values_14) <- common_colnames
    #  # colnames(p_values_15) <- common_colnames
    #  # complete_pvalues_all <- rbind(p_values_12,p_values_13,p_values_14,p_values_15)
    # 
    #  
    #  ##TODO: make it more professional
    #  ggplot(roc_data12,aes(x=fdp,y=tpr)) + geom_line()
    #  
    #  
    #  num_wrong_local_pep_Proline_wo_FAIMS <- c("ARSRTPPSAPSQSR_3&11",
    #  "ARSRTPPSAPSQSR_3&13",
    #  "ATSLPSLDTPGELR_2",
    #  "FSDQAGPAIPTSNSYSK_14",
    #  "HTDDEMTGYVATR_12",
    #  "HTDDEMTGYVATR_7",
    #  "LMTGDTYTAHAGAK_8",
    #  "LPLTRSHNNFVAILDLPEGEHQYK_4",
    #  "NGSLKPGSSHR_9",
    #  "RLSSTSLASGHSVR_3&5",
    #  "RLSSTSLASGHSVR_4&5",
    #  "RLSSTSLASGHSVR_4&6",
    #  "RLSSTSLASGHSVR_5",
    #  "RLSSTSLASGHSVR_5&6",
    #  "RLSSTSLASGHSVR_6",
    #  "SFGSPNRAYTHQVVTR_1&10",
    #  "SNSTSSMSSGLPEQDR_1",
    #  "STGDPQGVIR_1",
    #  "STVASMMHR_2",
    #  "TAGTSFMMTPYVVTR_14",
    #  "TGMGSGSAGKEGGPFK_1",
    #  "TVSTSSQPEENVDR_4",
    #  "VIEDNEYTAR_8",
    #  "VSPSPTTYR_6",
    #  "VSPSPTTYR_7",
    #  "YATPQVIQAPGPR_3")
    #  
    #  num_wrong_local_pept_MQ_score_40 <- c("ARSRTPPSAPSQSR_3&13",
    #    "DIYSTDYYR_8",
    #    "ESKSSPRPTAEK_2",
    #    "FSDQAGPAIPTSNSYSK_14",
    #    "IQPAGNTSPR_7",
    #    "LPLTRSHNNFVAILDLPEGEHQYK_4",
    #    "LSYYEYDFER_3",
    #    "RLSSFVTK_7",
    #    "RLSSTSLASGHSVR_3&5",
    #    "RLSSTSLASGHSVR_3&6",
    #    "RLSSTSLASGHSVR_5",
    #    "RSMSPFRGPK_2",
    #    "SFGSPNRAYTHQVVTR_1&10",
    #    "STVASMMHR_2",
    #    "TAGTSFMMTPYVVTR_14",
    #    "TVSTSSQPEENVDR_4&5")
    #  
    #  png_width <- 1800
    #  png_height <- 1500
    #  
    #  venn.diagram(
    #    x = list(num_wrong_local_pep_Proline_wo_FAIMS, num_wrong_local_pept_MQ_score_40),
    #    #category.names = c("Proline_wo_FAIMS", "MQ_wo_FAIMS"),
    #    #fill = c("#00AFBB", "#D16103"),  # Specify custom colors for the sets
    #    alpha = 0.5,  # Set transparency level
    #    fontfamily = "sans",  # Specify font family
    #    fontface = "bold",  # Specify font face
    #    fontcolor = "black",  # Specify font color
    #    fontsize = 12,  # Specify font size
    #    cat.fontfamily = "sans",  # Specify category label font family
    #    cat.fontface = "bold",  # Specify category label font face
    #    cat.fontsize = 12,  # Specify category label font size
    #    cat.cex = 1.2,  # Specify category label expansion factor
    #    cex = 1,  # Specify overall expansion factor for the diagram
    #    filename = "venn.png",  # Optional: specify a filename to save the plot as an image
    #    width = png_width,
    #    height = png_height
    #  )
    #  
    #  
    #  
    #  
    #  merging_all_pvalues <- final_imputed_data %>% 
    #    select(col_sel_stat_analysis) %>%
    #    left_join(p_values_12,by="common_col_for_merging") %>%
    #    left_join(p_values_13,by="common_col_for_merging") %>%
    #    left_join(p_values_14,by="common_col_for_merging") %>%
    #    left_join(p_values_15,by="common_col_for_merging")
    #  
    #  merging_all_pvalues <- merging_all_pvalues %>%
    #    mutate(Pool=replace_na(Pool,"unexpected"))
    #  
    #  
    #  
    #  
    #  
    #    
    #  # 
    #  # 
    #  # 
    #  # 
    #  # 
    #  # #adjusted_p_values <- as.data.frame(p.adjust(p_values_A1_vs_A5, method = "BH", n = length(p_values_A1_vs_A5)))
    #  # adjusted_p_values <- lapply(stat_analysis, function(x) p.adjust (x, method = "BH", n = length(x)))
    #  # p_values_12[,"adj_pvalue12"] <- p.adjust(p_values_12$p_values_12, method = "BH")
    #  # 
    #  # all_adjusted_p_values <- cbind(unlist(adjusted_p_values[["p_values_12"]]),
    #  #                                unlist(adjusted_p_values[["p_values_13"]]),
    #  #                                unlist(adjusted_p_values[["p_values_14"]]),
    #  #                                unlist(adjusted_p_values[["p_values_15"]]))
    #  # colnames(all_adjusted_p_values) <- c("adj_pvalues_1_2","adj_pvalues_1_3","adj_pvalues_1_4","adj_pvalues_1_5")
    #  # 
    #  # 
    #  
    #  ##### p_values_for_all_ratio changed it later to adjusted!!!!
    #  df_volcano <- final_imputed_data %>% 
    #    select(!contains("mean_abundance")) %>%
    #    mutate(as.data.frame(p_values_for_all_ratio)) %>%
    #    #mutate(across(everything(), ~replace(., is.infinite(.), NA)))
    #    rowwise() %>%
    #    filter(across(where(is.numeric), ~!is.infinite(.)))
    #  
    #  
    #  pvalues_df <- merging_all_pvalues %>%
    #    pivot_longer(cols = starts_with("p_values"), names_to = "pvalues_col", values_to = "pvalues_value") %>%
    #    select(common_col_for_merging,Pool,isomericity,pvalues_col,pvalues_value)
    #  
    #  # Pivot the columns containing "Ratio"
    #  ratio_df <- merging_all_pvalues %>%
    #    pivot_longer(cols = ends_with("_Ratio"), names_to = "Ratio_col", values_to = "Ratio_value") %>%
    #    select(Ratio_col,Ratio_value)
    #  
    #  # Pivot the columns containing "exp_FC"
    #  exp_FC_df <- merging_all_pvalues %>%
    #    pivot_longer(cols = starts_with("exp_FC"), names_to = "exp_FC_col", values_to = "exp_FC_value") %>%
    #    select(exp_FC_col,exp_FC_value)
    #  
    #  volcano_final <- bind_cols(pvalues_df,ratio_df,exp_FC_df)
    # 
    #  actual_ratio <- c(2,10,20,100)
    #  actual_ratio <- data.frame(Ratio_col=unique(volcano_final$Ratio_col),log2(c(2,10,20,100)))
    #  
    #  log10_p_thresholds <- data.frame(ratio_col=unique(volcano_final$Ratio_col), -log10(p_thresholds))
    #  
    #  
    #  #### DETERMINE THRESHOLD OF EACH P VALUE COMPARISON BASED on 
    #  ####    WHERE WE CAN SEE THE CHANGE FROM TRUE to FALSE
    #  
    #  ## THIS WILL BE ADDED AS AN EXTRA LINE FOR EACH COMPARISION 
    #  ##     WHEN I MERGED ALL VOLCANO PLOTS INTO ONE
    #  
    #  volcano_final1 <- volcano_final %>% select(-Ratio_value) %>%
    #  mutate(isomericity = ifelse(is.na(isomericity), "False Positive", isomericity)) %>%
    #    unite(Pool_new, Pool, isomericity,sep = "_",remove = FALSE) %>%
    #    unite('new_col_coloring',Pool_new,Ratio_col,sep = "_",remove = FALSE) %>%
    #    mutate(new_col_coloring = if_else(grepl("ISO-REF", new_col_coloring), "ISO-REF", new_col_coloring)) %>%
    #    mutate(new_col_coloring = if_else(grepl("unexpected", new_col_coloring), "unexpected", new_col_coloring))
    #  
    #  #write.table(volcano_final1,file = paste0(file_path,"/","exp2wo_FAIMS_DDA_Proline_volcano_plot_data.tsv"),sep = "\t",row.names = F)
    #  
    #  ggplot(volcano_final1, aes(x = log2(exp_FC_value), y = -log10(pvalues_value))) +
    #    geom_point(aes(shape = Pool_new, color = new_col_coloring), size = 2.5) +
    #    #geom_line(aes(color = new_col_coloring), size = 1) +  # Add color aesthetic to geom_line()
    #    scale_color_manual(values = c("ISO-REF" = "#000000", "unexpected" = "#999999",
    #                                  "Others_multi_A1-A2_Ratio" = "#CC79A7",
    #                                  "Others_mono_A1-A2_Ratio" = "#CC79A7",
    #                                  "Others_multi_A1-A3_Ratio" = "#E69F00",
    #                                  "Others_mono_A1-A3_Ratio" = "#E69F00",
    #                                  "Others_multi_A1-A4_Ratio" = "#56B4E9",
    #                                  "Others_mono_A1-A4_Ratio" = "#56B4E9",
    #                                  "Others_multi_A1-A5_Ratio" = "#009E73",
    #                                  "Others_mono_A1-A5_Ratio" = "#009E73"),
    #                       labels = c('Non-variant', 'Variant non-isomeric A1 vs A2',
    #                                  'Variant non-isomeric A1 vs A3',
    #                                  'Variant non-isomeric A1 vs A4',
    #                                  'Variant non-isomeric A1 vs A5',
    #                                  'Variant isomeric A1 vs A2',
    #                                  'Variant isomeric A1 vs A3',
    #                                  'Variant isomeric A1 vs A4',
    #                                  'Variant isomeric A1 vs A5',
    #                                  'Unexpected')) +
    #    scale_shape_manual(values = c(16, 15, 12, 17),
    #                       labels = c('Non-variant', 'Variant isomeric', 'Variant non-isomeric', 'Unexpected')) +
    #    scale_y_continuous(limits = c(0, 7.2), breaks = seq(0, 7.2, by = 0.8)) +
    #    scale_x_continuous(limits = c(-7,7)) +
    #    #scale_y_continuous(breaks = seq(0, max(-log10(volcano_final1$pvalues_value)), length.out = 21)) +
    #    theme_bw() +
    #    theme(legend.text = element_text(size = 15),
    #          axis.title.x = element_text(size = 15),
    #          axis.title.y = element_text(size = 15),
    #          plot.title = element_text(size = 30),
    #          legend.title = element_text(size = 15),
    #          axis.text.x = element_text(size = 15),
    #          axis.title = element_text(size = 15),
    #          axis.text.y = element_text(size = 15)) +
    #    expand_limits(x = 0, y = 0) +
    #    geom_vline(data = actual_ratio, aes(xintercept = actual_ratio$log2.c.2..10..20..100..),color=c("#CC79A7","#E69F00","#56B4E9","#009E73"), size = 1, show.legend = FALSE) +
    #    geom_hline(data = log10_p_thresholds, aes(yintercept = log10_p_thresholds$X.log10.p_thresholds.),color=c("#CC79A7","#E69F00","#56B4E9","#009E73"), size = 1, linetype = 2, show.legend = FALSE)+ 
    #    labs(title = "Experiment 2 - DDA no FAIMS processed by Proline", color = "Classes", shape="Type")
    #  
    #    #coord_fixed()
    #  
    #  final_imputed_data %>% 
    #    pivot_longer(cols = starts_with("mean_"),
    #                 names_to = "Abundance_col",
    #                 values_to = "Abundance_value") %>%
    #  ggplot(aes(x = Abundance_col, y = common_col_for_merging, fill = Abundance_value)) +
    #    geom_tile(color = "white") +
    #    scale_fill_viridis() +
    #    labs(x = "Condition", y = "Peptide", fill = "Abundance") +
    #    ggtitle("Heatmap")
    #  
    #  
    #  
    #  ggplot(df_volcano, aes(x=df_volcano$`exp_FC_A1/A4`, y=-log10(pvalues_1_4),color=Pool)) + 
    #    geom_point(aes(Pool)) + geom_hline(yintercept = pvalue_threshold, size=1) + 
    #    geom_vline(xintercept = actual_ratio[3], size=1)
    #  #### Distribution of experimental quantitative ratios with two different pools
    #  
    # 
    #  
    #  library(ggplot2)
    #  
    #  ## This can be used as an object name
    #  df_for_figure <- final_imputed_data %>% 
    #    select(contains(col_sel_stat_analysis)) %>%
    #    pivot_longer(cols = col_sel_stat_analysis[3:11], names_to = "Theo_Exp", values_to = "Ratios",values_drop_na = T) %>%
    #    #select(!contains(c("common_col_for_merging"))) %>%
    #    filter(grepl("exp_FC",Theo_Exp)) %>%
    #    filter_at(vars(Pool), all_vars(!is.na(.))) %>%
    #    filter_at(vars(Ratios), all_vars(!is.infinite(.)))
    #  
    #  p1 <- gg_density(data_set = df_for_figure, 
    #                   x_df = df_for_figure$Ratios,
    #                   fill_df = df_for_figure$Theo_Exp,
    #                   color_df = df_for_figure$Pool,
    #                   header="Distribution of experimental quantitative ratios with two different pools",
    #                   facet_df = "Theo_Exp",
    #                   x_lab = "log10(Ratios)",
    #                   color_lab= "Pool",
    #                   fill_lab = "Ratios")
    #  
    #  
    #  #### Distribution of mean abundance of every sample with two different pools
    #  
    #  figure_with_mean_abundance <- final_imputed_data %>%
    #    select(starts_with("mean_") | contains(c("Pool", "common_col_for_merging"))) %>%
    #    pivot_longer(cols = starts_with("mean"),
    #                 names_to = "Theo_Exp",
    #                 values_to = "Mean_abundance",
    #                 values_drop_na = T) %>%
    #    filter_at(vars(Pool), all_vars(!is.na(.)))
    #  
    #  
    #  p2 <- gg_density(data_set = figure_with_mean_abundance, 
    #                   x_df = figure_with_mean_abundance$Mean_abundance,
    #                   fill_df = figure_with_mean_abundance$Theo_Exp,
    #                   color_df = figure_with_mean_abundance$Pool,
    #                   header="Distribution of mean abundance of every sample with two different pools",
    #                   facet_df = "Theo_Exp",
    #                   x_lab = "log10(Ratios)",
    #                   color_lab= "Pool",
    #                   fill_lab = "Ratios")
    #  
    #  # final_imputed_normalized_data %>%
    #  #   select(starts_with("mean_") | contains(c("Pool", "common_col_for_merging"))) %>%  
    #  #   filter_at(vars(Pool), all_vars(!is.na(.))) %>%
    #  
    #  #### MANUAL PLOTTING   
    #  #
    #  #   ggplot(aes(x =log10(mean_abundances_A1) , y = log10(mean_abundances_A5), color=Pool)) +
    #  #   geom_point()+
    #  #   #facet_wrap(vars(df_for_figure$common_col_for_merging))  + 
    #  #   geom_smooth(formula = y ~ x,method = "loess", colour = "green", fill = "green") +
    #  #   theme_minimal() +
    #  #   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
    #  #         axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
    #  #         plot.title = element_text(size=20),
    #  #         legend.title=element_text(size=15),
    #  #         axis.text=element_text(size=15),
    #  #         axis.title=element_text(size=15)
    #  #   )
    #  #### PLOTTING WITH FUNCTION
    #  #
    #  # gg_density_mean_abun <- function(data_set = final_imputed_normalized_data,
    #  #                                  x_df = final_imputed_normalized_data$mean_abundances_A1,
    #  #                                  y_df = final_imputed_normalized_data$mean_abundances_A5,
    #  #                                  color_df =final_imputed_normalized_data$Pool,
    #  #                                  #header,
    #  #                                  facet_df= "common_col_for_merging")
    #  
    #  ### BOX-PLOT: Experimental Quantity Ratio of Synthetic Peptides  
    #  
    #  p3 <- gg_boxplt_exp_ratio(data_set = df_for_figure, 
    #                            x_df = df_for_figure$Theo_Exp,
    #                            y_df = df_for_figure$Ratios,
    #                            fill_df = df_for_figure$Pool,
    #                            header="Experimental Quantity Ratio of Synthetic Peptides",
    #                            x_lab="Sample Names",
    #                            y_lab="Abundance Ratios",
    #                            fill_lab = "Pool")
    #  
    #  ### HALF-BOX-PLOT & HALF-SCATTER-PLOT: Experimental Quantity Ratio of Synthetic Peptides  
    #  library(gghalves)
    #  
    #  p4 <- gg_half_boxplt_exp_ratio(data_set = df_for_figure, 
    #                                 x_df = df_for_figure$Theo_Exp,
    #                                 y_df = df_for_figure$Ratios,
    #                                 fill_df = df_for_figure$Pool,
    #                                 header="Experimental Quantity Ratio of Synthetic Peptides",
    #                                 x_lab="Sample Names",
    #                                 y_lab="Abundance Ratios",
    #                                 fill_lab = "Pool")
    #  
    #  
    #  
    #  ### VIOLIN-PLOT: Experimental Quantity Ratio of Synthetic Peptides   
    #  
    #  ### TODO: fix y scaling without trimming 
    #  p5 <- gg_violin_exp_ratio(data_set = df_for_figure, 
    #                            x_df = df_for_figure$Theo_Exp,
    #                            y_df = df_for_figure$Ratios,
    #                            fill_df = df_for_figure$Pool,
    #                            header="Experimental Quantity Ratio of Synthetic Peptides",
    #                            x_lab="Sample Names",
    #                            y_lab="Abundance Ratios",
    #                            fill_lab = "Pool",
    #                            trim=TRUE)
    #  
    #  ### SAVE ALL PLOT AUTOMATICALLY (without giving a custom file name)
    #  sapply(1:5,function(x) ggsave(filename = paste0("p",x,".tiff"),
    #                                width = 50, height = 40, 
    #                                path = file_path,
    #                                units = "cm",
    #                                get(paste0("p",x)),
    #                                device = "tiff", #".svg"
    #  ))
    #  
    #  ### SAVE PLOTS SEPARATELY (with custom file name)
    #  # ggsave(filename = "Experimental_Quantity_Ratio_of_Synthetic_Peptides",
    #  #        path = "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/",
    #  #        plot = p5,
    #  #        device = "tiff")
    #  
    #  
    #  # library(limma)
    #  # sample_size <- 5
    #  # # Do t-test
    #  # # The code below does t-test for each row. Because of that, multiple test correction (like BH, Bonferoni)
    #  # 
    #  # ## 1 ## TIDY DATA FOR T-TEST
    #  # df_for_ttest<- final_imputed_normalized_data %>%
    #  #   select(starts_with("mean_") | contains(c("Pool", "common_col_for_merging"))) %>%  
    #  #   filter_at(vars(Pool), all_vars(!is.na(.)))
    #  # 
    #  # ## 2 ## FUNCTION FOR T-TEST ALL CONC. ACROSS A1
    #  # ttest_func <- function(A1,rest){
    #  #   t.test(A1,rest,alternative = "two.sided", var.equal = TRUE)}
    #  # 
    #  # p_values_for_all_ratio <- data.frame(1:dim(df_for_ttest)[1])
    #  # 
    #  # ## 3 ## AUTOMIZED T-TEST ALL CONC. ACROSS A1
    #  # for(k in 2:sample_size){
    #  #   
    #  #   assign(paste0("t_test_res_A1_to_A",k) , 
    #  #          lapply(df_for_ttest$mean_abundances_A1, ttest_func,
    #  #                 rest=df_for_ttest[,k]))
    #  #   
    #  #   assign(paste0("p_values_for_ratio",k),
    #  #          lapply(get(paste0("t_test_res_A1_to_A",k)), function (x) x[c('p.value')]))
    #  #   
    #  #   p_values_for_all_ratio <- cbind(p_values_for_all_ratio,
    #  #                                   as.data.frame(unlist(get(paste0("p_values_for_ratio",k)))))
    #  # }
    #  # 
    #  # 
    #  # 
    #  # #adjusted_p_values <- as.data.frame(p.adjust(p_values_A1_vs_A5, method = "BH", n = length(p_values_A1_vs_A5)))
    #  # adjusted_p_values <- lapply(p_values_for_all_ratio[,-1], function(x) p.adjust (x, method = "BH", n = length(x)))
    #  # 
    #  # 
    #  # 
    #  # 
    #  # 
    #  # #plot(correct_identifed_peps$mean_abundances_log10_A1,correct_identifed_peps$mean_abundances_log10_A5, pch = 16, col = "blue")
    #  # #abline(h = mean(na.omit(correct_identifed_peps$mean_abundances_log10_A1)) - mean(na.omit(correct_identifed_peps$mean_abundances_log10_A5)), col = "red")
    #  # 
    #  # boxplot(correct_identifed_peps$mean_abundances_log10_A1,correct_identifed_peps$mean_abundances_log10_A5)
    #  # ggplot(correct_identifed_peps, aes(x = mean_abundances_log10_A1, y = mean_abundances_log10_A5)) +
    #  #   geom_point(aes(color = "red", size = 5))# +
    #  # #scale_color_discrete(name = "") +
    #  # #scale_size_discrete(name = "Size")
    #  # # 
    #  # # Take -log10() of results
    #  # ggplot(log_10_filtered_abundances_rowMeans_A1_A5,aes(x=log_10_filtered_abundances_rowMeans_A1_A5$mean_abundances_log10_A1,
    #  #                                                      y =log_10_filtered_abundances_rowMeans_A1_A5$mean_abundances_log10_A5,
    #  # )) +
    #  #   geom_point()+
    #  #   #facet_wrap(vars(df_for_figure$common_col_for_merging))  + 
    #  #   #geom_smooth(method = "lm", colour = "green", fill = "green") +
    #  #   theme_light()
    #  # 
    #  # library(gginference)
    #  # ggttest(t.test(na.omit(correct_identifed_peps$mean_abundances_log10_A1),na.omit(correct_identifed_peps$mean_abundances_log10_A5), alternative = "two.sided", var.equal = TRUE))
    #  # 
    #  # 
    #  # 
    #  # df_for_figure_exp <- correct_identifed_peps[,c(11,74:77)]
    #  # df_for_figure <- correct_identifed_peps[,c(11:15,74:77)]
    #  # melt_df_for_figure_exp <- melt(df_for_figure_exp)
    #  # 
    #  # df_for_figure_theo <- correct_identifed_peps[,c(11:15)]
    #  # melt_df_for_figure_theo <- melt(df_for_figure_theo)
    #  # 
    #  # p1 <- ggplot(melt_df_for_figure_theo,aes(x =melt_df_for_figure_theo$variable , y =log2(melt_df_for_figure_theo$value),fill = melt_df_for_figure_theo$Pool)  ) +
    #  #   geom_boxplot() +
    #  #   theme_light() + #scale_y_continuous(limits = c(-60, 60),breaks = seq(-60, 60, by = 20)) +
    #  #   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
    #  #         axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
    #  #         plot.title = element_text(size=20),
    #  #         legend.title=element_text(size=15),
    #  #         axis.text=element_text(size=15),
    #  #         axis.title=element_text(size=15)
    #  #   ) +   stat_boxplot(geom = "errorbar") + ggtitle("Theoretical Quantity Ratio of Synthetic Peptides") +
    #  #   scale_x_discrete(labels=c("A1/A2","A1/A2","A1/A4","A1/A5")) +
    #  #   labs(x="Sample Names",y="Abundance Ratios", color="Pool Names")+
    #  #   scale_fill_brewer(palette="Set1")
    #  # 
    #  # p2 <-ggplot(melt_df_for_figure_exp,aes(x =melt_df_for_figure_exp$variable , y =log2(melt_df_for_figure_exp$value),fill = melt_df_for_figure_exp$Pool)  ) +
    #  #   geom_boxplot() +
    #  #   theme_light() +
    #  #   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
    #  #         axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
    #  #         plot.title = element_text(size=20),
    #  #         legend.title=element_text(size=15),
    #  #         axis.text=element_text(size=15),
    #  #         axis.title=element_text(size=15)
    #  #   ) +   stat_boxplot(geom = "errorbar") + ggtitle("Experimental  Quantity Ratio of Synthetic Peptides") +
    #  #   scale_x_discrete(labels=c("A1/A2","A1/A2","A1/A4","A1/A5")) +
    #  #   labs(x="Sample Names",y="Abundance Ratios", color="Pool Names") +
    #  #   scale_fill_brewer(palette="Set1")
    #  # 
    #  # p3 <-ggplot(melt_df_for_figure_exp,aes(x =melt_df_for_figure_exp$variable , y =melt_df_for_figure_exp$value,color = melt_df_for_figure_exp$Pool)  ) +
    #  #   geom_point() +
    #  #   theme_light() +
    #  #   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
    #  #         axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
    #  #         plot.title = element_text(size=20),
    #  #         legend.title=element_text(size=15),
    #  #         axis.text=element_text(size=15),
    #  #         axis.title=element_text(size=15)
    #  #   ) +   ggtitle("Experimental Quantity Ratio of Synthetic Peptides") +
    #  #   scale_x_discrete(labels=c("A1/A2","A1/A3","A1/A4","A1/A5")) +
    #  #   labs(x="Sample Names",y="Abundance Ratios", color="Pool Names") +
    #  #   scale_color_brewer(palette="Set1")
    #  # 
    #  # p4 <-ggplot(melt_df_for_figure_theo,aes(x =melt_df_for_figure_theo$variable , y =melt_df_for_figure_theo$value,color = melt_df_for_figure_theo$Pool)  ) +
    #  #   geom_point() +
    #  #   theme_light() +
    #  #   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
    #  #         axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
    #  #         plot.title = element_text(size=20),
    #  #         legend.title=element_text(size=15),
    #  #         axis.text=element_text(size=15),
    #  #         axis.title=element_text(size=15),
    #  #   ) +   ggtitle("Theoretical Quantity Ratio of Synthetic Peptides") +
    #  #   scale_x_discrete(labels=c("A1/A2","A1/A3","A1/A4","A1/A5")) +
    #  #   labs(x="Sample Names",y="Abundance Ratios", color="Pool Names") +
    #  #   scale_color_brewer(palette="Set1")
    #  # 
    #  # 
    #  # library(ggpubr)
    #  # ggarrange(p1, p2, common.legend = TRUE, legend="right")
    #  # ggarrange(p4, p3, common.legend = TRUE, legend="right")
    #  # 
    #  # 
    #  # iso_exp <- df_for_figure_exp %>% filter(grepl("ISO-REF",Pool))
    #  # melt_iso_theo <- melt(iso_exp)
    #  # ggplot(melt_iso_theo,aes(x =melt_iso_theo$variable , y =log2(melt_iso_theo$value),color = melt_iso_theo$Pool)  ) +
    #  #   geom_violin()  + #geom_point() +
    #  #   theme_light() +
    #  #   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
    #  #         axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
    #  #         plot.title = element_text(size=20),
    #  #         legend.title=element_text(size=15),
    #  #         axis.text=element_text(size=15),
    #  #         axis.title=element_text(size=15)
    #  #   ) +   ggtitle("Experimental Quantity Ratio of Synthetic Peptides") +
    #  #   scale_x_discrete(labels=c("A1/A2","A1/A3","A1/A4","A1/A5")) +
    #  #   labs(x="Sample Names",y="Abundance Ratios", color="Pool Names") +
    #  #   scale_color_brewer(palette="Set1")
    #  # 
    #  # 
    #  # other_exp <- df_for_figure_exp %>% filter(grepl("Others",Pool))
    #  # melt_other_exp <- melt(other_exp)
    #  # ggplot(melt_other_exp,aes(x =melt_other_exp$variable , y =log2(melt_other_exp$value),color = melt_other_exp$Pool)  ) +
    #  #   geom_violin() + #geom_point() +
    #  #   theme_light() +
    #  #   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
    #  #         axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
    #  #         plot.title = element_text(size=20),
    #  #         legend.title=element_text(size=15),
    #  #         axis.text=element_text(size=15),
    #  #         axis.title=element_text(size=15)
    #  #   ) +   ggtitle("Experimental Quantity Ratio of Synthetic Peptides") +
    #  #   scale_x_discrete(labels=c("A1/A2","A1/A3","A1/A4","A1/A5")) +
    #  #   labs(x="Sample Names",y="Abundance Ratios", color="Pool Names") +
    #  #   scale_color_brewer(palette="Set1")
    #  
    #  
    
  
  
  
}

