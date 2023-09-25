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

final_proline_pep_quant_analysis_syn <- function(file_path,
                                                 file_name,
                                                 sheet_name,
                                                 theo_file_path,
                                                 theo_file_name,
                                                 sheet_theo_name,
                                                 background_species,
                                                 selected_spcies,
                                                 exp_id,
                                                 loc_filter_opt,
                                                 loc_filter,
                                                 exp_design,
                                                 fdr_threshold,
                                                 acquisiton_type,
                                                 software_name,
                                                 test_type,
                                                 num_reps,
                                                 actual_ratio,
                                                 subtitle){
    ### Source code was taken from here: https://rdrr.io/github/singjc/mstools/src/R/getModificationPosition.R
    ## The code was modified based on what I want and based on software input tyoe (DIANN-and MaxQuant)
    #source("D:/dev/Desktop_copy/PHD/data_analysis/scripts/getModificationPositionMQ_func_edit_v1_1634.R")
    source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/get_modification_func/getModificationPosition_general_change_condition_current_mod_sequence_MQ_Spectronaut.R")
    source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/ggplot/ggplot_functions.R")
    source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/roc_curve/roc_curve_generation_proline_edit.R")
    
    
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
    
    pep_list_w_theo <- read.xlsx(paste0(theo_file_path, theo_file_name), sheet = sheet_theo_name)
    pep_list_w_theo_quant <- pep_list_w_theo[,-1]
    
    common_col_theo_quant <- as.data.frame(paste(pep_list_w_theo_quant$Phosphopeptide.sequence,
                                                 pep_list_w_theo_quant$modified.position.in.peptide, sep = "_"))
    colnames(common_col_theo_quant) <- "pep_with_pos"
    pep_list_w_theo_quant_new <- cbind(common_col_theo_quant,pep_list_w_theo_quant)
    
    
    #################################################
    pep_list_w_theo_unique <- pep_list_w_theo %>% 
      select(Phospopeptide.sequence,Pool) %>%
      distinct(Phospopeptide.sequence,.keep_all = TRUE) %>%
      rename(Sequence = Phospopeptide.sequence) %>%
      rename(Pool_for_seq_merge=Pool)

    ecoli_seq <- quant_peptides %>% 
      select(Sequence,Modifications,Proteins) %>%
      filter(!grepl(selected_spcies, Proteins) & !grepl("CON__", Proteins)) %>%
      #pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
      #filter(grepl(background_species,Proteins)) %>%
      #distinct(Sequence,.keep_all = T) %>%
      mutate(species=background_species) 
    
    ecoli_seq_dist <- quant_peptides %>% 
      select(Sequence,Modifications,Proteins) %>%
      filter(!grepl(selected_spcies, Proteins) & !grepl("CON__", Proteins)) %>%
      #pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
      #filter(grepl(background_species,Proteins)) %>%
      distinct(Sequence,.keep_all = T) %>%
      mutate(species=background_species) 
    
    
    all_seq <- quant_peptides %>% 
      filter(grepl(selected_spcies,Proteins) & 
               grepl("Phospho",Modifications)) %>%
      mutate(species=selected_spcies) %>%
      full_join(pep_list_w_theo_unique,by="Sequence") %>%
      mutate_at("Pool_for_seq_merge", ~replace_na(.,"Unexpected")) %>%
      mutate(Pool_for_seq_merge= ifelse(is.na(species),"missing",Pool_for_seq_merge)) %>%
      filter(!grepl("Unexpected",Pool_for_seq_merge)) %>%
      bind_rows(ecoli_seq) %>% 
      mutate(Pool_for_seq_merge= ifelse(is.na(Pool_for_seq_merge),background_species,Pool_for_seq_merge))
    
    all_seq_syn <- all_seq %>%
      select(Sequence,Modifications,Proteins,Pool_for_seq_merge) %>% 
      distinct(Sequence, .keep_all = T) %>%
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
    
    write.table(all_seq_syn, file=paste0(file_path,"outputs_with_new_script/Experiment2",software_name,"_number_of_unique_sequence_for_each_species.txt"),sep = "\t",col.names = T,row.names = F)
    #################################################
    
    
    ### IMPUTATION
    
    # Calculate 1 percent quantile of each sample
    
    imputed_values <- quant_peptides %>% 
        group_by(Experiment) %>% 
        summarise(first_quantile=quantile(Intensity,probs=0.01,na.rm=TRUE))
    
    imputed_values_vec <- as.vector(imputed_values$first_quantile)
    # # Impute missing values
    # 
    # abundances_for_impute <- quant_peptides %>% select(Experiment,Intensity) %>%
    #     group_by(Experiment) %>%
    #     mutate(imputed_intensity=case_when(grepl("A1-R1",Experiment) ~imputed_values_vec[1],
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
    # 
    ### LOCALIZATION THRESHOLD ###
    # quant_phospho <- quant_peptides %>%
    #   #select(-Experiment,Intensity) %>% 
    #   #mutate(abundances_for_impute) %>% 
    #   #filter(grepl(9606,`Taxonomy IDs`))%>%
    #   filter(grepl(selected_spcies,Proteins)) %>% 
    #   filter(grepl("Phospho",Modifications)) %>% drop_na(`Modified sequence`)
    # 
    quant_phospho <- all_seq %>% 
      mutate(
        extracted_values = sapply(str_extract_all(`Phospho (STY) Probabilities`, "\\(\\d+(\\.\\d+)?\\)"), function(x) {
          values <- as.numeric(str_extract_all(x, "\\d+(\\.\\d+)?"))
          if (length(values) == 0 || all(is.na(values))) NA else max(values)
        })
        
        #select(Sequence,Modifications,`Modified sequence`,Proteins,`Phospho (STY) Probabilities`) %>%
        #drop_na(`Phospho (STY) Probabilities`) %>%
        # mutate(
        #   extracted_values = sapply(str_extract_all(`Phospho (STY) Probabilities`, "\\(\\d+\\.\\d+\\)"), function(x) {
        #     values <- as.numeric(str_extract_all(x, "\\d+")) #.\\d+
        #     if (length(values) == 0) NA else max(values)
        #   })
      )%>% relocate(extracted_values,.after = Sequence) %>% #drop_na(extracted_values) %>%
      filter(grepl(selected_spcies,Proteins) & grepl("Phospho",Modifications)) %>%
      drop_na(`Modified sequence`)
   

    df2 <- apply(X = as.data.frame(quant_phospho[,"Modified sequence"]),1,function(x){getModificationPosition_general(mod_seq = x,software_name = "MQ_214")})
    
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
        mutate(pep_with_pos = paste(value.y,value.x, sep = "_")) %>%
        select(!value.y) %>%
        rename("mods_id" = "name.x",
               "phospho_positions" ="value.x",
               "all_mods_with_mod_type" = "mods.x")
    
    
    #getModificationPosition_MQ(mod_seq = id_syn_pep$`Modified sequence`,F)
    
    final_results_with_common_col <- mutate(result_with_common_col,quant_phospho)
    
    
    barplt_df <- final_results_with_common_col %>%
        select(Sequence,Experiment,Intensity,pep_with_pos,Proteins,extracted_values) %>%
        separate(Experiment, into = c("Exp_id","Sample_id","Rep_id"),sep = "-",remove = F) %>%
        #mutate(sample_rep_id_seq = paste(pep_with_pos, Sample_id,Rep_id, sep = "_")) %>%
        group_by(pep_with_pos,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
        slice(which.max(Intensity)) %>%
        ungroup()
    
    barplt_phospho_seq <- final_results_with_common_col %>%
      select(Sequence,Experiment,Intensity,Proteins,extracted_values) %>%
      separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F) %>%
      mutate(sample_rep_id_seq = paste(Sequence, Sample_id,Rep_id, sep = "_")) %>%
      group_by(sample_rep_id_seq,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
      slice(which.max(Intensity)) %>%
      ungroup() %>% 
      mutate(Software_name=software_name) %>%
      mutate(Acquisition_type=acquisiton_type)
    
    write.table(barplt_phospho_seq, file = paste0(file_path,"Number_of_human_phospho_sequences_",
                                                  software_name,"_Experiment",exp_id,".txt"),
                sep = "\t",row.names = F)
    
    ####### ADDITIONAL PLOT TO DISPLAY MISSING and UNEXPECTED PEPTIDES ########
    site_prob <- barplt_df %>% 
      select(extracted_values, pep_with_pos,Experiment) %>% 
      pivot_wider(names_from = "Experiment",
                  values_from = "extracted_values") %>% 
      relocate(exp_design) %>%
      rowwise() %>%
      mutate(max_value = max(c_across(all_of(exp_design)), na.rm = TRUE)) %>%
      select(!exp_design)
    #mutate(max_value=ifelse(!is.numeric(max_value),NA,max_value))
    #filter(!grepl(-Inf,max_value))
    
     df_merge_syn <- barplt_df %>%
      select(pep_with_pos,Experiment,Intensity, Proteins) %>% 
      pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
      relocate(exp_design) %>%
      full_join(site_prob, by="pep_with_pos") %>%
      full_join(pep_list_w_theo_quant_new,by="pep_with_pos") %>% 
      mutate_at("Pool", ~replace_na(.,"Unexpected")) %>%
      mutate(Pool= ifelse(is.na(Proteins),"missing",Pool)) %>%
      select(pep_with_pos,starts_with(exp_design),Pool) %>%
      mutate(soft_name=software_name, ion_mobility=acquisiton_type)
    
    p11 <- gg_barplt_id_pep_count(data_set = df_merge_syn,
                                  x_df = df_merge_syn$Pool,
                                  fill_df = df_merge_syn$Pool,
                                  ymax = 20000,
                                  header = "Total number of quantified phospho-site across each sample",
                                  caption_lab = "NA values are removed. \n No filtering based on localization",
                                  x_lab = "Sample id",
                                  fill_lab =  "Sample id",
                                  y_lab = "Number of identified peptides",
                                  subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
    
    write.table(df_merge_syn,file = paste0(file_path,"outputs_with_new_script/Count_of_missing_unexpected_correct_phospho-sites_",
                                           software_name,"_Experiment",exp_id,".txt"),
                sep = "\t",row.names = F)
    #############################################################################
### COMPARED TO pd AND proline, MQ DOES NOT HAVE ANY QUERY THAT CONTAINS BACKGROUND SPECIES.
    ## THUS, WE FILTERED ECOLI DATA USING NOT SELECTING HUMAN AND CONTAMINANTS
    ## EXAMPLE ID WAS CONFRIMED ON UNIPROT = P0AAX3 -> YBIJ_ECOLI
    
    barplt_df_ecoli <- quant_peptides %>% 
        filter(!grepl(selected_spcies, Proteins) & !grepl("CON__", Proteins)) %>%
        select(Sequence,Experiment,Intensity,Proteins) %>%
        separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F) %>%
        mutate(sample_rep_id_seq = paste(Sequence, Sample_id,Rep_id, sep = "_")) %>%
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
   
    if(loc_filter_opt == TRUE){
      
      barplt_df_wide <- barplt_df %>%  ## If you select "charge" column, it will bring multiple rows for one seq
        select(Sequence,Experiment,Intensity,pep_with_pos,Proteins) %>%
        pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
        relocate(exp_design,.after = where(is.character)) %>%
        left_join(site_prob, by="pep_with_pos") %>%
        mutate(species=selected_spcies) %>%
        filter(max_value >= loc_filter) %>%
        relocate(pep_with_pos, .after = Sequence) 
      #relocate(exp_design,.after = "Proteins")
    }else{
      barplt_df_wide <- barplt_df %>%  ## If you select "charge" column, it will bring multiple rows for one seq
        select(Sequence,Experiment,Intensity,pep_with_pos,Proteins) %>%
        pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
        left_join(site_prob, by="pep_with_pos") %>%
        mutate(species=selected_spcies) %>%
        #filter(max_value >= loc_filter) %>%
        relocate(pep_with_pos, .after = Sequence) 
    }
   
    barplt_df_ecoli_wide <- barplt_df_ecoli %>% 
        select(Sequence,Experiment, Proteins,Intensity) %>%
        pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
        mutate(species=background_species) %>% relocate(exp_design,.after = "Proteins")
    extract_df <- barplt_df_wide %>% full_join(pep_list_w_theo_quant_new,by="pep_with_pos")
    
    write.table(extract_df,file=paste0(file_path,"outputs_with_new_script/Maxquant_transposed_data_merge_with_correct_pep_list.txt"),sep = "\t",col.names = T,row.names = F)
    
    # Nothing is changed
    filtered_abundances<-barplt_df_wide[rowSums(!is.na(select(barplt_df_wide,starts_with(exp_design))))>0,]
    filtered_abundances_ecoli <-barplt_df_ecoli_wide[rowSums(!is.na(select(barplt_df_ecoli_wide,starts_with(exp_design))))>0,]
    
    
    library(kableExtra)
    
    na_phospho_selected <- apply(X = is.na(filtered_abundances %>% select(Sequence, Proteins,starts_with("E2"))), MARGIN = 2, FUN = sum)
    na_ecoli <- apply(X = is.na(filtered_abundances_ecoli %>% select(Sequence,Proteins,starts_with("E2"))), MARGIN = 2, FUN = sum)
    
    na_table <- bind_rows(na_phospho_selected,na_ecoli)
    na_table$species <- c(selected_spcies,background_species)
    na_table$total <- c(dim(filtered_abundances)[1],dim(filtered_abundances_ecoli)[1])
    
    na_table %>% select(-Sequence) %>%
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
    
    ## REMOVE SEQUENCE COLUMN AFTER NA TABLE
    filtered_abundances <- filtered_abundances %>% select(!Sequence)
    
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
    
    quant_peptides_ECOLI_density_plot <- filtered_abundances_ecoli %>%
        select(!starts_with("E")) %>%
        bind_cols(abundances_ecoli_rowMeans) %>%
        tibble() %>%
        rename_with(~ paste0("mean abundance",1:5), matches("^row")) %>%
        pivot_longer(cols = starts_with("mean"), 
                     values_to = "Intensity",
                     names_to = "sample_ids",
                     values_drop_na = T) %>%
        mutate(sample_id_seq = paste(Sequence, sample_ids, sep = "_"))
    
    ecoli_density_plot<- quant_peptides_ECOLI_density_plot %>%
        select(contains(c("sample_ids","Intensity","species"))) 
    
    
    colnames(abundances_rowMeans) <- paste0("mean abundance",1:sample_size)
    
    #### MEAN ABUNDANCE RATIO WITH  DENSITY PLOT ####
    ### BEFORE IMPUTATION ###
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
        #separate(accession, into = c("prot_id","species","position"),sep = "_",remove = F)
    
    
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
        rename_with(~ paste0("pep_with_pos"), matches("^Seq")) %>%
        bind_rows(filtered_abundances) 
    
    
    
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
    
    final_imputed_data_syn <- final_imputed_data %>% filter(grepl(selected_spcies,species))
    
    final_imputed_data_ecoli <- final_imputed_data  %>% filter(!grepl(selected_spcies, species))
    
    df_merge <- final_imputed_data_syn %>%
        left_join(pep_list_w_theo_quant_new,by="pep_with_pos") %>% 
        mutate_at("Pool", ~replace_na(.,"Unexpected")) %>%
        bind_rows(final_imputed_data_ecoli) %>%
        mutate_at("Pool", ~replace_na(.,background_species)) 
    
    #write.table(final_imputed_data, file = "final_imputed_normalized_data_PAL _T_cell_Exp3_( 5 conc 3reps)_NoFAIMS_DDA_with_cont_230206_2023-02-07_0947.txt",sep = "\t",row.names = F)
    
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
    
    
    ymax <- max(-log10(merge_stat_df_final$P.Value)) + 0.5
    y_decrement <- 0.45
    
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
        geom_label(data = point_count_y_axis, aes(x = log2(actual_ratio_val), y = y_pos,fill=new_col_coloring, label = n),size=10, colour="white",show.legend = FALSE) 
    
    merge_stat_df_final_text <- merge_stat_df_final %>% mutate(soft_name=paste0(software_name)) %>% mutate(acq_type=paste0(acquisiton_type))
    
    write.table(merge_stat_df_final_text,file = paste0(file_path,"volcano_plot_",software_name,"_",acquisiton_type,".txt"),sep = 
                  "\t",col.names = T,row.names = F)
    df_roc <- merge_stat_df_final %>%
      select(pep_with_pos, Pool,P.Value)
    #filter(!grepl("Unexpected",Pool))
    
    ### ROC analysis custom func
    df_roc <- df_roc[order(df_roc$P.Value),]
    
    df_roc_func <- compute_roc_curve(df=df_roc, flag = "Others",expected = (4*141))
    
    p14 <- ggplot(df_roc_func, aes(y=as.numeric(tpr), x = as.numeric(fdp))) +
      geom_path(size=1.5) +  #scale_x_reverse() + 
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
           subtitle = paste(subtitle,"including unexpected"), color="Pool Type") + xlim(c(0,100)) + ylim(c(0,100))
    
    write.table(df_roc_func, file = paste0(file_path,"outputs_with_new_script/custom_Roc_analysis_",exp_id,"_",software_name,"_",".txt"),sep = "\t",row.names = F)
    
    #### ROC Analysis using pROC 
    
    df_roc$variant <- ifelse(df_roc$Pool == "Others", TRUE, FALSE)
    #df_roc$non_var <- ifelse(df_roc$Pool == "ISO-REF", TRUE, FALSE)
    
    library(pROC)
    # Calculate ROC curve for raw p-values
    roc_raw_variant <- roc(df_roc$variant, df_roc$P.Value)
    #tpr_and_fpr_variant  <- cbind(roc_raw_variant$sensitivities,
     #                             roc_raw_variant$specificities,
    #                              "Variant Pool")
    
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
      geom_path(size=1.5)  + theme_bw() + #+  scale_x_reverse()
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
    
    sapply(1:14,function(x) ggsave(filename = paste0("p",x,".tiff"),
                                   width = 50, height = 45, 
                                   path = paste0(file_path,"/outputs_with_new_script/"),
                                   units = "cm",
                                   get(paste0("p",x)),
                                   device = "tiff", #".svg"
    ))
    
    
    
}


  
  
  
  
  
  