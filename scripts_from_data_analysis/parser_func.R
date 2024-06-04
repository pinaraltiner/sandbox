#
#library(PhosR)
library(stringr)
library(dplyr)
library(data.table)
library(openxlsx)
library(ggplot2)
library(tidyr)
library(readr)
library(purrr)
library(tibble)
# pep_quant_analysis_bio <- function(file_path,
#                                    file_name,
#                                    sheet_name,
#                                    selected_species,
#                                    acquisiton_type,
#                                    exp_id,
#                                    background_species,
#                                    numerator,
#                                    software_name,
#                                    test_type,
#                                    exp_design,
#                                    num_reps){
  
  pep_quant_parser <- function(file_path,
                               file_name,
                               sheet_name,
                               selected_species,
                               background_species,
                               exp_design,
                               num_reps,
                               numerator,
                               create_impute_vals,
                               software_name,
                               output_dir_name
  ){
    source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/ggplot/ggplot_functions.R")
    
    #quant_peptides <- read.xlsx("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_analysis/experiment_3/Proline_data_analysis/quant_pep_06022023/PAL _T_cell_Exp3_( 5 conc 3reps)_NoFAIMS_DDA_with_cont_230206_2023-02-07_0947.xlsx",sheet = "Best PSM from protein sets")
    sample_size <- length(exp_design) / num_reps
    sample_names <- paste0("A",1:sample_size)
    comparisons <- NULL
    
    for (i in 1:sample_size){
      tmp <- paste0(sample_names[numerator], "/",sample_names[i])
      comparisons[i] <- tmp
      rm(tmp)
    }
    #comparisons <- comparisons[-1]
    comparisons <- comparisons[-numerator]
    
    intended_dir <-paste0(file_path,output_dir_name)
    
    if(dir.exists(intended_dir)){
      new_path <- intended_dir
      
    }else{
      dir.create(intended_dir)
      new_path <- list.dirs(intended_dir)
      
    }
    
    if (software_name == "Proline"){
      
      quant_peptides <- read.xlsx(paste0(file_path,file_name), sheet = sheet_name,sep.names = " ")
      
      abundances_for_impute_names <- quant_peptides %>% 
        select(starts_with("abundance"))  ## spectrum_title remove it because it was not make it as rownames (has duplicates)
      
      ordered_abun_cols <- names(abundances_for_impute_names)[order(names(abundances_for_impute_names), decreasing = FALSE)]
      
      abundances_for_impute <- quant_peptides %>%
        select(ordered_abun_cols) %>%
        rename_with(~exp_design,matches("^abundance"))
      
      quant_peptides_cor_abun <- quant_peptides %>% 
        relocate(ordered_abun_cols) %>%
        rename_with(~exp_design,matches("^abundance"))
      
      if(create_impute_vals==TRUE){
        # Calculate 5 percent quantile of each sample
        impute_values <- apply(abundances_for_impute, 2 , quantile , probs = 0.05 , na.rm = TRUE )
        impute_values <- as.numeric(impute_vals$x)
        write.table(impute_values,file=paste0(new_path,"/impute_values",software_name,".txt"),sep = "\t",row.names = F)
        
      }
      
      quant_peptides_refined <- quant_peptides_cor_abun %>%
        select(sequence,
               modifications,
               #ptm_protein_positions,
               ptm_score,
               #charge,
               #spectrum_title,
               accession,
               contains(exp_design),
        ) %>%
        rename(Sequence=sequence) %>%
        rename(Modifications=modifications) %>%
        separate(accession, into = c("Protein","species"),remove = T,sep="_") 
        
      quant_peptides_phospho <- quant_peptides_refined %>%
        filter(grepl(selected_species,species) & !grepl("CON__",species) & grepl("Phospho",Modifications)) %>%
        #mutate(species=selected_species) %>%
        mutate(phospho_pos = proline_phospho_pos_extraction(Modifications)) %>%
        mutate(pep_with_pos=paste0(Sequence,"_",phospho_pos)) %>%
        select(!c(phospho_pos,Sequence)) 
      
      quant_peptides_phospho <- quant_peptides_phospho %>% mutate(id=1:nrow(quant_peptides_phospho)) 
      
      long_quant_peptides_phospho <- quant_peptides_phospho %>% 
        select(pep_with_pos, contains(exp_design),id) %>%
        pivot_longer(cols = starts_with("E3"), 
                     values_to = "Intensity",
                     names_to = "Experiment",
                     values_drop_na = F) %>%
        group_by(pep_with_pos,Experiment) %>% 
        slice(which.max(Intensity)) %>%
        ungroup()
      
      quant_peptides_phospho_sel_cols <- quant_peptides_phospho %>% 
        select(Modifications,ptm_score,Protein,id,pep_with_pos,species)
      
      filtered_quant_phospho <- long_quant_peptides_phospho %>%
        pivot_wider(names_from = Experiment,values_from = Intensity) %>%
        relocate(exp_design) %>%
        left_join(quant_peptides_phospho_sel_cols,by=c("pep_with_pos","id"))
        
      quant_peptides_ecoli <- quant_peptides_refined %>% 
        filter(grepl(background_species,species) & !grepl("CON__",species)) %>%
        #mutate(species=background_species) %>%
        mutate(pep_with_pos=Sequence)
      
      quant_peptides_ecoli <- quant_peptides_ecoli %>%
        mutate(id=1:nrow(quant_peptides_ecoli)) 
      
      quant_peptides_ecoli_sel_cols <- quant_peptides_ecoli %>% 
        select(Modifications,ptm_score,Protein,id,pep_with_pos,species)
      
      long_quant_peptides_ECOLI <- quant_peptides_ecoli %>%
        select(pep_with_pos, contains(exp_design),id) %>%
        pivot_longer(cols = starts_with("E3"), 
                     values_to = "Intensity",
                     names_to = "Experiment",
                     values_drop_na = F) %>%
        separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "_",remove = F) %>%
        mutate(sample_rep_id_seq = paste(pep_with_pos, Sample_id,Rep_id, sep = "_")) %>%
        group_by(sample_rep_id_seq,Experiment) %>% 
        slice(which.max(Intensity)) %>%
        ungroup()
   
      filtered_quant_peptides_ECOLI <- long_quant_peptides_ECOLI %>%
        select(pep_with_pos, Experiment,Intensity, id) %>%
        pivot_wider(names_from = Experiment,values_from = Intensity) %>%
        left_join(quant_peptides_ecoli_sel_cols,by=c("pep_with_pos","id"))
      
      final_df <- filtered_quant_phospho %>% bind_rows(filtered_quant_peptides_ECOLI) %>%  relocate(exp_design,.after = id)
      
    }else if (software_name == "Proteome Discoverer"){
    
      quant_peptides <- read.table(paste0(file_path,file_name), sep = "\t", header = T)
      
      quant_peptides_refined <- quant_peptides %>%
        select(Sequence, Modifications,
               Modification.Pattern,
               Protein.Accessions,
               Master.Protein.Descriptions,
               Number.of.PSMs,
               contains("Abundances.Normalized."),
               ) %>%
        rename_with(~ exp_design, starts_with("Abundances.Normalized"))
  
      if(create_impute_vals==TRUE){
        
        abundances_for_impute <- quant_peptides %>% 
          select(starts_with("Abundances.Normalized")) %>% #### BEFORE RENAME IT BE SURED THAT COLUMNS ARE THE SAME ORDER AS EXP_DESIGN
          rename_with(~ exp_design, starts_with("Abundances.Normalized"))
        
        # Calculate 5 percent quantile of each sample
        impute_values <- apply(abundances_for_impute, 2 , quantile , probs = 0.05 , na.rm = TRUE )
        write.table(impute_values,file=paste0(new_path,"/impute_values",software_name,".txt"),sep = "\t",row.names = F)
      }
      
      quant_peptides_phospho <- quant_peptides_refined %>%
        filter(grepl(selected_species,Master.Protein.Descriptions) & !grepl("CON__",Master.Protein.Descriptions)) %>%
        filter(grepl("Phospho",Modifications)) %>%
        filter(!grepl("positions not distinguishable", Modification.Pattern)) %>%
        mutate(species=selected_species) %>%
        rowwise() %>%
        mutate(results = list(extract_phospho_numbers(Modifications)),
               phospho_pos = results[[1]],
               phospho_score = results[[2]]) %>%
        select(!results) %>%
        mutate(ptm_score=as.numeric(phospho_score)) %>%
        mutate(pep_with_pos=paste0(Sequence,"_",phospho_pos)) %>% 
        select(!phospho_pos) %>%
        filter(!grepl("NA",phospho_score)) %>%
        ungroup()
    
      quant_peptides_phospho <- quant_peptides_phospho %>% bind_cols(id=1:nrow(quant_peptides_phospho))
      
      long_quant_peptides_phospho <- quant_peptides_phospho %>% 
        select(pep_with_pos, contains(exp_design),id) %>%
        pivot_longer(cols = starts_with("E3"), 
                     values_to = "Intensity",
                     names_to = "Experiment",
                     values_drop_na = F) %>%
        group_by(pep_with_pos,Experiment) %>% 
        slice(which.max(Intensity)) %>%
        ungroup()
      
      quant_peptides_phospho_sel_cols <- quant_peptides_phospho %>% 
        select(Modifications,
               ptm_score,
               Master.Protein.Descriptions,
               Protein.Accessions,
               id,pep_with_pos,species)
      
      filtered_quant_phospho <- long_quant_peptides_phospho %>%
        pivot_wider(names_from = Experiment,values_from = Intensity) %>%
        relocate(exp_design) %>%
        left_join(quant_peptides_phospho_sel_cols,by=c("pep_with_pos","id"))
      
      quant_peptides_ecoli <- quant_peptides_refined %>% 
        filter(grepl("Escherichia coli",Master.Protein.Descriptions) & !grepl("CON__",Master.Protein.Descriptions)) %>%
        mutate(species=background_species) %>%
        mutate(pep_with_pos=Sequence) 
      
      quant_peptides_ecoli <- quant_peptides_ecoli %>% bind_cols(id=1:nrow(quant_peptides_ecoli))
      
      quant_peptides_ecoli_sel_cols <- quant_peptides_ecoli %>%
        select(Modifications,Master.Protein.Descriptions,Protein.Accessions,id,pep_with_pos,species)
      
      long_quant_peptides_ECOLI <- quant_peptides_ecoli %>%
        select(pep_with_pos, contains(exp_design),id) %>%
        pivot_longer(cols = starts_with("E3"), 
                     values_to = "Intensity",
                     names_to = "Experiment",
                     values_drop_na = F) %>%
        separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "_",remove = F) %>%
        mutate(sample_rep_id_seq = paste(pep_with_pos, Sample_id,Rep_id, sep = "_")) %>%
        group_by(sample_rep_id_seq,Experiment) %>% 
        slice(which.max(Intensity)) %>%
        ungroup()
      
      filtered_quant_peptides_ECOLI <- long_quant_peptides_ECOLI %>%
        select(pep_with_pos, Experiment,Intensity, id) %>%
        pivot_wider(names_from = Experiment,values_from = Intensity) %>%
        left_join(quant_peptides_ecoli_sel_cols,by=c("pep_with_pos","id"))
      
      final_df <- filtered_quant_phospho %>% bind_rows(filtered_quant_peptides_ECOLI) %>% 
        relocate(exp_design,.after = id) %>%
        rename(Protein=Protein.Accessions)
      
      }else if (software_name == "MaxQuant"){
      
      source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/get_modification_func/getModificationPosition_general_change_condition_current_mod_sequence_MQ_Spectronaut.R")
      
      quant_peptides <- read_tsv(paste0(file_path,file_name))
      
      quant_peptides_refined <- quant_peptides %>%
        select(Sequence,
               `Modified sequence`,
               `Phospho (STY) Probabilities`,
               Modifications,
               Proteins,
               id,
               `Protein group IDs`,
               `Taxonomy IDs`,
               Experiment,
               Intensity)
      
      if(create_impute_vals ==TRUE){
        imputed_values <- quant_peptides_refined %>% 
          group_by(Experiment) %>% 
          summarise(fifth_quantile=quantile(Intensity,probs=0.05,na.rm=TRUE))
        
        imputed_values_vec <- as.numeric(imputed_values$fifth_quantile)
        
        write.table(imputed_values_vec,file=paste0(new_path,"/impute_values",software_name,".txt"),sep = "\t",row.names = F)
      }
      
      quant_phospho <- quant_peptides_refined %>% 
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
        filter(grepl(10090,`Taxonomy IDs`) & !grepl("CON__",Proteins) & grepl("Phospho",`Modified sequence`)) %>%
        drop_na(`Modified sequence`)
      
      ## PHOSPHO-POSITION EXTRACTION ##
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
               "all_mods_with_mod_type" = "mods.x") %>%
        mutate(species=selected_species)
      
      quant_phospho_complete <- mutate(result_with_common_col,quant_phospho) %>%
        group_by(pep_with_pos,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
        slice(which.max(Intensity)) %>%
        ungroup()
        
      df_reversible_loc_val <- quant_phospho_complete %>%
        select(pep_with_pos,Experiment,id, `Protein group IDs`,extracted_values) %>% ## ADD COLUMN YOU WANT) 
        mutate(unique_row_name = paste(id,`Protein group IDs`,sep="@")) %>%
        select(-id,-`Protein group IDs`) %>%
        pivot_wider(names_from = Experiment,values_from = c(unique_row_name,extracted_values)) %>%
        select(!starts_with("unique_row_name")) %>% #:can be used to combine the object that you use
        rowwise() %>%
        mutate(ptm_score = max(c_across(where(is.numeric)), na.rm = TRUE)) %>%
        select(ptm_score)
        
      df_reversible_proteins <- quant_phospho_complete %>%
        select(Experiment,pep_with_pos,Proteins,id, `Protein group IDs`) %>%
        mutate(unique_row_name = paste(id,`Protein group IDs`,sep="@")) %>%
        select(-id,-`Protein group IDs`) %>%
        pivot_wider(names_from = "Experiment",values_from = c(unique_row_name,Proteins)) %>%
        select(!starts_with("unique_row_name")) %>%
        #mutate(comb_prots = apply(across(all_of(paste0("Proteins_",exp_design))), 1, paste, collapse = "@"))
        rowwise() %>%
        mutate(comb_prots = paste(na.omit(c_across(all_of(paste0("Proteins_",exp_design)))), collapse = "@")) %>%
        ungroup() %>%
        select(comb_prots)
      
      df_reversible_mods <- quant_phospho_complete %>%
        select(Experiment,pep_with_pos,Modifications,id,`Protein group IDs`) %>%
        mutate(unique_row_name = paste(id,`Protein group IDs`,sep="@")) %>%
        select(-id,-`Protein group IDs`) %>%
        pivot_wider(names_from = "Experiment",values_from = c(unique_row_name,Modifications)) %>%
        select(!starts_with("unique_row_name")) %>%
        #mutate(comb_prots = apply(across(all_of(paste0("Proteins_",exp_design))), 1, paste, collapse = "@"))
        rowwise() %>%
        mutate(comb_mod = paste(na.omit(c_across(all_of(paste0("Modifications_",exp_design)))), collapse = "@")) %>%
        ungroup() %>%
        select(comb_mod)
      
       ### FINAL DF FROM SELECTED SPECIES 
      pepwithpos_pivot_int<- quant_phospho_complete %>%
        select(Experiment,Intensity,pep_with_pos) %>%
        pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
        bind_cols(df_reversible_loc_val,df_reversible_proteins,df_reversible_mods) %>%
        filter(!is.infinite(ptm_score)) %>%
        mutate(species=selected_species)
      
      quant_peptides_ecoli <- quant_peptides_refined %>%
        filter(grepl(83333,`Taxonomy IDs`) &!grepl("CON__",Proteins)) %>%
        mutate(species = background_species) %>%
        mutate(pep_with_pos=Sequence) %>%
        separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F) %>%
        mutate(sample_rep_id_seq = paste(pep_with_pos, Sample_id,Rep_id, sep = "-")) %>%
        group_by(sample_rep_id_seq,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
        slice(which.max(Intensity)) %>%
        ungroup()
      
      df_reversible_proteinsECOLI <- quant_peptides_ecoli %>%
        select(Experiment,pep_with_pos,Proteins,id, `Protein group IDs`) %>%
        mutate(unique_row_name = paste(id,`Protein group IDs`,sep="@")) %>%
        select(-id,-`Protein group IDs`) %>%
        pivot_wider(names_from = "Experiment",values_from = c(unique_row_name,Proteins)) %>%
        select(!starts_with("unique_row_name")) %>%
        #mutate(comb_prots = apply(across(all_of(paste0("Proteins_",exp_design))), 1, paste, collapse = "@"))
        rowwise() %>%
        mutate(comb_prots = paste(na.omit(c_across(all_of(paste0("Proteins_",exp_design)))), collapse = "@")) %>%
        ungroup() %>%
        select(comb_prots)
      
      df_reversible_modsECOLI <- quant_peptides_ecoli %>%
        select(Experiment,pep_with_pos,Modifications,id,`Protein group IDs`) %>%
        mutate(unique_row_name = paste(id,`Protein group IDs`,sep="@")) %>%
        select(-id,-`Protein group IDs`) %>%
        pivot_wider(names_from = "Experiment",values_from = c(unique_row_name,Modifications)) %>%
        select(!starts_with("unique_row_name")) %>%
        #mutate(comb_prots = apply(across(all_of(paste0("Proteins_",exp_design))), 1, paste, collapse = "@"))
        rowwise() %>%
        mutate(comb_mod = paste(na.omit(c_across(all_of(paste0("Modifications_",exp_design)))), collapse = "@")) %>%
        ungroup() %>%
        select(comb_mod)
      
      ### FINAL DF FROM BACKGROUND SPECIES
      pepwithpos_pivot_intECOLI<- quant_peptides_ecoli %>%
        select(Experiment,Intensity,pep_with_pos) %>%
        pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
        relocate(exp_design) %>%
        bind_cols(df_reversible_proteinsECOLI,df_reversible_modsECOLI) %>%
        mutate(species=background_species)

      
      final_df <- pepwithpos_pivot_int %>% bind_rows(pepwithpos_pivot_intECOLI) %>% 
        relocate(exp_design,.after = ptm_score) %>% 
        rename(Protein=comb_prots) %>%
        rename(Modifications=comb_mod)
        
        
      #### df_reversible is a reversible object to retrieve any features from the very beginning,
      ## It can be done by using this object thanks to unique_row_name column (in the quant_phospho_complete)
      #### ONLY THING THAT SHOULD BE DONE TO USE PIVOT_WIDER() WITH "df_reversible" THEN USE left_join()
      ## unique_row_name is the combination of id and protein group id column from the evidence.txt
      ################################################################################
      
      # df_reversible <- quant_phospho_complete %>%
      #   select(pep_with_pos,Experiment,id, `Protein group IDs`,extracted_values) %>% ## ADD COLUMN YOU WANT) 
      #   mutate(unique_row_name = paste(id,`Protein group IDs`,sep="@")) %>%
      #   select(-id,-`Protein group IDs`) %>%
      #   pivot_wider(names_from = Experiment,values_from = c(unique_row_name,extracted_values)) %>%
      #   select(!starts_with("unique_row_name")) #%>%
      #bind_cols() : can be used to combine the object that you use
  
      # df_reversible_ecoli <- quant_peptides_ecoli %>% 
      #   select(Sequence,Experiment, Intensity,id, `Protein group IDs`) %>%
      #   mutate(unique_row_name = paste(id,`Protein group IDs`,sep="@")) %>%
      #   select(-id,-`Protein group IDs`) %>%
      #   pivot_wider(names_from = Experiment,values_from = c(unique_row_name,Intensity))
      # 
    
    }else if (software_name=="Spectronaut"){
      
      quant_peptides <- read_tsv(paste0(file_path,file_name))
      
      quant_peptides_with_cond <- quant_peptides %>% 
        mutate(Experiment=paste0(R.Condition,"-R",R.Replicate)) %>%
        rename("Intensity"= "EG.TotalQuantity (Settings)") %>%
        mutate(Modifications= `EG.PrecursorId`)
      
      quant_peptides_refined <- quant_peptides_with_cond %>%
        select(PEP.GroupingKey,
               Experiment,
               Intensity,
               `EG.PTMProbabilities [Phospho (STY)]`,
               `EG.PTMPositions [Phospho (STY)]`,
               `EG.PrecursorId`,
               PG.ProteinLabel,
               Modifications,
               Intensity) %>%
        mutate(Modifications=EG.PrecursorId)
      
      if(create_impute_vals ==TRUE){
        
        imputed_values <- quant_peptides_with_cond  %>%
        group_by(Experiment) %>% 
        summarise(fifth_quantile=quantile(Intensity,probs=0.05,na.rm=TRUE))
        
        imputed_values_vec <- as.numeric(imputed_values$fifth_quantile)
        write.table(imputed_values_vec,file=paste0(new_path,"/impute_values",software_name,".txt"),sep = "\t",row.names = F)
      }
      
      quant_phospho <- quant_peptides_refined %>% 
        #select(PEP.GroupingKey, EG.PrecursorId, PG.ProteinLabel,Ex) %>%
        filter(grepl(selected_species,PG.ProteinLabel) & !grepl("CON__",PG.ProteinLabel)) %>%
        filter(grepl("Phospho",EG.PrecursorId)) %>%
        mutate(species=selected_species) 
      
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
    
      df <- t(as.data.frame(lapply(max_prob_and_pos, unlist)))
      # Rename the columns as needed
      colnames(df) <- c("ptm_score","ptm_position","ptm_type")
      rownames(df) <- NULL
      
      quant_phospho_complete <- quant_phospho %>% 
        bind_cols(df) %>% 
        #filter(!grepl("undistinguishable",undistinguishable)) %>%
        #separate(position1, into = c("phospho_score","phospho_pos"),sep = "_") %>%
        mutate(pep_with_pos=paste(PEP.GroupingKey,ptm_position,sep = "_")) %>%
        mutate(ptm_score=as.numeric(ptm_score)) %>%
        group_by(pep_with_pos,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
        slice(which.max(Intensity)) %>%
        ungroup()
    
      df_reversible_loc_val <- quant_phospho_complete %>%
        mutate(id=1:nrow(quant_phospho_complete)) %>%
        select(pep_with_pos,Experiment,id,PG.ProteinLabel,ptm_score) %>% ## ADD COLUMN YOU WANT) 
        mutate(unique_row_name = paste(id,PG.ProteinLabel,sep="@")) %>%
        select(-id,-PG.ProteinLabel) %>%
        pivot_wider(names_from = Experiment,values_from = c(unique_row_name,ptm_score)) %>%
        select(!starts_with("unique_row_name")) %>% #:can be used to combine the object that you use
        rowwise() %>%
        mutate(ptm_score = max(c_across(where(is.numeric)), na.rm = TRUE)) %>%
        select(ptm_score)
      
      
      df_reversible_mods <- quant_phospho_complete %>%
        mutate(id=1:nrow(quant_phospho_complete)) %>%
        select(Experiment,pep_with_pos,PG.ProteinLabel,id,Modifications) %>%
        mutate(unique_row_name = paste(id,PG.ProteinLabel,sep="@")) %>%
        select(-id) %>%
        pivot_wider(names_from = "Experiment",values_from = c(unique_row_name,Modifications)) %>%
        select(!starts_with("unique_row_name")) %>%
        #mutate(comb_prots = apply(across(all_of(paste0("Proteins_",exp_design))), 1, paste, collapse = "@"))
        rowwise() %>%
        mutate(comb_mods = paste(na.omit(c_across(all_of(paste0("Modifications_",exp_design)))), collapse = "@")) %>%
        ungroup() %>%
        select(comb_mods)
      
      df_reversible_proteins <- quant_phospho_complete %>%
        mutate(id=1:nrow(quant_phospho_complete)) %>%
        select(Experiment,pep_with_pos,PG.ProteinLabel,id) %>%
        mutate(unique_row_name = paste(id,PG.ProteinLabel,sep="@")) %>%
        select(-id) %>%
        pivot_wider(names_from = "Experiment",values_from = c(unique_row_name,PG.ProteinLabel)) %>%
        select(!starts_with("unique_row_name")) %>%
        #mutate(comb_prots = apply(across(all_of(paste0("Proteins_",exp_design))), 1, paste, collapse = "@"))
        rowwise() %>%
        mutate(comb_prots = paste(na.omit(c_across(all_of(paste0("PG.ProteinLabel_",exp_design)))), collapse = "@")) %>%
        ungroup() %>%
        select(comb_prots)
      
      ### FINAL DF FROM SELECTED SPECIES 
      pepwithpos_pivot_int<- quant_phospho_complete %>%
        select(Experiment,Intensity,pep_with_pos) %>%
        pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
        relocate(exp_design) %>%
        bind_cols(df_reversible_loc_val,df_reversible_proteins,df_reversible_mods) %>%
        filter(!is.infinite(ptm_score)) %>%
        mutate(species=selected_species)
      
      quant_peptides_ecoli <- quant_peptides_refined %>%
        filter(grepl(background_species,PG.ProteinLabel) &!grepl("CON__",PG.ProteinLabel)) %>%
        mutate(species = background_species) %>%
        mutate(pep_with_pos=PEP.GroupingKey) %>%
        separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F) %>%
        mutate(sample_rep_id_seq = paste(pep_with_pos, Sample_id,Rep_id, sep = "-")) %>%
        group_by(sample_rep_id_seq,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
        slice(which.max(Intensity)) %>%
        ungroup()
      
      df_reversible_proteinsECOLI <- quant_peptides_ecoli %>%
        mutate(id=1:nrow(quant_peptides_ecoli)) %>%
        select(Experiment,pep_with_pos,PG.ProteinLabel,id) %>%
        mutate(unique_row_name = paste(id,PG.ProteinLabel,sep="@")) %>%
        select(-id) %>%
        pivot_wider(names_from = "Experiment",values_from = c(unique_row_name,PG.ProteinLabel)) %>%
        select(!starts_with("unique_row_name")) %>%
        #mutate(comb_prots = apply(across(all_of(paste0("Proteins_",exp_design))), 1, paste, collapse = "@"))
        rowwise() %>%
        mutate(comb_prots = paste(na.omit(c_across(all_of(paste0("PG.ProteinLabel_",exp_design)))), collapse = "@")) %>%
        ungroup() %>%
        select(comb_prots)
      
      df_reversible_proteins_modsECOLI <-quant_peptides_ecoli %>% mutate(id=1:nrow(quant_peptides_ecoli)) %>%
        select(Experiment,pep_with_pos,Modifications,id,PG.ProteinLabel) %>%
        mutate(unique_row_name = paste(id,Modifications,sep="@")) %>%
        select(-id) %>%
        pivot_wider(names_from = "Experiment",values_from = c(unique_row_name,Modifications)) %>%
        select(!starts_with("unique_row_name")) %>%
        #mutate(comb_prots = apply(across(all_of(paste0("Proteins_",exp_design))), 1, paste, collapse = "@"))
        rowwise() %>%
        mutate(comb_mods = paste(na.omit(c_across(all_of(paste0("Modifications_",exp_design)))), collapse = "@")) %>%
        ungroup() %>%
        select(comb_mods)
      
      ### FINAL DF FROM BACKGROUND SPECIES
      pepwithpos_pivot_intECOLI<- quant_peptides_ecoli %>%
        select(Experiment,Intensity,pep_with_pos) %>%
        pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
        bind_cols(df_reversible_proteinsECOLI,df_reversible_proteins_modsECOLI) %>%
        mutate(species=background_species)
      
      final_df <- pepwithpos_pivot_int %>% bind_rows(pepwithpos_pivot_intECOLI) %>% 
        relocate(exp_design,.after = ptm_score) %>% 
        rename(Modifications=comb_mods) %>%
        rename(Protein=comb_prots)
      
    }else if(software_name=="DIANN"){ ### IF THIS IS DIANN, mapping file has to be given as parameter in this parser!!
      library(purrr)
      source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/get_modification_func/getModificationPosition_func_for_all_mods.R")
      
      quant_peptides <- read_tsv(paste0(file_path,file_name),col_names = T)
      mapping_file <- read.table(mapping,sep = "\t",header = T)
      exp_design <- mapping_file$Experiment
      
      quant_peptides_with_cond <- quant_peptides %>% 
        left_join(mapping_file,by="Run") %>%
        rename("Sequence"="Stripped.Sequence") %>%
        rename("Intensity"= "PG.Normalised") 
      
      if(create_impute_vals ==TRUE){
        imputed_values <- quant_peptides_with_cond  %>%
          group_by(Experiment) %>% 
          summarise(fifth_quantile=quantile(Intensity,probs=0.05,na.rm=TRUE))
        
        #imputed_values_vec <- as.vector(imputed_values$first_quantile)
        write.table(impute_values,file=paste0(new_path,"/impute_values",software_name,".txt"),sep = "\t",row.names = F)
      }
      
      quant_peptides_refined <- quant_peptides_with_cond %>%
        select(Sequence,
               Experiment,
               Intensity,
               Modified.Sequence,
               PTM.Site.Confidence,
               Protein.Names,
               Experiment,
               Intensity)
      
      
      quant_phospho <- quant_peptides_refined %>% 
        filter(grepl(selected_species, Protein.Names) & !grepl("CON__",Protein.Names)) %>%
        filter(grepl("UniMod:21",Modified.Sequence))
      
      df2 <- apply(quant_phospho[,"Modified.Sequence"],1,getModificationPosition_)
      
      #apply(X = as.data.frame(quant_phospho[,"Modified.Sequence"]),1,function(x){getModificationPosition_general(mod_seq = x,software_name = )})
      
      test <- map_dfr(df2, enframe)
      
      ## OUTPUT FORMAT DOES NOT SUITABLE FOR DISTINGUSING BTW MODS 
      results1 <- map_dfr(df2, ~ enframe(.x)) %>%
        filter(grepl("modification",name) | grepl("pep_seq", name)) %>%
        mutate(value = map_chr(value, str_c, collapse="&")) %>%
        mutate(mods=case_when(grepl("(UniMod:21",fixed = T,name) ~ "phospho",
                              grepl("(UniMod:35",fixed = T,name) ~ "Oxidation",
                              grepl("(UniMod:1",fixed = T,name) ~ "N-term_Acetyl",
                              TRUE ~ ""))
      
      ## RE-SHAPING DATA TO HAVE ONE SEQ with ALL MODS in ONE ROW
      
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
        full_join(filter(results_with_index, name == "pep_seq"), by = "id") %>% 
        ungroup() %>%
        ## SELECTING USEFUL COLUMNS
        select(c(name.x,value.x,mods.x,value.y))
      
      result_with_common_col <- reshaped_results %>%
        filter(grepl("modification_(UniMod:21)",fixed = T,name.x)) %>%
        mutate(pep_with_pos = paste(value.y,value.x, sep = "_")) %>%
        select(!value.y) %>%
        rename("mods_unimod_id" = "name.x",
               "phospho_positions" ="value.x",
               "all_mods_with_mod_type" = "mods.x")
      
      
      quant_phospho_complete <- mutate(result_with_common_col,quant_phospho) %>%
        group_by(pep_with_pos,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
        slice(which.max(Intensity)) %>%
        ungroup()
      
      df_reversible_loc_val <- quant_phospho_complete %>%
        mutate(id=1:nrow(quant_phospho_complete)) %>%
        select(pep_with_pos,Experiment,id, Protein.Names,PTM.Site.Confidence) %>% ## ADD COLUMN YOU WANT) 
        mutate(unique_row_name = paste(id,Protein.Names,sep="@")) %>%
        select(-id,-Protein.Names) %>%
        pivot_wider(names_from = Experiment,values_from = c(unique_row_name,PTM.Site.Confidence)) %>%
        select(!starts_with("unique_row_name")) %>% #:can be used to combine the object that you use
        rowwise() %>%
        mutate(ptm_score = max(c_across(where(is.numeric)), na.rm = TRUE)) %>%
        select(ptm_score)
      
      df_reversible_proteins <- quant_phospho_complete %>%
        mutate(id=1:nrow(quant_phospho_complete)) %>%
        select(Experiment,pep_with_pos,Protein.Names,id) %>%
        mutate(unique_row_name = paste(id,Protein.Names,sep="@")) %>%
        select(-id) %>%
        pivot_wider(names_from = "Experiment",values_from = c(unique_row_name,Protein.Names)) %>%
        select(!starts_with("unique_row_name")) %>%
        #mutate(comb_prots = apply(across(all_of(paste0("Proteins_",exp_design))), 1, paste, collapse = "@"))
        rowwise() %>%
        mutate(comb_prots = paste(na.omit(c_across(all_of(paste0("Protein.Names_",exp_design)))), collapse = "@")) %>%
        ungroup() %>%
        select(comb_prots)
      
      ### FINAL DF FROM SELECTED SPECIES 
      pepwithpos_pivot_int<- quant_phospho_complete %>%
        select(Experiment,Intensity,pep_with_pos) %>%
        pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
        bind_cols(df_reversible_loc_val,df_reversible_proteins) %>%
        filter(!is.infinite(ptm_score)) %>%
        mutate(species=selected_species)
      #######TODO: FIX THE PROBLEM HERE: FILTERING OF BACKGROUND SPECIES (every Prot. Name is NA except SELECTED sPECIES)
      ####### RIGHT NOW, THERE IS NO WAY TO DISTINGUISH ECOLI & CONT. and RETRIEVE THEIR PROT_NAMEs
      quant_peptides_ecoli <- quant_peptides_refined %>%
        filter(!grepl(selected_species,Protein.Names)) %>% # & !grepl("CON__",Protein.Names)
        mutate(species = background_species) %>%
        mutate(pep_with_pos=Sequence) %>%
        separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F) %>%
        mutate(sample_rep_id_seq = paste(pep_with_pos, Sample_id,Rep_id, sep = "-")) %>%
        group_by(sample_rep_id_seq,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
        slice(which.max(Intensity)) %>%
        ungroup()
      
      df_reversible_proteinsECOLI <- quant_peptides_ecoli %>%
        select(Experiment,pep_with_pos,Protein.Names,id) %>%
        mutate(unique_row_name = paste(id,`Protein group IDs`,sep="@")) %>%
        select(-id) %>%
        pivot_wider(names_from = "Experiment",values_from = c(unique_row_name,Protein.Names)) %>%
        select(!starts_with("unique_row_name")) %>%
        #mutate(comb_prots = apply(across(all_of(paste0("Proteins_",exp_design))), 1, paste, collapse = "@"))
        rowwise() %>%
        mutate(comb_prots = paste(na.omit(c_across(all_of(paste0("Protein.Names_",exp_design)))), collapse = "@")) %>%
        ungroup() %>%
        select(comb_prots)
      ### FINAL DF FROM BACKGROUND SPECIES
      pepwithpos_pivot_intECOLI<- quant_peptides_ecoli %>%
        select(Experiment,Intensity,pep_with_pos) %>%
        relocate(exp_design) %>%
        pivot_wider(names_from = "Experiment",values_from = "Intensity") %>%
        bind_cols(df_reversible_proteinsECOLI) %>%
        mutate(species=background_species)
      
      final_df <- pepwithpos_pivot_int %>% bind_rows(pepwithpos_pivot_intECOLI)
      
    }
    write.table(final_df,file=paste0(new_path,"/refined_input",software_name,".txt"),sep = "\t",row.names = F)
  }
  