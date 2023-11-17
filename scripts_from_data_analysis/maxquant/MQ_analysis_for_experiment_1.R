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
library(gtools)
###############################################
## MAXQUANT DATA ANALYSIS

### Source code was taken from here: https://rdrr.io/github/singjc/mstools/src/R/getModificationPosition.R
## The code was modified based on what I want and based on software input tyoe (DIANN-and MaxQuant)
#source("D:/dev/Desktop_copy/PHD/data_analysis/scripts/getModificationPositionMQ_func.R")
#source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/get_modification_func/getModificationPositionMQ_func_edit.R")
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/get_modification_func/getModificationPosition_general_change_condition_current_mod_sequence_MQ_Spectronaut.R")
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/ggplot/ggplot_functions.R")
#all_files <- list.files(paste0(file_path,all_dirs[i],"/"),pattern = ".xlsx")

#comb_result <- NULL
#comb_ecoli <- NULL
final_MQ_pep_quant_analysis_exp1 <- function(file_path,
                                             curr_dir,
                                             #file_name,
                                             theo_file_path,
                                             theo_file_name,
                                             sheet_theo_name,
                                             selected_species,
                                             exp_id,
                                             #exp_design,
                                             background_species,
                                             acquisiton_type,
                                             software_name){
  
 
    assign(paste0("tmp"), read.table(file=paste0(file_path,curr_dir,"/evidence.txt"),sep = "\t",header = T))
    
    phospho_tmp <- tmp %>% 
      filter(grepl("Phospho",Modified.sequence) & grepl(selected_species,Proteins)) %>%
      drop_na(Modified.sequence) %>%
      mutate(curr_dir) %>%
      mutate(acq_type=acquisiton_type) %>%
      mutate(soft_name=software_name) %>%
      separate(Experiment,into =c("Exp_id","Sample_id","tmp","inj_id"),sep = "-",remove = F) %>%
      mutate(new_col=paste(Exp_id,Sample_id,sep = "_"))
    #mutate(phospho_pos = proline_phospho_pos_extraction(Modifications)) %>%
    #mutate(pep_with_pos = paste(Sequence,phospho_pos,sep = "_")) %>%
    #distinct(pep_with_pos,.keep_all = TRUE)
    
    #assign(paste0(all_dirs[i]), phospho_tmp)
    #assign(paste0("comb_result"),bind_rows(comb_result, phospho_tmp))
    
    if (background_species %in% c("ECOLI","Escherichia coli")){
      
      if (any(curr_dir %in% c("exp1_with_FAIMS_with_Ecoli_OXPAL230421"))) {

        ecoli_tmp <- tmp %>%
          filter(!grepl(selected_species,Proteins) & !grepl("CON_", Proteins)) %>%
          mutate(Experiment = ifelse(Raw.file == "OXPAL230421_08_-45","E1-M2-coli-inj1",Experiment)
          ) %>%
          mutate(curr_dir) %>%
          mutate(acq_type=acquisiton_type) %>%
          mutate(soft_name=software_name) %>%
          group_by(Sequence,Experiment) %>%
          slice(which.max(Intensity)) %>%
          ungroup() %>%
          #distinct(Sequence,.keep_all = TRUE) %>%
          separate(Experiment,into =c("Exp_id","Sample_id","tmp","inj_id"),sep = "-",remove = F) %>%
          mutate(new_col=paste(Exp_id,Sample_id,sep = "_"))

      }else{
        
        ecoli_tmp <- tmp %>%
          filter(!grepl(selected_species,Proteins) & !grepl("CON_", Proteins)) %>%
          #mutate(Experiment = ifelse(Raw.file == "OXPAL230421_08_-45","E1-M2-coli-inj1",Experiment)
          #) %>%
          mutate(curr_dir) %>%
          mutate(acq_type=acquisiton_type) %>%
          mutate(soft_name=software_name) %>%
          group_by(Sequence,Experiment) %>%
          slice(which.max(Intensity)) %>%
          ungroup() %>%
          #distinct(Sequence,.keep_all = TRUE) %>%
          separate(Experiment,into =c("Exp_id","Sample_id","tmp","inj_id"),sep = "-",remove = F) %>%
          mutate(new_col=paste(Exp_id,Sample_id,sep = "_"))

      }
      
      plot6 <- gg_barplt_id_pep_count(data_set = ecoli_tmp,
                                      x_df = ecoli_tmp$Raw.file,
                                      fill_df = ecoli_tmp$new_col,
                                      ymax = 20000,
                                      header = paste("Total number of identified synthetic phospho-sequences across each pool",sep=" "),
                                      caption_lab = "Wrong Sequences were removed for the futher analysis.",
                                      x_lab = "Sample id",
                                      fill_lab =  "Sample id",
                                      y_lab = "Number of identified Sequence",
                                      subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type,curr_dir," data processed by ", software_name)) +
        scale_fill_brewer(palette = "Dark2") + theme(axis.text.x = element_text(angle = 90))
      
      write.table(ecoli_tmp,file=paste0(file_path,curr_dir,"/Exp",exp_id,
                                        "id_ecoli_seq.tsv"),
                  sep = "\t",col.names = T,row.names = F)
      
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
      mutate(pep_with_pos = paste(value.y,value.x, sep = "_")) %>%
      select(!value.y) %>%
      rename("mods_id" = "name.x",
             "phospho_positions" ="value.x",
             "all_mods_with_mod_type" = "mods.x")
    
    
    #getModificationPosition_MQ(mod_seq = id_syn_pep$`Modified sequence`,F)
    
    final_results_with_common_col <- phospho_tmp %>%
      mutate(result_with_common_col) %>%
      mutate(extracted_values = sapply(str_extract_all(Phospho..STY..Probabilities,
                                                       "\\(\\d+(\\.\\d+)?\\)"), function(x) {
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
      )%>% relocate(extracted_values,.after = Sequence)
        #drop_na(extracted_values) %>%

    
    
    ## READ THEO LIST
    pep_list_w_theo <- read.xlsx(paste0(theo_file_path, theo_file_name), sheet = sheet_theo_name)
    pep_list_w_theo_quant <- pep_list_w_theo[,-1]
    common_col_theo_quant <- as.data.frame(paste(pep_list_w_theo_quant$Phosphopeptide.sequence,
                                                 pep_list_w_theo_quant$modified.position.in.peptide, sep = "_"))
    
    ## ADD COMMON COLUMN TO MERGE WITH EXP. DATA
    colnames(common_col_theo_quant) <- "pep_with_pos"
    pep_list_w_theo_quant_new <- cbind(common_col_theo_quant,pep_list_w_theo_quant)
    
    
    ## EXTRACTION OF UNIQUE SEQUENCES
    pep_list_w_theo_unique <- pep_list_w_theo %>% 
      distinct(Phospopeptide.sequence,.keep_all = TRUE) %>%
      rename(Sequence = Phospopeptide.sequence) %>% 
      select(Sequence,pool_id) %>%
      mutate(Pool_for_seq_merge="Correct")#rename(Pool_for_seq_merge=pool_id)
    
    ## MERGING WITH SEQUENCE TO ELIMINATE WRONG SEQUENCES
    if(background_species %in% c("ECOLI","Escherichia coli")){
      
      comb_result_seq <- final_results_with_common_col %>% 
        select(Sequence, pep_with_pos,extracted_values,
               Phospho..STY..Probabilities,
               new_col,inj_id,
               Experiment,Raw.file,
               `curr_dir`, MS.MS.scan.numbers,
                Retention.time,Intensity, #Marked.as
               Proteins) %>%
        full_join(pep_list_w_theo_unique,by="Sequence") %>%
        mutate_at("Pool_for_seq_merge", ~replace_na(.,"Unexpected")) %>%
        mutate(Pool_for_seq_merge= ifelse(is.na(Raw.file),"missing",Pool_for_seq_merge)) %>%
        #filter(!grepl("Unexpected",Pool_for_seq_merge))
        mutate(acq_type=acquisiton_type) %>%
        mutate(soft_name=software_name) %>%
        #separate(`all_dirs[i]`,into = c("tmp","injs","tmp2"),sep = "_",remove = F) %>%
        #select(!c(tmp,tmp2)) %>%
        #separate(injs,into = c("tmp","expid","sample_id","Ecoli","inj"),sep = "-") %>%
        #mutate(new_col=paste(expid,sample_id,Ecoli,sep = "_")) %>%
        mutate(new_col_inj = paste(new_col,inj_id,sep = "_")) %>%
        mutate(new_col_inj= ifelse(is.na(new_col_inj),"missing",new_col_inj))
      
      
      comb_result_seq_dist <- comb_result_seq %>% 
        group_by(Sequence,new_col) %>% 
        distinct(Sequence,.keep_all = TRUE)%>% 
        ungroup() 
      
      write.table(comb_result_seq,file=paste0(file_path,curr_dir,"/Exp",exp_id,
                                              "merge_theo_list_id_phospho_sequences.tsv",sep = "_"),
                  sep = "\t",col.names = T,row.names = F)
      
      plot1 <- gg_barplt_id_pep_count(data_set = comb_result_seq_dist,
                                      x_df = comb_result_seq_dist$Experiment, #new_col
                                      fill_df = comb_result_seq_dist$Pool_for_seq_merge,
                                      ymax = 250,
                                      header = paste("Total number of identified synthetic phospho-sequences across each pool",sep=" "),
                                      caption_lab = "Wrong Sequences were removed for the futher analysis.",
                                      x_lab = "Sample id",
                                      fill_lab =  "Sample id",
                                      y_lab = "Number of identified Sequence",
                                      subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) +
        theme(axis.text.x = element_text(angle = 90))
      
    }else{
      
      comb_result_seq <- final_results_with_common_col %>% 
        select(Sequence, pep_with_pos,
               Phospho..STY..Probabilities,extracted_values,
               new_col,inj_id,
               Experiment,Raw.file,
               curr_dir, MS.MS.scan.numbers,
               Retention.time,Intensity, #Marked.as
               Proteins) %>%
        full_join(pep_list_w_theo_unique,by="Sequence") %>%
        mutate_at("Pool_for_seq_merge", ~replace_na(.,"Unexpected")) %>%
        mutate(Pool_for_seq_merge= ifelse(is.na(Raw.file),"missing",Pool_for_seq_merge)) %>%
        #filter(!grepl("Unexpected",Pool_for_seq_merge))
        mutate(acq_type=acquisiton_type) %>%
        mutate(soft_name=software_name)
        #separate(`all_dirs[i]`,into = c("tmp","injs","tmp2"),sep = "_",remove = F) %>%
        #select(!c(tmp,tmp2)) %>%
        #separate(injs,into = c("tmp","expid","sample_id","Ecoli","inj"),sep = "-") %>%
        #mutate(new_col=paste(expid,sample_id,Ecoli,sep = "_")) %>%
       
      
      comb_result_seq_dist <- comb_result_seq %>% 
        group_by(Sequence,new_col) %>% 
        distinct(Sequence,.keep_all = TRUE)%>% 
        ungroup() 
      
      write.table(comb_result_seq,file=paste0(file_path,curr_dir,"/Exp",exp_id,
                                              "merge_theo_list_id_phospho_sequences.tsv",sep = "_"),
                  sep = "\t",col.names = T,row.names = F)
      
      plot1 <- gg_barplt_id_pep_count(data_set = comb_result_seq_dist,
                                      x_df = comb_result_seq_dist$Raw.file ,
                                      fill_df = comb_result_seq_dist$Pool_for_seq_merge ,
                                      ymax = 250,
                                      header = paste("Total number of identified synthetic phospho-sequences across each pool",sep=" "),
                                      caption_lab = "Wrong Sequences were removed for the futher analysis.",
                                      x_lab = "Sample id",
                                      fill_lab =  "Sample id",
                                      y_lab = "Number of identified Sequence",
                                      subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) + 
        scale_fill_brewer(palette = "Paired") + theme(axis.text.x = element_text(angle = 90))
      
    }
      comb_result_pep <- comb_result_seq %>% 
        filter(!grepl("Unexpected",Pool_for_seq_merge) & !grepl("missing",Pool_for_seq_merge)) %>%
        select(!Pool_for_seq_merge) %>%
        full_join(pep_list_w_theo_quant_new,by="pep_with_pos") %>% ## IF FULL_JOIN IS USED, 
        mutate(Pool_for_pep_merge=NA) %>%
        mutate(Pool_for_pep_merge=ifelse(is.na(Pool),"Wrong Localization","Correct")) %>%  ## MISSING PEPTIDE INFO COULD BE OBTAINED
        mutate(Pool_for_pep_merge= ifelse(is.na(Raw.file),"missing",Pool_for_pep_merge)) %>%
        mutate(Pool=ifelse(is.na(Pool),"Unexpected",Pool)) %>%
        mutate(Pool= ifelse(is.na(Raw.file),"missing",Pool)) %>%
        mutate(isomericity=ifelse(is.na(isomericity),"Unexpected",isomericity)) %>%
        mutate(isomericity=ifelse(is.na(Raw.file),"missing",isomericity))# %>%
        #drop_na(pool_id) ## ELIMINATION OF UNEXPECTED PEPTIDES
        #mutate(acq_type=acquisiton_type) %>%
        #mutate(soft_name=software_name)
      
      write.table(comb_result_pep ,file=paste0(file_path,curr_dir,"/Exp",exp_id,
                                               "merge_theo_list_id_phospho_sites.tsv",sep = "_"),
                                              
                  sep = "\t",col.names = T,row.names = F)
      
      comb_result_pep_max_int  <- comb_result_pep %>% 
        group_by(pep_with_pos,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
        slice(which.max(Intensity)) %>% ## ELIMINATE MULTIPLE CHARGES
        ungroup() 
      
      plot5 <- comb_result_pep %>% distinct(pep_with_pos,.keep_all = T) %>%
        ggplot( aes(x=Pool,fill=isomericity)) +
        geom_bar() + theme_bw() +
        geom_text(color="black",aes(label=after_stat(count)),
                  show.legend = F,stat = "count",size=8,
                  position = position_stack(vjust = 0.5)) + 
        scale_fill_brewer(palette = "Dark2") +
        labs(title=paste("Number of Identified phospho-sites based on pools of Experiment 2",
                         acquisiton_type),x="Sample id",y="Number of identified Sequence"
             ) +
        theme_minimal() +
        theme(legend.text = element_text(size=30), 
              axis.title.x = element_text(size=30),
              axis.title.y = element_text(size=30),
              plot.title = element_text(size=35),
              plot.subtitle = element_text(size = 25),
              plot.caption = element_text(size = 25),
              legend.title=element_text(size=30),
              axis.text.x = element_text(size=20,angle = 90),
              axis.text.y = element_text(size = 30),
              axis.title=element_text(size=30))
      
      
      
      plot2 <- ggplot(comb_result_pep_max_int, aes(x=Raw.file,fill=Pool_for_pep_merge)) +
        geom_bar() + theme_bw() +
        geom_text(color="black",aes(label=after_stat(count)),
                  show.legend = F,stat = "count",size=8,
                  position = position_stack(vjust = 0.5)) + 
        scale_fill_brewer(palette = "Dark2") +
        labs(title=paste("Number of Identified phospho-sites for each pool \n","Experiment",
                         exp_id, acquisiton_type),x="Sample id",caption = "Pools were missing",y="Number of identified Sequence") + 
        theme_minimal() +
        theme(legend.text = element_text(size=30), 
              axis.title.x = element_text(size=30),
              axis.title.y = element_text(size=30),
              plot.title = element_text(size=35),
              plot.subtitle = element_text(size = 25),
              plot.caption = element_text(size = 25),
              legend.title=element_text(size=30),
              axis.text.x = element_text(size=20,angle = 90),
              axis.text.y = element_text(size = 30),
              axis.title=element_text(size=30))
     
      ##################### SITE CONFIDENCE #####################
      threshold <- seq(0,1,length=100)
      final_df <- NULL
      
      for (i in threshold){
        # assign(paste0("tmp",i), comb_result_pep %>% 
        #   filter(max_value_column > threshold[i]) %>%
        #   count(Pool_for_pep_merge,pool_id.x))
        df <- comb_result_pep_max_int %>% 
          
          filter(extracted_values > i) %>%
          count(Pool_for_pep_merge) %>% mutate(threshold_val = i)
        
        final_df <-bind_rows(final_df,df)
        rm(df)
      }
      
      write.table(final_df, file=paste0(file_path,curr_dir,"/Experiment",exp_id,software_name,"_number_of_sites_with_scores.tsv"),sep = "\t",col.names = T,row.names = F)
    
      wrong_df <- final_df %>% filter(grepl("Wrong",Pool_for_pep_merge))
      correct_df <- final_df %>% filter(!grepl("Wrong",Pool_for_pep_merge))
      
      plot4 <- ggplot(correct_df,aes(y=n, x=threshold_val)) + geom_line(size=2) +
        scale_fill_brewer(palette = "Dark2") +
        labs(title=paste("Number of Correctly Identified phospho-sites at various thresholds \n","Experiment",
                         exp_id, acquisiton_type),x="Sample id",y="Number of identified Sequence",
             subtitle = acquisiton_type) + 
        theme_minimal() +
        theme(legend.text = element_text(size=30), 
              axis.title.x = element_text(size=30),
              axis.title.y = element_text(size=30),
              plot.title = element_text(size=35),
              plot.subtitle = element_text(size = 25),
              plot.caption = element_text(size = 25),
              legend.title=element_text(size=30),
              axis.text.x = element_text(size=20,angle = 90),
              axis.text.y = element_text(size = 30),
              axis.title=element_text(size=30)) #+ 
      #scale_x_continuous(limits = c(0, 100),breaks = seq(from = 0, to = 100, by = 10)) +  # Set the ticks for the x-axis
      # scale_y_continuous(limits = c(0, 180),breaks = seq(from = 0, to = 180, by = 10))
      
      #################### ADJACENT ASSESSMENT ##################
      elements <- c("S", "T", "Y")
      
      test <- permutations(n = length(elements), r = 1, v = elements, repeats.allowed = TRUE)
      test <- as.data.frame(test)
      test[,"comb"] <- as.data.frame(paste0(test$V1))
      
      test2 <- permutations(n = length(elements), r = 2, v = elements, repeats.allowed = TRUE)
      test2 <- as.data.frame(test2)
      test2[,"comb"] <- as.data.frame(paste0(test2$V1,test2$V2))
      #colnames(test2)[3] <- "comb"
      
      test3 <- permutations(n = length(elements), r = 3, v = elements, repeats.allowed = TRUE)
      test3 <- as.data.frame(test3)
      test3[,"comb"] <- as.data.frame(paste0(test3$V1,test3$V2,test3$V3))
      #colnames(test3)[4] <- "comb"
      
      possibilites <-test2 %>% bind_rows(test3) %>% select(comb)
      possibilites <- possibilites$comb
      
      sing_STY <-test$comb
      
      comb_pep_adj_dist_df_neg <- comb_result_pep_max_int %>% mutate(
        Contains_STY = rowSums(sapply(possibilites, function(query) grepl(query, Sequence))) > 0,
        STY_count = rowSums(sapply(possibilites, function(query) str_count(Sequence, query)))) %>% 
        mutate(STY_count = ifelse(Contains_STY == FALSE,(rowSums(sapply(sing_STY, function(query) str_count(Sequence, query)))*-1),STY_count))
      #mutate(stat=ifelse(STY_count < 0,"negative","positive"))
      
      comb_pep_adj_dist_df_pos <- comb_result_pep_max_int %>% mutate(
        Contains_STY = rowSums(sapply(possibilites, function(query) grepl(query, Sequence))) > 0,
        STY_count = rowSums(sapply(possibilites, function(query) str_count(Sequence, query)))) %>% 
        mutate(STY_count = ifelse(Contains_STY == FALSE,rowSums(sapply(sing_STY, function(query) str_count(Sequence, query))),STY_count))
      #mutate(stat=ifelse(STY_count < 0,"negative","positive"))
      
      
      gg_comb_pep_adj_dist_df <- comb_pep_adj_dist_df_neg %>% 
        group_by(Contains_STY) %>% 
        count(STY_count) %>%
        mutate(ifelse(Contains_STY==FALSE,(-1*(n)),n)) %>%
        mutate(ifelse(STY_count < 0, abs(STY_count),STY_count))
      
      plot8 <- ggplot(gg_comb_pep_adj_dist_df, aes(y= `ifelse(Contains_STY == FALSE, (-1 * (n)), n)`,x=`ifelse(STY_count < 0, abs(STY_count), STY_count)`,fill=Contains_STY)) +
        geom_col() +  geom_text(aes(label = n), vjust = -0.5,size=10) +
        labs(x = "Contains_STY", y = "STY_count",title = "Total number adjacent and separate STY counts",caption = paste(software_name,"Experiment",exp_id,acquisiton_type,sep = " ")) +
        theme_minimal() +
        scale_fill_brewer(palette = "Dark2")+ 
        theme(legend.text = element_text(size=30), 
              axis.title.x = element_text(size=30),
              axis.title.y = element_text(size=30),
              plot.title = element_text(size=35),
              plot.subtitle = element_text(size = 25),
              plot.caption = element_text(size = 25),
              legend.title=element_text(size=30),
              axis.text.x = element_text(size=20),
              axis.text.y = element_text(size = 30),
              axis.title=element_text(size=30)) +
        scale_x_continuous(limits = c(0,7),breaks = seq(from =0, to = 7, by = 1))
      
      
      plot7 <- ggplot(comb_pep_adj_dist_df_neg, aes(x= STY_count,fill=Contains_STY)) +
        geom_bar() +
        geom_text(aes(label=after_stat(count)),stat = "count", position=position_dodge(width=0.9), vjust=-0.25,size=10)+
        labs(x = "Contains_STY", y = "STY_count",title = "Total number of adjacent and non-adjacent STY residues",caption = paste(software_name,"Experiment",exp_id,acquisiton_type,sep = " ")) +
        theme_minimal() +
        scale_fill_brewer(palette = "Dark2")+
        theme(legend.text = element_text(size=30), 
              axis.title.x = element_text(size=30),
              axis.title.y = element_text(size=30),
              plot.title = element_text(size=35),
              plot.subtitle = element_text(size = 25),
              plot.caption = element_text(size = 25),
              legend.title=element_text(size=30),
              axis.text.x = element_text(size=20),
              axis.text.y = element_text(size = 30),
              axis.title=element_text(size=30)) +
        scale_x_continuous(limits = c(-6,7),breaks = seq(from =-6, to = 7, by = 1))
      
      ##### ptmRS SCORE DISTRIBUTION BASED ON ADJACENT DISTANCES ##### 
      # 
      # comb_pep_adj_dist_prop <- comb_pep_adj_dist_df_pos %>%
      #   mutate(STY_prop = as.numeric(sapply(ptmRS.Best.Site.Probabilities, function(x) {
      #     str_extract(x, "\\d+(?:\\.\\d+)?(?=;|$)") #%>%
      #     #ifelse((x<10), str_extract(x, "\\d+\\.\\d{2,3}|\\d+"), .)
      #   })))
      # 
      plot9 <- ggplot(comb_pep_adj_dist_df_pos, aes(x= extracted_values,fill=Contains_STY)) +
        geom_density( alpha=0.8) +
        labs(x = "Contains_STY", y = "STY_count",title = "Distributions of the number of phosphopeptides holding adjacent and non-adjacent STY residues",
             caption = paste(software_name,"Experiment",exp_id,acquisiton_type,sep = " ")) +
        theme_minimal() +
        scale_fill_brewer(palette = "Dark2")+
        theme(legend.text = element_text(size=30), 
              axis.title.x = element_text(size=30),
              axis.title.y = element_text(size=30),
              plot.title = element_text(size=35),
              plot.subtitle = element_text(size = 25),
              plot.caption = element_text(size = 25),
              legend.title=element_text(size=30),
              axis.text.x = element_text(size=20),
              axis.text.y = element_text(size = 30),
              axis.title=element_text(size=30))
    
    
  comb_result_pep_cor <- comb_result_pep_max_int %>% 
    filter(grepl("Correct",Pool_for_pep_merge))
  
  if(background_species == "Escherichia coli"){
    
    plot3 <- ggplot(comb_result_pep_cor, aes(x=new_col_inj,fill=pool_id.y)) +
      geom_bar() + theme_bw() +
      geom_text(color="black",aes(label=after_stat(count)),
                show.legend = F,stat = "count",size=8,
                position = position_stack(vjust = 0.5)) + 
      scale_fill_brewer(palette = "Dark2") +
      labs(title=paste("Distribution of Correctly Identified phospho-sites for each pool \n","Experiment",
                       exp_id, acquisiton_type),x="Sample id",y="Number of identified Sequence") + 
      theme_minimal() +
      theme(legend.text = element_text(size=30), 
            axis.title.x = element_text(size=30),
            axis.title.y = element_text(size=30),
            plot.title = element_text(size=35),
            plot.subtitle = element_text(size = 25),
            plot.caption = element_text(size = 25),
            legend.title=element_text(size=30),
            axis.text.x = element_text(size=20,angle = 90),
            axis.text.y = element_text(size = 30),
            axis.title=element_text(size=30))
    
  }else{
    
    plot3 <- ggplot(comb_result_pep_cor, aes(x=comb_result_pep_cor$Raw.file,fill=pool_id.y)) +
      geom_bar() + theme_bw() +
      geom_text(color="black",aes(label=after_stat(count)),
                show.legend = F,stat = "count",size=8,
                position = position_stack(vjust = 0.5)) + 
      scale_fill_brewer(palette = "Dark2") +
      labs(title=paste("Distribution of Correctly Identified phospho-sites for each pool \n","Experiment",
                       exp_id, acquisiton_type),x="Sample id",y="Number of identified Sequence") + 
      theme_minimal() +
      theme(legend.text = element_text(size=30), 
            axis.title.x = element_text(size=30),
            axis.title.y = element_text(size=30),
            plot.title = element_text(size=35),
            plot.subtitle = element_text(size = 25),
            plot.caption = element_text(size = 25),
            legend.title=element_text(size=30),
            axis.text.x = element_text(size=20,angle = 90),
            axis.text.y = element_text(size = 30),
            axis.title=element_text(size=30))
    
    
  }
    
  
  # rm(tmp,phospho_tmp)
  # if (background_species %in% c("ECOLI","Escherichia coli")){
  #   rm(ecoli_tmp)
  # }

  
  plt_obj <- ls(pattern="plot")
  plot_obj <- plt_obj[!is.na(plt_obj)]
  sapply(1:length(plot_obj),function(x) ggsave(filename = paste0("p",x,".tiff"),
                                               width = 60, height = 45, 
                                               path = paste0(file_path,curr_dir,"/outputs_with_new_script/"),
                                               units = "cm",
                                               get(plot_obj[x]),
                                               device = "tiff", #".svg"
  ))

    
  
}
  
  # exp1_wo_faims_wo_Ecoli <- read_tsv("evidence.txt")
  # #exp1_wo_faims_with_Ecoli <- read_tsv("evidence.txt")
  # 
  # # ## NECESSARY FOR COUNTING NUM OF IDENTIFIED ECOLIPEPTIDES
  # # id_ecoli_pep <- exp2_wo_faims %>% 
  # #   filter(grepl(83333,`Taxonomy IDs`)) %>% 
  # #   group_by(Sequence) %>% 
  # #   summarise(num_id_ecoli_pep = n(), mods=unique(Modifications))
  # 
  # 
  # ## NECESSARY FOR COUNTING NUM OF IDENTIFIED PHOSPHO PEPTIDES
  # id_syn_phospho_pep <- exp1_wo_faims_wo_Ecoli %>% 
  #   filter(grepl(9606,`Taxonomy IDs`)) %>% 
  #   filter(grepl("Phospho",Modifications)) 
  # 
  # 
  # 
  # #filtered_abundances<-final_results_with_common_col[rowSums(!is.na(select(quant_peptides_with_all,starts_with("abundance_"))))>0,]
  # 
  # 
  # 
  # sheet_theo_name <- "final_pep_list_with_pos"
  # theo_path = "D:/dev/Desktop_copy/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/"
  # theo_file = "pep_list_ordered_in_pool_id_with_pos.xlsx"
  # pep_list_w_theo <- read.xlsx(paste0(theo_path, theo_file), sheet = sheet_theo_name)
  # #pep_list_w_theo <- pep_list_w_theo[,-1]
  # 
  # 
  # common_col_theo <- as.data.frame(paste(pep_list_w_theo$Sequence,
  #                                        pep_list_w_theo$modified.position.in.peptide, sep = "_"))
  # colnames(common_col_theo) <- "common_col_merging"
  # pep_list_w_theo_new <- cbind(common_col_theo,pep_list_w_theo)
  # 
  # #####
  # # Combine data frame (we will continue with this for further step)
  # df_merge_phosphosite <- final_results_with_common_col %>% 
  #   left_join(pep_list_w_theo_new,by="common_col_merging")
  # #####
  # 
  # mapping_path <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_1/MaxQuant/"
  # mapping_file <- "exp1_wo_FAIMS_wo_Ecoli_mapping.txt"
  # ## MAPPING_PATH
  # raw_files_order <- as.data.frame(read_tsv(paste0(mapping_path,mapping_file)))
  # 
  # 
  # ### EXPERIMENTAL DATA SORTING BASED ON POOL ID
  # for (i in 1:dim(raw_files_order)[1]){
  #   
  #   query <- paste0("_", strsplit(raw_files_order[i,"raw_file"], "_")[[1]][2])
  #   assign(paste0(raw_files_order[i,"pool_id"],query), final_results_with_common_col %>% filter(grepl(query,`Raw file`)))
  #   
  #   rm(query)
  #   
  # }
  # pool_size <- 8
  # ## TO MERGE RAWS FILES ARE ASSOCIATED TO THE SAME POOL_ID
  # for (j in 1:pool_size){
  #   
  #   object_list <- mget(ls(pattern=paste0("pool",j)))
  #   assign(paste0("pool",j,"_all"), do.call(rbind, object_list))
  #   
  #   rm(object_list)
  #   
  #   # new_col_name <- "common_col_merging"
  #   # pool_j_all <- get(paste0("pool", j, "_all"))
  #   # pool_j_all[[new_col_name]] <- paste(pool_j_all$Stripped.Sequence, pool_j_all$value, sep = "_")
  #   #assign(paste0("pool", j, "_all"), pool_j_all)
  #   
  #   ## WITHOUT INSIDE FOR LOOP
  #   #pool1_all[,"common_col_merging"] <- paste(pool1_all$Stripped.Sequence,
  #   #pool1_all$value, sep="_")
  #   
  #   tmp <- pep_list_w_theo_new %>%
  #     filter(grepl(paste0("pool",j), pool_id)) %>%
  #     ## DON'T WORRY ABOUT WARNING MESSAGE,
  #     ## DUPLICATES WILL BE ELIMINATED IN THE NEXT STEP
  #     full_join(get(paste0("pool",j,"_all")),by="Sequence") ## TODO:TRY left_join() without filtering duplicates
  #   
  #   # ## REMOVE MULTIPLE PHOSPHOSITE THAT ARE RELATED TO THE SAME PEPTIDE
  #   assign(paste0("pool",j,"merged"),tmp %>%
  #            filter(duplicated(Sequence) == FALSE))
  #   
  #   ## UNEXPECTEDLY IDENTIFIED - FALSE POSITIVES
  #   assign(paste0("unexp_id_pool",j),get(paste0("pool",j,"merged")) %>% 
  #            filter(is.na(`Raw file`) == FALSE & is.na(Well.position) == TRUE) %>%
  #            mutate(type=paste0("unexpected_pool",j)))
  #   
  #   ## CORRECTLY IDENTIFIED - TRUE POSTIVIES
  #   assign(paste0("correct_id_pool",j),get(paste0("pool",j,"merged")) %>% 
  #            filter(is.na(`Raw file`) == FALSE & is.na(Well.position) == FALSE)%>%
  #            mutate(type=paste0("correct_pool",j)))
  #   
  #   ### MISSED IDENTIFIED - FALSE NEGATIVES
  #   assign(paste0("missed_id_pool",j),get(paste0("pool",j,"merged")) %>% 
  #            filter(is.na(`Raw file`) == TRUE & is.na(Well.position) == FALSE)%>%
  #            mutate(type=paste0("missed_pool",j)))
  #   
  # }
  # 
  # type_of_id <- c("correct_id_pool","missed_id_pool","unexp_id_pool")  
  # 
  # for (k in 1:length(type_of_id)){
  #   
  #   assign(paste0("final_",type_of_id[k],"list"),mget(ls(pattern=paste0(type_of_id[k]))))
  #   assign(paste0("ffinal_",type_of_id[k]), do.call(rbind, get(paste0("final_",type_of_id[k],"list"))))
  #   #assign(paste0("ffinal_",type_of_id[k]),cbind(get(paste0("ffinal_",type_of_id[k])),paste0(type_of_id[k])))
  #   
  # }
  # 
  # 
  # final_all_id_list <- mget(ls(pattern=paste0("ffinal")))
  # final_all_id <- do.call(rbind, final_all_id_list)
  # 
  # 
  # output_path <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_1/MaxQuant/exp1_wo_FAIMS_wo_Ecoli_OXPAL230121/modified_version_of_evidence/"
  # output_file <- "modified_processable_version_evidence.txt"
  # write_tsv(final_all_id,paste0(output_path,output_file))
  # 
  # library(viridis)
  # library(ggplot2)
  # final_all_id  %>% 
  #   group_by(type) %>% 
  #   summarise(n()) %>%
  #   separate(type,c('type_name','pool_id')) %>%
  #   ggplot( aes(fill=type_name,x=pool_id,y=`n()`)) + geom_bar(stat="identity",position = "dodge")+
  #   geom_text(aes(label=`n()`),colour = "gray", size = 4,position = position_dodge(0.9),vjust=1.6) +
  #   theme_minimal() +
  #   theme(legend.text = element_text(size=15), 
  #         axis.title.x = element_text(size = 15),
  #         axis.title.y = element_text(size = 15),
  #         plot.title = element_text(size=20),
  #         legend.title=element_text(size=15),
  #         axis.text=element_text(size=15),
  #         axis.title=element_text(size=15)) + 
  #   scale_fill_viridis(discrete = TRUE,option = "viridis")
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # intensity_match <- final_results_with_common_col %>%
  #   select(common_col_merging,Intensity,Experiment) %>% 
  #   tibble() %>%
  #   mutate(row = row_number()) %>%
  #   pivot_wider(names_from = "Experiment",
  #               values_from = "Intensity")
  # 
  # pep_list_with_pool_id <- read.xlsx(paste0("D:/dev/Desktop_copy/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/",
  #                                           "pep_list_ordered_in_pool_id_with_pos.xlsx"), sheet = "final_pep_list_with_pos")
  # 
  # pep_list_with_pool_id[,"common_col_merging"] <- paste(pep_list_with_pool_id$sequence,
  #                                                       pep_list_with_pool_id$modified.position.in.peptide,
  #                                                       sep="_")
  # 
  # 
  # 
  # 
  # 
  # 
  
  
  

  