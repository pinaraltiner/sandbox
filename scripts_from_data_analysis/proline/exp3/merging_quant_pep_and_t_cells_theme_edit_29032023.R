#
#library(PhosR)
library(stringr)
library(dplyr)
library(data.table)
library(openxlsx)
library(ggplot2)
library(tidyr)
library(readr)
library(patchwork)

final_pep_quant_analysis_bio <- function(file_path,
                                             file_name,
                                             sheet_name,
                                         actual_ratio,
                                         loc_filter_opt,
                                         loc_filter,
                                             selected_species,
                                             acquisiton_type,
                                             exp_id,
                                             background_species,
                                             numerator,
                                             software_name,
                                             test_type,
                                             exp_design,
                                             num_reps,
                                         subtitle){
  
  source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/ggplot/ggplot_functions.R")
  source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/parser_func.R")
  source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/roc_curve/new_roc_curve_generation_with_custom_threshold.R")
  
  pep_quant_parser(file_path = file_path,
                   file_name = file_name,
                   num_reps=3,
                   numerator=1,
                   sheet_name=sheet_name,
                   selected_species = selected_species,
                   background_species = background_species,
                   exp_design = experimental_design,
                   create_impute_vals = TRUE,
                   software_name = software_name,output_dir_name = "refined_output")


  output_dir_name = "refined_output"
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
  
  
  quant_peptides <- read_tsv(paste0(file_path,output_dir_name,"/refined_input",software_name,".txt"))
  impute_vals <- read_tsv(paste0(file_path,output_dir_name,"/impute_values",software_name,".txt"))
  
  quant_peptides_ECOLI <- quant_peptides %>% 
      filter(grepl(background_species,species)) #%>%
      #relocate(exp_design,.after = id)
  
  
  all_seq_dist <- quant_peptides %>%
    distinct(pep_with_pos, .keep_all = T) %>%
    #bind_rows(ecoli_seq_dist) %>%
    mutate(acq_type=acquisiton_type) %>%
    mutate(soft_name=software_name)
  
  plot13 <- gg_barplt_id_pep_count(data_set = all_seq_dist,
                                   x_df = all_seq_dist$species,
                                   fill_df = all_seq_dist$species,
                                   ymax = nrow(all_seq_dist),
                                   size_num=10,
                                   header = paste("Total number of identified unique phosphorylated", selected_species,"peptides \n and", background_species,"sequences across each sample",sep=" "),
                                   caption_lab = "NA values are removed.",
                                   x_lab = "Sample id",
                                   fill_lab =  "Sample id",
                                   y_lab = "Number of identified sequences",
                                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
  write.table(all_seq_dist, file=paste0(file_path,
                                        "Experiment",
                                        exp_id,
                                        software_name,
                                        "_number_of_unique_sequence&peptides",
                                        background_species,selected_species,
                                        ".txt"),sep = "\t",col.names = T,row.names = F)
  #################################################  

  df_id_pep <- quant_peptides %>% 
    filter(grepl(selected_species,species) & 
             grepl("Phospho",Modifications)) %>%
    tibble() %>%
    pivot_longer(cols = starts_with("E3"), 
                 values_to = "Intensity",
                 names_to = "Experiment",
                 values_drop_na = T) %>%
    separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F)
  ## SAME STRATEGIES ABOVE (3rd) WAS APPLIED TO BACKGROUND AS WELL
  df_id_pep_ecoli <- quant_peptides_ECOLI %>% 

    pivot_longer(cols = starts_with("E3"), 
                 values_to = "Intensity",
                 names_to = "Experiment",
                 values_drop_na = T) %>%
    separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F)
    
  
  
  plot1 <- gg_barplt_id_pep_count(data_set = df_id_pep,
                               x_df = df_id_pep$Sample_id,
                               fill_df = df_id_pep$Rep_id,
                               ymax = 20000,
                               size_num = 8,
                               header = "Total number of quantified phospho-peptides across each sample",
                               subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name),
                               x_lab = "Sample id",
                               fill_lab =  "Sample id",
                               y_lab = "Number of identified peptides",
                               caption_lab = "Before applying localization filtering \n After removing multiple charages.")
  
  
  plot2 <- gg_barplt_id_pep_count(data_set = df_id_pep_ecoli,
                               x_df = df_id_pep_ecoli$Sample_id,
                               fill_df = df_id_pep_ecoli$Rep_id,
                               ymax = 20000,
                               size_num = 8,
                               header = "Total number of quantified Ecoli sequences across each sample",
                               caption_lab = "There were no localization filtering applied.\n After removing multiple charages.",
                               x_lab = "Sample id",
                               fill_lab =  "Sample id",
                               y_lab = "Number of identified peptides",
                               subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
                              
  if(loc_filter_opt == TRUE){
    
    quant_phospho_peptides <- quant_peptides %>%
      filter(grepl(selected_species,species) & 
               grepl("Phospho",Modifications)) %>%
      rowwise() %>%
      filter(as.numeric(ptm_score) >= loc_filter) #%>%
      #relocate(exp_design,.after = id)
      #mutate(phospho_pos = proline_phospho_pos_extraction(Modifications))

  }else{
    quant_phospho_peptides <- quant_peptides %>%
      filter(grepl(selected_species,species) & 
               grepl("Phospho",Modifications)) #%>%

      #relocate(exp_design,.after = id)
      #filter(as.numeric(ptm_score) >= loc_filter)  %>%  
      #mutate(phospho_pos = proline_phospho_pos_extraction(Modifications)) %>%
  }
  
  ########## ########## ########## ########## ########## ########## ########## ##########  
  
 ### DATA FILTERING ###
  filtered_abundances <-  filter_NA(df = quant_phospho_peptides,samp_names = sample_names,num_allowed_NA = 1,num_expected_nonNA = 3)
  filtered_abundances_ecoli <-  filter_NA(df = quant_peptides_ECOLI,samp_names = sample_names,num_allowed_NA = 5,num_expected_nonNA = 1)
  
  filtered_abundances_ecoli_bfr_impt <- filtered_abundances_ecoli
  filtered_abundances_bfr_impt <- filtered_abundances
  ############# PROTEIN COUNT #####################
  barplt_prot_ecoli <- filtered_abundances_ecoli %>% 
    pivot_longer(cols = starts_with("E3"), 
                 values_to = "Intensity",
                 names_to = "Experiment",
                 values_drop_na = F) %>%
    separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F) %>%
    #mutate(sample_rep_id_seq = paste(Sequence, Sample_id,Rep_id, sep = "_")) %>%
    group_by(Protein,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    slice(which.max(Intensity)) %>%
    ungroup() %>%
    mutate(acquisiton_type=acquisiton_type)
  
  plot16 <- gg_barplt_id_pep_count(data_set = barplt_prot_ecoli,
                                   x_df = barplt_prot_ecoli$Sample_id,
                                   fill_df = barplt_prot_ecoli$Rep_id,
                                   ymax = 20500,
                                   size_num = 10,
                                   header = "Total number of quantified Ecoli proteins across each sample",
                                   caption_lab = "After applying data filtering",
                                   x_lab = "Sample id",
                                   fill_lab =  "Sample id",
                                   y_lab = "Number of identified proteins",
                                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))

  barplt_prot_phospho <- filtered_abundances %>% 
    pivot_longer(cols = starts_with("E3"), 
                 values_to = "Intensity",
                 names_to = "Experiment",
                 values_drop_na = F) %>%
    separate(Experiment, into = c("Exp_id","Sample_id", "Rep_id"), sep = "-",remove = F) %>%
    #mutate(sample_rep_id_seq = paste(Sequence, Sample_id,Rep_id, sep = "_")) %>%
    group_by(Protein,Experiment) %>% ## sample_rep_id_seq allowed us to keep one sequence for each sample
    slice(which.max(Intensity)) %>%
    ungroup() %>%
    mutate(acquisiton_type=acquisiton_type)
  
  plot17 <- gg_barplt_id_pep_count(data_set = barplt_prot_phospho,
                                   x_df = barplt_prot_phospho$Sample_id,
                                   fill_df = barplt_prot_phospho$Rep_id,
                                   ymax = 20500,
                                   size_num = 10,
                                   header = "Total number of quantified mouse proteins across each sample",
                                   caption_lab = "After applying data filtering",
                                   x_lab = "Sample id",
                                   fill_lab =  "Sample id",
                                   y_lab = "Number of identified proteins",
                                   subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name))
  
  barplt_prot_all <- barplt_prot_ecoli %>% bind_rows(barplt_prot_phospho) 
  write.table(barplt_prot_all, file=paste0(file_path,output_dir_name,"/Exp",exp_id,software_name,"_number_of_quantified_proteins.txt"),sep = "\t",col.names = T,row.names = F)
  
  # ## THIS RESHAPING IS ONLY FOR ELIMINATION OF MULTIPLE PHOSPHO-SITES and ECOLI PEPTIDES
  #  ## ELIMINATION STEP IS NOT NECESSARY FOR ECOLI, 1st STRATEGY can be used only (this will decrease lines of code)
  # barplt_df_wide <- barplt_df %>% 
  #   select(pep_with_pos, accession, sample_ids,intensity,spectrum_title) %>%
  #   pivot_wider(names_from = "sample_ids",values_from = "intensity") %>%
  #   separate(spectrum_title,into = c(paste0("tmp",1:6),"spec_tit"),sep = ";",remove = FALSE) %>%
  #   select(!c(paste0("tmp",1:6))) %>%
  #   separate(spec_tit,into = c("tmp","raw_file"),sep = ":") %>%
  #   select(!tmp) # %>%
  #   #separate(raw_file_path, into = c(paste0("tmp",1:7),"raw_file"),sep = "/") %>%
  #   #select(!c(paste0("tmp",1:7))) %>%
  #   #mutate(raw_file=str_remove(raw_file, "\"")) %>%
  #   #mutate(raw_file=str_remove(raw_file, ".raw"))
  # 
  # barplt_df_ecoli_wide <- barplt_df_ecoli %>% 
  #   select(sequence, accession, sample_ids,intensity,spectrum_title) %>%
  #   pivot_wider(names_from = "sample_ids",values_from = "intensity") %>%
  #   separate(spectrum_title,into = c(paste0("tmp",1:6),"spec_tit"),sep = ";",remove = FALSE) %>%
  #   select(!c(paste0("tmp",1:6))) %>%
  #   separate(spec_tit,into = c("tmp","raw_file"),sep = ":") %>%
  #   select(!tmp)
  
  filtered_abundances[,"na_val"] <- apply(X = !is.na(select(filtered_abundances,contains(exp_design))), MARGIN = 1, FUN = sum)
  
  phospho_completeness <- filtered_abundances %>% 
    count(na_val) %>%
    mutate(data_complete=((n/dim(filtered_abundances)[1])*100)) %>% mutate(species=selected_species)
  
  filtered_abundances_ecoli[,"na_val"] <- apply(X = !is.na(select(filtered_abundances_ecoli,contains(exp_design))), MARGIN = 1, FUN = sum)
  
  completeness <- filtered_abundances_ecoli %>% 
    count(na_val) %>%
    mutate(data_complete=((n/dim(filtered_abundances_ecoli)[1])*100)) %>%
    mutate(species=background_species) %>%
    bind_rows(phospho_completeness)
  
  write.table(completeness,file = paste0(file_path,output_dir_name,"/data_completeness",acquisiton_type,software_name,".txt"))
  
  plot23 <- ggplot(completeness, aes(x=na_val,y=data_complete,color=species)) + geom_point(size=2.5) +
    geom_line(size=2)+
    scale_x_reverse(limits=c(18,0),breaks=seq(0, 18, by = 2)) +
    ## If you look for is.na() in apply function, you should use the one below:
    #### scale_x_continuous(limits=c(0,18),breaks=seq(0, 18, by = 2)) +
    scale_y_continuous(limits = c(0,100), breaks = seq(from =0, to=100,by=10)) +
    theme_minimal() +
    theme(legend.text = element_text(size=30), 
          axis.title.x = element_text(size=30),
          axis.title.y = element_text(size=30),
          plot.title = element_text(size=35),
          plot.subtitle = element_text(size = 20),
          legend.title=element_text(size=30),
          axis.text=element_text(size=30),
          axis.title=element_text(size=30),
          strip.text.x = element_text(
            size = 15
          )
    ) + scale_color_manual(values = c("#3182bd","#a6bddb"))+
    ggtitle(label = paste("Data completeness of",background_species,"and",selected_species)) +
    labs(x="n Sample", y="% of peptides",subtitle = paste("Experiment",exp_id,software_name,acquisiton_type,"\n",file_name))
  
  ## DENSITY PLOT OF BEFORE IMPUTATION 

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
  
  colnames(abundances_rowMeans) <- paste0("Mean Abundance of A",1:sample_size)
  
  colnames(abundances_ecoli_rowMeans) <- paste0("Mean Abundance of A",1:sample_size)
  
  quant_peptides_ECOLI_density_plt <- filtered_abundances_ecoli %>%
    select(!starts_with("E")) %>%
    bind_cols(abundances_ecoli_rowMeans) %>%
    tibble() %>%
    pivot_longer(cols = starts_with("mean"), 
                 values_to = "Intensity",
                 names_to = "Experiment",
                 values_drop_na = T)# %>%
    #mutate(sample_id_seq = paste(pep_with_pos, Experiment, sep = "_"))
  
  ecoli_density_plt<- quant_peptides_ECOLI_density_plt %>%
    select(c(Experiment,Intensity,species))
  
  #### MEAN ABUNDANCE RATIO WITH  DENSITY PLOT ####
      ### BEFORE IMPUTATION ###
  quant_phospho_density_plt <- filtered_abundances %>%
    select(!starts_with("E")) %>%
    bind_cols(abundances_rowMeans) %>% 
    #rename_with(~ paste0("mean_abun",1:5), matches("^row")) %>%
    tibble() %>% #mutate(pep_with_pos = sequence) %>% ###  At this stage, no need for phospho-position#   
    pivot_longer(cols = starts_with("mean"), 
                 values_to = "Intensity",
                 names_to = "Experiment",
                 values_drop_na = T)# %>%
    #mutate(sample_id_seq = paste(pep_with_pos, Experiment, sep = "_"))
  
  density_df <-quant_phospho_density_plt %>%
    select(c(Experiment,Intensity,species))%>%
    bind_rows(ecoli_density_plt)
    
  plot3 <- gg_density(data_set = density_df, 
                      x_df = density_df$Intensity,
                      fill_df = density_df$species,
                      color_df = NULL,
                      header="Distribution of mean abundance of every sample before imputation",
                      facet_df = "Experiment",
                      x_lab = "log10(intensities)",
                      color_lab= "",
                      fill_lab = "species",
                      subtitle_txt = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name)) +scale_fill_brewer(palette = "Set1")

  ### RATIO SUPPRESION ASSESSMENT ### 
  
  plot15 <- gg_raincloud(data_set = density_df,
               x_df = density_df$Experiment,
               y_df = density_df$Intensity,
               fill_df = density_df$species,
               header = "Distribution of mean abundance of every sample before imputation",
               x_lab = "Sample Names",
               y_lab = " Density of log10(Mean Abundance)",
               fill_lab = "Sample Names",
               caption_lab = "",
               subtitle_txt = "") + scale_fill_brewer(palette = "Set1")

  ###################################
  
  is.imputed_df_syn <- filtered_abundances %>%  pivot_longer(cols = starts_with(exp_design),
                                                        names_to = "Experiment",
                                                        values_to = "Intensity",
                                                        values_drop_na = F) %>%
    mutate(is.imputed=FALSE) %>%
    mutate(is.imputed=ifelse(is.na(Intensity), TRUE,is.imputed)) 

  is.imputed_df_ecoli <- filtered_abundances_ecoli %>%  pivot_longer(cols = starts_with(exp_design),
                                                                names_to = "Experiment",
                                                                values_to = "Intensity",
                                                                values_drop_na = F) %>%
    mutate(is.imputed=FALSE) %>%
    mutate(is.imputed=ifelse(is.na(Intensity), TRUE,is.imputed)) 
  
  impute_values <- as.numeric(impute_vals$x)
  # Impute missing values
  for (j in 1:length(impute_values)){
    # Number NA
    #num_NA <- length(abundances_for_impute_all[,j+2][is.na(abundances_for_impute_all[,j+2])])
    
    filtered_abundances[,j+2][is.na(filtered_abundances[,j+2])] <- impute_values[j]
    filtered_abundances_ecoli[,j+2][is.na(filtered_abundances_ecoli[,j+2])] <- impute_values[j]
    
    #abundances_for_impute_all[,j+2][is.na(abundances_for_impute_all)[,j+2]] <- impute_values[j]
    # After imputation number of imputed values
    #num_imp <-length(abundances_for_impute_all[,j+2][(abundances_for_impute_all[,j+2]==impute_values[j])])
    
    # This is verification of imputation is done successfully
    # Because we expect to see that number of imputed values should be the same amount as number of NA
    #print(setequal(num_NA,num_imp))
    #print(num_NA)
    #print(num_imp)
  }
  
  abundances_all_aft_imputation <- filtered_abundances %>%
    #rename_with(~ paste0("pep_with_pos"), matches("^seq")) %>%
    bind_rows(filtered_abundances_ecoli) 
  
  ### DISTRIBUTION OF IMPUTED VALUES ACROSS non-NA values
  is.imputed_df <- is.imputed_df_ecoli %>% 
    bind_rows(is.imputed_df_syn)
  
  imputed_dataset <- abundances_all_aft_imputation %>% 
    pivot_longer(cols = starts_with(exp_design),
                 names_to = "Experiment",
                 values_to = "Intensity",
                 values_drop_na = F) %>%
    left_join(is.imputed_df,by=c("pep_with_pos","Experiment")) %>%
    mutate(Intensity.y=ifelse(is.na(Intensity.y),0,Intensity.y))
  
  write.table(imputed_dataset,file=paste0(file_path,output_dir_name,"/imputed_dataset",software_name,".txt"),sep = "\t",row.names = F)
  
  p21  <- imputed_dataset %>% filter(grepl(selected_species,species.x)) %>% 
    ggplot(aes(x=log2(Intensity.x),fill=is.imputed)) + geom_histogram(bins = 30) +
    theme_minimal() + scale_fill_brewer(palette = "Set1",direction = -1) +
    theme(legend.text = element_text(size = 45), #aspect.ratio=6.5/11, 
          axis.title.x = element_text(size = 45),
          axis.title.y = element_text(size = 45),
          plot.title = element_text(size = 55),
          legend.title = element_text(size = 45),
          axis.text.x = element_text(size = 45),
          axis.title = element_text(size = 45),
          axis.text.y = element_text(size = 45),
          plot.subtitle = element_text(size = 45)) +
    labs( y= "Count of Intensity", x="Intenisity",
          title = paste("Distribution of imputed values \n",selected_species), 
          subtitle = paste('Experiment 2 ', acquisiton_type, " data processed by ", software_name))
  
  p22 <- imputed_dataset %>% filter(grepl(background_species,species.x)) %>% 
    ggplot(aes(x=log2(Intensity.x),fill=is.imputed)) + geom_histogram(bins = 30) +
    theme_minimal() + scale_fill_brewer(palette = "Set1",direction = -1) +
    theme(legend.text = element_text(size = 45), #aspect.ratio=6.5/11, 
          axis.title.x = element_text(size = 45),
          axis.title.y = element_text(size = 45),
          plot.title = element_text(size = 55),
          legend.title = element_text(size = 45),
          axis.text.x = element_text(size = 45),
          axis.title = element_text(size = 45),
          axis.text.y = element_text(size = 45),
          plot.subtitle = element_text(size = 45)) +
    labs( y= "Count of Intensity", x="Intenisity",
          title = paste(background_species)) 
  #subtitle = paste('Experiment 2 ', acquisiton_type, " data processed by ", software_name))
  
  plot22 <- p21/p22
  
  ## ADDITIONAL IMPUTATION METHOD with MICE()
  # library(tidyverse)
  # library(tidyr)
  # library(mice)
  # quant_phospho_peptides <- as.data.frame(quant_phospho_peptides)
  # rownames(quant_phospho_peptides) <- paste0(quant_phospho_peptides$pep_with_pos,quant_phospho_peptides$species,"@",1:nrow(quant_phospho_peptides))
  # 
  # intensities <- quant_phospho_peptides %>%
  #   select(starts_with(exp_design))
  # 
  # #intensities_short <- intensities[1:10,]
  # barplt_df_ecoli_wide <- as.data.frame(barplt_df_ecoli_wide)
  # rownames(barplt_df_ecoli_wide) <- paste0(barplt_df_ecoli_wide$sequence,barplt_df_ecoli_wide$species,"@",1:nrow(barplt_df_ecoli_wide))
  # 
  # intensitiesECOLI <- barplt_df_ecoli_wide %>%
  #   select(starts_with(exp_design))
  # 
  # imp_intensities <-  mice(intensities,m=5,maxit=50,meth='cart',seed=500)
  # imp_intensitiesECOLI <- mice(intensitiesECOLI,m=5,maxit=50,meth='cart',seed=500)
  # 
  # completeData <- complete(imp_intensities,2)
  # completeDataECOLI <- complete(imp_intensitiesECOLI,2)
  # #pattern <- md.pattern(select(quant_phospho_peptides,starts_with(exp_design)))
  # #library(VIM)
  # #aggr_plot <- aggr(intensities, col=c('navyblue','red'),
  #                   #numbers=TRUE, sortVars=TRUE,
  #                   #labels=names(intensities), cex.axis=.7,
  #                   #gap=3, ylab=c("Histogram of missing data","Pattern"))
  # 
  # completeDataECOLIs <-completeDataECOLI %>%
  #   mutate(tmp=rownames(completeDataECOLI)) %>%
  #   separate(tmp,into=c("pep_with_pos","indx"),sep="@") %>%
  #   select(!indx) %>%mutate(species="ECOLI")
  #   
  # abundances_all_aft_imputation <- completeData %>% 
  #   mutate(tmp=rownames(completeData)) %>%
  #   separate(tmp, into = c("pep_with_pos","species"),sep = "@") %>%
  #   separate(pep_with_pos, into = c("pep","pos","species"),sep = "_") %>%
  #   mutate(pep_with_pos=paste0(pep,"_",pos)) %>%
  #   select(!c(pep,pos)) %>%
  #   bind_rows(completeDataECOLIs)
  # 
  # 
  # # UNNECESSARY TO KEEP DATA BEFORE IMPUTATION
  # abundances_all_aft_imputation <- barplt_df_ecoli_wide %>%
  #   rename_with(~ paste0("pep_with_pos"), matches("^seq")) %>%
  #   bind_rows(barplt_df_wide) 
    #select(sequence, accession) %>% #spectrum_title
    #bind_cols(abundances_for_before_impt)
  
  # Take log10 
  #log_10_abundances_only_for_impute <- log10(abundances_only_for_impute)
  #colnames(log_10_abundances_only_for_impute) <- paste0(colnames(abundances_only_for_impute),"_log10")
  
  ## No need is right now. I will cont with "log_10_abundances_only_for_impute" for rowMeans and FC
  #log_10_filtered_abundances <- cbind(filtered_abundances,log_10_abundances_only_for_impute)
  
  # Take mean of triplicates of each sample 
  # Ask sample_size additional parameter
  
  
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
  for(An in 1:cols){
    filtered_abundances_rowMeans[,paste0("exp_FC_A1/A",An)] <- filtered_abundances_rowMeans[,1]/filtered_abundances_rowMeans[,An]
    
  }
  rmv_col <- paste0("exp_FC_A",numerator,"/A",numerator)
  filtered_abundances_rowMeans <- filtered_abundances_rowMeans %>% select(!rmv_col)
  
  # To calculate all binary combination in the data frame
  #mat <- do.call(cbind, lapply(cols, function(xj) 
  #  sapply(cols, function(xi) (filtered_abundances_rowMeans[, xj]/(filtered_abundances_rowMeans[, xj])))))
  #colnames(mat) <-  outer(names(filtered_abundances_rowMeans), names(filtered_abundances_rowMeans), paste0)
  
  final_imputed_data <- cbind(abundances_all_aft_imputation, filtered_abundances_rowMeans,filtered_abundances_log10,filtered_abundances_log10_rowMeans) #filtered_abundances
  final_imputed_data_ecoli <- final_imputed_data  %>% filter(!grepl(selected_species, species))
  
  write.table(final_imputed_data, file = paste0(file_path,output_dir_name,"/final_imputed&filtered_data",software_name,".txt"),sep = "\t",row.names = F)
  
  df_mean_ab_after_impt <- final_imputed_data %>% 
    select(contains("aft_imp") | contains("species")) %>% # contains("accession")
    tibble() %>% 
    #separate(accession, into = c("uniprot_id", "species"), remove = F) %>%
    pivot_longer(cols = contains("aft_imp"),
                 names_to = "Mean_abundance",
                 values_to = "values")
    
  df_FC_ratio_after_impt <- final_imputed_data %>% 
    select(starts_with("exp_")| contains("species")) %>%# contains("accession")
    #separate(accession, into = c("uniprot_id", "species"), remove = F) %>%
    tibble() %>% 
    pivot_longer(cols = starts_with("exp_"),
                 names_to = "exp_FC",
                 values_to = "values") %>%
    mutate(actual_ratio_val= case_when(grepl(comparisons[1],exp_FC) ~ actual_ratio[1],
                                       grepl(comparisons[2],exp_FC) ~actual_ratio[2],
                                       grepl(comparisons[3],exp_FC) ~actual_ratio[3],
                                       grepl(comparisons[4],exp_FC) ~actual_ratio[4])) %>%
    mutate(log2_act_val=log2(actual_ratio_val)) %>%
    mutate(log2_exp_val=log2(values)) 
    
  df_FC_ratio_after_impt_mouse <- df_FC_ratio_after_impt %>%
    filter(grepl(selected_species,species))
  
  median_val <- df_FC_ratio_after_impt_mouse %>% group_by(exp_FC) %>%
    summarise(exp_median=median(log2_exp_val))
  
  #test <- df_FC_ratio_after_impt %>% group_by(exp_FC) %>% summarise(min_val=min(log2_exp_val),max_val=max(log2_exp_val))
  
  plot21 <-gg_quant_ratio_acc(data_set = df_FC_ratio_after_impt_mouse,
                              x_df = df_FC_ratio_after_impt_mouse$log2_act_val,
                              y_df = df_FC_ratio_after_impt_mouse$log2_exp_val,
                              color_df = df_FC_ratio_after_impt_mouse$exp_FC,
                              median_col = "red",
                              header="Quantitative Ratio Assessment",
                              x_lab="log2(Actual Ratio)",
                              y_lab="log2(Experimental Ratio",
                              color_lab="Comparisons",
                              subtitle_txt=paste("Experiment", exp_id,"data acquired from", acquisiton_type,"processed by ",software_name)
  ) #scale_color_manual(values = c("#08519c","#3182bd","#6baed6","#a6bddb"))
  #+ geom_errorbar(data = test,aes(ymin = test$min_val, ymax=  test$max_val,color=test$exp_FC), width=0.5) 
  
  df_FC_ratio_after_impt_final <- df_FC_ratio_after_impt %>% mutate(acquisition=acquisiton_type) %>% mutate(software_name=software_name)
  
  write.table(df_FC_ratio_after_impt_final,file = paste0(file_path,output_dir_name,"/df_FC_ratio_after_impt",software_name,".txt"),sep = "\t",row.names =F )
  #
  
  plot4 <- gg_density(data_set = df_mean_ab_after_impt, 
                   x_df = df_mean_ab_after_impt$values,
                   fill_df = df_mean_ab_after_impt$species,
                   color_df = NULL,
                   header="Distribution of mean abundance of every sample after imputation",
                   facet_df = "Mean_abundance",
                   x_lab = "log10(values)",
                   color_lab= "",
                   fill_lab = "Sample Names",
                   subtitle_txt = "") +scale_fill_brewer(palette = "Set1")
  
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
  
  ########################
  ############################## ############################
  ############      ############    CV ASSESSMENT BEFORE & AFTER IMPUTATION      ############      ############  
  ### BEFEORE  ### 
  
  df_merge_before_imp <- filtered_abundances_bfr_impt %>%
    bind_rows(filtered_abundances_ecoli_bfr_impt)
  
  
  abundances_all_before_impt <- NULL
  
  for (k in 1:(sample_size)){
    assign(paste0("df_merge_before_imp_A",k),as.data.frame(rowMeans(df_merge_before_imp  %>% select(contains(paste0("A",k))))))
    
    abundances_all_before_impt<- bind_cols(abundances_all_before_impt,get(paste0("df_merge_before_imp_A",k)))
  }
  colnames(abundances_all_before_impt) <- paste0("mean_abundance_A",1:(sample_size))
  
  df_merge_before_imp_all <- df_merge_before_imp %>% bind_cols(abundances_all_before_impt)
  
  df_before_impt_CV <- CV_calculator(data = df_merge_before_imp_all,
                                     abundance_col = "mean_abundance_A",
                                     iter = 5,replicate = num_reps,
                                     exp_id=3,
                                     acquisiton_type = acquisiton_type,
                                     software_name = software_name)
  
  write.table(df_before_impt_CV,file = paste0(file_path,output_dir_name,"/CV_calculation_before_imputation",acquisiton_type,software_name,".txt"))
  
  ############      ############                     ############      ############  
  ### AFTER  ### 
  
  df_merge_aft_impt <- final_imputed_data %>%  #mutate(Pool=ifelse(Pool=="Diluted",paste0(Pool,"_",isomericity),Pool)) %>%
    select(!contains("log10") &!contains("exp_FC"))
  
  df_aft_impt_CV <- CV_calculator(data = df_merge_aft_impt,abundance_col = "mean_abundances_aft_imp_A",
                                  iter = 5,replicate = num_reps,
                                  exp_id = 3,
                                  acquisiton_type = acquisiton_type,
                                  software_name = software_name)
  
  write.table(df_aft_impt_CV,file = paste0(file_path,output_dir_name,"/CV_calculation_after_imputation",acquisiton_type,software_name,".txt"))
  
  library(ggridges)
  library(paletteer) 
  
  
  p18 <- df_before_impt_CV %>%
    ggplot(aes(x=CV_values,y=df_before_impt_CV$species,
               color=CV_samples
    )) + 
    geom_density_ridges(scale = 0.9,
                        jittered_points = TRUE,
                        position = position_points_jitter(width = 0.05, height = 0),
                        point_shape = '*',
                        point_size = 3,
                        point_alpha = 1,
                        alpha = 0,linewidth=1.5)+ 
    #scale_fill_brewer(palette ="PuOr")+
    scale_color_manual(values = c("#08519c","#3182bd","#6baed6","#a6bddb","#d0d1e6"))+
    #scale_color_manual(values = c("#7fc97f","#beaed4","#fb8072","#fdc086","#386cb0"))+
    #scale_color_manual(values = c("#6A51A3","#6A51A3","#FD8D3C","#FD8D3C"))+ #"grey68"
    scale_linetype_manual(values = c("solid","dashed","solid","dashed"))+
    theme_minimal() +
    scale_y_discrete(expand = expand_scale(mult = c(0, 0)))+
    #facet_wrap(~Pool) +
    theme(legend.text = element_text(size = 45), #aspect.ratio=6.5/11, 
          axis.title.x = element_text(size = 45),
          axis.title.y = element_text(size = 45),
          plot.title = element_text(size = 55),
          legend.title = element_text(size = 45),
          axis.text.x = element_text(size = 45),
          axis.title = element_text(size = 45),
          axis.text.y = element_text(size = 45),
          plot.subtitle = element_text(size = 45)
    ) +
    labs( y= "Intensity", x="Percentage of CVs",
          #title = paste("Experiment 2 Percentage of CVs", acquisiton_type, " data processed by ", software_name), 
          subtitle = paste("Before imputation")) +
    #xlim(0,100)
    xlim(0,round(max(df_before_impt_CV$CV_values))) 
  # scale_y_continuous(
  # limits = c(0,100),
  # breaks = seq(0, 100,10))
  #####################
  p19 <- df_aft_impt_CV %>%
    ggplot(aes(x=CV_values,y=species,
               color=CV_samples
    )) + 
    geom_density_ridges(scale = 0.9,
                        jittered_points = TRUE,
                        position = position_points_jitter(width = 0.05, height = 0),
                        point_shape = '*',
                        point_size = 3,
                        point_alpha = 1,
                        alpha = 0,linewidth=1.5)+ 
    #scale_fill_brewer(palette ="PuOr")+
    scale_color_manual(values = c("#08519c","#3182bd","#6baed6","#a6bddb","#d0d1e6"))+
    #scale_color_manual(values = c("#6A51A3","#6A51A3","#FD8D3C","#FD8D3C"))+ #"grey68"
    scale_linetype_manual(values = c("solid","dashed","solid","dashed"))+
    theme_minimal() +
    scale_y_discrete(expand = expand_scale(mult = c(0, 0)))+
    #facet_wrap(~Pool) +
    theme(legend.text = element_text(size = 45), #aspect.ratio=6.5/11, 
          axis.title.x = element_text(size = 45),
          axis.title.y = element_text(size = 45),
          plot.title = element_text(size = 55),
          legend.title = element_text(size = 45),
          axis.text.x = element_text(size = 45),
          axis.title = element_text(size = 45),
          axis.text.y = element_text(size = 45),
          plot.subtitle = element_text(size = 45)) +
    labs( y= "Intensity", x="Percentage of CVs",
          #title = paste("Experiment 2 Percentage of CVs", acquisiton_type, " data processed by ", software_name), 
          subtitle = paste("After imputation")) +
    xlim(0,round(max(df_before_impt_CV$CV_values)))
  #scale_y_continuous(
  #limits = c(0,max(df_aft_impt_CV$CV_values)), 
  #breaks = seq(0, max(df_aft_impt_CV$CV_values),10))
  
  plot18 <- p18/p19
  #############
  p20 <- df_before_impt_CV %>%
    ggplot(aes(x=CV_samples,y=df_before_impt_CV$CV_values,
               color=CV_samples,fill=species
    )) + 
    geom_boxplot(alpha=0.65)+
    facet_wrap(~species) +
    scale_fill_brewer(palette ="Set1")+
    scale_color_manual(values = c("#08519c","#3182bd","#6baed6","#a6bddb","#d0d1e6"))+
   
    theme_minimal() +
    
    #facet_wrap(~Pool) +
    theme(legend.text = element_text(size = 45), #aspect.ratio=6.5/11, 
          axis.title.x = element_text(size = 45),
          axis.title.y = element_text(size = 45),
          plot.title = element_text(size = 55),
          legend.title = element_text(size = 45),
          axis.text.x = element_text(size = 45),
          axis.title = element_text(size = 45),
          axis.text.y = element_text(size = 45),
          plot.subtitle = element_text(size = 45),
          strip.text = element_text(size=30)
    ) +
    labs( y= "Percentage of CV (%)",x="Sample ID",
          title = paste("Experiment 3 Percentage of CVs", acquisiton_type, " data processed by ", software_name), 
          subtitle = paste("Before imputation")) +
    #xlim(0,100)
    ylim(0,round(max(df_before_impt_CV$CV_values))) 
  # scale_y_continuous(
  # limits = c(0,100),
  # breaks = seq(0, 100,10))
  #####################
  p21 <- df_aft_impt_CV %>%
    ggplot(aes(x=CV_samples,y=df_aft_impt_CV$CV_values,
               color=CV_samples,fill=species
    )) + 
    geom_boxplot(alpha=0.65)+
    facet_wrap(~species) +
    scale_fill_brewer(palette ="Set1")+
    scale_color_manual(values = c("#08519c","#3182bd","#6baed6","#a6bddb","#d0d1e6"))+
    
    theme_minimal() +
    
    #facet_wrap(~Pool) +
    theme(legend.text = element_text(size = 45), #aspect.ratio=6.5/11, 
          axis.title.x = element_text(size = 45),
          axis.title.y = element_text(size = 45),
          plot.title = element_text(size = 55),
          legend.title = element_text(size = 45),
          axis.text.x = element_text(size = 45),
          axis.title = element_text(size = 45),
          axis.text.y = element_text(size = 45),
          plot.subtitle = element_text(size = 45),
          strip.text = element_text(size=30)
    ) +
    labs(y= "Percentage of CV (%)", x="Sample ID",
          title = paste("Experiment 3 Percentage of CVs", acquisiton_type, " data processed by ", software_name), 
          subtitle = paste("After imputation")) +
    #xlim(0,100)
    ylim(0,round(max(df_aft_impt_CV$CV_values))) 
  
  plot19 <- p20/p21
###########
  plot5 <- gg_density(data_set = df_FC_ratio_after_impt,
                   x_df = df_FC_ratio_after_impt$values,
                   fill_df = df_FC_ratio_after_impt$exp_FC,
                   color_df = df_FC_ratio_after_impt$species,
                   header="Distribution of Fold change Ratio of every sample after imputation",
                   facet_df = "exp_FC",
                   x_lab = "log10(values)",
                   color_lab= "",
                   fill_lab = "Sample Names",
                   subtitle_txt = "") + scale_color_brewer(palette = "Set1")
  # 

  ### BOX-PLOT: Experimental Quantity Ratio of Phospho Peptides  
  
  plot6 <- gg_boxplt_exp_ratio(data_set = df_FC_ratio_after_impt, 
                            x_df = df_FC_ratio_after_impt$exp_FC,
                            y_df = df_FC_ratio_after_impt$values,
                            fill_df = df_FC_ratio_after_impt$species,
                            header="Experimental Quantity Ratio of T-cell Phospho Peptides",
                            x_lab="Sample Names",
                            y_lab="Abundance Ratios",
                            fill_lab = "Sample Names",
                            subtitle_txt = "") +scale_fill_brewer(palette = "Set1")
  
  
  ### HALF-BOX-PLOT & HALF-SCATTER-PLOT: Experimental Quantity Ratio of Synthetic Peptides  
  library(gghalves)
  
  plot7 <- gg_half_boxplt_exp_ratio(data_set = df_FC_ratio_after_impt, 
                                 x_df = df_FC_ratio_after_impt$exp_FC,
                                 y_df = df_FC_ratio_after_impt$values,
                                 fill_df = df_FC_ratio_after_impt$species,
                                 header="Experimental Quantity Ratio of T-cell Phospho Peptides with Background",
                                 x_lab="Sample Names",
                                 y_lab="Abundance Ratios",
                                 fill_lab = "Sample Names",
                                 subtitle_txt = "") +scale_fill_brewer(palette = "Set1")
  
  ### VIOLIN-PLOT: Experimental Quantity Ratio of Synthetic Peptides   
  
  ### TODO: fix y scaling without trimming 
  plot8 <- gg_violin_exp_ratio(data_set = df_FC_ratio_after_impt, 
                            x_df = df_FC_ratio_after_impt$exp_FC,
                            y_df = df_FC_ratio_after_impt$values,
                            fill_df = df_FC_ratio_after_impt$species,
                            header="Experimental Quantity Ratio of T-cell Phospho Peptides with Background",
                            x_lab="Sample Names",
                            y_lab="Abundance Ratios",
                            fill_lab = "Sample Names",
                            trim=TRUE,
                            subtitle_txt = "") +scale_fill_brewer(palette = "Set1")
  
  
  
  # ############################ # Do t-test # ############################
  # 
  # # The code below does t-test for each row. Because of that, multiple test correction (like BH, Bonferoni)
  # p_values_for_all_ratio <- NULL
  # for(k in 2:sample_size){
  #   
  #   assign(paste0("t_test_res_A1_to_A",k) , lapply(na.omit(correct_identifed_peps$mean_abundances_A1),
  #                                                  na.omit(correct_identifed_peps[,paste0("mean_abundances_A",k)], 
  #                                                          function(x,y) t.test(x,y,alternative = "two.sided", var.equal = TRUE))))
  #   assign(paste0("p_values_A1_vs_A",k), cbind(p_values_for_all_ratio, get(paste0("t_test_res_A1_to_A",k))$p.value))
  #   
  # }
  # 
  # #adjusted_p_values <- as.data.frame(p.adjust(p_values_A1_vs_A5, method = "BH", n = length(p_values_A1_vs_A5)))
  # adjusted_p_values <- lapply(p_values_for_all_ratio, function(x) p.adjust (x, method = "BH", n = length(x)))
  # 
  # #plot(correct_identifed_peps$mean_abundances_log10_A1,correct_identifed_peps$mean_abundances_log10_A5, pch = 16, col = "blue")
  # #abline(h = mean(na.omit(correct_identifed_peps$mean_abundances_log10_A1)) - mean(na.omit(correct_identifed_peps$mean_abundances_log10_A5)), col = "red")
  # 
  # boxplot(correct_identifed_peps$mean_abundances_log10_A1,correct_identifed_peps$mean_abundances_log10_A5)
  # ggplot(correct_identifed_peps, aes(x = mean_abundances_log10_A1, y = mean_abundances_log10_A5)) +
  #   geom_point(aes(color = "red", size = 5))# +
  # #scale_color_discrete(name = "") +
  # #scale_size_discrete(name = "Size")
  # # 
  # # Take -log10() of results
  # ggplot(log_10_filtered_abundances_rowMeans_A1_A5,aes(x=log_10_filtered_abundances_rowMeans_A1_A5$mean_abundances_log10_A1,
  #                                                      y =log_10_filtered_abundances_rowMeans_A1_A5$mean_abundances_log10_A5,
  # )) +
  #   geom_point()+
  #   #facet_wrap(vars(df_for_figure$common_col_for_merging))  + 
  #   #geom_smooth(method = "lm", colour = "green", fill = "green") +
  #   theme_light()
  # 
  # library(gginference)
  # ggttest(t.test(na.omit(correct_identifed_peps$mean_abundances_log10_A1),na.omit(correct_identifed_peps$mean_abundances_log10_A5), alternative = "two.sided", var.equal = TRUE))
  # 
  # 
  # 
  # 
  
  ttest_func <- function(x, y) {
    # if (sum(!is.na(x)) < 2 | sum(!is.na(y)) < 2) {
    #   return(NA)
    # }else{
    #   
    # }
    t.test(x, y,var.equal=TRUE,alternative = c("two.sided"))$p.value
  }
  
  wilcox_func <- function(x, y) {
    # if (sum(!is.na(x)) < 2 | sum(!is.na(y)) < 2) {
    #   return(NA)
    # }else{
    #   
    # }
    wilcox.test(as.numeric(x), as.numeric(y),alternative = c("two.sided"))$p.value
  }
  ## TODO: ADD LIMMA
  stat_analysis <- final_imputed_data %>%
    select(pep_with_pos,species, starts_with("log10_") | starts_with("mean_log10_") | starts_with("exp_FC")) #spectrum_title,accession
  
  rownames(stat_analysis) <- paste0(stat_analysis$pep_with_pos,"@",stat_analysis$species,"@",(1:nrow(stat_analysis))) #stat_analysis$spectrum_title,"@"
  # if(test_type== "t.test" | test_type== "wilcoxon"){
  #     library(multtest)
  #   all_pvalues <- NULL
  #   all_adjust_pval <- NULL
    # for (i in 2:sample_size){
    #   p_values_tmp <- NULL
    #   
    #   for(j in 1:dim(stat_analysis)[1]){
    #     
    #     if(test_type=="t.test"){
    #       p_values_tmp[j] <- ttest_func(select(stat_analysis,contains("A1_") & contains("log10_"))[j,],     ### FOR DIFFERENT KIND OF EXP SETUP, 
    #                                     select(stat_analysis,contains(paste0("A",i,"_")) & contains("log10_"))[j,] ) ## It should be defined as an input.
    #       
    #     }else if(test_type=="wilcoxon"){
    #       p_values_tmp[j] <- wilcox.test(select(stat_analysis,contains("A1_") & contains("log10_"))[j,],     ### FOR DIFFERENT KIND OF EXP SETUP, 
    #                                      select(stat_analysis,contains(paste0("A",i,"_")) & contains("log10_"))[j,])
    #       ## It should be defined as an input.
    #     }
    #     
    #   }
    #   p_values_tmp <- as.data.frame(p_values_tmp)
    #   colnames(p_values_tmp) <- paste0("pvalues_A1","/","A",i)
    #   all_pvalues <- bind_cols(all_pvalues,p_values_tmp)
    #   rownames(all_pvalues) <- row.names(stat_analysis)
    #   rm(p_values_tmp)
    #   
    #   #### p-value adjust at a given FDR ####
    #   #procs<-c("Bonferroni","Holm","Hochberg","SidakSS","SidakSD","BH","BY","ABH","TSBH")
    #   adjust_pval_tmp <-  mt.rawp2adjp(all_pvalues[,paste0("pvalues_A1","/","A",i)],
    #                                    proc="BH", alpha=0.05)
    #   qval <- data.frame(adjust_pval_tmp$adjp,adjust_pval_tmp$index)[order(adjust_pval_tmp$index),2]
    #   qval <- as.data.frame(qval)
    #   colnames(qval) <- paste0("adjust_pval_A1","/","A",i)
    #   all_adjust_pval <- bind_cols(all_adjust_pval,qval)
    #   rownames(all_adjust_pval) <- row.names(stat_analysis)
    #   rm(adjust_pval_tmp,qval)
    #   
    # }
    if (test_type == "t.test" || test_type == "wilcoxon") {
      library(multtest)
        # Pre-allocate memory for results
        
        num_iterations <- sample_size - 1
        all_pvalues <- matrix(NA, nrow(stat_analysis), num_iterations)
        all_adjust_pval <- matrix(NA, nrow(stat_analysis), num_iterations)
        
        # Loop through iterations using lapply
        for (i in 2:sample_size) {
          p_values_tmp <- lapply(1:dim(stat_analysis)[1], function(j) {
            if (test_type == "t.test") {
              ttest_func(select(stat_analysis, contains(paste0("A",numerator,"-")) & contains("log10_"))[j,],
                         select(stat_analysis, contains(paste0("A", i, "-")) & contains("log10_"))[j,])
            } else if (test_type == "wilcoxon") {
              wilcox.test(select(stat_analysis, contains(paste0("A",numerator,"-")) & contains("log10_"))[j,],
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
        
        all_pvalues <- as.data.frame(all_pvalues)
        all_adjust_pval <- as.data.frame(all_adjust_pval)
        
        colnames(all_pvalues) <- paste0("pvalues_A1/A", 2:sample_size)
        colnames(all_adjust_pval) <- paste0("adjust_pval_A1/A", 2:sample_size)
        
        rownames(all_pvalues) <- row.names(stat_analysis)
        rownames(all_adjust_pval) <- row.names(stat_analysis)
    
    all_pvalues_common_col <- stat_analysis %>% 
      select(pep_with_pos, species) %>% #spectrum_title, accession
      bind_cols(all_pvalues) %>% 
      pivot_longer(cols = starts_with("pvalues_"), values_to = "pvalues", names_to ="p_ratios") %>%
      separate(p_ratios, into = c("tmp","ratio"),sep = "_") %>%
      select(!tmp) %>%
      mutate(common_col = paste(pep_with_pos,species,ratio,1:((sample_size-1)*nrow(stat_analysis)),sep="@")) #spectrum_title accession
    
    ## THE BEST WAY TO DO is this:
    merge_stat_df <- final_imputed_data %>%
      select(pep_with_pos, species, starts_with("exp_FC")) %>% #accession
      pivot_longer(cols = starts_with("exp_FC"), values_to = "fold_change_values", names_to ="fold_change_ratios") %>%
      separate(fold_change_ratios, into = c("tmp","tmp1","ratio"),sep = "_") %>%
      select(!c(tmp,tmp1)) %>%
      mutate(common_col = paste(pep_with_pos,species,1:((sample_size-1)*nrow(final_imputed_data)),sep="@")) %>% #accession
      bind_cols(all_pvalues_common_col$ratio,all_pvalues_common_col$pvalues) %>%
      rename_with(.col =6 , ~"ratio1") %>%
      rename_with(.col=7, ~ "pvalues") #%>%
      #separate(accession, into = c("prot_id","species"),sep = "_")
    
    
    
    actual_ratio_col <- merge_stat_df %>%
      select(ratio1) %>% distinct() %>%
      mutate(actual_ratio_val= case_when(grepl(comparisons[1],ratio1) ~ actual_ratio[1],
                                         grepl(comparisons[2],ratio1) ~actual_ratio[2],
                                         grepl(comparisons[3],ratio1) ~actual_ratio[3],
                                         grepl(comparisons[4],ratio1) ~actual_ratio[4]))
    col_vline <- rep("#000000",4)
    actual_ratio_col <- as.data.frame(actual_ratio_col %>% bind_cols(col_vline)) %>% rename(col=3)
    
    merge_stat_df_final <- merge_stat_df
    
    plot9_t_test <- ggplot(merge_stat_df_final ,aes(x =log2(fold_change_values), y = -log10(merge_stat_df$pvalues))) +
      geom_point(size = 9,aes(shape=species,color=species)) + #, aes(shape=merge_stat_df_final$species)
      facet_wrap(~ratio) +
      #scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 0.8)) +
      #scale_x_continuous(limits = c(-3,3),breaks = seq(-3, 3, by = 0.8)) +
      scale_color_brewer(palette = "Set1") +
      #scale_y_continuous(breaks = seq(0, max(-log10(volcano_final1$pvalues_value)), length.out = 21)) +
      theme_bw() +
      theme(legend.text = element_text(size = 45),
            axis.title.x = element_text(size = 45),
            axis.title.y = element_text(size = 45),
            plot.title = element_text(size = 55),
            legend.title = element_text(size = 45),
            axis.text.x = element_text(size = 45),
            axis.title = element_text(size = 45),
            axis.text.y = element_text(size = 45),
            plot.subtitle = element_text(size = 45)) +
      labs(title =  paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name), subtitle = "T-test was used")+
      scale_x_continuous(limits = c(-5, 10),breaks = seq(from = -5, to = 10, by = 2)) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#006D2C",size=2) +
      geom_vline(data = actual_ratio_col, aes(xintercept = log2(actual_ratio_val), show.legend = FALSE),size=2) 
    #expand_limits(x = 0, y = 0) +
    #geom_vline(data = actual_ratio, aes(xintercept = actual_ratio$X.1....log2.c.2..10..20..100..., size = 1, show.legend = FALSE)) + #color=c("#CC79A7","#E69F00","#56B4E9","#009E73")
    #geom_hline(data = log10_p_thresholds, aes(yintercept = log10_p_thresholds$X.log10.p_thresholds.),color=c("#CC79A7","#E69F00","#56B4E9","#009E73"), size = 1, linetype = 2, show.legend = FALSE)+ 
      
  }else if(test_type=="limma"){
    library(limma)

      #### ONE COMPARISION AT A TIME WITH SEPARATE DESIGN MATRIX ####
      
      design_matrix <- model.matrix(~factor(c(rep(2,num_reps),rep(1,num_reps))))
      merge_stat_df <-NULL
      for ( i in 2:sample_size){
        # Change only the colname iteratively makes fit to every comparison
        colnames(design_matrix) <- c("Intercept", paste0("A1-A",i))
        #print(colnames(design_matrix))
        # Col selection for each comparison
        assign(paste0("df_A1vsA",i),stat_analysis %>% select(1:2 | contains("log10_E3_A1") | contains(paste0("log10_E3_A",i))))
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
      
      ## MULTIPLE PIVOTING IN ONE DATAFRAME DOES NOT WORK - RATHER THAN THAT, BIND_ROWS() WAS USED ABOVE.
      
      merge_stat_df_final <- merge_stat_df %>%  #
        rename_with(.col =1 , ~"mult_col") %>%
        rename_with(.col=8, ~ "ratio") %>%
        separate(mult_col, into = c("pep_with_pos","species","id"),sep = "@") #%>% #"spectrum_title" accession
        #separate(accession, into = c("uniprot_id","species"),sep = "_")
      
      plot9_limma <- gg_volcano(data_set = merge_stat_df_final,
                   x_df = merge_stat_df_final$logFC,  
                   y_df = -log10(merge_stat_df_final$P.Value),
                   color_df = merge_stat_df_final$species,
                   facet_df = "ratios",
                   header = paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name), ## Specifications of the header will be asked as an input.
                   color_lab = "Species",
                   subtitle_txt = "Limma was used",
                   x_lab = "log2(fold_change_values)",
                   y_lab = "-log10(p_values)") +
        scale_x_continuous(limits = c(-5, 10),breaks = seq(from = -5, to = 10, by = 2)) +
        scale_color_brewer(palette = "Set1") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#006D2C",size=2) #+
        #geom_vline(data = actual_ratio_col, aes(xintercept = log2(actual_ratio_val), show.legend = FALSE),size=2) 
      
      #select(accession, spectrum_title, contains("P.value") | contains("logFC")) %>%
      #pivot_longer(cols = contains("P.value"), values_to = "pvalues", names_to ="p_ratios") %>%
      #pivot_longer(cols = contains("logFC"), values_to = "fold_change_values", names_to ="fold_change_ratios")
      
      # Filter the results based on the FDR threshold
      # Add column that defined values either above or lower the threshold as Boolean 
      
      #desired_fdr_threshold <- 0.005
      #merge_stat_df$sig <- merge_stat_df$adj.P.Val < desired_fdr_threshold
      #significant_results_A1_vs_A2 <- subset(results_A1_vs_A2, adj.P.Val <= desired_fdr_threshold)
      
  }else{
      print("Statistical test could not be assessed. Check the input files!")
    }
  # plot24 <- ggplot(merge_stat_df_final ,aes(x =merge_stat_df_final$logFC, y = -log10(merge_stat_df_final$P.Value), color=merge_stat_df_final$species)) +
  #   geom_point(size = 2,aes(shape=merge_stat_df_final$species)) + #, aes(shape=merge_stat_df_final$species)
  #   facet_wrap(~ratio) +
  #   #scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 0.8)) +
  #   #scale_x_continuous(limits = c(-3,3),breaks = seq(-3, 3, by = 0.8)) +
  #   scale_color_brewer(palette = "Set1") +
  #   #scale_y_continuous(breaks = seq(0, max(-log10(volcano_final1$pvalues_value)), length.out = 21)) +
  #   theme_bw() +
  #   theme(legend.text = element_text(size = 15),
  #         axis.title.x = element_text(size = 15),
  #         axis.title.y = element_text(size = 15),
  #         plot.title = element_text(size = 30),
  #         legend.title = element_text(size = 15),
  #         axis.text.x = element_text(size = 15),
  #         axis.title = element_text(size = 15),
  #         axis.text.y = element_text(size = 15)) +
  #   labs(title =  paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name), subtitle = "Limma was used")
  # 
  merge_stat_df_final_text <- merge_stat_df_final %>% 
    mutate(soft_name=paste0(software_name)) %>%
    mutate(acq_type=paste0(acquisiton_type))
  
  write.table(merge_stat_df_final_text,file = paste0(file_path,output_dir_name,"/volcano_plot_",software_name,"_",acquisiton_type,".txt"),sep = 
                "\t",col.names = T,row.names = F)
  
  df_roc <- merge_stat_df_final %>%
    select(pep_with_pos,species,contains("Value"))
  
  colnames(df_roc)[4] <- "P.Value"
  ### ROC analysis custom func
  df_roc_order <- df_roc[order(df_roc$P.Value),]
  size_variying_pep_size <- df_roc_order %>% count(species)
  
  df_roc_func <- compute_roc_curve_exp3(df=df_roc_order, flag = "MOUSE",expected =length(comparisons)*dim(filtered_abundances)[1])
  
  # Define a function to compute TPR and FDR for a given threshold
  # compute_metrics <- function(threshold) {
  #   filtered_data <- df_roc_order %>% filter(pvalue < threshold)
  #   tp <- sum(filtered_data$flag_present)
  #   fp <- sum(filtered_data$flag_absent)
  #   
  #   fdr <- (fp / (fp + tp)) * 100
  #   tpr <- (tp / expected) * 100
  #   
  #   return(c(threshold, tpr, fdr))
  # }
  # 
  # # Use sapply to apply the function over all thresholds
  # results <- sapply(thresholds_num, compute_metrics)
  # 
  
  write.table(df_roc_func, file = paste0(file_path,output_dir_name,"/new_custom_Roc_analysis_",exp_id,"_",software_name,"_",".txt"),sep = "\t",row.names = F)
  #write.table(df_roc_func_filt, file = paste0(file_path,output_dir_name,"/new_custom_Roc_analysis_",exp_id,"_",software_name,"_","filtered.txt"),sep = "\t",row.names = F)
  
  
  plot14 <- ggplot(df_roc_func, aes(y=tpr, x = fdr)) +
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
         subtitle = "", color="Pool Type")
  
  
  #### ROC Analysis using pROC 
  
  df_roc$variant <- ifelse(df_roc$species == selected_species, TRUE, FALSE)
  #df_roc$non_var <- ifelse(df_roc$Pool == "ISO-REF", TRUE, FALSE)
  
  library(pROC)
  # Calculate ROC curve for raw p-values
  roc_raw_variant <- roc(df_roc$variant, df_roc$P.Value)
  #tpr_and_fpr_variant  <- cbind(roc_raw_variant$sensitivities,
  #                              roc_raw_variant$specificities,
  #                              "Variant Pool")
  
  fpr <- as.data.frame(1 - roc_raw_variant$specificities)
  tpr_and_fpr_variant  <- cbind(roc_raw_variant$sensitivities,
                                fpr,#roc_raw_variant$specificities,
                                selected_species) 
   
  roc_plt_df <- as.data.frame(tpr_and_fpr_variant) %>%
    bind_cols(software_name)
    #bind_rows(as.data.frame(tpr_and_fpr_non_var)) %>% 

  
  colnames(roc_plt_df) <- c("sensitivity", "fpr","species","Software_name")
  
  
  plot10 <- roc_plt_df %>% 
    ggplot( aes(y=as.numeric(sensitivity), x = as.numeric(fpr), color=species)) +
    geom_path(size=1.5) + theme_bw() + #scale_x_reverse() 
    theme(legend.text = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          plot.title = element_text(size = 25),
          legend.title = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title = element_text(size = 20),
          axis.text.y = element_text(size = 20)) +
    scale_color_brewer(palette = "Dark2") +
    labs(y="True Positive Rate \n (Sensitivity)", x="False Discovery Rate",
         title =  paste("Experiment - ", exp_id, acquisiton_type, " data processed by ", software_name), 
         subtitle ="", color="Pool Type")
  
  write.table(roc_plt_df, file = paste0(file_path,output_dir_name,"/pRoc_analysis_",exp_id,"_",software_name,"_",".txt"),sep = "\t",row.names = F)
  
  # sapply(1:14,function(x) ggsave(filename = paste0("p",x,".tiff"),
  #                                width = 60, height = 45, 
  #                                path = paste0(file_path,"/outputs_with_new_script/"),
  #                                units = "cm",
  #                                get(paste0("p",x)),
  #                                device = "tiff", #".svg"
  # ))
  
  plt_obj <- ls(pattern="plot")
  plt_obj <- plt_obj[!is.na(plt_obj)]
  sapply(1:length(plt_obj),function(x) ggsave(filename = paste0("p",x,".png"),
                                              width = 90, height = 60, 
                                              path = paste0(file_path,"/refined_output/"),
                                              units = "cm",
                                              get(plt_obj[x]),
                                              device = "png", #".svg"
  ))
  
  
    
   
  }

  
  
  
  



