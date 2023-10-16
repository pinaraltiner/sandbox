#library(stringr)
library(dplyr)
#library(data.table)
library(openxlsx)
library(tidyr)
library(ggplot2)
#library(tidyverse)  

theo_file_path= "D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/"
theo_file_name="Synthetic peptides list_theo_conc_corrected_isomericity.xlsx"#"Synthetic peptides list_theo_conc_corrected_pool_id_iso_count_final.xlsx"
sheet_theo_name = "ISO-ref and OTHER with FC"
pep_list_w_theo <- read.xlsx(paste0(theo_file_path, theo_file_name), sheet = sheet_theo_name)
pep_list_w_theo_quant <- pep_list_w_theo[,-1]

common_col_theo_quant <- as.data.frame(paste(pep_list_w_theo_quant$Phosphopeptide.sequence,
                                             pep_list_w_theo_quant$modified.position.in.peptide, sep = "_"))

colnames(common_col_theo_quant) <- "pep_with_pos"
pep_list_w_theo_quant_new <- cbind(common_col_theo_quant,pep_list_w_theo_quant)

test <- pep_list_w_theo_quant_new %>% 
    count(Phosphopeptide.sequence) %>% left_join(pep_list_w_theo_quant_new,by="Phosphopeptide.sequence")

p<-  ggplot(test, aes(x=n, color=`#.of.sites.in.peptide`,fill=Pool)) +
  geom_bar(position = "dodge") + scale_fill_manual(values = c("#FF0000","#800080"),
                                                   labels = c("Non-variant = 35", "Variant = 141")) +
  #geom_segment(aes(x = 3, y = 70, xend = 2, yend = 40),size=1,
               #arrow = arrow(length = unit(0.5, "cm")),show.legend = FALSE, inherit.aes = FALSE) +
  theme_minimal() +
  theme(legend.text = element_text(size=30), 
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        plot.title = element_text(size=40),
        #plot.subtitle = element_text(size = 20),
        plot.caption = element_text(size = 25),
        legend.title=element_text(size=30),
        axis.text=element_text(size=30),
        axis.title=element_text(size=30),
        plot.subtitle = element_text(size = 30)) +
  geom_text(aes(label=after_stat(count)),
            stat = "count", position=position_dodge(width=0.9), vjust=-0.25,size=12) +
  #annotate(geom="text", 
           #label="\n TVSTSSQPEENVDR has two isomers (3 and 3&5) \n and are collected in the variant pool. ",
          # x=4.2, y=85,
          # color="black",size=11) + 
  labs(x="Number of phospho-sites", y="Phospho-site counts", 
       title = "Number of phospho-sites in the synthetic phospho-peptide library")


ggsave(filename = paste0("p2",".tiff"),
       width = 50, height = 45, 
       path = paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/combine_figures/"),
       units = "cm",
       get(paste0("p")),
       device = "tiff", #".svg"
)

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######                    
                 
                
file_paths <- c(paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/",
                       c("MQ_214_optimize_noFAIMS/","MQ_214_optimize_withFAIMS/","MQ_timstof/"),
                       c("outputs_with_new_script/")),
                
                paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/PD_data_analysis/", 
                       c("target_decoy_no_FAIMS/","percolator_no_FAIMS/",
                         "target_decoy_with_FAIMS/","percolator_with_FAIMS/"),
                      c("outputs_with_new_script/")),
                
                paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/",
                       c("exp2_re_injection/","with_FAIMS/"), 
                       c("outputs_with_new_script/")),
                "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_timsdata/outputs_with_new_script/",
                "D:/dev/Pinar/PHD/wet_lab_experiments/DIA_data_analysis/experiment_2/Spectronaut/correct_norm/outputs_with_new_script/")

software_names <- c(rep("MaxQuant v2.1.4",times=3),rep("Proteome Discoverer",times=4),rep("Proline",times=3),"Spectronaut")
acquisiton_types <- c("DDA Exploris no FAIMS","DDA Exploris with FAIMS", "DDA TIMS-TOF","target decoy no FAIMS","percolator no FAIMS",
                      "target decoy with FAIMS","percolator with FAIMS","DDA Exploris no FAIMS","DDA Exploris with FAIMS", "DDA TIMS-TOF","DIA Exploris no FAIMS")


acquisiton_types2<- c("DDA Exploris no FAIMS&target decoy",
                      "DDA Exploris with FAIMS&target decoy",
                      "DDA TIMS-TOF&target decoy",
                      "DDA Exploris no FAIMS&target decoy",
                      "DDA Exploris no FAIMS&percolator",
                      "DDA Exploris with FAIMS&target decoy",
                      "DDA Exploris with FAIMS&percolator",
                      "DDA Exploris no FAIMS&target decoy",
                      "DDA Exploris with FAIMS&target decoy",
                      "DDA TIMS-TOF&target decoy",
                      "DIA Exploris no FAIMS&directDIA")

###### Roc_analysis_2_Software_name_.txt and new_custom_Roc_analysis_2_Software_name_.txt######  pRoc_analysis_2_Proline_.txt
all_results <- NULL
for( i in 1:length(acquisiton_types)){
  assign(paste0(software_names[i],"_",acquisiton_types[i]),read.table(file = paste0(file_paths[i], paste0("new_custom_Roc_analysis_2","_",software_names[i],"_.txt")), sep = "\t",header = T))
  
  assign(paste0(software_names[i],"_",acquisiton_types[i]), 
         get(paste0(software_names[i],"_",acquisiton_types[i])) %>% #custom_Roc_analysis_2 
           mutate(new_col=paste(acquisiton_types2[i],software_names[i],sep = ".")) %>%
           mutate(soft_name=software_names[i])) 
  all_results <- bind_rows(all_results,get(paste0(software_names[i],"_",acquisiton_types[i])))
}

#custom_Roc_analysis_2_MaxQuant v2.1.4_

#var <- c("Variant")#,"Non-variant")
for (i in 1: length(software_names)){
  
  assign(paste0("roc_",software_names[i]), all_results %>% filter(grepl(acquisiton_types2[i],new_col)) %>% #filter(grepl(software_names[i],Software_name)) %>%
    ggplot( aes(y=as.numeric(tpr), x = as.numeric(fdp), color=soft_name)) +
    geom_path(size=1.5) + #scale_x_reverse() + 
  xlim(seq(0:100,by=10))+ylim(seq(0:100,by=10)) +
    #scale_alpha_manual(values = c(0.5,1)) +
    #facet_wrap(~new_col) + theme_bw() +
    theme_bw() +
    theme(legend.text = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text = element_text(size = 25),
          plot.title = element_text(size = 25),
          legend.title = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          plot.subtitle = element_text(size = 20)) +
    scale_color_brewer(palette = "Dark2") +
    labs(y="True Positive Rate \n (Sensitivity)", x="False Discovery Rate \n (FDR)",
         title =  paste("Experiment - ", 2, " ROC analysis in dataset 2"), 
         subtitle = paste("Including unexpected \n", acquisiton_types2[i]), color="Software Name"))


  
  ggsave(filename = paste0("custom_roc",acquisiton_types2[i],".tiff"),
         width = 58, height = 48,
         path = paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/combine_figures/"),
         units = "cm",
         get(paste0("roc_",software_names[i])),
         device = "tiff", #".svg"
  )

}

all_results_new_cat <- all_results %>% 
  filter(grepl("DDA Exploris no FAIMS&target decoy.MaxQuant v2.1.4",new_col) |
           grepl("DDA Exploris no FAIMS&target decoy.Proteome Discoverer",new_col) | 
           grepl("DDA Exploris no FAIMS&percolator.Proteome Discoverer",new_col) |
           grepl("DDA Exploris no FAIMS&target decoy.Proline",new_col) |
           grepl("DIA Exploris no FAIMS&directDIA.Spectronaut",new_col)) %>%
  separate(new_col, into = c("acq_type","second"),sep = "&") 

all_results_new_cat$new_cat <- factor(all_results_new_cat$second, levels = c("target decoy.MaxQuant v2.1.4",
                                                             "target decoy.Proline",
                                                             "target decoy.Proteome Discoverer",
                                                             "percolator.Proteome Discoverer",
                                                             "directDIA.Spectronaut"))

iter_val <- unique(all_results$new_col)
points_roc_curve <- NULL
for (i in 1:length(iter_val)){
  
  df_iter_val <- all_results %>% filter(grepl(iter_val[i],new_col)) %>% filter(threshold >=0.05)

  tmp <- df_iter_val$threshold[1]
  tpr <- df_iter_val$tpr[i]
  fdr <- df_iter_val$fdr[i]
  points_roc_curve$value[i] <- tmp 
  
  points_roc_curve$id[i] <- iter_val[i] 
  points_roc_curve$tpr[i] <- tpr
  points_roc_curve$fdr[i] <- fdr
}
points_roc_curve <- as.data.frame(points_roc_curve)
point_df_nofaims <- points_roc_curve %>% filter(grepl("no FAIMS",id))

roc_dda <- all_results_new_cat %>%
  ggplot(aes(y=all_results_new_cat$tpr, x = all_results_new_cat$fdr, color=all_results_new_cat$new_cat)) + #y=tpr, x=fdp #y=sensitivity, x=fpr
  geom_path(size=2.5) +
  #scale_x_reverse() + 
  #scale_x_continuous(breaks = seq(from= 0, to =(max(all_results$fpr)+20),by=5)) +
  #scale_y_continuous(breaks = seq(from =0,to=(max(all_results$sensitivity)+20),by=5)) +
  #xlim(c(0,100))+ylim(c(0,100)) +
  #scale_alpha_manual(values = c(0.5,1)) +
  #facet_wrap(~new_col) + theme_bw() +
  theme_bw() +
  theme(legend.text = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        strip.text = element_text(size = 30),
        plot.title = element_text(size = 35),
        legend.title = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.title = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        plot.subtitle = element_text(size = 30)) +
  #xlim(c(0,100)) + ylim(c(0,100)) +
  geom_point(data=point_df_nofaims, aes(x=point_df_nofaims$fdr,y=point_df_nofaims$tpr),color="black",size=5)+
  #scale_x_continuous(breaks = seq(from =0, to=100, by=5)) +
  scale_y_continuous(breaks = seq(from =0, to=115, by=5)) +
  scale_color_brewer(palette = "Set1") +
  labs(y="True Positive Rate \n (Sensitivity)", x="False Discovery Rate \n (FDR)", #"False Positive Rate\n (1-specificity)", #
       title =  paste("Experiment - ", 2, " ROC analysis in dataset 2"), 
       subtitle = paste("Including unexpected \n", "DDA Exploris no FAIMS"), color="Software Name")

ggsave(filename = paste0("TESTnew_custom_roc_dda_with_all_pvalue_threshold.tiff"),
       width = 58, height = 48,
       path = paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/combine_figures/after_eliminating_wrong_seq/"),
       units = "cm",
       get(paste0("roc_dda")),
       device = "tiff", #".svg"
)
rm(all_results_new_cat)

all_results_new_cat<- all_results %>% 
  filter(grepl("DDA Exploris with FAIMS&target decoy.MaxQuant v2.1.4",new_col) |
           grepl("DDA Exploris with FAIMS&target decoy.Proteome Discoverer",new_col) | 
           grepl("DDA Exploris with FAIMS&percolator.Proteome Discoverer",new_col) |
           grepl("DDA Exploris with FAIMS&target decoy.Proline",new_col) |
           grepl("DIA Exploris no FAIMS&directDIA.Spectronaut",new_col)) %>%
  separate(new_col, into = c("acq_type","second"),sep = "&")

all_results_new_cat$new_cat <- factor(all_results_new_cat$second, levels = c("target decoy.MaxQuant v2.1.4",
                                                                             "target decoy.Proline",
                                                                             "target decoy.Proteome Discoverer",
                                                                             "percolator.Proteome Discoverer",
                                                                             "directDIA.Spectronaut"))

point_df_withfaims <- points_roc_curve %>% filter(grepl("with FAIMS",id) | grepl("DIA",id))

roc_dda_faims <-  all_results_new_cat%>%
  ggplot( aes(y=tpr, x = fdr, color=new_cat)) + #y=tpr, x=fdp
  geom_path(size=2.5) + #scale_x_reverse() + 
  #scale_x_continuous(breaks = seq(from= 0, to =(max(all_results$tpr)+20),by=5)) +
  #scale_y_continuous(breaks = seq(from =0,to=(max(all_results$tpr)+20),by=5)) +
  #scale_alpha_manual(values = c(0.5,1)) +
  #facet_wrap(~new_col) + theme_bw() +
  #scale_x_continuous(breaks = seq(from =0, to=100, by=5)) +
  scale_y_continuous(breaks = seq(from =0, to=115, by=5)) +
  theme_bw() +
  theme(legend.text = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        strip.text = element_text(size = 30),
        plot.title = element_text(size = 35),
        legend.title = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.title = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        plot.subtitle = element_text(size = 30)) +
  geom_point(data=point_df_withfaims, aes(x=point_df_withfaims$fdr,y=point_df_withfaims$tpr),color="black",size=5)+
  scale_color_brewer(palette = "Set1") +
  labs(y="True Positive Rate \n (Sensitivity)", x="False Discovery Rate \n (FDR)", #"False Positive Rate\n (1-specificity)", #
       title =  paste("Experiment - ", 2, " ROC analysis in dataset 2"), 
       subtitle = paste("Including unexpected \n", "DDA Exploris with FAIMS"), color="Software Name")

ggsave(filename = paste0("TESTnew_custom_roc_dda_withfaims_with_all_pvalue_threshold.tiff"),
       width = 58, height = 48,
       path = paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/combine_figures/after_eliminating_wrong_seq/"),
       units = "cm",
       get(paste0("roc_dda_faims")),
       device = "tiff", #".svg"
)

rm(all_results_new_cat)
all_results_new_cat <-  all_results %>% 
  filter(grepl("DDA TIMS-TOF&target decoy.MaxQuant v2.1.4",new_col) |
           grepl("DDA TIMS-TOF&target decoy.Proline",new_col) |
           grepl("DIA Exploris no FAIMS&directDIA.Spectronaut",new_col)) %>%
  separate(new_col, into = c("acq_type","second"),sep = "&")

all_results_new_cat$new_cat <- factor(all_results_new_cat$second, levels = c("target decoy.MaxQuant v2.1.4",
                                                                             "target decoy.Proline",
                                                                             "target decoy.Proteome Discoverer",
                                                                             "percolator.Proteome Discoverer",
                                                                             "directDIA.Spectronaut"))
point_df_tims <- points_roc_curve %>% filter(grepl("TOF",id) | grepl("DIA",id))

dda_tims <- all_results_new_cat %>%
  ggplot( aes(y=tpr, x = fdr, color=new_cat)) + #y=tpr, x=fdp
  geom_path(size=2.5) + #scale_x_reverse() + 
  #scale_x_continuous(breaks = seq(from= 0, to =(max(all_results$tpr)+20),by=5)) +
  #scale_y_continuous(breaks = seq(from =0,to=(max(all_results$tpr)+20),by=5)) +
  #scale_alpha_manual(values = c(0.5,1)) +
  #facet_wrap(~new_col) + theme_bw() +
  theme_bw() +
  theme(legend.text = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        strip.text = element_text(size = 30),
        plot.title = element_text(size = 35),
        legend.title = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.title = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        plot.subtitle = element_text(size = 30)) +
  geom_point(data=point_df_tims, aes(x=point_df_tims$fdr,y=point_df_tims$tpr),color="black",size=5)+
  scale_color_manual(values = c("#E41A1C", "#377EB8","#FF7F00")) +
  #scale_x_continuous(breaks = seq(from =0, to=100, by=5)) +
  scale_y_continuous(breaks = seq(from =0, to=115, by=5)) +
  #scale_color_brewer(palette = "Set1") +
  labs(y="True Positive Rate \n (Sensitivity)", x="False Discovery Rate \n (FDR)", #"False Positive Rate\n (1-specificity)", #
       title =  paste("Experiment - ", 2, " ROC analysis in dataset 2"), 
       subtitle = paste("Including unexpected \n", "DDA TIMSTOF"), color="Software Name")

ggsave(filename = paste0("TESTnew_custom_roc_dda_TIMS_with_all_pvalue_threshold.tiff"),
       width = 58, height = 48,
       path = paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/combine_figures/after_eliminating_wrong_seq/"),
       units = "cm",
       get(paste0("dda_tims")),
       device = "tiff", #".svg"
)

######## custom_Roc_analysis_2_Software_name_.txt with/out unexpected######  
source("D:/dev/Pinar/PHD/sandbox/benchmarking_scripts/scripts_from_data_analysis/roc_curve/roc_curve_generation_proline_edit.R")

#volcano_plot_Proteome Discoverer_DDA Exploris no FAIMS.txt
vp.Spectronaut.noFaims.DIA <- vp.Spectronaut.noFaims.DIA %>% 
  mutate(soft_name="Spectronaut") %>% 
  mutate(acq_type="DIA Exploris no FAIMS")

#iter <- ls()
iter <- c("vp_Proline_DDA.Exploris.no.FAIMS",
          "vp_Proline_DDA.Exploris.with.FAIMS",
          "vp_Proline_DDA.TIMS.TOF",
          "vp_Proteome.Discoverer.TD.noFAIMS",
          "vp_Proteome.Discoverer.P.noFAIMS",
          "vp_Proteome.Discoverer.TD.Faims",
          "vp_Proteome.Discoverer.P.Faims",
          "vp_MaxQuant.v2.1.4.no.FAIMS",
          "vp_MaxQuant.v2.1.4.with.FAIMS",
          "vp_MaxQuant.v2.1.4.TIMS.TOF",
          "vp.Spectronaut.noFaims.DIA")

acquisiton_types2<- c("DDA Exploris no FAIMS&target decoy",
                      "DDA Exploris with FAIMS&target decoy",
                      "DDA TIMS-TOF&target decoy",
                      "DDA Exploris no FAIMS&target decoy",
                      "DDA Exploris no FAIMS&percolator",
                      "DDA Exploris with FAIMS&target decoy",
                      "DDA Exploris with FAIMS&percolator",
                      "DDA Exploris no FAIMS&target decoy",
                      "DDA Exploris with FAIMS&target decoy",
                      "DDA TIMS-TOF&target decoy",
                      "DIA Exploris no FAIMS")

software_names <- c(rep("MaxQuant v2.1.4",times=3),rep("Proteome Discoverer",times=4),rep("Proline",times=3),"Spectronaut")

all_results <- NULL
for( i in 1:length(iter)){
  
 tmp <- get(paste0(iter[i])) %>% 
   filter(!grepl("Unexpected",Pool)) %>% 
   select(Pool,P.Value,soft_name,acq_type) %>%
   arrange(P.Value)
 
  tmp2 <- compute_roc_curve(df=tmp,flag = "Others",expected = (141*4))
  tmp3 <- tmp2 %>% mutate(soft_name=tmp$soft_name[1]) %>% mutate(acq_type=acquisiton_types2[i])
  all_results <- bind_rows(all_results,tmp3)
  rm(tmp,tmp2,tmp3)
}

for (i in 1:length(software_names)){
  
  assign(paste0("new_roc_",software_names[i]), all_results %>% filter(grepl(software_names[i],soft_name)) %>% #filter(grepl(software_names[i],Software_name)) %>%
         ggplot( aes(y=as.numeric(tpr), x = as.numeric(fdp), color=acq_type)) +
         geom_path(size=1.5) + #scale_x_reverse() + 
         xlim(c(0,100))+ylim(c(0,100)) +
         #scale_alpha_manual(values = c(0.5,1)) +
         #facet_wrap(~new_col) + theme_bw() +
         theme(legend.text = element_text(size = 20),
               axis.title.x = element_text(size = 20),
               axis.title.y = element_text(size = 20),
               strip.text = element_text(size = 25),
               plot.title = element_text(size = 25),
               legend.title = element_text(size = 20),
               axis.text.x = element_text(size = 20),
               axis.title = element_text(size = 20),
               axis.text.y = element_text(size = 20),
               plot.subtitle = element_text(size = 20)) +
         scale_color_brewer(palette = "Dark2") +
         labs(y="True Positive Rate \n (Sensitivity)", x="False Discovery Rate \n (FDR)",
              title =  paste("Experiment - ", 2, software_names[i]," ROC analysis in dataset 2"), 
              subtitle = paste("Unexpected are excluded \n", software_names[i]), color="Software Name"))
  ggsave(filename = paste0("new_roc_",software_names[i],".tiff"),
         width = 58, height = 48,
         path = paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/combine_figures/"),
         units = "cm",
         get(paste0("new_roc_",software_names[i])),
         device = "tiff", #".svg"
  )

}

###### Count_of_missing_unexpected_correct_phospho-sites_with_all_col_Proteome Discoverer_Experiment2.txt ###### 
all_results <- NULL
for( i in 1:length(acquisiton_types)){
    assign(paste0(software_names[i],"_",acquisiton_types[i]),read.table(file = paste0(file_paths[i], paste0("Count_of_missing_unexpected_correct_phospho-sites_",software_names[i],"_Experiment2.txt")), sep = "\t",header = T))
    
    assign(paste0(software_names[i],"_",acquisiton_types[i]), get(paste0(software_names[i],"_",acquisiton_types[i])) %>% mutate(new_col=acquisiton_types2[i]))
    all_results <- bind_rows(all_results,get(paste0(software_names[i],"_",acquisiton_types[i])))
}

all_results_rmv_missing <- all_results %>% filter(!grepl("missing",Pool))
Pool_type <- c("ISO-REF","Others","Unexpected") #missing
Pool_type_name <- c("Non-variant Pool","Variant Pool", "Unexpected") #"missing"
for (i in 1:length(Pool_type)){
  
  assign(paste0("plot_count_",Pool_type[i]), all_results %>% select(!starts_with('E2')) %>% 
    mutate(combine =paste0(ion_mobility," ",new_col)) %>%
    group_by(soft_name, combine) %>% 
    count(Pool)  %>%
    filter(grepl(Pool_type[i],Pool)) %>%
    ggplot(aes(x=soft_name,y=n,fill=combine)) + 
    geom_col(position = "dodge") +
    scale_fill_brewer(palette = 2,direction=-1) +
    theme_minimal() +
    theme(legend.text = element_text(size=25), 
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25),
          plot.title = element_text(size=30),
          plot.subtitle = element_text(size = 20),
          plot.caption = element_text(size = 20),
          legend.title=element_text(size=25),
          axis.text=element_text(size=25),
          axis.title=element_text(size=25)) +
    ggtitle("Number of Phospho-site found in the dataset 2") + 
    scale_y_continuous(breaks = seq(from=0, to=200,by=500)) +
    geom_text(aes(label=n),stat = "identity", position=position_dodge(width=0.9), vjust=-0.25,size=6) +
    labs(x="Software name",y="Number of Phospho-isomer", fill= "Acquisiton and Validation Method", subtitle = Pool_type_name[i]))
  
  ggsave(filename = paste0("plot_count_",Pool_type[i],".tiff"),
         width = 50, height = 45, 
         path = paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/"),
         units = "cm",
         get(paste0("plot_count_",Pool_type[i])),
         device = "tiff", #".svg"
  )
  
}
display.brewer.pal(n = 9, name = 'Purples')
brewer.pal(n, 'Purples')
brewer.pal(n = 5, name='Purples')


new_acquisition_type_for_theo <- c("DDA Exploris no FAIMS&percolator",
                                   "DDA Exploris no FAIMS&target decoy",#""DDA TIMS-TOF&target decoy",
                                   "DIA Exploris no FAIMS&directDIA") 



theo_pep_all <- NULL
for ( i in 1:length(unique(acquisiton_types2))){
  assign(paste0("theo_pep",i ),pep_list_w_theo_quant_new %>%
           select(pep_with_pos, Pool) %>% mutate(soft_name="Theoretical Count") %>% 
           mutate(new_col=unique(new_acquisition_type_for_theo)[i]))
  theo_pep_all <- bind_rows(theo_pep_all,get(paste0("theo_pep",i)))    
}

all_results_com <- all_results_rmv_missing %>% select(!starts_with('E2')) %>% bind_rows(theo_pep_all)

all_results_com$new_cat <- factor(all_results_com$Pool, levels = c("Unexpected", "Others","ISO-REF" ))

all_results_com$new_facet <- factor(all_results_com$new_col, levels = c("DDA Exploris no FAIMS&percolator",
                                                                "DDA Exploris no FAIMS&target decoy",
                                                                "DIA Exploris no FAIMS&directDIA",
                                                                "DDA Exploris with FAIMS&percolator",
                                                                "DDA Exploris with FAIMS&target decoy",
                                                                "DDA TIMS-TOF&target decoy"
                                                                ))
all_results_com$soft_cat <- factor(all_results_com$soft_name, levels = c("MaxQuant v2.1.4","Proline","Proteome Discoverer","Spectronaut","Theoretical Count"))


phospho_site_count_dataset2 <- all_results_com %>% 
  filter(grepl( "DDA Exploris no FAIMS&target decoy",new_facet) & grepl("Proteome Discoverer",soft_name) | grepl("DIA",new_facet)) %>%
      ### DDA all soft
  #filter(grepl("target decoy",new_facet) | grepl('TIMS-TOF',new_facet)) %>%
         ## PERCO vs TD
  #filter(grepl("Proteome Discoverer",soft_name) & grepl("DDA Exploris no FAIMS",new_facet) | grepl("Theoretical Count", soft_name) & !grepl("DIA", new_facet)) %>%
  
  ggplot(aes(x=soft_cat,fill=new_cat)) + #soft_name
  geom_bar(position = "stack",width = 0.5) +
  facet_wrap(~new_facet,labeller = labeller(new_facet = label_wrap_gen(width = 60)))  +
  ylim(c(0,240)) +
  #scale_fill_brewer(palette = "Set2")+
  scale_fill_manual(values = c("#B3B3B3","#6A51A3","#FB6A4A"),labels=c("Wrong Localization","Fixed Pool","Spiked Pool")) + 
  theme_bw() +
  theme(legend.text = element_text(size=25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        plot.title = element_text(size=30),
        plot.subtitle = element_text(size = 20),
        strip.text = element_text(size=20),
        plot.caption = element_text(size = 20),
        legend.title=element_text(size=25),
        axis.text.x=element_text(angle = 90,size=25),
        axis.text.y=element_text(size=25),
        axis.title=element_text(size=25)) +
  ggtitle("Number of Phospho-site found in the dataset 2") + 
  #scale_y_continuous(breaks = seq(from=0, to=200,by=500)) +
  geom_text(color="black",aes(label=after_stat(count)),show.legend = F,stat = "count",size=13,position = position_stack(vjust = 0.5))+ #position=position_dodge(width=0.9), vjust=-0.25,size=6) +
  labs(x="Software Name",y="Number of Phospho-isomer", fill= "Pool Types", subtitle = "") #+
  #scale_x_discrete(label=c("Non-variant Pool","Variant Pool","Unexpected","Missing"))

  
  ggsave(filename = paste0("comparison_perco_td_nofaims_pd.tiff"),
         width =50, height =45, # 70, height = 60, 
         path = paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/combine_figures/after_eliminating_wrong_seq/"),
         units = "cm",
         phospho_site_count_dataset2,
         device = "tiff", #".svg"
  )
  
  
  ####### ##### ###### ###### ###### ###### ###### ###### ###### 

  
for( i in 1:length(acquisiton_types)){
    
assign(paste0("p",i,"_",software_names[i]), 
       all_results %>% filter(grepl("Variant",Pool_type)) %>% filter(grepl(software_names[i],Software_name)) %>%
           ggplot( aes(y=as.numeric(sensitivity), x = as.numeric(specificity), color=new_col)) +
           geom_path(size=1.5) +scale_x_reverse() + theme_bw() +
         theme(legend.text = element_text(size = 30),
               axis.title.x = element_text(size = 30),
               axis.title.y = element_text(size = 30),
               plot.title = element_text(size = 35),
               legend.title = element_text(size = 30),
               axis.text.x = element_text(angle = 90,size = 30),
               axis.title = element_text(size = 30),
               axis.text.y = element_text(size = 30),
               strip.text = element_text(size=25)) +
           scale_color_brewer(palette = "Set1") +
           labs(y="True Positive Rate \n (Sensitivity)", x="False Positive Rate \n (Specificity)",
                title =  paste("Experiment - ", 2, " data processed by ", software_names[i]), 
                subtitle = paste(subtitle), color="Data Acquisition Methods"))
    
    ggsave(filename = paste0("p",i,"_",software_names[i],".tiff"),
           width = 50, height = 45, 
           path = paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/"),
           units = "cm",
           get(paste0("p",i,"_",software_names[i])),
           device = "tiff", #".svg"
    )
}
    

all_results1 <- NULL
for( i in 1:length(acquisiton_types2)){
    assign(paste0(software_names[i],"_",acquisiton_types[i]), get(paste0(software_names[i],"_",acquisiton_types[i])) %>% mutate(new_col=acquisiton_types2[i]))
    all_results1 <- bind_rows(all_results1,get(paste0(software_names[i],"_",acquisiton_types2[i])))
}

for (i in 1:3){
    
    all_results2 <- all_results %>% filter(grepl("Variant",Pool_type)) %>% 
        #filter(grepl(acquisiton_types[i],new_col)) %>%
        separate(new_col,into = c("instrument_mod","validation_mod"),sep = "&") %>%
        mutate(soft_val_type=paste(Software_name,validation_mod,sep = "_")) %>%
        filter(grepl(unique(instrument_mod)[i],instrument_mod))
    
    instrument_mod_tit <- unique(all_results2$instrument_mod)
    
    assign(paste0("p",i), ggplot(all_results2, aes(y=as.numeric(sensitivity), x = as.numeric(specificity), color=soft_val_type ,linetype=Software_name)) +
               geom_path(size=1.5) +scale_x_reverse() + theme_bw() +
               #scale_linetype_manual(values = c("solid" = "solid", "dashed" = "dashed")) +
             theme(legend.text = element_text(size = 30),
                   axis.title.x = element_text(size = 30),
                   axis.title.y = element_text(size = 30),
                   plot.title = element_text(size = 35),
                   legend.title = element_text(size = 30),
                   axis.text.x = element_text(angle = 90,size = 30),
                   axis.title = element_text(size = 30),
                   axis.text.y = element_text(size = 30),
                   strip.text = element_text(size=25)) +
               scale_color_brewer(palette = "Set1") +
               labs(y="True Positive Rate \n (Sensitivity)", x="False Positive Rate \n (Specificity)",
                    title =  paste("Experiment - ", 2, " data obtained from ", instrument_mod_tit), 
                    subtitle = paste(subtitle), color="Software Names"))
    
    
    
    ggsave(filename = paste0("p",i,"_",instrument_mod_tit,".tiff"),
           width = 50, height = 68, 
           path = paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/roc_analysis"),
           units = "cm",
           get(paste0("p",i)),
           device = "tiff", #".svg"
    )
    rm(all_results2,instrument_mod_tit)
    
}

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 

#### Experiment2 SOFTWARE_NAME _number_of_unique_sequence_for_each_species.txt
all_results <- NULL
for( i in 1:length(software_names)){
  assign(paste0(software_names[i],"_",acquisiton_types[i]),read.table(file = paste0(file_paths[i],paste0("Experiment2",software_names[i],"_number_of_unique_sequence_for_each_species.txt")), sep = "\t",header = T))
  
  #assign(paste0(software_names[i],"_",acquisiton_types[i]), get(paste0(software_names[i],"_",acquisiton_types[i])) %>% mutate(new_col=acquisiton_types[i]))
  all_results <- bind_rows(all_results,get(paste0(software_names[i],"_",acquisiton_types[i])))
}



all_results1 <- NULL
for( i in 1:length(software_names)){
  
  assign(paste0(software_names[i],"_",acquisiton_types[i]), get(paste0(software_names[i],"_",acquisiton_types[i])) %>% mutate(new_col=acquisiton_types2[i]))
  
  tmp <- get(paste0(software_names[i],"_",acquisiton_types[i])) %>% select(1,Pool_for_seq_merge,acq_type,soft_name,new_col) %>% rename(phospho_seq=1)
  
  all_results1 <- bind_rows(all_results1,tmp)
  rm(tmp)
}

# colnames(`Proteome Discoverer_percolator no FAIMS`)[5] <- "Phosphopeptide.sequence"
# colnames(`Proteome Discoverer_target decoy no FAIMS`)[5] <- "Phosphopeptide.sequence"
# colnames(`Proteome Discoverer_percolator with FAIMS`)[5] <- "Phosphopeptide.sequence"
# colnames(`Proteome Discoverer_percolator with FAIMS`)[5] <- "Phosphopeptide.sequence"
# colnames(`Proteome Discoverer_target decoy with FAIMS`)[5] <- "Phosphopeptide.sequence"

plot_data <- all_results1 %>% 
  filter(!grepl("ECOLI",Pool_for_seq_merge) & !grepl("Escherichia coli",Pool_for_seq_merge)) %>% 
  separate(new_col,into = c("acq_type","val_type"),sep = "&") %>%
  mutate(new_col= paste(acq_type,val_type,sep = " "))
                   #filter(grepl(acquisiton_types2[i],new_col)) %>%
plot_data$new_cat <- factor(plot_data$Pool_for_seq_merge, levels = c("ISO-REF", "Others", "missing"))
plot_data$new_facet <- factor(plot_data$new_col, levels = c("DDA Exploris no FAIMS percolator",
                                                                        "DDA Exploris no FAIMS target decoy",
                                                                        "DIA Exploris no FAIMS directDIA",
                                                                        "DDA Exploris with FAIMS percolator",
                                                                        "DDA Exploris with FAIMS target decoy",
                                                                        "DDA TIMS-TOF target decoy"
))

phospho_count_dataset1 <- ggplot(plot_data, aes(x =new_cat , fill=soft_name)) +
         geom_bar(position = "stack") +
         #geom_text(aes(label=after_stat(count)),stat = "count", position=position_dodge(width=0.9), vjust=-0.25,size=10) +
         facet_wrap(~new_facet,
                    labeller = labeller(new_facet = label_wrap_gen(width = 25))) + 
         #scale_alpha_manual(values=c(0.25, 0.50,0.75, 1)) +
         theme_bw() +
         #scale_linetype_manual(values = c("solid" = "solid", "dashed" = "dashed")) +
         theme(legend.text = element_text(size = 30),
               axis.title.x = element_text(size = 30),
               axis.title.y = element_text(size = 30),
               plot.title = element_text(size = 35),
               legend.title = element_text(size = 30),
               axis.text.x = element_text(angle = 90,size = 30),
               axis.title = element_text(size = 30),
               axis.text.y = element_text(size = 30),
               strip.text = element_text(size=25)) +
          scale_fill_brewer(palette = "Set2")+
         #scale_fill_manual(values = c("#E41A1C","#984EA3","#F781BF")) +
         labs(y="Count of Phospho-sequence", x="Pool name",
              title =  paste("Number of identified phospho sequences across all dataset 2"), 
              subtitle = "", fill="Pool Type",alpha="Software Name") +
  scale_x_discrete(label=c("Non-variant Pool","Variant Pool","Unexpected","Missing")) +
  geom_text(color="black",aes(label=after_stat(count)),show.legend = F,stat = "count",size=14,position = position_stack(vjust = 0.5))+ #position=position_dodge(width=0.9), vjust=-0.25,size=6) +
  labs(x="Pool Types",y="Number of Phospho-isomer", fill= "Pool Name", subtitle = "") #+



ggsave(filename = paste0("phospho_seq_count_dataset2.tiff"),
       width = 58, height = 68,
       path = paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/combine_figures/after_eliminating_wrong_seq/"),
       units = "cm",
       phospho_count_dataset1,
       device = "tiff", #".svg"
)




