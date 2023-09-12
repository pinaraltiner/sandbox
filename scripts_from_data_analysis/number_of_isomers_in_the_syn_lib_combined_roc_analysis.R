

test <- pep_list_w_theo_quant_new %>% 
    count(Phosphopeptide.sequence) %>% left_join(pep_list_w_theo_quant_new,by="Phosphopeptide.sequence")

                   p<-  ggplot(test, aes(x=n, color=`#.of.sites.in.peptide`,fill=Pool)) +
                    geom_bar(position = "dodge") + scale_fill_brewer(palette = "Dark2",labels = c("Non-variant = 35", "Variant = 144")) +
                    geom_segment(aes(x = 3, y = 80, xend = 2, yend = 40),size=1,
                                     arrow = arrow(length = unit(0.5, "cm")),show.legend = FALSE, inherit.aes = FALSE) +
                    theme_minimal() +
                    theme(legend.text = element_text(size=25), 
                          axis.title.x = element_text(size = 25),
                          axis.title.y = element_text(size = 25),
                          plot.title = element_text(size=35),
                          plot.subtitle = element_text(size = 20),
                          plot.caption = element_text(size = 20),
                          legend.title=element_text(size=25),
                          axis.text=element_text(size=25),
                          axis.title=element_text(size=25)) +
                    geom_text(aes(label=after_stat(count)),
                    stat = "count", position=position_dodge(width=0.9), vjust=-0.25,size=12) +
                    annotate(geom="text", 
                             label="\n TVSTSSQPEENVDR has two isomers (3 and 3&5) \n and are collected in the variant pool. ",
                             x=4.2, y=85,
                             color="black",size=11) + 
                       labs(x="Number of phospho-sites", y="Phospho-site counts", 
                            title = "Number of phospho-sites in the synthetic phospho-peptide library")
                
                   
                   ggsave(filename = paste0("p",".tiff"),
                          width = 50, height = 45, 
                          path = paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/"),
                          units = "cm",
                          get(paste0("p")),
                          device = "tiff", #".svg"
                   )

                   
                   
                   
file_paths <- c(paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/",c("MQ_214_optimize_noFAIMS/","MQ_214_optimize_withFAIMS/","MQ_timstof/")),
                paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/PD_data_analysis/",
                       c("target_decoy_no_FAIMS/","percolator_no_FAIMS/",
                         "target_decoy_with_FAIMS/","percolator_with_FAIMS/")),
                paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/", c("exp2_re_injection/","with_FAIMS/")),
                "D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_timsdata/")

software_names <- c(rep("MaxQuant v2.1.4",times=3),rep("Proteome Discoverer",times=4),rep("Proline",times=3))
acquisiton_types <- c("DDA Exploris no FAIMS","DDA Exploris with FAIMS", "DDA TIMS-TOF","target decoy no FAIMS/","percolator no FAIMS/",
                      "target decoy with FAIMS/","percolator with FAIMS/","DDA Exploris no FAIMS","DDA Exploris with FAIMS", "DDA TIMS-TOF")

all_results <- NULL
for( i in 1:length(acquisiton_types)){
    assign(paste0(software_names[i],"_",acquisiton_types[i]),read.table(file = paste0(file_paths[i], paste0("Roc_analysis_2","_",software_names[i],"_.txt")), sep = "\t",header = T))
    
    assign(paste0(software_names[i],"_",acquisiton_types[i]), get(paste0(software_names[i],"_",acquisiton_types[i])) %>% mutate(new_col=acquisiton_types[i]))
    all_results <- bind_rows(all_results,get(paste0(software_names[i],"_",acquisiton_types[i])))
}


for( i in 1:length(acquisiton_types)){
    
assign(paste0("p",i,"_",software_names[i]), 
       all_results %>% filter(grepl("Variant",Pool_type)) %>% filter(grepl(software_names[i],Software_name)) %>%
           ggplot( aes(y=as.numeric(sensitivity), x = as.numeric(specificity), color=new_col)) +
           geom_path(size=1.5) +scale_x_reverse() + theme_bw() +
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
    

acquisiton_types2<- c("DDA Exploris no FAIMS@target decoy",
                      "DDA Exploris with FAIMS@target decoy",
                      "DDA TIMS-TOF@target decoy",
                      "DDA Exploris no FAIMS@target decoy",
                      "DDA Exploris no FAIMS@percolator",
                      "DDA Exploris with FAIMS@target decoy",
                      "DDA Exploris with FAIMS@percolator",
                      "DDA Exploris no FAIMS@target decoy",
                      "DDA Exploris with FAIMS@target decoy",
                      "DDA TIMS-TOF@target decoy")

all_results1 <- NULL
for( i in 1:length(acquisiton_types2)){
    assign(paste0(software_names[i],"_",acquisiton_types[i]), get(paste0(software_names[i],"_",acquisiton_types[i])) %>% mutate(new_col=acquisiton_types2[i]))
    all_results1 <- bind_rows(all_results1,get(paste0(software_names[i],"_",acquisiton_types2[i])))
}

for (i in 1:3){
    
    all_results2 <- all_results %>% filter(grepl("Variant",Pool_type)) %>% 
        #filter(grepl(acquisiton_types[i],new_col)) %>%
        separate(new_col,into = c("instrument_mod","validation_mod"),sep = "@") %>%
        mutate(soft_val_type=paste(Software_name,validation_mod,sep = "_")) %>%
        filter(grepl(unique(instrument_mod)[i],instrument_mod))
    
    instrument_mod_tit <- unique(all_results2$instrument_mod)
    
    assign(paste0("p",i), ggplot(all_results2, aes(y=as.numeric(sensitivity), x = as.numeric(specificity), color=soft_val_type ,linetype=Software_name)) +
               geom_path(size=1.5) +scale_x_reverse() + theme_bw() +
               #scale_linetype_manual(values = c("solid" = "solid", "dashed" = "dashed")) +
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
                    title =  paste("Experiment - ", 2, " data obtained from ", instrument_mod_tit), 
                    subtitle = paste(subtitle), color="Software Names"))
    
    
    
    ggsave(filename = paste0("p",i,"_",instrument_mod_tit,".tiff"),
           width = 50, height = 45, 
           path = paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/roc_analysis"),
           units = "cm",
           get(paste0("p",i)),
           device = "tiff", #".svg"
    )
    rm(all_results2,instrument_mod_tit)
    
}
