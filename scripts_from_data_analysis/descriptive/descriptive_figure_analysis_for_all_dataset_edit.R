#
library(readr)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(viridis)
library(gridExtra)
setwd("D:/dev/Desktop_copy/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/")
setwd("D:/dev/Desktop_copy/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/list_of_all_collected_rawfiles/")
list_table <- read.delim(file="shortcuts_sequence_of_sample_injection_13042023.txt",sep = "\t",header = T)

list_table <- readxl::read_xlsx("list_of_all_collected_rawfiles/shortcuts_sequence_of_sample_injection_with_filenames_1_26052023.xlsx",
                                sheet = "shortcuts_sequence_of_sample_in")

sample_injection <- list_table %>% 
  tibble() %>%
  pivot_longer(!Shortcuts_Sample_Names,names_to = "parameters", values_to = "sample_injection") 

sample_inj <- ggplot(sample_injection,
                       aes(x=factor(parameters,
                                    levels=c("DDA.with.FAIMS",
                                             "DDA.without.FAIMS",
                                             "DIA.with.FAIMS",
                                             "DIA.without.FAIMS",
                                             "DIA.with.TIMS_ToF",
                                             "DDA.with.TIMS_Tof")),
                           y=Shortcuts_Sample_Names,
                           fill=as.factor(sample_injection))) + 
  geom_tile(color="gray") +
  theme_bw() + 
  theme(legend.text = element_text(size=30), 
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        plot.title = element_text(size=30),
        legend.title=element_text(size=30),
        axis.text.x=element_text(size=30),
        axis.text.y=element_blank(),
        axis.title=element_text(size=30),
        axis.line = element_blank(),
        strip.text.x = element_text(size = 30)
  ) + labs(y="") + labs(y="",x="parameters") +
  #scale_y_continuous(breaks = seq(from=0,to=40,by=5)) +
  #scale_fill_viridis(discrete = TRUE,option = "rocket",direction = -1)
  scale_fill_manual(values = c("#fcfdbf","#e75263"))


#list_table_data_analysis <- read.delim(file="sequence_of_data_analysis_representation_charc_13042023.txt",sep = "\t",header = T)


# list_table_data_analysis <- readxl::read_xlsx("list_of_all_collected_rawfiles/shortcuts_sequence_of_sample_injection_with_filenames_1_26052023.xlsx",
#                                               sheet = "data_analysis")

list_table_data_analysis <- read.delim(file="sequence_of_data_analysis_representation_charc_29062023.txt",sep = "\t",header = T)

analyzed_data <- list_table_data_analysis %>%
  tibble() %>%
  pivot_longer(cols = Proline:Spectronaut, names_to = "Software_name",values_to = "count") %>%
  mutate(exp_type= case_when(grepl("Coli",Shortcuts_Sample_Names,fixed = T) ~ "exp1#2",
                             grepl("E1",Shortcuts_Sample_Names,fixed = T) ~ "exp1#1",
                             grepl("E2",Shortcuts_Sample_Names,fixed = T) ~ "exp2",
                             TRUE ~ "exp3")) %>% 
  drop_na() %>% 
  mutate(count = factor(count)) 

group_col <- analyzed_data %>% 
  ggplot( aes(x="",y=factor(Shortcuts_Sample_Names,levels = c("E1-M1","E1-M2","E1-M3","E1-M4",
                                                              "E1-M5",
                                                              "E1-M6",
                                                              "E1-M7",
                                                              "E1-M8",
                                                              "E1-M1-Coli",
                                                              "E1-M2-Coli",
                                                              "E1-M3-Coli",
                                                              "E1-M4-Coli",
                                                              "E1-M5-Coli",
                                                              "E1-M6-Coli",
                                                              "E1-M7-Coli",
                                                              "E1-M8-Coli",
                                                              "E2-A1-R1","E2-A1-R2","E2-A1-R3","E2-A2-R1","E2-A2-R2","E2-A2-R3","E2-A3-R1","E2-A3-R2","E2-A3-R3","E2-A4-R1",
                                                              "E2-A4-R2","E2-A4-R3","E2-A5-R1","E2-A5-R2","E2-A5-R3","E3-A1-R1","E3-A1-R2","E3-A1-R3","E3-A2-R1","E3-A2-R2","E3-A2-R3","E3-A3-R1","E3-A3-R2",
                                                              "E3-A3-R3","E3-A4-R1","E3-A4-R2","E3-A4-R3","E3-A5-R1","E3-A5-R2","E3-A5-R3"),exclude = NA),
              fill=exp_type)) + geom_tile() + labs(y='Sample Shorcuts',x="d",fill="") + theme_minimal() + theme(legend.position = "none",axis.text.y=element_text(size=15))




exp_summary <- ggplot(analyzed_data,aes(x=Software_name,
             y=factor(Shortcuts_Sample_Names,levels = c("E1-M1",
                                                        "E1-M2",
                                                        "E1-M3",
                                                        "E1-M4",
                                                        "E1-M5",
                                                        "E1-M6",
                                                        "E1-M7",
                                                        "E1-M8",
                                                        "E1-M1-Coli",
                                                        "E1-M2-Coli",
                                                        "E1-M3-Coli",
                                                        "E1-M4-Coli",
                                                        "E1-M5-Coli",
                                                        "E1-M6-Coli",
                                                        "E1-M7-Coli",
                                                        "E1-M8-Coli",
                                                        "E2-A1-R1","E2-A1-R2","E2-A1-R3","E2-A2-R1","E2-A2-R2","E2-A2-R3","E2-A3-R1","E2-A3-R2","E2-A3-R3","E2-A4-R1",
                                                        "E2-A4-R2","E2-A4-R3","E2-A5-R1","E2-A5-R2","E2-A5-R3","E3-A1-R1","E3-A1-R2","E3-A1-R3","E3-A2-R1","E3-A2-R2","E3-A2-R3","E3-A3-R1","E3-A3-R2",
                                                        "E3-A3-R3","E3-A4-R1","E3-A4-R2","E3-A4-R3","E3-A5-R1","E3-A5-R2","E3-A5-R3"),exclude = NA),
             fill=count)) +
             #color=Shortcuts_Sample.Names,
             #pattern=as.factor(count))) + 
  geom_tile(na.rm = T,color="gray") +
  
  #scale_pattern_discrete(choices = count) +
  facet_wrap(~factor(parameters, levels=c("DDA with FAIMS",
                                           "DDA without FAIMS",
                                           "DIA with FAIMS",
                                           "DIA without FAIMS",
                                           "DIA TIMSToF",
                                           "DDA TIMSToF")),nrow = 1,scales = 'free_x') +
  #facet_grid(exp_type~.) +
  theme_bw() + 
  theme(legend.text = element_text(size=30), 
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        plot.title = element_text(size=30),
        legend.title=element_text(size=30),
        axis.text.x=element_text(size=20),
        axis.text.y=element_blank(),
        axis.title=element_text(size=30),
        axis.line = element_blank(),
        strip.text.x = element_text(size = 30)
  ) + labs(y="") +
  #scale_y_continuous(breaks = seq(from=0,to=40,by=5)) +
  #scale_fill_viridis(discrete = TRUE,option = "rocket",direction = -1)
  scale_fill_manual(values = c("#fcfdbf","#e75263"))
  
  

final <- grid.arrange(group_col, exp_summary, ncol = 2, widths = c(0.5, 4))

final1 <- grid.arrange(group_col, sample_inj, ncol = 2, widths = c(0.5, 4))
ggsave(filename=paste0("with_groups_sequence_of_sample_inj_representation.tiff"), 
       plot = final1 ,path = "D:/dev/Desktop_copy/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/list_of_all_collected_rawfiles/",
       scale = 1,device= "tiff",width = 2800,height = 1500, units = "px",dpi = 100)


ggsave(filename=paste0("with_groups_sequence_of_data_analysis_representation_current.tiff"), 
       plot = final ,path = "D:/dev/Desktop_copy/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/list_of_all_collected_rawfiles/",
       scale = 1,device= "tiff",width = 4000,height = 1500, units = "px",dpi = 100)


final_data <- analyzed_data %>% group_by(Shortcuts_Sample.Names,parameters) %>%
  left_join(analyzed_data,sample_injection,by=c("parameters","Shortcuts_Sample.Names"))
