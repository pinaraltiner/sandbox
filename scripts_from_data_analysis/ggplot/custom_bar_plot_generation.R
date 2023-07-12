library(reshape2)
library(ggplot2)
library(dplyr)
#library(patchwork)
library(ggpubr)
#SINCE YOU USE ACTUAL NUM OF PEP FOR EACH SAMPLE, YOU DON'T NEED TO ADD THIS COL AS ADDITIONAL GEOM_BAR()

data <- read.csv("D:/dev/Desktop_copy/PHD/wet_lab_experiments/QC_result/eyers_syn_pep_num_id_pep_with_custom_fasta.csv", sep = ",", header = T)
data <- data[-c(27:30),-c(3,4,8)]
for (i in 1:8){
  assign(paste0("data",i), filter(data, grepl(paste0("Pool",i,"_"), data$Pool_name)))
  assign(paste0("melt_data",i), melt(get(paste0("data",i))))
  p <- ggplot(data=get(paste0("melt_data",i)), aes(x=Pool_name, y=value, fill=variable)) + #reorder( +colnames_data)
    #geom_bar(stat="identity", color="black", width= 0.50) + 
    geom_col(width= 0.50,color="black",position = position_stack(reverse = TRUE)) +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=45, hjust=0.9)) + coord_flip() +
    scale_y_continuous(limits = c(0, 100)) + geom_text(aes(label=value), size=8,position=position_stack(reverse = TRUE,vjust = 0.5)) #position=position_dodge(width=0.9), vjust=-0.25)
  
  p <- p +  geom_bar(data= actual_num_pep[1,], width= 0.50, aes(x='Actual num of pep', y=value, fill=variable),
                     stat = "identity", color="black", position=position_dodge2()) + 
    scale_y_continuous(limits = c(0, 100)) + scale_fill_manual(values=c("blue","#CC79A7","#009E73", "#D55E00")) 
  
  assign(paste0("p",i),p)
  
}

#change color of actual number of pep and change legend title

ggarrange(p1, p2, p3, p4, p5,p6,p7,p8, ncol=4, nrow=2, common.legend = TRUE, legend="right")

p1+p2+p3+p4/p5+p6+p7+p8

colnames_data <- c("Num_correctly_identified_pep","Num_unidentified_pep","Num_unexpectedly_identified_pep")

#p <- p + geom_text(aes(label=ifelse(actual_num_pep[1,1] == "Pool1_BCC_62.5fmol_w_Ecoli_syn_fasta",
#                                    actual_num_pep[1,3] , aes(label=value), 
#                                    size=8,position=position_stack(reverse = TRUE,vjust = 0.5)), position = position_dodge(width = .9), vjust = 0))
  #data=actual_num_pep[1,],aes(label=actual_num_pep[1,1]), size=8)



data <- data[1:9,-c(1,4,8)]


actual_num_pep <- melt_data[1:9,]
p <- ggplot(data=melt_data[-c(1:9),], aes(x=Pool_name, y=value, fill=variable)) +
  geom_bar(stat="identity", color="black", width= 0.50) + 
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=0.9)) + coord_flip()

p <- p +  geom_bar(data= actual_num_pep, width= 0.50, aes(x=paste0(Pool_name,"1"), y=value, fill=variable),
                   stat = "identity", color="black", position=position_dodge2()) + scale_x_discrete(labels=c('Actual num of pep'))

p1


p <- ggplot(data=get(paste0("melt_data",1)), aes(x=Pool_name, y=value, fill=variable)) + #reorder( +colnames_data)
  #geom_bar(stat="identity", color="black", width= 0.50) + 
  geom_col(width= 0.50,color="black",position = position_stack(reverse = TRUE)) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=0.9)) + coord_flip() +
  scale_y_continuous(limits = c(0, 70)) + geom_text(aes(label=value), size=8,position=position_stack(reverse = TRUE,vjust = 0.5))#position=position_dodge(width=0.9), vjust=-0.25)

p <- p +  geom_bar(data= actual_num_pep[1,], width= 0.50, aes(x=paste0(Pool_name,"1"), y=value, fill=variable),
                   stat = "identity", color="black", position=position_dodge2()) + scale_y_continuous(limits = c(0, 70)) + scale_x_discrete(labels=c('Actual num of pep'))


