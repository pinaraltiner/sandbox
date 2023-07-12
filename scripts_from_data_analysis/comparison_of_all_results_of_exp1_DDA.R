#
library(stringr)
library(dplyr)
library(data.table)
library(openxlsx)
library(ggplot2)

a <- setwd("D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_1/Proline/processing")
list_files <- list.files(a)
list_files <- list_files[list_files != c("comparison_of_all_results_of_exp1_DDA.R","results")]

final_ecoli_count <- NULL

fasta_syn_pep <- read.xlsx("D:/dev/Desktop_copy/PHD/wet_lab_experiments/Synthetic peptides list_theo_conc_added_080122.xlsx", sheet = "ISO-ref and OTHER with FC")
common_col_fasta_syn_pep <- as.data.frame(paste(fasta_syn_pep$Phospopeptide.sequence, fasta_syn_pep$modified.position.in.peptide, sep = "_"))
colnames(common_col_fasta_syn_pep) <- "common_col_for_merging"
common_col_fasta_syn_pep_new <- cbind(common_col_fasta_syn_pep,fasta_syn_pep)



for (i in 1:length(list_files)){
  assign(paste0(list_files[i],i),read.xlsx(paste0(list_files[i]), sheet = "Best PSM from protein sets"))
  assign(paste0(list_files[i],i,"filtered_ecoli"),get(paste0(list_files[i],i)) %>%
           filter(grepl("_ECOLI",accession)) %>%
           filter(!grepl("CON__",accession)))
  assign(paste0("ecoli_dim_",i), dim(get(paste0(list_files[i],i,"filtered_ecoli"))))
  
  ## where final_ecoli_count is used
  final_ecoli_count <- rbind(final_ecoli_count,get(paste0("ecoli_dim_",i)))

  
  assign(paste0(list_files[i],i),read.xlsx(paste0(list_files[i]), sheet = "Best PSM from protein sets"))
  assign(paste0(list_files[i],i,"filtered"),get(paste0(list_files[i],i)) %>% 
  filter(grepl("_HUMAN",accession)) %>% 
  filter(grepl("Phospho",modifications)) %>% 
  filter(!grepl("CON__",accession)))

  
  phospho_ptm_pos <- lapply(get(paste0(list_files[i],i,"filtered"))[,"modifications"], function(each_ptm_protein_positions) {
    
    ptm_list <- as.list(strsplit(each_ptm_protein_positions,"; ", fixed=TRUE)[[1]]) # Split ptm_protein_position depending on ";"
    ptm_list <- ptm_list[grepl("Phospho", ptm_list, fixed = TRUE)] # Extract only which contains "Phospho"
    phospho_positions <- lapply(ptm_list, function(ptm) { 
      #sub('Phospho \\(([A-Z]\\d+)\\)', "\\d+", ptm) #then, remove "Phospho" and remain only positions
      as.character(str_extract(ptm, "\\d+"))
      
    })
    
    phospho_positions_as_str <- paste(phospho_positions, collapse="&") #combine each position with "|"
    
  })
  
  phospho_ptm_pos_df <- t(as.data.frame(phospho_ptm_pos))
  rownames(phospho_ptm_pos_df) <- 1:length(phospho_ptm_pos_df)
  
  assign(paste0(list_files[i],i,"filtered","modifications"),cbind(get(paste0(list_files[i],i,"filtered")),paste(get(paste0(list_files[i],i,"filtered"))[,"sequence"], as.data.frame(phospho_ptm_pos_df)[,1], sep = "_")))
  rm(phospho_ptm_pos_df)
  
}

## TODO: automatize here
timsdata_E1_M1_M8wo_coli_pool_5414 <- `20230301-1825_TTP_002911_ONJ_TR_Toulouse_E1_M1_M8_pool_5414_DA6.0.434_1.mgf_2023-03-14_1447.xlsx1filteredmodifications`
timsdata_E1_M1_M8wcoli_pool_5415 <- `20230301-1913_TTP_002912_ONJ_TR_Toulouse_E1_M1_M8coli_pool_5415_DA6.0.434_1.mgf_2023-03-14_1446.xlsx2filteredmodifications`
timsdata_E1_M1_M8wo_coli_pool_5419_CE <- `20230302-1428_TTP_002916_ONJ_TR_Toulouse_E1_M1_M8_pool_CE_5419_DA6.0.434_1.mgf_2023-03-14_1447.xlsx3filteredmodifications`
timsdata_E1_M1_M8wcoli_pool_5420_CE <- `20230302-1517_TTP_002917_ONJ_TR_Toulouse_E1_M1_M8coli_pool_CE_5420_DA6.0.434_1.mgf_2023-03-14_1446.xlsx4filteredmodifications`


exp1_w_Ecoli_noFAIMS_DDA <- `PAL _Phosphopeptides exp1_w_Ecoli_noFAIMS_DDA_230208_2023-03-14_1424.xlsx5filteredmodifications`
exp1_wo_Ecoli_noFAIMS_DDA <- `PAL _Phosphopeptides exp1_wo_Ecoli_noFAIMS_DDA_230208_2023-03-14_1424.xlsx6filteredmodifications`
exp1_wo_Ecoli_FAIMS_DDA <- `PAL_Synthetic phosphopeptide_Mix1_8 pure_with_FAIMS_62.5fmol _from PDmgf_synthetic fasta_2023-03-14_1434.xlsx7filteredmodifications`
exp1_w_Ecoli_FAIMS_DDA <- `PAL_Synthetic phosphopeptide_Mix1_8_with_Ecoli_with FAIMS_62.5fmol_from_PDmgf_synthetic_fasta_2023-03-14_1432.xlsx8filteredmodifications`
#

rm(list = ls()[grepl(".xlsx", ls())])

obj_list <-  ls()[grepl("coli", ls())]


for (k in 1:length(obj_list)){

  ### DOES NOT WORK !!!! 
  ## https://stackoverflow.com/questions/15303972/r-get-function-error
  #assign("common_col_for_merging",colnames(get(obj_list[k]))[35])
  
  ## INSTEAD, tmp has been created just to change colname
  tmp <- get(paste0(obj_list[k]))
  colnames(tmp)[35] <- "common_col_for_merging"
  assign(paste0(obj_list[k]),tmp)
  
  
  assign(paste0("combine_",obj_list[k]),left_join(get(obj_list[k]),common_col_fasta_syn_pep,by="common_col_for_merging"))
  assign(paste0(obj_list[k],"merge"), merge(get(obj_list[k]),common_col_fasta_syn_pep,by.y = "common_col_for_merging", by.x = "common_col_for_merging"))
  
  write.table(get(paste0("combine_",obj_list[k])),file =paste0(obj_list[k],"combine",".tsv"), sep = "\t", row.names = F )
  write.table(get(paste0(obj_list[k],"merge")),file =paste0(obj_list[k],"merge",".tsv"), sep = "\t", row.names = F )
  
}



                            

final_table <- data.frame(c("timsdata_E1_M1_M8wo_coli_pool_5414",
                 "timsdata_E1_M1_M8wcoli_pool_5415",
                 "timsdata_E1_M1_M8wo_coli_pool_5419_CE",
                 "timsdata_E1_M1_M8wcoli_pool_5420_CE",
                 "exp1_w_Ecoli_noFAIMS_DDA",
                 "exp1_wo_Ecoli_noFAIMS_DDA",
                 "exp1_wo_Ecoli_FAIMS_DDA",
                 "exp1_w_Ecoli_FAIMS_DDA"),
                 #c(117,	111,	98,	105,	99,	142,	44,	92),
                 #c(6,	29,	19,	13,	42,	43,	17,	46),
                 #c(133,	140,	117,	118,	141,	185,	61,	138),
                 final_ecoli_count[,1])

#colnames(final_table) <- c("method","TP","missing","sum","E.coli peps")

colnames(final_table) <- c("method","E.coli peps")



library(ggplot2)
library(reshape2)
melt_final_table <- melt(final_table)

ggplot(data=melt_final_table, aes(x=value,y=method,fill=variable)) +
  geom_bar(stat = "identity", position = "dodge") + geom_label(aes(label=value)) + 
  theme_light() + #scale_y_continuous(limits = c(-60, 60),breaks = seq(-60, 60, by = 20)) +
  theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
        axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
        plot.title = element_text(size=20),
        legend.title=element_text(size=15),
        axis.text=element_text(size=15),
        axis.text.x = element_text(size=15, angle = 90),
        axis.title=element_text(size=15)
  ) + ggtitle("Comparison of num. identified Syn. Phospho-peptides btw instruments") + 
  scale_y_discrete(labels(c(phospho_pep_label,ecoli_label)))

annotate("segment", x = c(0, 3.5), xend = c(3.6, 3.8), y = c(14, 14), yend = c(15, 15))+
  coord_cartesian(clip = "off", ylim = c(15, 25)) 

phospho_pep_label <- seq(from=0,to=150,by=10)

ecoli_label <- seq(from=10000, to=22000,by=1000)

