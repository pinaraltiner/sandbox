#
library(stringr)
library(dplyr)
library(data.table)
library(openxlsx)
library(ggplot2)

a <- setwd("D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_1/Proline_timstof_no_Ecoli/")
list_files <- list.files(a)[-c(1,2,11,12)]

#list_files <- list_files[!(list_files %in% c("comparison_of_all_results_of_exp1_DDA.R","results"))]


theo_file_path <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/"
theo_file_name <- "Synthetic peptides list_theo_conc_added_pool_id.xlsx"

#theo_file_name <- "Synthetic peptides list_theo_conc_added_pool_id_iso_count.xlsx"  ### TODO CHANGE me
sheet_theo_name <- "ISO-ref and OTHER with FC"


pep_list_w_theo_quant <- read.xlsx(paste0(theo_file_path, theo_file_name), sheet = sheet_theo_name)
pep_list_w_theo_quant <- pep_list_w_theo_quant[,-1]

common_col_theo_quant <- as.data.frame(paste(pep_list_w_theo_quant$Phosphopeptide.sequence,
                                             pep_list_w_theo_quant$modified.position.in.peptide, sep = "_"))
pool_id_seq <- as.data.frame(paste(pep_list_w_theo_quant$Phosphopeptide.sequence,
                                             pep_list_w_theo_quant$X12, sep = "_"))
colnames(common_col_theo_quant) <- "common_col_for_merging"
colnames(pool_id_seq) <- "pool_id_seq"
pep_list_w_theo_quant_new <- cbind(common_col_theo_quant,pep_list_w_theo_quant,pool_id_seq)


#final_ecoli_count <- NULL

for (i in 1:length(list_files)){
  assign(paste0("pool_",i),read.xlsx(paste0(list_files[i]), sheet = "Best PSM from protein sets"))
  # assign(paste0("pool_",i,"_coli"),get(paste0("pool_",i)) %>%
  #          filter(grepl("ECOLI",accession)) %>%
  #          filter(!grepl("CON__",accession))) %>%
  #          filter(!grepl("HUMAN",accession))

  assign(paste0("pool_",i,"_phospho"),get(paste0("pool_",i)) %>%
           filter(grepl("HUMAN",accession)) %>%
           filter(grepl("Phospho",modifications)))
           #filter(!grepl("CON__",accession))) %>%
           #filter(!grepl("ECOLI",accession))

  assign(paste0("pool_",i,"_phospho"), cbind(get(paste0("pool_",i,"_phospho")),paste0("pool",i)))
  #assign(paste0("ecoli_dim_",i), cbind(dim(get(paste0("pool_",i,"_coli")))[1],paste0("pool",i)))
  
  #final_ecoli_count <- rbind(final_ecoli_count,get(paste0("ecoli_dim_",i)))
}
  ## where final_ecoli_count is used
  

  final_phospho_df <- rbind(pool_1_phospho,pool_2_phospho,
                            pool_3_phospho,pool_4_phospho,
                            pool_5_phospho,pool_6_phospho,
                            pool_7_phospho,pool_8_phospho)
  
  colnames(final_phospho_df)[35] <- "pool_id"
  
  
  phospho_ptm_pos <- lapply(final_phospho_df$modifications, function(each_ptm_protein_positions) {
    
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
  
  final_phospho_df_common_col <- final_phospho_df %>%
    mutate(common_col_for_merging = paste(sequence,phospho_ptm_pos_df, sep = "_")) %>%
    mutate(pool_id_seq=paste(sequence,pool_id,sep = "_"))


library(ggplot2)
  # 
  # colnames(final_ecoli_count) <- c("num_id_Ecoli_peptides", "Pool")
  # 
  # ggplot(as.data.frame(final_ecoli_count), aes(y=num_id_Ecoli_peptides,x=Pool,fill=Pool)) +
  #   geom_col() + geom_label(aes(label=num_id_Ecoli_peptides)) + 
  #   theme_light() + #scale_y_continuous(limits = c(-60, 60),breaks = seq(-60, 60, by = 20)) +
  #   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
  #         axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
  #         plot.title = element_text(size=20),
  #         legend.title=element_text(size=15),
  #         axis.text=element_text(size=15),
  #         axis.text.x = element_text(size=15, angle = 90),
  #         axis.title=element_text(size=15)
  #   ) + ggtitle("Number of identified E.coli peptides in TIMS-TOF Pro")


  
colnames(pep_list_w_theo_quant_new)[2] <- "sequence"
library(tidyr)


df_rmv_redunc <- final_phospho_df_common_col %>% 
  separate(spectrum_title, into = c("a", "b", "c","d","e","f","raw_files"), sep = ";") %>%
  select(!a:f)  

df_merge_localization <- df_rmv_redunc %>%
  #separate(raw_files, into = c("a","raw_file"), sep = ":") %>%
  group_by(common_col_for_merging,raw_files) %>%
  slice(which.max(ptm_score)) %>%
  ungroup() %>%
  full_join(pep_list_w_theo_quant_new,by="common_col_for_merging")
  
  

df_merge_sequence <- df_rmv_redunc %>%
  #separate(raw_files, into = c("a","raw_file"), sep = ":") %>%
  group_by(sequence,raw_files) %>%
  slice(which.max(ptm_score)) %>%
  ungroup() %>%
  full_join(pep_list_w_theo_quant_new, by="pool_id_seq")

#write.table(df_merge_localization,file=paste0(a,"/","PAL_Exp1_no_Ecoli_DDA_TIMSdata_localization_merge.tsv"),sep = "\t",row.names = F,col.names = T)
# 
#write.table(df_merge_sequence,file=paste0(a,"/","PAL_Exp1_no_Ecoli_DDA_TIMSdata_sequence_merge.tsv"),sep = "\t",row.names = F,col.names = T)


df_merge_sequence %>% select(pool_id_seq) %>% 
  separate(pool_id_seq, into = c("sequence","pool_id"), sep = "_") %>%
  count(pool_id) %>%
  ggplot(aes(x=pool_id,y=n,fill=pool_id)) +
  geom_col() + geom_label(aes(label=n)) + 
  theme_light() + #scale_y_continuous(limits = c(-60, 60),breaks = seq(-60, 60, by = 20)) +
  theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
        axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
        plot.title = element_text(size=20),
        legend.title=element_text(size=15),
        axis.text=element_text(size=15),
        axis.text.x = element_text(size=15, angle = 90),
        axis.title=element_text(size=15)
  ) + ggtitle("Comparison of num. identified Syn. Phospho-peptides btw instruments")



df_merge_localization %>% select(common_col_for_merging,
                                         pool_id) %>% 
  count(pool_id) %>% drop_na() %>%
  ggplot(aes(x=pool_id,y=n,fill=pool_id)) +
  geom_col() + geom_label(aes(label=n)) + 
  theme_light() + #scale_y_continuous(limits = c(-60, 60),breaks = seq(-60, 60, by = 20)) +
  theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
        axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
        plot.title = element_text(size=20),
        legend.title=element_text(size=15),
        axis.text=element_text(size=15),
        axis.text.x = element_text(size=15, angle = 90),
        axis.title=element_text(size=15)
  ) + ggtitle("Comparison of num. localized Syn. Phospho-peptides btw instruments")

