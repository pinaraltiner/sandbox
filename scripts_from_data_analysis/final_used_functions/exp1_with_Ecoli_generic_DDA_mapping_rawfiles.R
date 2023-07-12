#
library(stringr)
library(dplyr)
library(data.table)
library(openxlsx)
library(ggplot2)

a <- setwd("D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_1/Proline/exp1_with_ecoli_no_FAIMS/")
all_list_files <- list.files(a)
list_files <- all_list_files[(grepl(".xlsx",all_list_files))]

#list_files <- list_files[!(list_files %in% c("comparison_of_all_results_of_exp1_DDA.R","results"))]


theo_file_path <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/"
theo_file_name <- "Synthetic peptides list_theo_conc_corrected_pool_id_iso_count_final.xlsx" 

#theo_file_name <- "Synthetic peptides list_theo_conc_added_pool_id_iso_count.xlsx"  ### TODO CHANGE me
sheet_theo_name <- "ISO-ref and OTHER with FC"

map_file_path <- "D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_1/"
map_file_name <- "Experiment1_with_Ecoli_noFAIMS_DDA_raw_files_pool_id.txt"
#theo_file_name <- "Synthetic peptides list_theo_conc_added_pool_id_iso_count.xlsx"  ### TODO CHANGE me
sheet_theo_name <- "ISO-ref and OTHER with FC"


pep_list_w_theo_quant <- read.xlsx(paste0(theo_file_path, theo_file_name), sheet = sheet_theo_name)
pep_list_w_theo_quant <- pep_list_w_theo_quant[,-1]

map_file_raw_to_pool <- read.table(paste0(map_file_path, map_file_name),sep="\t",header = TRUE)


common_col_theo_quant <- as.data.frame(paste(pep_list_w_theo_quant$Phosphopeptide.sequence,
                                             pep_list_w_theo_quant$modified.position.in.peptide, sep = "_"))

colnames(common_col_theo_quant) <- "common_col_for_merging"

#pool_id_seq <- as.data.frame(paste(pep_list_w_theo_quant$Phosphopeptide.sequence,
#pep_list_w_theo_quant$pool_id, sep = "_"))

#colnames(pool_id_seq) <- "pool_id_seq"


pep_list_w_theo_quant_new <- cbind(common_col_theo_quant,pep_list_w_theo_quant)



pep_list_theo <- pep_list_w_theo_quant_new %>% 
  left_join(map_file_raw_to_pool, by="pool_id",multiple="all") %>%
  mutate(rawfile_seq=paste(Phosphopeptide.sequence,raw_file,sep = "."))


exp1_noFAIMS_with_Ecoli <- data.frame(pool_id = rep(1:8, each = 3),
                                      inj_id = rep(c("inj1","inj2","inj3"),8))
                             
final_ecoli_count <- NULL

for (i in 1:length(list_files)){
  assign(paste0("pool_",i),read.xlsx(paste0(list_files[i]), sheet = "Best PSM from protein sets"))
  assign(paste0("pool_",i,"_coli"),get(paste0("pool_",i)) %>%
           filter(grepl("ECOLI",accession)) %>%
           filter(!grepl("CON__",accession))) %>%
           filter(!grepl("HUMAN",accession))

  assign(paste0("pool_",i,"_phospho"),get(paste0("pool_",i)) %>%
           filter(grepl("HUMAN",accession)) %>%
           filter(grepl("Phospho",modifications)) %>%
           filter(!grepl("CON__",accession))) %>%
           filter(!grepl("ECOLI",accession))

  assign(paste0("pool_",i,"_phospho"), cbind(get(paste0("pool_",i,"_phospho")),paste0("pool",exp1_noFAIMS_with_Ecoli[i,1])))
  assign(paste0("ecoli_dim_",i), cbind(dim(get(paste0("pool_",i,"_coli")))[1],paste0("pool",exp1_noFAIMS_with_Ecoli[i,1])))
  
  final_ecoli_count <- rbind(final_ecoli_count,get(paste0("ecoli_dim_",i)))
}
  ## where final_ecoli_count is used
final_ecoli_count_df <- cbind(final_ecoli_count, exp1_noFAIMS_with_Ecoli[,-1])


### Collect all separated phospho data into one object
object_names <- paste0("pool_", 1:24, "_phospho")

object_list <- mget(object_names)

# Combine the objects into a single data frame
final_phospho_df <- do.call(rbind, object_list)
  
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
    mutate(common_col_for_merging = paste(sequence,phospho_ptm_pos_df, sep = "_"))


library(ggplot2)
  
  colnames(final_ecoli_count_df) <- c("num_id_Ecoli_peptides", "Pool","inj_id")
  
  as.data.frame(final_ecoli_count_df) %>% tibble() %>% mutate(comb_pool_inj = paste(Pool,inj_id, sep="_")) %>%
    ggplot(aes(y=num_id_Ecoli_peptides,
             x=comb_pool_inj, fill=inj_id)) +
    geom_col(position = "dodge",width = 0.98) + 
    geom_text(
      aes(x = comb_pool_inj, y = num_id_Ecoli_peptides, label = num_id_Ecoli_peptides),
      position = position_dodge(width = 1),
      vjust = -0.5, size = 5) + 
    scale_x_discrete(label=rep(c("pool1", "pool2","pool3",   ## Since ordering of numbers 
                                 "pool4","pool5", "pool6",    # based on injection does not work,
                                 "pool7", "pool8"),each=3)) +  # First combined pool and 
    theme_light() +                                            # inj_id and changed the x axis labels manually  
    theme(legend.text = element_text(size=20), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
          axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20),
          plot.title = element_text(size=30),
          legend.title=element_text(size=20),
          axis.text=element_text(size=20),
          axis.text.x = element_text(size=20, angle = 90),
          axis.title=element_text(size=20),
          plot.subtitle = element_text(size = 20)
    ) + ggtitle("Number of identified E.coli peptides in Exploris no FAIMS") +
    labs(subtitle = "Experiment 1 DDA processed by Proline", 
         x="Pool ID", y = "Number of E.coli peptides") +
    guides(fill = guide_legend(title = "Number of injection")) 
    
      


  
colnames(pep_list_theo)[2] <- "sequence"
library(tidyr)


# df_rmv_redunc <- final_phospho_df_common_col %>% 
#   separate(spectrum_title, into = c("a", "b", "c","d","e","f","raw_files"), sep = ";") %>%
#   select(!a:f) %>% 
#   separate(raw_files, into = c('tmp', 'raw_file'), sep = ":") %>%
#   mutate(common_col_for_merging = paste(sequence,phospho_ptm_pos_df, sep = "_")) %>%
#   mutate(pool_id_seq=paste(sequence,pool_id,raw_file,sep = "_"))
#   


df_rmv_redunc <- final_phospho_df_common_col %>% 
  separate(spectrum_title, into = c("a", "b", "c","d","e","f","spec_title"), sep = ";") %>%
  select(!a:f) %>%
  separate(spec_title, into = c("tmp","raw_files"), sep = ":")



#df_merge_localization <- df_rmv_redunc %>%
  # #separate(raw_files, into = c("a","raw_file"), sep = ":") %>%
  # group_by(common_col_for_merging,raw_files) %>%
  # slice(which.max(ptm_score)) %>%
  # ungroup() %>%
  # full_join(pep_list_w_theo_quant_new,by="common_col_for_merging")
  
df_merge_localization <- df_rmv_redunc %>%
  #separate(raw_files, into = c("a","raw_file"), sep = ":") %>%
  group_by(common_col_for_merging,raw_files) %>%
  slice(which.max(ptm_score)) %>%
  ungroup() %>%
  full_join(pep_list_theo,by="common_col_for_merging",multiple="all")



# df_merge_sequence <- df_rmv_redunc %>%
#   #separate(raw_files, into = c("a","raw_file"), sep = ":") %>%
#   group_by(sequence,raw_files) %>%
#   slice(which.max(ptm_score)) %>%
#   ungroup() %>%
#   full_join(pep_list_w_theo_quant_new, by="pool_id_seq")

df_merge_sequence <- df_rmv_redunc %>%
  mutate(rawfile_seq = paste(sequence,raw_files,sep=".")) %>%
  #separate(raw_files, into = c("a","raw_file"), sep = ":") %>%
  group_by(sequence,raw_files) %>%
  slice(which.max(ptm_score)) %>%
  ungroup() %>%
  full_join(pep_list_theo, by="rawfile_seq",multiple="all")


df_merge_localization_type <- df_merge_localization %>% 
  mutate(type=case_when(is.na(Protein) ~ "unexpected",
                        is.na(raw_files) ~ "missing",
                        TRUE ~ "correct"))


df_merge_sequence_type <- df_merge_sequence %>% 
  mutate(type=case_when(is.na(Protein) ~ "unexpected",
                        is.na(raw_files) ~ "missing",
                        TRUE ~ "correct"))


#write.table(df_merge_localization,file=paste0(a,"/","PAL_Exp1_with_coli_DDA_Exploris_npFAIMS_localization_merge_with_rawfiles.tsv"),sep = "\t",row.names = F,col.names = T)
# 
#write.table(df_merge_sequence,file=paste0(a,"/","PAL_Exp1_with_coli_DDA_Exploris_npFAIMS_sequence_merge_with_rawfiles.tsv"),sep = "\t",row.names = F,col.names = T)


df_merge_localization_type %>%
  select(common_col_for_merging,pool_id.x,
         raw_files,type) %>% 
  count(raw_files,type) %>% #drop_na() %>%
  ggplot(aes(x=raw_files,y=n,fill=type)) + 
  geom_col(position = "dodge",width = 0.98) + 
  geom_text(
    aes(x = raw_files, y = n, label = n, group = type),
    position = position_dodge(width = 1),
    vjust = -0.5, size = 5
  )+
  theme_light() + #scale_y_continuous(limits = c(-60, 60),breaks = seq(-60, 60, by = 20)) +
  theme(legend.text = element_text(size=20), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size=20),
        legend.title=element_text(size=20),
        axis.text=element_text(size=20),
        axis.text.x = element_text(size=20, angle = 90),
        axis.title=element_text(size=30)
  ) + ggtitle("Comparison of num. localized Syn. Phospho-peptides") +
  labs(subtitle = "Experiment 1 with E.coli no FAIMS DDA proceesed by Proline")


df_merge_localization_type %>%
  select(common_col_for_merging,
         raw_file,type) %>% 
  count(raw_file,type) %>% #drop_na() %>%
  ggplot(aes(x=raw_file,y=n,fill=type)) + 
  geom_col(position = "dodge",width = 0.98) + 
  geom_text(
    aes(x = raw_file, y = n, label = n, group = type),
    position = position_dodge(width = 1),
    vjust = -0.5, size = 5
  )+
  theme_light() + #scale_y_continuous(limits = c(-60, 60),breaks = seq(-60, 60, by = 20)) +
  theme(legend.text = element_text(size=20), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size=20),
        legend.title=element_text(size=20),
        axis.text=element_text(size=20),
        axis.text.x = element_text(size=20, angle = 90),
        axis.title=element_text(size=30)
  ) + ggtitle("Comparison of num. localized Syn. Phospho-peptides") +
  labs(subtitle = "Experiment 1 with E.coli no FAIMS DDA proceesed by Proline")


df_merge_sequence_type %>%
  separate(rawfile_seq, into = c("sequence","raw_files"), sep = "\\.") %>%
  count(raw_files,type) %>% #drop_na() %>%
  ggplot(aes(x=raw_files,y=n,fill=type)) + 
  geom_col(position = "dodge",width = 0.98) + geom_text(
    aes(x = raw_files, y = n, label = n, group = type),
    position = position_dodge(width = 1),
    vjust = -0.5, size = 5
  )+
  theme_light() + #scale_y_continuous(limits = c(-60, 60),breaks = seq(-60, 60, by = 20)) +
  theme(legend.text = element_text(size=20), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
        axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20),
        plot.title = element_text(size=20),
        legend.title=element_text(size=20),
        axis.text=element_text(size=20),
        axis.text.x = element_text(size=20, angle = 90),
        axis.title=element_text(size=30)
  )  + ggtitle("Comparison of num. identified Syn. Phospho-peptides") +
  labs(subtitle = "Experiment 1 no E.coli no FAIMS DDA proceesed by Proline")




df_merge_sequence_type %>% select(rawfile_seq,ptm_score,type) %>% 
  separate(rawfile_seq, into = c("sequence","raw_files"), sep = "\\.") %>%
  drop_na(ptm_score) %>%
  #count(pool_id) %>%
  ggplot(aes(x= ptm_score)) +
  geom_density(aes(fill=type),alpha=0.5) + 
  theme_light() + #scale_y_continuous(limits = c(-60, 60),breaks = seq(-60, 60, by = 20)) +
  theme(legend.text = element_text(size=20), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
        axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20),
        plot.title = element_text(size=30),
        legend.title=element_text(size=20),
        axis.text=element_text(size=20),
        axis.text.x = element_text(size=20, angle = 90,hjust=0.95,vjust=0.2),
        axis.title=element_text(size=30),
        strip.text.x = element_text(size = 20),
        plot.subtitle = element_text(size = 20)
  ) + ggtitle("Distribution of PTM-score of identified Syn. Phospho-peptides across each run") +
  labs(subtitle = "Experiment 1 with E.coli no FAIMS DDA proceesed by Proline") +
  scale_y_continuous(limits = c(0,6))+
  facet_wrap(~raw_files,ncol= 8)



###############################################################################

## OLDER VERSION OF THESE PLOTS WERE FOUND BELOW!!!!


df_merge_sequence_type %>% select(rawfile_seq,ptm_score,type,isomeric_count) %>% 
  separate(rawfile_seq, into = c("sequence","raw_files"), sep = "\\.") %>% 
  drop_na(ptm_score) %>%
  #count(pool_id) %>%
  ggplot(aes(y=raw_files,x= sequence)) +
  geom_point(aes(size=ptm_score,color=raw_files)) + 
  theme_light() + #scale_y_continuous(limits = c(-60, 60),breaks = seq(-60, 60, by = 20)) +
  theme(legend.text = element_text(size=20), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
        axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20),
        plot.title = element_text(size=20),
        legend.title=element_text(size=20),
        axis.text=element_text(size=20),
        axis.text.x = element_text(size=15, angle = 90,hjust=0.95,vjust=0.2),
        axis.title=element_text(size=30),
        plot.subtitle = element_text(size = 20)
  ) + ggtitle("Distribution of identified Syn. Phospho-peptides across each pool") +
  labs(subtitle = "Experiment 1 with E.coli no FAIMS DDA proceesed by Proline") 
#facet_wrap(~type)




df_merge_sequence_type %>% select(rawfile_seq,ptm_score,type,isomeric_count) %>% 
  separate(rawfile_seq, into = c("sequence","raw_files"), sep = "\\.") %>% 
  mutate(seq_iso_count = paste(sequence,isomeric_count,sep="_")) %>%
  drop_na(ptm_score) %>%
  #count(pool_id) %>%
  ggplot(aes(x=seq_iso_count,y= ptm_score)) +
  geom_col(aes(fill=type)) + 
  theme_light() + #scale_y_continuous(limits = c(-60, 60),breaks = seq(-60, 60, by = 20)) +
  theme(legend.text = element_text(size=20), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
        axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20),
        plot.title = element_text(size=20),
        legend.title=element_text(size=20),
        axis.text=element_text(size=20),
        axis.text.x = element_text(size=15, angle = 90,hjust=0.95,vjust=0.2),
        axis.title=element_text(size=30),
        strip.text.x = element_text(size = 20),
        plot.subtitle = element_text(size = 20)
  ) + ggtitle("Distribution of identified Syn. Phospho-peptides across each pool") +
  labs(subtitle = "Experiment 1 with E.coli no FAIMS DDA proceesed by Proline") +
  facet_wrap(~raw_files,ncol= 1)






df_merge_localization_type %>% select(common_col_for_merging,pool_id.x,ptm_score,type) %>% 
  drop_na(ptm_score) %>%
  #count(pool_id) %>%
  ggplot(aes(y=common_col_for_merging,x=pool_id.x)) +
  geom_point(aes(size=ptm_score,color=pool_id.x)) + 
  theme_light() + #scale_y_continuous(limits = c(-60, 60),breaks = seq(-60, 60, by = 20)) +
  theme(legend.text = element_text(size=20), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
        axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20),
        plot.title = element_text(size=20),
        legend.title=element_text(size=20),
        axis.text=element_text(size=20),
        axis.text.x = element_text(size=15, angle = 90,hjust=0.95,vjust=0.2),
        axis.title=element_text(size=30),
        plot.subtitle = element_text(size = 20)
  ) + ggtitle("Distribution of localized Syn. Phospho-peptides across each pool") +
  labs(subtitle = "Experiment 1 with E.coli no FAIMS DDA proceesed by Proline") +
  facet_wrap(~type)



df_merge_sequence %>% select(rawfile_seq) %>% 
  separate(rawfile_seq, into = c("sequence","rawfiles"), sep = "\\.") %>%
  count(rawfiles) %>%
  ggplot(aes(x=rawfiles,y=n,fill=rawfiles)) +
  geom_col() + geom_label(aes(label=n)) + 
  theme_light() + #scale_y_continuous(limits = c(-60, 60),breaks = seq(-60, 60, by = 20)) +
  theme(legend.text = element_text(size=20), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
        axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20),
        plot.title = element_text(size=20),
        legend.title=element_text(size=20),
        axis.text=element_text(size=20),
        axis.text.x = element_text(size=20, angle = 90),
        axis.title=element_text(size=30)
  ) + ggtitle("Comparison of num. identified Syn. Phospho-peptides btw instruments")

