#
library(openxlsx)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

proline_phospho_pos_extraction <- function(df){
  # Extraction of phospho positions from quant peptides object
  phospho_ptm_pos <- lapply(df, function(each_ptm_protein_positions) {
    
    ptm_list <- as.list(strsplit(each_ptm_protein_positions,"; ", fixed=TRUE)[[1]]) # Split ptm_protein_position depending on ";"
    ptm_list <- ptm_list[grepl("Phospho", ptm_list, fixed = TRUE)] # Extract only which contains "Phospho"
    phospho_positions <- lapply(ptm_list, function(ptm) { 
      #sub('Phospho \\(([A-Z]\\d+)\\)', "\\d+", ptm) #then, remove "Phospho" and remain only positions
      as.character(str_extract(ptm, "\\d+"))
      
    })
    
    phospho_positions_as_str <- paste(phospho_positions, collapse="&") #combine each position with "|"
    
  })
  # Data conversion 
  phospho_ptm_pos_df <- t(as.data.frame(phospho_ptm_pos))
  rownames(phospho_ptm_pos_df) <- 1:length(phospho_ptm_pos_df)
  
  return(phospho_ptm_pos_df)
  
}

theo_file_path= "D:/dev/Pinar/PHD/wet_lab_experiments/Eyers_syn_peptides_experiment/"
theo_file_name="Synthetic peptides list_theo_conc_corrected_isomericity.xlsx"#"Synthetic peptides list_theo_conc_corrected_pool_id_iso_count_final.xlsx"
sheet_theo_name = "ISO-ref and OTHER with FC"
pep_list_w_theo <- read.xlsx(paste0(theo_file_path, theo_file_name), sheet = sheet_theo_name)
pep_list_w_theo_quant <- pep_list_w_theo[,-1]

common_col_theo_quant <- as.data.frame(paste(pep_list_w_theo_quant$Phosphopeptide.sequence,
                                             pep_list_w_theo_quant$modified.position.in.peptide, sep = "_"))

colnames(common_col_theo_quant) <- "pep_with_pos"
pep_list_w_theo_quant_new <- cbind(common_col_theo_quant,pep_list_w_theo_quant)


file_paths <- paste0("D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/",c("exp2_re_injection/","with_FAIMS/"))
file_names <- c("PAL _Phosphopeptides exp2 ( 5 conc 3reps) DDA_230117 with Design_2023-06-07_1003.xlsx",
               "PAL _Exp2_( 5 conc 3reps)_withFAIMS_DDA_25052023_noDesign - correct_2023-06-13_1514.xlsx")

file_ids <- c("nofaims","faims")

for (i in 1:length(file_paths)){
  assign(file_ids[i], read.xlsx(paste0(file_paths[i],file_names[i]),sheet =  "Best PSM from protein sets"))
}

faims_cols <- faims %>% select(sequence,
                         modifications,
                         psm_score,
                         calculated_mass,
                         spectrum_title,
                         accession,
                         is_validated_for_quanti,
                         rt) %>%
  filter(grepl("HUMAN",accession)& !grepl("CON__",accession)) %>%
  filter(grepl("Phospho",modifications)) %>%
  mutate(acq_type="faims") %>%
  separate(spectrum_title,into = c("first_cycle",
                                   "last_cycle",
                                   "first_scan",
                                   "last_scan",
                                   "first_time",
                                   "last_time",
                                   "raw_file",
                                   "tmp"),sep = ";") %>%
  select(!c(first_cycle,
            last_cycle,
            last_scan,
            first_time,
            last_time,
            tmp)) %>%
  separate(first_scan,into = c("tmp","first_scan"),sep = ":") %>% select(!tmp) %>%
  separate(raw_file,into = c("tmp","raw_file"),sep = ":") %>% select(!tmp) 

faims_cols <- faims_cols %>%
  mutate(pep_with_pos = paste(sequence,proline_phospho_pos_extraction(faims_cols$modifications),sep="_"))


nofaims_cols <- nofaims %>% select(sequence,
                               modifications,
                               psm_score,
                               calculated_mass,
                               spectrum_title,
                               accession,
                               is_validated_for_quanti,
                               rt) %>%
  filter(grepl("HUMAN",accession) & !grepl("CON__",accession)) %>% mutate(acq_type="no_faims") %>% 
  filter(grepl("Phospho",modifications)) %>%
  separate(spectrum_title,into = c("first_cycle",
                                   "last_cycle",
                                   "first_scan",
                                   "last_scan",
                                   "first_time",
                                   "last_time",
                                   "raw_file",
                                   "tmp"),sep = ";") %>%
  select(!c(first_cycle,
            last_cycle,
            last_scan,
            first_time,
            last_time,
            tmp)) %>%
  separate(first_scan,into = c("tmp","first_scan"),sep = ":") %>% select(!tmp) %>%
  separate(raw_file,into = c("tmp","raw_file"),sep = ":") %>% select(!tmp) 

nofaims_cols <- nofaims_cols %>%
  mutate(pep_with_pos = paste(sequence,proline_phospho_pos_extraction(nofaims_cols$modifications),sep="_"))


###############################################################################
merge_df_rowwise <- bind_rows(nofaims_cols,faims_cols)
merge_df_rowwise_phospho <- merge_df_rowwise %>% filter(grepl("Phospho",modifications))
###############################################################################

colnames(faims_cols) <- c(paste0("faims_",colnames(faims_cols)[-c(10,11)]),colnames(faims_cols)[c(10,11)])
colnames(nofaims_cols) <- c(paste0("nofaims_",colnames(nofaims_cols)[-c(10,11)]),colnames(nofaims_cols)[c(10,11)])


merge_df <- nofaims_cols %>%
  full_join(faims_cols,by = "pep_with_pos") %>%
  left_join(pep_list_w_theo_quant_new,by="pep_with_pos") %>%
  mutate(Pool= ifelse(is.na(Pool),"Unexpected",Pool))

for(i in 1:nrow(merge_df)){
  if((is.na(merge_df$acq_type.x)[i] ==  FALSE) & (is.na(merge_df$acq_type.y)[i] == FALSE)){
    
    if(merge_df$nofaims_sequence[i] == merge_df$faims_sequence[i]){
      
      merge_df$result[i] <- "Position found in both"
      
    }else{
      
      merge_df$result[i] <- "Different Localization."
      
    }
  }else{
    
    merge_df$result[i] <- paste(merge_df$acq_type.x[i],merge_df$acq_type.y[i],sep = "&")
    
  }
}
p1 <- merge_df %>% group_by(result) %>% count(result) %>%
  bind_cols(ypos = c(298,150,250)) %>%
  ggplot(aes(x="", y=n, fill=result)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_void() + labs(title = "Relation between with/out FAIMS based on the identified phospho-peptides ")+
  geom_text(aes(y = ypos, label = n), color = "white", size=15) +
  scale_fill_brewer(palette="Set1")


p2 <- merge_df %>% group_by(result,Pool) %>% 
  count(result) %>% ggplot(aes(x="",y=n,fill=Pool)) +
  geom_bar(stat = "identity",width = 1) +
  facet_wrap(~result) +
  theme_bw() +
  geom_text(aes(label=(n),size=75),position = position_stack(vjust = 0.5),
            show.legend = FALSE) 
  
only_with_faims <- merge_df %>% filter(grepl("NA&faims",result))
only_no_faims <- merge_df %>% filter(grepl("no_faims&NA",result))



p3 <- ggplot(data=merge_df_rowwise,
       aes(x=psm_score,
           fill=acq_type)) +
  geom_density(alpha=0.6) + theme_bw() +
  scale_fill_brewer(palette = "Set1") + 
  labs(title = "Density Plot of PSM score data from DDA with/out FAIMS")

p4 <- ggplot(data=merge_df_rowwise,
       aes(x=rt,
           fill=acq_type)) +
  geom_density(alpha=0.6) + theme_bw() +
  scale_fill_brewer(palette = "Set1") + 
  labs(title = "Density Plot of Retention Time data from DDA with/out FAIMS")

p5 <- ggplot(data=merge_df_rowwise,
       aes(x=is_validated_for_quanti,
           fill=acq_type)) +
  geom_bar(alpha=0.6,position = "dodge") + theme_bw() +
  scale_fill_brewer(palette = "Set1") + 
  geom_text(aes(label=after_stat(count)),
            stat = "count", position=position_dodge(width=0.9), vjust=-0.25,size=10) +
  labs(title = "Bar Plot of Quantified peptides data from DDA with/out FAIMS")

sapply(1:5,function(x) ggsave(filename = paste0("p",x,".tiff"),
                               width = 50, height = 40, 
                               path ="D:/dev/Pinar/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/comparison_dda_with_without_faims/",
                               units = "cm",
                               get(paste0("p",x)),
                               device = "tiff", #".svg"
))


