

df <- read.delim("D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/comparision_DDA_noFAIMS_mq_proline_pd/comparison_wrong_localized_pep_MQ_Proline_with_out_FAIMS/number_id_pep_and_seq_comparison_MQ_both_version_Proline_PD_with_and_without_FAIMS.txt")
df_noFAIMS <- df[,-3]
df_no_ecoli <- df_noFAIMS[-2,]
colnames(df_no_ecoli)[c(3,5)] <- c("MaxQuant v2.1.4","MaxQuant v1.6.3.4") 
df_long <- df_no_ecoli %>% pivot_longer(values_to = "number_of_pep",cols = -type)

ggplot(df_long, aes(x = type, y = number_of_pep,fill=name)) +
  geom_col(position="dodge") +
  theme_bw() +
  geom_text(aes(label=number_of_pep), position=position_dodge(width=0.9), vjust=-0.25)+
  theme(legend.text = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 30),
        legend.title = element_text(size = 15),
        axis.text.x = element_text(size = 15,angle = 90),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 15)) +
  expand_limits(x = 0, y = 0) +
  labs(title = "Experiment 2 - DDA no FAIMS - Software Comparison",x="Category" ,fill = "Software Names", y="Number of Peptides")



roc_curve_old_v214 <- read.delim("D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/MQ_data_analysis/exp2_wo_FAIMS/exp2_wo_FAIMSwith_MBR/roc_curve_data_exp2_noFAIMS_MQ_min_mode_score_40_consider_isoref_fp.txt") 
#roc_curve_combine_proline_mq_score_40 <- read.delim("D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/Proline_data_analysis/exp2_re_injection/exp2_no_FAIMS_proline_roc_curve_data_consider_isoref_fp.tsv")

#roc_curve_pd_wo_FAIMS <- read.delim("D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/PD_data_analysis/target_decoy_no_FAIMS/new/ROC_curve_data_PD_target_decoy_noFAIMS.txt")
roc_curve_MQ_1634 <- read.delim("D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/MQ_data_analysis/exp2_wo_FAIMS/rerun_using_1634/roc_curve_data_exp2_noFAIMS_MQ_version1634_rename.txt")
roc_curve_MQ_214_optimized<- read.delim("D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/MQ_214_optimize_param_noFAIMS/roc_curve_data_exp2_noFAIMS_MQ_data_optimized_param_v214_consider_isoref_fp.txt")

#roc_curve_proline <- roc_curve_combine_proline_mq_score_40 %>% filter(!grepl("MaxQuant",Software))


roc_mq_1634 <- roc_curve_data_exp2_noFAIMS_MQ_min_mode_score_40_1634_consider_isoref_fp[,-6]
colnames(roc_mq_1634) <- colnames(exp2_no_FAIMS_proline_roc_curve_data_consider_isoref_fp)

roc_mq_214 <- roc_curve_data_exp2_noFAIMS_MQ_min_mode_score_40_consider_isoref_fp[,-6]
colnames(roc_mq_214) <- colnames(exp2_no_FAIMS_proline_roc_curve_data_consider_isoref_fp)

roc_curve_MQ_214_optimized$type <-"v.2.1.4_optimized"
roc_curve_old_v214$type <-"v.2.1.4_custom"
roc_curve_MQ_1634$type <-"v.1.6.3.4"

final_roc_df <- rbind(roc_curve_MQ_214_optimized[,c("fdp","tpr","type","Comparison")],
                      roc_curve_old_v214[,c("fdp","tpr","type","Comparison")],
                      roc_curve_MQ_1634[,c("fdp","tpr","type","Comparison")])
                      #ROC_curve_data_PD_target_decoy_noFAIMS_consider_isoref_false)

write.table(final_roc_df,file="ROC_curve_comparison_btw_MQversions.tsv",sep = "\t",row.names = F)
library(plotROC)
rocplot <- ggplot(df, aes(m = predictions, d = labels))+ geom_roc(n.cuts=20,labels=FALSE)
rocplot + style_roc(theme = theme_grey) + geom_rocci(fill="pink") 


line_types <- c("solid","dotted","dotdash","dashed" ) #
##TODO: make it more professional
# ggplot(final_roc_df,aes(x=fdp,y=tpr)) + 
#   #geom_line(size=1.2,aes(color=Software)) +
#   geom_line(n.cuts=20,labels=FALSE) +
#   scale_color_manual(values = c("#CC79A7", "#E69F00", "#56B4E9", "#009E73"))+
#                      #labels=c('A1 vs A2','A1 vs A3','A1 vs A4','A1 vs A5')) +
#   scale_linetype_manual(values = line_types)+
#   theme_bw() +
#   theme(legend.text = element_text(size=15), 
#         axis.title.x = element_text(size = 15),
#         axis.title.y = element_text(size = 15),
#         plot.title = element_text(size=30),
#         legend.title=element_text(size=15),
#         axis.text.x=element_text(size=15),
#         axis.title=element_text(size=15),
#         axis.text.y = element_text(size = 15)) + 
#   expand_limits(x = 0, y = 0) +
#   scale_y_continuous(limits = c(0,100)) +
#   #facet_wrap(~Ratio_col,scales = "free_x") +
#   labs(title = "Experiment 2 - DDA without FAIMS - \n Comparision of MaxQuant version",
#        y="True Positive Rate \n (Sensitivity)",
#        x="False Positive Rate \n  (1-specificity)")
# 


ggplot(final_roc_df,aes(x=fdp,y=tpr,color=Comparison,linetype=type)) + geom_line(size=1) +
  scale_color_manual(values = c("#CC79A7", "#E69F00", "#56B4E9", "#009E73"))+
  scale_linetype_manual(values = line_types)+
  theme_bw() +
  theme(legend.text = element_text(size=15), 
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size=30),
        legend.title=element_text(size=15),
        axis.text.x=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.y = element_text(size = 15)) + 
  expand_limits(x = 0, y = 0) +
  scale_y_continuous(limits = c(0,100)) +
  labs(title = "Experiment 2 - DDA without FAIMS - \n Comparision of MaxQuant version",
                                              y="True Positive Rate \n (Sensitivity)",
                                              x="False Positive Rate \n  (1-specificity)")


roc_curve_data_exp2_noFAIMS_MQ_min_mode_score_40_1634_consider_isoref_fp <- read.delim("D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/MQ_data_analysis/exp2_wo_FAIMS/rerun_using_1634/roc_curve_data_exp2_noFAIMS_MQ_min_mode_score_40_1634_consider_isoref_fp.txt")
ROC_curve_data_PD_target_decoy_noFAIMS_consider_isoref_false <- read.delim("D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/PD_data_analysis/target_decoy_no_FAIMS/new/ROC_curve_data_PD_target_decoy_noFAIMS_consider_isoref_false.txt")
exp2_no_FAIMS_proline_roc_curve_data_consider_isoref_fp <- read.delim("D:/dev/Desktop_copy/PHD/wet_lab_experiments/DDA_data_analysis/experiment_2/Proline_data_analysis/exp2_re_injection/exp2_no_FAIMS_proline_roc_curve_data_consider_isoref_fp.tsv")








