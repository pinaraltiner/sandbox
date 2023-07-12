library(tidyverse)
library(limma)
library(QFeatures)
library(msqrob2)
library(plotly)

setwd("D:/dev/Desktop_copy/PHD/wet_lab_experiments/PD_data_analysis/")
data <- read.table("Multiconsensus_TCellsExp3_230224_PeptideGroups.txt", sep = "\t", header = T)
data_abundace <- data[,64:78]
colnames(data_abundace) <- c("Abundance.F27.Sample.Ecoli.A.1.500.Tcells.phosphopeptides", "Abundance.F32.Sample.Ecoli.A.2.500.Tcells.phosphopeptides",
                             "Abundance.F37.Sample.Ecoli.A.3.500.Tcells.phosphopeptides", "Abundance.F26.Sample.Ecoli.B.1.500.Tcells.phosphopeptides" ,
                             "Abundance.F31.Sample.Ecoli.B.2.500.Tcells.phosphopeptides", "Abundance.F36.Sample.Ecoli.B.3.500.Tcells.phosphopeptides" ,
                             "Abundance.F25.Sample.Ecoli.C.1.500.Tcells.phosphopeptides", "Abundance.F30.Sample.Ecoli.C.2.500.Tcells.phosphopeptides" ,
                             "Abundance.F35.Sample.Ecoli.C.3.500.Tcells.phosphopeptides", "Abundance.F24.Sample.Ecoli.D.1.500.Tcells.phosphopeptides" ,
                             "Abundance.F29.Sample.Ecoli.D.2.500.Tcells.phosphopeptides", "Abundance.F34.Sample.Ecoli.D.3.500.Tcells.phosphopeptides" ,
                             "Abundance.F23.Sample.Ecoli.E.1.500.Tcells.phosphopeptides", "Abundance.F28.Sample.Ecoli.E.2.500.Tcells.phosphopeptides" ,
                             "Abundance.F33.Sample.Ecoli.E.3.500.Tcells.phosphopeptides")

rmv_na <- na.omit(data_abundace)
# #c("Abundance_F27_Sample_E_Coli_A1_1_500_Tcells_phosphopeptides", "Abundance_F32_Sample_E_Coli_A1_2_500_Tcells_phosphopeptides",
# "Abundance_F37_Sample_E_Coli_A1_3_500_Tcells_phosphopeptides", "Abundance_F26_Sample_E_Coli_A2_1_500_Tcells_phosphopeptides" ,
# "Abundance_F31_Sample_E_Coli_A2_2_500_Tcells_phosphopeptides", "Abundance_F36_Sample_E_Coli_A2_3_500_Tcells_phosphopeptides" ,
# "Abundance_F25_Sample_E_Coli_A3_1_500_Tcells_phosphopeptides", "Abundance_F30_Sample_E_Coli_A3_2_500_Tcells_phosphopeptides" ,
# "Abundance_F35_Sample_E_Coli_A3_3_500_Tcells_phosphopeptides", "Abundance_F24_Sample_E_Coli_A4_1_500_Tcells_phosphopeptides" ,
# "Abundance_F29_Sample_E_Coli_A4_2_500_Tcells_phosphopeptides", "Abundance_F34_Sample_E_Coli_A4_3_500_Tcells_phosphopeptides" ,
# "Abundance_F23_Sample_E_Coli_A5_1_500_Tcells_phosphopeptides", "Abundance_F28_Sample_E_Coli_A5_2_500_Tcells_phosphopeptides" ,
# "Abundance_F33_Sample_E_Coli_A5_3_500_Tcells_phosphopeptides")

pe <- readQFeatures(
  table = rmv_na,
  #fnames = 1,
  ecol = 1:15,
  name = "peptideRaw", sep="\t")
colnames(pe)

cond <- which(
   strsplit(colnames(pe)[[1]][4], split = "")[[1]] == "B") # find where condition is stored
# #24
# colData(pe)$condition <- substr(colnames(pe), cond, cond) %>%
#   unlist %>%  
#   as.factor
# 

colData(pe)$rep <- rep(rep(paste0("rep",1:3)),5) %>% as.factor
colData(pe)$condition <- substr(colnames(pe), cond, cond) %>%
     unlist %>%  
     as.factor
colData(pe)$spikeConcentration <- rep(c(A = 170, B = 85, C = 17, D = 8.5, E = 1.7),each = 3)





rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)
pe <- zeroIsNA(pe, "peptideRaw") # convert 0 to NA

pe <- logTransform(pe, base = 2, i = "peptideRaw", name = "peptideLog")


#Drop peptides that were only identified in one sample
# this is the threshold -> if we see one peptide in the at least
# two samples we will keep them. 

 #### DOES NOT WORK !!! ####
pe_test <- filterFeatures(pe, ~nNonZero >= 2)
nrow(pe[["peptideLog"]])



pe <- normalize(pe, 
                i = "peptideLog", 
                name = "peptideNorm", 
                method = "center.median")


pe[["peptideNorm"]] %>% 
  assay %>%
  as.data.frame() %>%
  gather(sample, intensity) %>% 
  mutate(condition = colData(pe)[sample,"condition"]) %>%
  ggplot(aes(x = intensity,group = sample,color = condition)) + 
  geom_density()


boxplot(assay(pe[["peptideNorm"]]),
        col = palette()[-1],
        main = "Peptide distribtutions after normalisation", ylab = "intensity"
)

pe[["peptideNorm"]] %>% 
  assay %>%
  limma::plotMDS(col = as.numeric(colData(pe)$condition))


#pe <- aggregateFeatures(pe,
                       # i = "peptideNorm",
                       # fcol = "Proteins",
                       # na.rm = TRUE,
                       # name = "pep")


pe <- msqrob(object = pe, i = "peptideNorm", formula = ~condition)

getCoef(rowData(pe[["peptideNorm"]])$msqrobModels[[1]])

# Explore the design of the model that we specified using the the package
library(ExploreModelMatrix)
VisualizeDesign(colData(pe),~condition)$plotlist[[1]]
### EXPLANALATION of THE FIGURE: 


L <- makeContrast(contc("(Intercept)",
                    "conditionB",
                    "conditionC",
                    "conditionD",
                    "conditionE"
                  )
)
pe <- hypothesisTest(object = pe, i = "peptideNorm", contrast = L, overwrite=T)

volcano <- ggplot(rowData(pe[["peptideNorm"]])$conditionA,
                  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)) +
  geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) + theme_minimal()
volcano




# Summarization
summaryPlot <- pe[["peptideNorm"]][
  rowData(pe[["peptideNorm"]])$Proteins == "P12081ups|SYHC_HUMAN_UPS",
  colData(pe)$lab=="lab2"&colData(pe)$condition %in% c("A","E")] %>%
  assay %>%
  as.data.frame %>%
  rownames_to_column(var = "peptide") %>%
  gather(sample, intensity, -peptide) %>% 
  mutate(condition = colData(pe)[sample,"condition"]) %>%
  ggplot(aes(x = peptide, y = intensity, color = sample, group = sample, label = condition), show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  geom_text(show.legend = FALSE) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlab("Peptide") + 
  ylab("Intensity (log2)")




