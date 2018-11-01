# dissect compiled and curated Output from Homer
# Author: Joern Pezoldt
# Date: 12.10.2018
# Type of script: Exploratory

########### Libraries
library(pheatmap)



############ PATH & Global Variables ----------------------------

sample_ID <- "SPF"
# path to binary matrix of TFs detected using homer via motifs identified
PATH_input_TF_homer_matrix <- paste("/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/Motif/Homer_Output/allHomerResults/size_500bp/", sample_ID, sep = "")

# path to regulons identified
cell_subsets_1 <- "Adventi"
cell_subsets_2 <- "nonAdventi"

PATH_input_regulons_scenic <- "/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/mLN_SPF"

############ Load tables ----------------------------------------
#scenic
cell_subsets_regulon_1 <- read.table(paste(PATH_input_regulons_scenic,"/",cell_subsets_1,"/AUCell_scores_core_",cell_subsets_1,".csv", sep = ""), dec = ".", sep = ",")
cell_subsets_regulon_2 <- read.table(paste(PATH_input_regulons_scenic,"/",cell_subsets_2,"/AUCell_scores_core_",cell_subsets_2,".csv", sep = ""), dec = ".", sep = ",")

#motif TFs ATAC
TFs_ATAC <- read.table(paste(PATH_input_TF_homer_matrix,"/",sample_ID,"_existence_matrix_known_homer_compile.txt", sep = ""), dec = ".", sep = "\t")

############ Extract TFs of interest ----------------------------
Unique_TFs_NoDEG_Open <- subset(TFs_ATAC, 
                                  Open_DEG_SPF_lowRNAseq == 0 &  
                                  Closed_DEG_SPF_lowRNAseq == 0 &
                                  No_DAR_UP_SPF_lowRNAseq == 0 &
                                  No_DAR_DOWN_SPF_lowRNAseq == 0 &
                                  NoDEG_Closed_SPF_lowRNAseq == 0 &
                                  NoDEG_Open_SPF_lowRNAseq == 1)

############ Overlap and heatmap -------------------------------
cell_subsets_regulon_TFs_1 <- subset(cell_subsets_regulon_1, rownames(Unique_TFs_NoDEG_Open) %in% GeneSymbol)
cell_subsets_regulon_TFs_2 <- subset(cell_subsets_regulon_2, rownames(Unique_TFs_NoDEG_Open) %in% GeneSymbol)

pheatmap(as.matrix(cell_subsets_regulon_TFs_1[,4:ncol(cell_subsets_regulon_TFs_1)]),
         scale = "none")

# Continue Here:
# 1) add cell identification system
