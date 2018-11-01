# Identify Clusters which express TF signatures
# Author: Joern Pezoldt
# Date: 31.10.2018
# Type of script: Exploratory


# Library
library("Seurat")
library("cowplot")
library("dplyr")
library("pheatmap")


#####
#PATHs & Global Variables
#####
sample_1 <- "mLN"
sample_2 <- "pLN"
PATH_input_scRNAseq <- "/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/seurat/LN_Commnesals/mLN_pLN_SPF"
PATH_input_TFs <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/Motif/Homer_Output/allHomerResults/size_500bp/SPF/SPF_existence_matrix_known_homer_compile.txt"
PATH_output <- "/home/pezoldt/NAS2/pezoldt/Analysis/01_Integrated/scRNA_TFs"

#####
#Load Data
#####
LN_minus <- readRDS(paste(PATH_input_scRNAseq,"/LN_minus_mLN_pLN_SPF.Rds",sep = ""))
#DEGs across all clusters
LN_minus.markers <- readRDS(paste(PATH_input_scRNAseq,"/mLN_pLN_SPF_DEG_all.Rds",sep = ""))
#DEGs calculated per cluster
mLN_pLN_DEGs_per_cluster_tobit <- read.delim(paste(PATH_input_scRNAseq,"/mLN_pLN_DEGs_per_cluster_tobit.csv",sep = ""),
                                             dec = ".", sep = ",")
#TFs binary per setting
TFs_binary <- read.delim(PATH_input_TFs)

#####
#Compare TF across subsets
#####
#Select TF list
TFs_of_interest <- rownames(subset(TFs_binary, 
                                    Closed_DEG_SPF_lowRNAseq == 0 &
                                    No_DAR_DOWN_SPF_lowRNAseq == 1 &
                                    No_DAR_UP_SPF_lowRNAseq == 0 &
                                    NoDEG_Closed_SPF_lowRNAseq == 0 &
                                    NoDEG_Open_SPF_lowRNAseq == 0 &
                                    Open_DEG_SPF_lowRNAseq  == 0))

l_DEGs_per_cluster <- split(mLN_pLN_DEGs_per_cluster_tobit, mLN_pLN_DEGs_per_cluster_tobit$cluster_i)

#Overlap of TFs with clusters
TFs_in_cluster <- lapply(l_DEGs_per_cluster, function(x){
  overlap_i <- subset(x, GeneSymbol %in% TFs_of_interest)
})


