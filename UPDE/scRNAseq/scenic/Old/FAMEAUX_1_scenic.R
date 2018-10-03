#!/usr/bin/env Rscript

# Autor: Joern Pezoldt (code sections provided by Maria Litovchenko)
# 26.08.2018
# Function:
#1) QC single scRNA-seq dataset
#2) Prep for GRNBoost Input

#Libraries
library(biomaRt)
library(data.table)
library(AUCell)
library(GENIE3)
library(RcisTarget)
library(SCENIC)
library(Seurat)
library(cowplot)
library(dplyr)
library(SingleCellExperiment)

#####
#Global variables
#####
sample_ID <- "mLN_SPF"
path_scenicOptions <- paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/",sample_ID,sep="")
path_SO_int <- paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/",sample_ID, "/int/",sep="")
path_SO_out <- paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/",sample_ID, "/out/",sep="")
#scenic params
smallest_pop <- 0.005

#####
#Load Objects
#####
#set path to seurat analyzed sample .rds
seurat_analyzed <- readRDS("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/seurat/Output_NC_2018/mLN_FSC_PCA_24_nVarGene_1500_Res_1.3.rds")
#seurat_analyzed <- readRDS("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/seurat/Output_NC_2018/pLN_FSC_PCA_18_nVarGene_1000_Res_0.9.rds")


#####
#Generate sce-object from seurat object
#####
###mLN
#annotation file
#head(LN_minus_mLN@meta.data)
#head(LN_minus_mLN@ident)
#annotation_mLN <- data.frame(cell_type1=LN_minus_mLN@ident, row.names = rownames(LN_minus_mLN@meta.data))
#readcount file with only the cells within the annotation file
#count_matrix <- as.matrix(LN_minus_mLN@raw.data)
#count_matrix_mLN <- count_matrix[,rownames(annotation_mLN)]
#mLN_sce <- SingleCellExperiment(assays = list(counts = as.matrix(count_matrix_mLN)), colData = annotation_mLN)

###pLN
head(seurat_analyzed@meta.data)
head(seurat_analyzed@ident)
annotation_sce <- data.frame(cell_type1=seurat_analyzed@ident, row.names = rownames(seurat_analyzed@meta.data))
#readcount file with only the cells within the annotation file
count_matrix <- as.matrix(seurat_analyzed@raw.data)
count_matrix_sce <- count_matrix[,rownames(annotation_sce)]
sce_seurat_derived <- SingleCellExperiment(assays = list(counts = as.matrix(count_matrix_sce)),
                                colData = annotation_sce)

#####
#Sample cells for GRNBoost/GENIE3
#####
#Sample number of SC cells
#set.seed(123)
#sampled_SCs <- sample(sample_1_NO_LECs_BECs, 5000)

#Filter for mitochondiral and duplets and LECs and BECs out
#sample_1_seurat_minus <- FilterCells(object = sample_1_seurat, subset.names = c("nGene", "percent.mito"), 
 #                                    low.thresholds = c(min_gene_number, -Inf), high.thresholds = c(nGene, mito_cutoff),
  #                                   cells.use = sample_1_NO_LECs_BECs)
#sample_1_seurat_minus <- FilterCells(object = sample_1_seurat_minus, subset.names = c("Ackr4","Pecam1"),
 #                                    low.thresholds = c(-Inf, -Inf), high.thresholds = c(1,1))
#sampled_cells <- rownames(sample_1_seurat_minus@meta.data)

#extract count matrix
exprMat <- counts(sce_seurat_derived)
dim(exprMat)
saveRDS(exprMat, paste(path_SO_int, "/exprMat.Rds",sep=""))
#Note: Use Seurat output for SCENIC

##################
#SCENIC
##################
#####
#Settings
#####
#set storage path

#Set Rcistarget location
org="mgi"
dbDir="/data/software/R-3.5.0/library/RcisTarget/db" # RcisTarget databases location
myDatasetTitle=paste(sample_ID, "non-endothelial SC")
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, datasetTitle=myDatasetTitle, nCores=12) 
#Set to mm10
scenicOptions@settings$dbs <- setNames("mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather", "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
saveRDS(path_SO_int, file=paste(path_SO_int,"scenicOptions.Rds",sep=""))

#####
#Store CellInformation
#####
cellInfo <- colData(sce_seurat_derived)
cellInfo$nGene <- colSums(exprMat>0)
cellInfo <- data.frame(cellInfo)
head(cellInfo)
saveRDS(cellInfo, paste(path_SO_int,"cellInfo.Rds",sep=""))

# Color to assign to the variables (same format as for NMF::aheatmap)
total_cell_types <- levels(cellInfo$cell_type1)
#Make a vector of colors to match to cell types
colors <- palette(rainbow(20))
set.seed(123)
colors_sample <- sample(colors, length(total_cell_types))
colVars <- list(CellType=setNames(colors_sample, 
                                  total_cell_types))
saveRDS(colVars, paste(path_SO_int,"colVars.Rds",sep=""))
#plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

#####
#QC
#####
#Summary tables
nCellsPerGene <- apply(exprMat, 1, function(x) sum(x>0))
nCountsPerGene <- apply(exprMat, 1, sum)
summary(nCellsPerGene)
summary(nCountsPerGene)
max(exprMat)
sum(exprMat>0) / sum(exprMat==0)

# 1) Eliminate low count cells
#Keep only cells that have sufficient UMIs
nUMIsPerCell <- apply(exprMat, 2, sum)
nGenesPerCell <- apply(exprMat, 2, function(x) sum(x>0))
hist(nGenesPerCell)
nGenesPerCell_Left <- names(nGenesPerCell)[which(nGenesPerCell > 1500)]
#Select the cells to keep (Threshold see above)
exprMat <- exprMat[,nGenesPerCell_Left]

# 2) Thresh for sufficiently represented genes
#Keep only cells that have sufficient UMIs
nUMIsPerCell <- apply(exprMat, 2, sum)
nGenesPerCell <- apply(exprMat, 2, function(x) sum(x>0))
hist(nGenesPerCell)
nGenesPerCell_Left <- names(nGenesPerCell)[which(nGenesPerCell > 1500)]

#Keep only genes that have sufficient read numbers
# Expressed with 3 UMIs in 1% of the cells
minReads <- 3*smallest_pop*ncol(exprMat)
genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > minReads)]
genesTotal <- names(nCountsPerGene)
length(genesTotal)
length(genesLeft_minReads)

#Make sure that genes relevant for small populations are not excluded (e.g. Perycytes)
# Currently set to 1%
minSamples <- ncol(exprMat)*smallest_pop
nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > minSamples)]
length(genesLeft_minCells)

#Only genes in RcisTarget will be used
# Check overlap
motifRankings <- importRankings(getDatabases(scenicOptions)[[1]]) # either one, they should have the same genes
genesInDatabase <- colnames(getRanking(motifRankings))
genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]
length(genesLeft_minCells_inDatabases)

# Check whether any relevant gene / potential gene of interest is missing:
interestingGenes <- c("Ccl9", "Nkx2-3", "Pdpn", "Il7")
interestingGenes[which(!interestingGenes %in% genesLeft_minCells_inDatabases)]

#Filter ExpressionMatrix for the genes present in database
genesKept <- genesLeft_minCells_inDatabases
#Select the genes to keep (Threshold see above)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered <- as.matrix(exprMat_filtered)
exprMat_filtered <- log2(exprMat_filtered+1)
saveRDS(exprMat_filtered, paste(path_SO_int, "exprMat_filtered_log.Rds",sep=""))

#Delete matrix
#rm(exprMat)
setwd(path_scenicOptions)
exportsForGRNBoost(exprMat_filtered, scenicOptions)
#####
#Correlation analysis
#####
#Detection of positive and negative regulations possible thus split dataset
corrMat <- cor(t(as.matrix(exprMat_filtered)), method="spearman")
# (Only the rows for TFs will be needed needed):
#allTFs <- getDbTfs(scenicOptions)
#corrMat <- corrMat[which(rownames(corrMat) %in% allTFs),]
saveRDS(corrMat, paste(path_SO_int,"1.2_corrMat.Rds",sep=""))

#####
#Run GRNBoost
#####
#Switch to Python
#Use -> GRNBoost_mLNSPF_500.py
#Observations: Multi-Threading critical, for high cell numbers (>1000) Ascersion Error

