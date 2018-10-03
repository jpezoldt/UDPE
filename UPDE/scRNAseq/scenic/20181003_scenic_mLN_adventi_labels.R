#!/usr/bin/env Rscript

# Autor: Joern Pezoldt (code sections provided by Maria Litovchenko)
# 14.08.2018
# Function:
#1) QC single scRNA-seq dataset
#2) Perform SCENIC with GENIE3

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
#Load and label mLN
#####
#####
#Load Objects
#####
LN_minus_mLN <- readRDS("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/seurat/Output_NC_2018/mLN_FSC_PCA_24_nVarGene_1500_Res_1.3.rds")
#LN_minus_pLN <- readRDS("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/seurat/Output_NC_2018/pLN_FSC_PCA_18_nVarGene_1000_Res_0.9.rds")


#####
#Parameters Used
#####
#Seurat----------------------------------------------------
#mLN
mLN_VarGens <- 1000
mLN_Resolution <- 1.3
mLN_Dims <- 24
mLN_Param <- as.character(c(mLN_VarGens,mLN_Resolution,mLN_Dims))
#Scenic---------------------------------------------------
#Smallest expected cell population
smallest_pop <- 0.01
#Path store RDS
cell_type = "Adventi"
path_scenic_rds <- paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/mLN_SPF/", cell_type, "/int", sep="")
path_scenic_rds_minus_int <- paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/mLN_SPF/", cell_type, "/", sep="")
#"/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/mLN_pilot"
#/data/pezoldt/Analysis/scRNAseq/scenic

#####
#Generate sce-object from seurat object
#####
###mLN
#annotation file
head(LN_minus_mLN@meta.data)
head(LN_minus_mLN@ident)
annotation_mLN <- data.frame(cell_type1=LN_minus_mLN@ident, row.names = rownames(LN_minus_mLN@meta.data))
#readcount file with only the cells within the annotation file
count_matrix <- as.matrix(LN_minus_mLN@raw.data)
count_matrix_mLN <- count_matrix[,rownames(annotation_mLN)]
mLN_sce <- SingleCellExperiment(assays = list(counts = as.matrix(count_matrix_mLN)), colData = annotation_mLN)
TSNEPlot(object = LN_minus_mLN, pt.size = 0.3, do.label = TRUE)

#subset SingleCellExperimentObject
# mLN adventitial cells
adventi_cells <- rownames(subset(annotation_mLN, cell_type1 %in% c("Cd34+Has1+","Cd34+Gdf10+","Cd34+Aldh1a2+",
                                                          "Cd34+Ackr3+","Cd34+Cd248+")))
annotation_mLN_adventi <- subset(annotation_mLN, cell_type1 %in% c("Cd34+Has1+","Cd34+Gdf10+","Cd34+Aldh1a2+",
                                                                   "Cd34+Ackr3+","Cd34+Cd248+"))
count_matrix_adventi <- count_matrix_mLN[,adventi_cells]
mLN_adventi_sce <- SingleCellExperiment(assays = list(counts = as.matrix(count_matrix_adventi)), colData = annotation_mLN_adventi)

#####
#Extract RCM for cell subsets
#####
#Grab CD34+ cluster IDs
Adventi_clusters <- c("Cd34+Gdf10+",
                     "Cd34+Has1+","Cd34+Aldh1a2+","Cd34+Ackr3+",
                     "Cd34+Cd248+")

Adventi_cells <- names(LN_minus_mLN@ident[LN_minus_mLN@ident %in% Adventi_clusters])


#extract count matrix
exprMat_Adventi <- as.matrix(LN_minus_mLN@raw.data[,Adventi_cells])
saveRDS(exprMat_Adventi, file=paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/mLN_SPF/Adventi/int/exprMat_Adventi.Rds", sep=""))
expr_Adventi <- readRDS(paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/mLN_SPF/Adventi/int/exprMat_Adventi.Rds", sep=""))

##################
#SCENIC
##################
#####
#Settings
#####
#set storage path
setwd(path_scenic_rds)
#Set Rcistarget location
org="mgi"
dbDir="/data/software/R-3.5.0/library/RcisTarget/db" # RcisTarget databases location
myDatasetTitle="mLN-SPF non-Endothelial SC Adventi"
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, datasetTitle=myDatasetTitle, nCores=12) 
#Set to mm10
scenicOptions@settings$dbs <- setNames("mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather", "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
saveRDS(scenicOptions, file=paste(path_scenic_rds,"/scenicOptions.Rds", sep=""))

#####
#Store CellInformation
#####
cellInfo <- colData(mLN_adventi_sce)
cellInfo$nGene <- colSums(exprMat_Adventi>0)
cellInfo <- data.frame(cellInfo)
head(cellInfo)
saveRDS(cellInfo, file="int/cellInfo.Rds")

# Color to assign to the variables (same format as for NMF::aheatmap)
# Adventi
colVars <- list(CellType=setNames(c("brown", "orange", "green", "hotpink", "deepskyblue"), 
                                  c("Cd34+Gdf10+","Cd34+Has1+","Cd34+Aldh1a2+","Cd34+Ackr3+","Cd34+Cd248+")))

# nonAdventi
#colVars <- list(CellType=setNames(c("forestgreen", "darkorange", "magenta4", "hotpink", "deepskyblue"), 
 #                                c("Ccl19highMadcam1+","Inmt+Cxcl12+","Inmt+","Il6+Cxcl1","Ccl19+Il7+")))
saveRDS(colVars, file="int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))


#####
#QC
#####
#Summary tables
exprMat <- exprMat_Adventi
nCellsPerGene <- apply(exprMat, 1, function(x) sum(x>0))
nCountsPerGene <- apply(exprMat, 1, sum)
summary(nCellsPerGene)
summary(nCountsPerGene)
max(exprMat)
sum(exprMat>0) / sum(exprMat==0)

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
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered <- as.matrix(exprMat_filtered)
exprMat_filtered <- log2(exprMat_filtered+1)

#Delete matrix
#rm(exprMat)
exportsForGRNBoost(exprMat_filtered, scenicOptions)
#####
#Correlation analysis
#####
#Detection of positive and negative regulations possible thus split dataset
corrMat <- cor(t(as.matrix(exprMat_filtered)), method="spearman")
# (Only the rows for TFs will be needed needed):
#allTFs <- getDbTfs(scenicOptions)
#corrMat <- corrMat[which(rownames(corrMat) %in% allTFs),]
saveRDS(corrMat, file=getIntName(scenicOptions, "corrMat"))

#####
#load respective "scenicOptions" file
#####
scenicOptions <- readRDS(file=paste(path_scenic_rds,"/scenicOptions.Rds",sep=""))
exprMat <- readRDS(file=paste(path_scenic_rds,"/exprMat_",cell_type,".Rds",sep=""))

#####
#Run GRNBoost
#####
#Switch to Python
#Use -> GRNBoost_mLNSPF_500.py
#Observations: Multi-Threading critical, for high cell numbers (>1000) Ascersion Error

#####
#GRNBoost Output processing
#####
#Load GRN boost
#Load GRN boost
ex_network <- read.delim(paste(path_scenic_rds,"/GRNBoost_linklist.tsv",sep=""), header=FALSE)
#Export to use for scenic
#saveRDS(d_GRNBoost, file="int/1.4_GENIE3_linkList.Rds")
colnames(ex_network) <- c("TF","Target","weight")
saveRDS(ex_network, file=paste(path_scenic_rds,"/1.4_GENIE3_linkList.Rds", sep=""))

#####
#Perform AUC, Rcistarget
#####
setwd(path_scenic_rds_minus_int)

#scenicOptions <- readRDS(file=paste(path_scenic_rds,"/scenicOptions.Rds",sep=""))
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 6
scenicOptions@settings$seed <- 123

#Observation: Functions have problem with n(core) > 4 (observed 22.08.2018)
runSCENIC_1_coexNetwork2modules(scenicOptions)
#Note: Warning: In runSCENIC_1_coexNetwork2modules(scenicOptions) :
#         The following TFs are missing from the correlation matrix: Adarb1, Bcl6b, Carf, Erg, Fli1, Foxc2, Gata2, Hey1, Hhex, Hoxb5, Hoxb8, Hoxd4, Hoxd8, Irx3, Lyl1, Mecom, Prox1, Sox17, Sox7, Tal1, Tbx1, Tbx3, Tcf15, Zfp407
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat)

# Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=c(10,20,30), perpl=c(10,20,30))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=c(10,20,30), perpl=c(10,20,30), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# Plot as pdf (individual files in int/): 
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE)

par(mfcol=c(2,5))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, cex=.5, varName="cell_type1")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

# Using only "high-confidence" regulons (normally similar)
par(mfcol=c(2,5))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, cex=.5) # varName="cell_type1")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

#chosen t-SNE can then be saved as default to use for plots
scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 20
scenicOptions@settings$defaultTsne$perpl <- 10
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#####
#TSNE Plots
#####
par(mfcol=c(3,5))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, cex=.5, varName="cell_type1")
plot.new(); 
legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))






