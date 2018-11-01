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
#Load Objects
#####
LN_minus_mLN <- readRDS("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/seurat/Output_NC_2018/mLN_FSC_PCA_24_nVarGene_1500_Res_1.3.rds")
#LN_minus_pLN <- readRDS("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/seurat/Output_NC_2018/pLN_FSC_PCA_18_nVarGene_1000_Res_0.9.rds")


#####
#Parameters Used
#####
#Scenic---------------------------------------------------
#Smallest expected cell population
smallest_pop <- 0.01
#Path store RDS
cell_type = "nonAdventi"
path_scenic_rds <- paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/mLN_SPF/", cell_type, "/int", sep="")
path_scenic_rds_minus_int <- paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/mLN_SPF/", cell_type, "/", sep="")


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
clusters_of_interest <- c("Ccl19highMadcam1+","Inmt+Cxcl12+","Inmt+",
                          "Il6+Cxcl1","Ccl19+Il7+")
#coloring for clusters of interest
clusters_coloring <- c("forestgreen", "darkorange", "magenta4", "hotpink", "deepskyblue")
cell_to_scenic <- rownames(subset(annotation_mLN, cell_type1 %in% clusters_of_interest))
annotation_subset <- subset(annotation_mLN, cell_type1 %in% clusters_of_interest)
count_matrix_sce <- count_matrix_mLN[,cell_to_scenic]
mLN_subsets_sce <- SingleCellExperiment(assays = list(counts = as.matrix(count_matrix_sce)), colData = annotation_subset)

#####
#Extract RCM for cell subsets
#####
#Grab cluster IDs of interest
cell_to_scenic <- names(LN_minus_mLN@ident[LN_minus_mLN@ident %in% clusters_of_interest])


#extract count matrix
exprMat <- as.matrix(LN_minus_mLN@raw.data[,cell_to_scenic])
saveRDS(count_matrix_sce, file=paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/mLN_SPF/Adventi/int/exprMat_", cell_type, ".Rds", sep = ""))
exprMat <- readRDS(paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/mLN_SPF/Adventi/int/exprMat_Adventi.Rds", sep=""))

##################
#SCENIC
##################
#####
#Settings
#####
#set storage path
setwd(path_scenic_rds_minus_int)
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
cellInfo <- colData(mLN_subsets_sce)
cellInfo$nGene <- colSums(count_matrix_sce>0)
cellInfo <- data.frame(cellInfo)
head(cellInfo)
saveRDS(cellInfo, file="int/cellInfo.Rds")

# nonAdventi
colVars <- list(CellType=setNames(clusters_coloring,
                                  clusters_of_interest))
saveRDS(colVars, file="int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))


#####
#QC
#####
#Summary tables
exprMat <- count_matrix_sce
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
motifRankings <- importRankings(getDatabases(scenicOptions)[[1]])
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
#Run GRNBoost
#####
#Switch to Python
#Use -> GRNBoost_mLNSPF_500.py
#Observations: Multi-Threading critical, for high cell numbers (>1000) Ascersion Error

#####
#load respective "scenicOptions" file
#####
scenicOptions <- readRDS(file=paste(path_scenic_rds,"/scenicOptions.Rds",sep=""))
exprMat <- readRDS(file=paste(path_scenic_rds,"/exprMat_",cell_type,".Rds",sep=""))

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

#Observation: Functions have problem with n(core) > 6
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

par(mfcol=c(3,5))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, cex=.5, varName="cell_type1")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

# Using only "high-confidence" regulons (normally similar)
par(mfcol=c(3,5))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, cex=.5, varName="cell_type1")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

#chosen t-SNE can then be saved as default to use for plots
scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 20
scenicOptions@settings$defaultTsne$perpl <- 10
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#################
#Explore dataset
#################
loadInt(scenicOptions)
logMat <- log2(exprMat+1) # Better if it is logged/normalized
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat) #default t-SNE
savedSelections <- shiny::runApp(aucellApp)
##

#####
#Regulon t-SNE
#####
par(mfrow=c(2,4))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat,
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Irf1","Arntl","Fos","Fosl1",
                                                                                                   "Rel","Junb","Jun","Egr1")],],
                        plots="Expression")

#####
#Density plot to detect most likely stable states
#####
library(KernSmooth)
library(RColorBrewer)
dens2d <- bkde2D(tSNE_scenic$Y, 4)$fhat
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=3, drawlabels=FALSE)
c("Hoxa5","Arnt","Ets1","Tcf3","Fosl2","Fosl1","Rela","Rel","Nfkb1","Bach1","Irf1","Arntl","Fos","Junb","Jun","Elk3","Egr1")


#####
#Genes included in a regulon
#####
regulons <- loadInt(scenicOptions, "regulons")
regulons[c("Irf1","Arntl","Fos","Fosl1",
           "Rel","Junb","Jun","Egr1")]

regulons <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))
tail(cbind(onlyNonDuplicatedExtended(names(regulons))))

#####
#Check genes included in regulon
#####
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
regulonTargetsInfo <- RcisTarget::addLogo(regulonTargetsInfo, motifCol="bestMotif")
regulonTargetsInfo$Genie3Weight <- signif(regulonTargetsInfo$Genie3Weight, 2)
regulonTargetsInfo_highConfAnnot <- subset(regulonTargetsInfo, highConfAnnot == "TRUE")

#Grab regulons per cells
#check output files
loadInt(scenicOptions)
motifEnrichment_full <- loadInt(scenicOptions, "motifEnrichment_full")    


#regulons actually detected with associated genes
regulons <- loadInt(scenicOptions, "regulons")

#Genes per regulon
aucell_regulons <- loadInt(scenicOptions, "aucell_regulons")

#binary matrix for which gene belongs to which regulon
regulon_incidence_matrix_binary <- loadInt(scenicOptions, "regulons_incidMat")
regulon_incidence_matrix_binary[1:5,1:5]

#Names of genes
aucell_rankings <- loadInt(scenicOptions, "aucell_rankings")

#Names of regulons
#This table contains the non-binarized scores
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

#thresholds for all regulons
aucell_thresholds <- loadInt(scenicOptions, "aucell_thresholds")

#Only if data has been binarized
#aucell_binary_full <- loadInt(scenicOptions, "aucell_binary_full")

