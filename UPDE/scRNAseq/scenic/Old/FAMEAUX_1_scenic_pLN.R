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
LN_minus_pLN <- readRDS("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/seurat/Output_NC_2018/pLN_FSC_PCA_18_nVarGene_1000_Res_0.9.rds")

#####
#Parameters Used
#####
#Seurat----------------------------------------------------
#mLN
mLN_VarGens <- 1000
mLN_Resolution <- 1.3
mLN_Dims <- 24
mLN_Param <- as.character(c(mLN_VarGens,mLN_Resolution,mLN_Dims))
#pLN
pLN_VarGens <- 1500
pLN_Resolution <- 0.9
pLN_Dims <- 18
pLN_Param <- as.character(c(pLN_VarGens,pLN_Resolution,pLN_Dims))

#Scenic---------------------------------------------------
#Smallest expected cell population
smallest_pop <- 0.01
#Path store RDS
path_scenic_rds <- "/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/pLN_SPF"


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

###pLN
head(LN_minus_pLN@meta.data)
head(LN_minus_pLN@ident)
annotation_pLN <- data.frame(cell_type1=LN_minus_pLN@ident, row.names = rownames(LN_minus_pLN@meta.data))
#readcount file with only the cells within the annotation file
count_matrix <- as.matrix(LN_minus_pLN@raw.data)
count_matrix_pLN <- count_matrix[,rownames(annotation_pLN)]
pLN_sce <- SingleCellExperiment(assays = list(counts = as.matrix(count_matrix_pLN)),
                                colData = annotation_pLN)


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
exprMat <- counts(pLN_sce)
dim(exprMat)
saveRDS(exprMat_filtered, paste(path_scenic_rds, "/int/exprMat.Rds",sep=""))
#Note: Use Seurat output for SCENIC

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
myDatasetTitle="pLN-SPF non-Endothelial SC"
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, datasetTitle=myDatasetTitle, nCores=12) 
#Set to mm10
scenicOptions@settings$dbs <- setNames("mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather", "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#####
#Store CellInformation
#####
#!!!Note: dataset dependent
cellInfo <- colData(pLN_sce)
cellInfo$nGene <- colSums(exprMat>0)
cellInfo <- data.frame(cellInfo)
head(cellInfo)
saveRDS(cellInfo, file="int/cellInfo.Rds")

# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(CellType=setNames(c("forestgreen", "darkorange", "magenta4", "hotpink", "red3", "skyblue", "darkblue","red","blue","black"), 
                                  c("Il6+Cxcl1+","Ccl19+","Ccl19+Il7+","Inmt+","Inmt+Gdf10+","pSC","Ccl19high","Cd34+","PvC","Cxcl9+")))
saveRDS(colVars, file="int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

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
saveRDS(exprMat_filtered, paste(path_scenic_rds, "/int/exprMat_filtered_log.Rds",sep=""))

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
#GRNBoost Output processing
#####
#Load GRN boost
ex_5000SC_network <- read.delim("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/pLN_SPF/int/GRNBoost_NC_2018_SC_network.tsv", header=FALSE)
#Note:
# Column #1 is the TFs
# Column #2 is the Genes
#Build Matrix
#Split according to TFsl
#l_network <- split(ex_5000SC_network, ex_5000SC_network$V1)
#Note:
# No score for every TF <-> Gene couple
#head(ex_5000SC_network)
#tail(ex_5000SC_network)
# Fill gaps with zeros
#Generate with:
# Colnames = Genes/V2
# Rowname = TF/V1
# Scores = Content/V3
#head((data.frame(l_network[[1]]$V3)))
#head(t(as.matrix(l_network[[1]]$V3, colnames = l_network[[1]]$V2)))

#For each TF generate a row with the genes as colnames
#l_of_rows_network <- lapply(l_network, function(x){
#  v_TF_i <- x$V3
#  names(v_TF_i) <- x$V2
#  d_TF_i <- as.data.frame(v_TF_i)
#  colnames(d_TF_i) <- x$V1[1]
#  as.data.frame(t(d_TF_i))
#})

#Grab TF names in the order of the rows (see above)
#l_TF_names <- lapply(l_of_rows_network, function(x){
#  TF_name_i <- rownames(x)
#  TF_name_i
#})

#make table
#d_GRNBoost <- rbindlist(l_of_rows_network, fill = TRUE)
#rownames(d_GRNBoost) <- l_TF_names
#replace NAs with Zeros
#d_GRNBoost[is.na(d_GRNBoost)] <- 0

#Export to use for scenic
#saveRDS(d_GRNBoost, file="int/1.4_GENIE3_linkList.Rds")
colnames(ex_5000SC_network) <- c("TF","Target","weight")
saveRDS(ex_5000SC_network, file="/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/pLN_SPF/int/1.4_GENIE3_linkList.Rds")

#####
#Perform AUC, Rcistarget
#####
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 4
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
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, cex=.5, varName="cell_type1")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

#chosen t-SNE can then be saved as default to use for plots
scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 20
scenicOptions@settings$defaultTsne$perpl <- 10
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#Optional: Binarize AUC 
#Note: Needs https://cran.r-project.org/web/packages/rbokeh/index.html
logMat <- exprMat_filtered # Better if it is logged/normalized
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
savedSelections <- shiny::runApp(aucellApp)

# Save the modified thresholds:
#newThresholds <- savedSelections$thresholds
#scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
#saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
#saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
# scenicOptions@settings$devType="png"
runSCENIC_4_aucell_binarize(scenicOptions)

#Explore datasets
#Note: package rbokeh needs to be installed
#library(rbokeh)
#logMat <- exprMat # Better if it is logged/normalized
#aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat) #default t-SNE
#savedSelections <- shiny::runApp(aucellApp)
#tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
#aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
#
# Show TF expression:
#par(mfrow=c(1,1))
#Note: not functional
#AUCell::AUCell_plotTSNE(tSNE_scenic$Y,
#                        exprMat, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Nkx2-3", "Cd34", "Nfkb1", "Maf")],],
#                        plots="Expression")

#Density plot to detect most likely stable states
#library(KernSmooth)
#library(RColorBrewer)
#dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
#image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
#contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)

#Show several regulons simultaneously
#par(mfrow=c(1,2))
#Note: only works with binarized data
#regulonNames <- c( "Dlx1","Sox9")
#cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)
#text(-30,-25, attr(cellCol,"red"), col="red", cex=.7, pos=4)
#text(-30,-25-4, attr(cellCol,"green"), col="green3", cex=.7, pos=4)

#regulonNames <- list(red=c("Sox10", "Sox8"),
 #                    green=c("Maf"),
  #                   blue=c( "Neurod2"))
#cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC")
#text(-30,-25, attr(cellCol,"red"), col="red", cex=.7, pos=4)
#text(-30,-25-4, attr(cellCol,"green"), col="green3", cex=.7, pos=4)
#text(-30,-25-8, attr(cellCol,"blue"), col="blue", cex=.7, pos=4)

#####
#Regulon targets
#####
#Only regulons with 10 genes or more are scored
#regulons <- loadInt(scenicOptions, "regulons")
#regulons[c("Arnt", "Atf7","Irf2","Ctcf","Foxp2","Nfkb1","Nfkb2")]

#Analyze Links of TF to genes and Motifs
#regulons <- loadInt(scenicOptions, "aucell_regulons")
#head(cbind(onlyNonDuplicatedExtended(names(regulons))))
#regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
#regulonTargetsInfo <- RcisTarget::addLogo(regulonTargetsInfo, motifCol="bestMotif")
#regulonTargetsInfo$Genie3Weight <- signif(regulonTargetsInfo$Genie3Weight, 2)
#colsToShow <- c("TF", "gene", "nMotifs", "bestMotif", "logo", "NES", "highConfAnnot", "Genie3Weight")
#DT::datatable(regulonTargetsInfo[TF=="Foxp2" & highConfAnnot==TRUE, colsToShow, with=F], escape=FALSE, filter="top")

#Generate Output files
# Note: Needs to be adapted
#motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
#motifEnrichment_selfMotifs_wGenes <- RcisTarget::addLogo(motifEnrichment_selfMotifs_wGenes)
#colsToShow <- c("motifDb", "logo", "NES", "geneSet", "TF_highConf") # "TF_lowConf", "enrichedGenes"
#DT::datatable(motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Dlx5", colsToShow, with=F], escape=FALSE, filter="top")

#Export Loom dataset
# DGEM (Digital gene expression matrix)
# (non-normalized counts)
# Make an SCE from data
#library(SingleCellExperiment)
#load("data/sceMouseBrain.RData")
#dgem <- counts(sceMouseBrain)
#head(colnames(dgem))  #should contain the Cell ID/name
# Export:
#scenicOptions@fileNames$output["loomFile",] <- "mouseBrain.loom"
#export2scope(scenicOptions, dgem)
