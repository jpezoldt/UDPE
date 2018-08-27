#!/usr/bin/env Rscript

# Autor: Joern Pezoldt (code sections provided by Maria Litovchenko)
# 26.08.2018
# Function:
#1) Run SCENIC with GRNBoost Outpu

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
#Global variable
#####
sample_ID <- "pLN_SPF"
#Smallest expected cell population
smallest_pop <- 0.01
#Paths
path_scenicOptions <- paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/",sample_ID,sep="")
path_SO_int <- paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/",sample_ID, "/int/",sep="")
path_SO_out <- paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/",sample_ID, "/out/",sep="")


#####
#GRNBoost Output processing
#####
#Load GRN boost
GRNBoost_link_list <- read.delim(paste(path_SO_int,"GRNBoost_linklist.tsv",sep=""), header=FALSE)
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
colnames(GRNBoost_link_list) <- c("TF","Target","weight")
saveRDS(GRNBoost_link_list, file=paste(path_SO_int,"1.4_GENIE3_linkList.Rds",sep=""))

#####
#Perform AUC, Rcistarget
#####
#Set Rcistarget location
org="mgi"
dbDir="/data/software/R-3.5.0/library/RcisTarget/db" # RcisTarget databases location
myDatasetTitle=paste(sample_ID, "non-endothelial SC")
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, datasetTitle=myDatasetTitle, nCores=12) 
#Set to mm10
scenicOptions@settings$dbs <- setNames("mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather", "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 4
scenicOptions@settings$seed <- 123

#Observation: Functions have problem with n(core) > 4 (observed 22.08.2018)
runSCENIC_1_coexNetwork2modules(scenicOptions)
#Note: Warning: In runSCENIC_1_coexNetwork2modules(scenicOptions) :
#         The following TFs are missing from the correlation matrix: Adarb1, Bcl6b, Carf, Erg, Fli1, Foxc2, Gata2, Hey1, Hhex, Hoxb5, Hoxb8, Hoxd4, Hoxd8, Irx3, Lyl1, Mecom, Prox1, Sox17, Sox7, Tal1, Tbx1, Tbx3, Tcf15, Zfp407
runSCENIC_2_createRegulons(scenicOptions)

#Load log normalized filtered expression matrix
exprMat <- readRDS(paste(path_SO_int, "exprMat_filtered_log.Rds",sep=""))
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
colVars <- readRDS(paste(path_SO_int,"colVars.Rds",sep=""))
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))


#chosen t-SNE can then be saved as default to use for plots
scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 20
scenicOptions@settings$defaultTsne$perpl <- 20
saveRDS(scenicOptions, file="int/scenicOptions.Rds")











#Optional: Binarize AUC 
#Note: Needs https://cran.r-project.org/web/packages/rbokeh/index.html
#logMat <- exprMat_filtered # Better if it is logged/normalized
#aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
#savedSelections <- shiny::runApp(aucellApp)

# Save the modified thresholds:
#newThresholds <- savedSelections$thresholds
#scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
#saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
#saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
# scenicOptions@settings$devType="png"
#runSCENIC_4_aucell_binarize(scenicOptions)

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
