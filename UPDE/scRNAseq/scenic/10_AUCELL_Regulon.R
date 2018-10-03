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

#####
#PATHs
#####
#Input required:
cell_type = "nonAdventi"
path_scenic_rds <- paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/mLN_SPF/", cell_type, "/int", sep="")
setwd(paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/mLN_SPF/", cell_type, sep=""))


#####
#load respective "scenicOptions" file
#####
scenicOptions <- readRDS(file=paste(path_scenic_rds,"/scenicOptions.Rds",sep=""))
exprMat <- readRDS(file=paste(path_scenic_rds,"/exprMat_",cell_type,".Rds",sep=""))

#####
#GRNBoost Output processing
#####
#Load GRN boost
ex_network <- read.delim(paste(path_scenic_rds,"/GRNBoost_linklist.tsv",sep=""), header=FALSE)
#Export to use for scenic
#saveRDS(d_GRNBoost, file="int/1.4_GENIE3_linkList.Rds")
colnames(ex_network) <- c("TF","Target","weight")
saveRDS(ex_network, file=paste(path_scenic_rds,"/1.4_GENIE3_linkList.Rds", sep=""))
#saveRDS(ex_network, file=getIntName(scenicOptions,"int/1.4_GENIE3_linkList.Rds"))

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

par(mfcol=c(3,3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, cex=.5)
#varName="CellType", 

# Using only "high-confidence" regulons (normally similar)
par(mfcol=c(3,3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, cex=.5)
#varName="CellType", 
