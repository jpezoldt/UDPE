#!/usr/bin/env Rscript

# Autor: Joern Pezoldt (code sections provided by Maria Litovchenko)
# 14.10.2018
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
library(rPython)

#####
#Load Seurat Objects
#####
PATH_seurat <- "/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/seurat/Output_NC_2018/mpLN_FSC_PCA_22_nVarGene_1500_Res_1.1.rds"

#####
#Global variables and PATHs
#####
#Scenic---------------------------------------------------
#Smallest expected cell population
smallest_pop <- 0.01
#Path store RDS
cell_type = "mSC"
organ = "pLN_mLN_SPF"

# PATH
# Note: The directory saved at "path_user_scripts" needs to be set for each new environment to enable access to python script
#       The script directory needs also be adapted in the python script
path_user_scripts <- "/home/pezoldt/NAS2/pezoldt/Scripts/scRNAseq/scenic/"
path_user_analysis <- "/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/"

# Generate PATHs for downsstream analsis
path_scenic_sample <- paste(path_user_analysis,organ,"/",cell_type,sep="")
path_scenic_sample_int <- paste(path_scenic_sample, "/int",sep="")
#create directory for storage
dir.create(path_scenic_sample_int, recursive = TRUE)

# Save data to python accessible file
python_access <- c(cell_type,organ,path_user_analysis)
write.table(python_access, paste(path_user_scripts,"name_space_for_python.txt", sep=""),sep="\t", row.names = FALSE)

#Setup SCENIC object
#organism
org="mgi"
#Loaction of database
dbDir="/data/software/R-3.5.0/library/RcisTarget/db" # RcisTarget databases location
#Name of dataset
myDatasetTitle=paste(organ, cell_type, sep=" ")
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, datasetTitle=myDatasetTitle, nCores=12) 
#Set to mm10
scenicOptions@settings$dbs <- setNames("mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather", "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
saveRDS(scenicOptions, file=paste(path_scenic_sample_int,"/scenicOptions.Rds", sep=""))

# Names defined by SEURAT for sce
# Note: Input required
# Print clusters available in Seurat object
seurat_cluster_check <- readRDS(PATH_seurat)
levels(seurat_cluster_check@ident)
clusters_of_interest <- c("mSC")

#coloring for clusters of interest
# Note: Input required
#clusters_coloring <- c("forestgreen")
#"red","black")
#,"yellow")


#########################
#Functions
#########################
#' makeSCEfromSEURAT
#' generates SCE object from SEURAT object including the annotation
#' @param PATH_seurat path to the file
#' @param clusters_of_interest character vector that contains the names of the clusters predefined by seurat
#' @return count_matrix object containt only the cells that align with the predefined annotation of in cluster_of_interest 

makeSCEfromSEURAT <- function(PATH_seurat, clusters_of_interest){
  
  #read seurat object
  object_seurat <- readRDS(PATH_seurat)
  #extract annotation
  annotation_seurat <- data.frame(cell_type1=object_seurat@ident, row.names = rownames(object_seurat@meta.data))
  
  #readcount file with only the cells within the annotation file
  count_matrix <- as.matrix(object_seurat@raw.data)
  count_matrix_seurat <- count_matrix[,rownames(annotation_seurat)]
  sce <- SingleCellExperiment(assays = list(counts = as.matrix(count_matrix_seurat)), colData = annotation_seurat)
  
  #subset sce object for the clusters of interest
  cell_to_scenic <- rownames(subset(annotation_seurat, cell_type1 %in% clusters_of_interest))
  annotation_subset <- subset(annotation_seurat, cell_type1 %in% clusters_of_interest)
  count_matrix_sce <- count_matrix[,cell_to_scenic]
  sce_from_seurat <- SingleCellExperiment(assays = list(counts = as.matrix(count_matrix_sce)), colData = annotation_subset)
  
  #Save color code for cells obtained from seurat distinguishing mLN and pLN cells
  cellInfo <- colData(sce_from_seurat)
  cellInfo$nGene <- colSums(count_matrix_sce>0)
  cellInfo <- data.frame(cellInfo)
  row_names <- rownames(cellInfo)
  l_cell_names <- strsplit(rownames(cellInfo), "_")
  l_cell_new_names <- lapply(l_cell_names, function(x){
    categorizer_i <- as.numeric(x[2])
    if(categorizer_i == 2){
      new_cell_name_i <- paste("pLN_",as.character(x[3]),"_",as.character(x[4]), sep = "")
      organ_i <- "_pLN"
    }
    if(categorizer_i == 1){
      new_cell_name_i <- paste("mLN_",as.character(x[3]),"_",as.character(x[4]), sep = "")
      organ_i <- "_mLN"
    }
    #return
    output_i <- c(new_cell_name_i,paste(cell_type,organ_i,sep=""))
    output_i
    
  })
  l_cell_new_name_column <- unlist(lapply(l_cell_new_names, function(x){x[1]}))
  l_cell_new_type_column <- unlist(lapply(l_cell_new_names, function(x){x[2]}))
  cellInfo <- data.frame(cell_type1 = l_cell_new_type_column, nGene = cellInfo$nGene)
  rownames(cellInfo) <- row_names
  saveRDS(cellInfo, file=paste(path_scenic_sample_int,"/cellInfo.Rds",sep=""))
  
  #Assign colors and cells
  clusters_of_interest <- c(as.character(cellInfo[1,1]), as.character(cellInfo[nrow(cellInfo),1]))
  clusters_coloring <- c("blue","green")
  colVars <- list(CellType=setNames(clusters_coloring,
                                    clusters_of_interest))
  saveRDS(colVars, file=paste(path_scenic_sample_int,"/colVars.Rds",sep = ""))
  #Grab cluster IDs of interest
  cell_to_scenic <- names(object_seurat@ident[object_seurat@ident %in% clusters_of_interest])
  #extract count matrix
  exprMat <- as.matrix(object_seurat@raw.data[,cell_to_scenic])
  saveRDS(count_matrix_sce, file=paste(path_scenic_sample_int,"/exprMat.Rds", sep = ""))
  
  #returned
  count_matrix_sce
}

# run function
count_matrix <- makeSCEfromSEURAT(PATH_seurat, clusters_of_interest)

#' makeGRNinput_CorrMat
#' generates SCE object from SEURAT object including the annotation
#' @param count_matrix sce object obtained from seurat
#' @param path_scenic_sample_int path to scenicOptions and files
#' @param smallest_pop smallest expeted cell population as percent of total
#' @return nothing, just stores CorrMat and GRNboost input at /int of scenic object

makeGRNinput_CorrMat <- function(count_matrix, smallest_pop,path_scenic_sample_int){
  #load scenicOptions
  scenicOptions <- readRDS(paste(path_scenic_sample_int,"/scenicOptions.Rds",sep=""))
  
  exprMat <- count_matrix
  nCellsPerGene <- apply(exprMat, 1, function(x) sum(x>0))
  nCountsPerGene <- apply(exprMat, 1, sum)
  
  #Keep only genes that have sufficient read numbers
  # Expressed with 3 UMIs in 1% of the cells
  minReads <- 3*smallest_pop*ncol(exprMat)
  genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > minReads)]
  genesTotal <- names(nCountsPerGene)
  
  #Make sure that genes relevant for small populations are not excluded (e.g. Perycytes)
  # Currently set to 1%
  minSamples <- ncol(exprMat)*smallest_pop
  nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
  genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > minSamples)]
  
  #Only genes in RcisTarget will be used
  # Check overlap
  motifRankings <- importRankings(getDatabases(scenicOptions)[[1]])
  genesInDatabase <- colnames(getRanking(motifRankings))
  genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]
  length(genesLeft_minCells_inDatabases)
  
  # Check whether any relevant gene / potential gene of interest is missing:
  #interestingGenes <- c("Ccl9", "Nkx2-3", "Pdpn", "Il7")
  #interestingGenes[which(!interestingGenes %in% genesLeft_minCells_inDatabases)]
  
  #Filter ExpressionMatrix for the genes present in database
  genesKept <- genesLeft_minCells_inDatabases
  exprMat_filtered <- exprMat[genesKept, ]
  exprMat_filtered <- as.matrix(exprMat_filtered)
  exprMat_filtered <- log2(exprMat_filtered+1)
  
  #Exprot files for GRNboost
  setwd(path_scenic_sample)
  exportsForGRNBoost(exprMat_filtered, scenicOptions)

  #Detection of positive and negative regulations possible thus split dataset
  print("Calculating correlation matrix")
  corrMat <- cor(t(as.matrix(exprMat_filtered)), method="spearman")
  # (Only the rows for TFs will be needed needed):
  #allTFs <- getDbTfs(scenicOptions)
  #corrMat <- corrMat[which(rownames(corrMat) %in% allTFs),]
  saveRDS(corrMat, file=getIntName(scenicOptions, "corrMat"))
  
}

# run function
makeGRNinput_CorrMat(count_matrix, smallest_pop,path_scenic_sample_int)

##### Run python script for GRNBoost
# Continue here: Something is not working, script runs in python
#python.load(paste(path_user_scripts, "FAMEAUX_GRNBoost_via_R.py", sep = ""))

#' scenic_from_GRNboost
#' calculate AUC, rcistarget object from SEURAT object including the annotation
#' @param path_scenic_sample_int path to scenicOptions and files
#' @return nothing scenic defined output

scenic_from_GRNboost <- function(path_scenic_sample_int){
  #load scenicOptions
  scenicOptions <- readRDS(file=paste(path_scenic_sample_int,"/scenicOptions.Rds",sep=""))
  #load exprMat
  exprMat <- readRDS(file=paste(path_scenic_sample_int,"/exprMat.Rds",sep=""))
  
  #####
  #GRNBoost Output processing
  #####
  #Load GRN boost
  ex_network <- read.delim(paste(path_scenic_sample_int,"/GRNBoost_linklist.tsv",sep=""), header=FALSE)
  #Export to use for scenic
  colnames(ex_network) <- c("TF","Target","weight")
  saveRDS(ex_network, file=paste(path_scenic_sample_int,"/1.4_GENIE3_linkList.Rds", sep=""))
  
  #####
  #Perform AUC, Rcistarget
  #####
  setwd(path_scenic_sample)
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
  "TRUE"
}

scenic_status <- scenic_from_GRNboost(path_scenic_sample_int)

#################
#Explore dataset
#################
scenicOptions <- readRDS(paste(path_scenic_sample_int, "/scenicOptions.Rds", sep = ""))
loadInt(scenicOptions)

#####
#Check clustering
#####
setwd(path_scenic_sample)
colVars <- readRDS(paste(path_scenic_sample_int, "/colVars.Rds", sep = ""))
cellInfo <- readRDS(paste(path_scenic_sample_int, "/cellInfo.Rds", sep = ""))



###Comparison tSNE maps
par(mfcol=c(3,5))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files(path_scenic_sample_int), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, cex=.5,
                         varName="cell_type1")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))


# Using only "high-confidence" regulons (normally similar)
par(mfcol=c(3,5))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, cex=.5, varName="cell_type1")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

#chosen t-SNE can then be saved as default to use for plots
scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 20
scenicOptions@settings$defaultTsne$perpl <- 20
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#####
#AUC scores per cell across TF networks
#####
#get cell names and annotation
cellInfo <- readRDS(paste(path_scenic_sample_int,"/cellInfo.Rds",sep=""))
#Names of regulons
#This table contains the non-binarized scores
matrix_AUC_scores <- as.data.frame(getAUC(loadInt(scenicOptions, "aucell_regulonAUC")))

#split into core and extended regulons
names_1 <- strsplit(rownames(matrix_AUC_scores), "_")
#Get number of genes for extended regulons
number_genes_extended <- unlist(sapply(names_1, function(x) {x[2]}))
number_genes_extended_1 <- strsplit(number_genes_extended, " ")
number_genes_extended_2 <- unlist(sapply(number_genes_extended_1, function(x) {x[2]}))
number_genes_extended_3 <- gsub("g","",number_genes_extended_2)
n_g_extended <- as.numeric(gsub("\\(|\\)","",number_genes_extended_3))

#all Names
names_non_extended <- unlist(sapply(names_1, function(x) {x[1]}))
names_non_extended_1 <- strsplit(names_non_extended, " ")
names_final <- unlist(sapply(names_non_extended_1, function(x) {x[1]}))

#Get numbers of genes for nonextended regulons
number_genes_nonextended <- unlist(sapply(names_non_extended_1, function(x) {x[2]}))
number_genes_nonextended_1 <- gsub("g","",number_genes_nonextended)
n_g_core <- as.numeric(gsub("\\(|\\)","",number_genes_nonextended_1))

#Write table
new_columns <- data.frame(GeneSymbol = names_final, n_GRN_extended = n_g_extended, n_GRN_core = n_g_core)
matrix_AUC_scores_final <- cbind(new_columns,matrix_AUC_scores)
matrix_AUC_scores_final_core <- matrix_AUC_scores_final[!is.na(matrix_AUC_scores_final$n_GRN_core),]
matrix_AUC_scores_final_extended <- matrix_AUC_scores_final[!is.na(matrix_AUC_scores_final$n_GRN_extended),]

write.table(matrix_AUC_scores_final, paste(path_scenic_sample,"/AUCell_scores_ext_core_",organ,"_",cell_type,".csv", sep=""), dec=".", sep=",")
write.table(matrix_AUC_scores_final_core, paste(path_scenic_sample,"/AUCell_scores_core_",organ,"_",cell_type,".csv", sep=""), dec=".", sep=",")
write.table(matrix_AUC_scores_final_extended, paste(path_scenic_sample,"/AUCell_scores_ext_",organ,"_",cell_type,".csv", sep=""), dec=".", sep=",")







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

#thresholds for all regulons
aucell_thresholds <- loadInt(scenicOptions, "aucell_thresholds")

#Only if data has been binarized
#aucell_binary_full <- loadInt(scenicOptions, "aucell_binary_full")

#####
#Extract RCM for cell subsets
#####
#Grab cluster IDs of interest
cell_to_scenic <- names(LN_minus_mLN@ident[LN_minus_mLN@ident %in% clusters_of_interest])
#extract count matrix
exprMat <- as.matrix(LN_minus_mLN@raw.data[,cell_to_scenic])
saveRDS(count_matrix_sce, file=paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/mLN_SPF/",cell_type,"/int/exprMat.Rds", sep = ""))
exprMat <- readRDS(paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/mLN_SPF/",cell_type,"/int/exprMat.Rds", sep = ""))






#########################
#Add ons
#########################
#####
#Genes included in a regulon
#####
regulons <- loadInt(scenicOptions, "regulons")
regulons[c("Irf7","Fos")]
#"Rel","Junb","Jun","Egr1"

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


#####
#Interactively check out the data
#####
loadInt(scenicOptions)
exprMat <- readRDS(paste(path_scenic_rds, "/exprMat.Rds", sep = ""))
logMat <- log2(exprMat+1) # Better if it is logged/normalized
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat) #default t-SNE
savedSelections <- shiny::runApp(aucellApp)

#####
#Plot AUC scores on t-SNE
#####
#par(bg = "black")
par(mfrow=c(1,2))

regulonNames <- c("Irf7", "Fos")
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)
text(-30,-25, attr(cellCol,"red"), col="red", cex=.7, pos=4)
text(-30,-25-4, attr(cellCol,"green"), col="green3", cex=.7, pos=4)

regulonNames <- list(red=c("Irf7"),
                     green=c("Fos"),
                     blue=c( "Irf1"))
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC")
text(-30,-25, attr(cellCol,"red"), col="red", cex=.7, pos=4)
text(-30,-25-4, attr(cellCol,"green"), col="green3", cex=.7, pos=4)
text(-30,-25-8, attr(cellCol,"blue"), col="blue", cex=.7, pos=4)

#####
#Regulon t-SNE
#####
par(mfrow=c(2,4))
tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat,
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Irf1","Arntl","Fos","Fosl1",
                                                                                                   "Rel","Junb","Jun","Egr1")],],
                        plots="Expression")

##