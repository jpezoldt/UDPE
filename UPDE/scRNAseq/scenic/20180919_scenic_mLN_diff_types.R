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
#Load and label mLN
#####
#Label mLN cells with 1_1 & 1_2
#Experiment 2
sample_1_2.data <- Read10X(data.dir = "/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/2018_NC_stroma/Exp_252/mLN_SPF/B6/filtered_gene_bc_matrices/mm10")
#sample_1_2.data <- Read10X(data.dir = "/Users/Pezoldt/PowerFolders/R/2017_252_scRNAseq/Data/mLN_SPF/B6/filtered_gene_bc_matrices/mm10")
sample_1_2.data@Dimnames[[2]] <- paste("sample_1_2_", c(1:length(sample_1_2.data@Dimnames[[2]])), sep = "")
#Experiment 1
sample_1_1.data <- read.csv("/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/2018_NC_stroma/Exp_244/mLNSPF.csv")
#sample_1_1.data <- read.csv("/Users/Pezoldt/PowerFolders/R/2017_244_scRNASeq/Data/mLNSPF.csv")
row.names(sample_1_1.data) <- sample_1_1.data$X
sample_1_1.data <- sample_1_1.data[,c(2:ncol(sample_1_1.data))]
colnames(sample_1_1.data) <- paste("sample_1_1_", c(1:ncol(sample_1_1.data)), sep = "")

#######
#Global variables
#######
#Seurat--------------------------------------------------
sample_1 <- "mLN"
#Important: If changing global variable Setup of pLN and mLN needs to be manually adjusted as the cluster IDs get changed
min_gene_number <- 250
mito_cutoff <- 0.045
nGene <- 3500
min_cells <- 20
resolution_single <- 1.2
#Parameters FSC only
pca_use <- 24
n_variable_genes <- 1500
resolution_use <- 1.3
#Scenic---------------------------------------------------
#Smallest expected cell population
smallest_pop <- 0.01
#Path store RDS
path_scenic_rds <- "/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/mLN_SPF/nonAdventi"
#"/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/mLN_pilot"
#/data/pezoldt/Analysis/scRNAseq/scenic

##############
#Setup sample_1
##############
#sample_1_seurat <- CreateSeuratObject(raw.data = sample_1.data, min.cells = min_cells, min.genes = min_gene_number)

#Merger
sample_1_1_seurat <- CreateSeuratObject(raw.data = sample_1_1.data, min.cells = min_cells, min.genes = min_gene_number)
sample_1_2_seurat <- CreateSeuratObject(raw.data = sample_1_2.data, min.cells = min_cells, min.genes = min_gene_number)
sample_1_seurat <- MergeSeurat(sample_1_1_seurat, sample_1_2_seurat, do.normalize = FALSE)

#Eliminate mitochondrial high cells and dublets
#Identify cells with mitochondrial readss
mito.genes <- grep(pattern = "^mt-", x = rownames(sample_1_seurat@data), value = TRUE)
percent.mito <- colSums(as.matrix(sample_1_seurat@raw.data[mito.genes, ]))/colSums(as.matrix(sample_1_seurat@raw.data))

# AddMetaData adds columns to object@data.info and good to store QC stats
sample_1_seurat <- AddMetaData(object = sample_1_seurat, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = sample_1_seurat, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(object = sample_1_seurat, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = sample_1_seurat, gene1 = "nUMI", gene2 = "nGene")

#Identify cells with high rpl reads
rpl.genes <- grep(pattern = "^Rpl", x = rownames(sample_1_seurat@data), value = TRUE)
percent.rpl <- colSums(as.matrix(sample_1_seurat@raw.data[rpl.genes, ]))/colSums(as.matrix(sample_1_seurat@raw.data))

# AddMetaData for RPL percentage
sample_1_seurat <- AddMetaData(object = sample_1_seurat, metadata = percent.rpl, col.name = "percent.rpl")


# We filter out cells according unique gene counts 
# and percentage mitochondrial reads
sample_1_seurat <- FilterCells(object = sample_1_seurat, subset.names = c("nGene", "percent.mito"), 
                               low.thresholds = c(min_gene_number, -Inf), high.thresholds = c(nGene, mito_cutoff))

#Normalize the expression and log-transform
sample_1_seurat <- NormalizeData(object = sample_1_seurat, normalization.method = "LogNormalize", 
                                 scale.factor = 10000)

#Find variable genes independent of expression level
sample_1_seurat <- FindVariableGenes(object = sample_1_seurat, mean.function = ExpMean, dispersion.function = LogVMR, 
                                     x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

#scale data
sample_1_seurat <- ScaleData(object = sample_1_seurat, vars.to.regress = c("nUMI", "percent.mito","percent.rpl"))

#perfrom PCA and print genes that define PCA
sample_1_seurat <- RunPCA(object = sample_1_seurat, pc.genes = sample_1_seurat@var.genes, do.print = TRUE, pcs.print = 1:5, 
                          genes.print = 5)

# Examine and visualize PCA results
PrintPCA(object = sample_1_seurat, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = sample_1_seurat, pcs.use = 1:6)
PCAPlot(object = sample_1_seurat, dim.1 = 1, dim.2 = 2)

#ProjectPCA scores each gene in the dataset 
sample_1_seurat <- ProjectPCA(object = sample_1_seurat, do.print = FALSE)
PCHeatmap(object = sample_1_seurat, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)


#Group the cells into clusters
sample_1_seurat <- FindClusters(object = sample_1_seurat, reduction.type = "pca", dims.use = 1:13, 
                                resolution = resolution_single, print.output = 0, save.SNN = TRUE,
                                force.recalc = TRUE)
PrintFindClustersParams(object = sample_1_seurat)

#perfrom Tsne
sample_1_seurat <- RunTSNE(object = sample_1_seurat, dims.use = 1:13, do.fast = TRUE)
TSNEPlot(object = sample_1_seurat, pt.size = 0.3, do.label = TRUE)

#####
#Eliminate LECs, BECs sample_1
#####
#Identify LECs and BECs
VlnPlot(object = sample_1_seurat, features.plot = c("Pecam1"))
sample_1_aExp <- AverageExpression(sample_1_seurat)
aExp_Pecam1 <- subset(sample_1_aExp, rownames(sample_1_aExp) == c("Pecam1"))
#Grab colnames/clusters that have Pecam1 expression < 1 on average per cluster
cluster_keep <- as.numeric(rownames(subset(as.data.frame(t(aExp_Pecam1)), Pecam1 < 1)))

#Get cells that are not within the cluster BEC and LEC cluster
all_cells <- sample_1_seurat@meta.data
all_cells_keep <- subset(all_cells, res.1.2 == 1000)

#attain all cells per cluster ID and concatenate the tables
for(i in 1:length(cluster_keep)){
  cluster_id_i <- cluster_keep[i]
  print(cluster_id_i)
  cells_cluster_i <- subset(all_cells, res.1.2 == cluster_id_i)
  print(nrow(cells_cluster_i))
  all_cells_keep <- rbind(all_cells_keep, cells_cluster_i)
}
sample_1_NO_LECs_BECs <- as.character(rownames(all_cells_keep))

#####
#Select LECs and BECs sample_1
#####
#Identify LECs and BECs
VlnPlot(object = sample_1_seurat, features.plot = c("Pecam1"))
sample_1_aExp <- AverageExpression(sample_1_seurat)
aExp_Pecam1 <- subset(sample_1_aExp, rownames(sample_1_aExp) == c("Pecam1"))
#Grab colnames/clusters that have Pecam1 expression < 1 on average per cluster
cluster_keep <- as.numeric(rownames(subset(as.data.frame(t(aExp_Pecam1)), Pecam1 > 1)))

#Get cells that are not within the cluster BEC and LEC cluster
all_cells <- sample_1_seurat@meta.data
all_cells_keep <- subset(all_cells, res.1.2 == 1000)

#attain all cells per cluster ID and concatenate the tables
for(i in 1:length(cluster_keep)){
  cluster_id_i <- cluster_keep[i]
  print(cluster_id_i)
  cells_cluster_i <- subset(all_cells, res.1.2 == cluster_id_i)
  print(nrow(cells_cluster_i))
  all_cells_keep <- rbind(all_cells_keep, cells_cluster_i)
}
sample_1_ONLY_LECs_BECs <- as.character(rownames(all_cells_keep))


##########################
#Analyse ONLY FSC
##########################
sample_1_1_seurat <- CreateSeuratObject(raw.data = sample_1_1.data, min.cells = min_cells, min.genes = min_gene_number)
sample_1_2_seurat <- CreateSeuratObject(raw.data = sample_1_2.data, min.cells = min_cells, min.genes = min_gene_number)
sample_1_seurat <- MergeSeurat(sample_1_1_seurat, sample_1_2_seurat, do.normalize = FALSE)
#Eliminate mitochondrial high cells and dublets
#Identify cells with mitochondrial readss
mito.genes <- grep(pattern = "^mt-", x = rownames(sample_1_seurat@data), value = TRUE)
percent.mito <- colSums(as.matrix(sample_1_seurat@raw.data[mito.genes, ]))/colSums(as.matrix(sample_1_seurat@raw.data))

# AddMetaData adds columns to object@data.info and good to store QC stats
sample_1_seurat <- AddMetaData(object = sample_1_seurat, metadata = percent.mito, col.name = "percent.mito")

#Identify cells with high rpl reads
rpl.genes <- grep(pattern = "^Rpl", x = rownames(sample_1_seurat@data), value = TRUE)
percent.rpl <- colSums(as.matrix(sample_1_seurat@raw.data[rpl.genes, ]))/colSums(as.matrix(sample_1_seurat@raw.data))

# AddMetaData for RPL percentage
sample_1_seurat <- AddMetaData(object = sample_1_seurat, metadata = percent.rpl, col.name = "percent.rpl")

#Sample number of SC cells
set.seed(123)
sampled_SCs <- sample(sample_1_NO_LECs_BECs, 5000)

#Filter for mitochondiral and duplets and LECs and BECs out
sample_1_seurat_minus <- FilterCells(object = sample_1_seurat, subset.names = c("nGene", "percent.mito"), 
                                     low.thresholds = c(min_gene_number, -Inf), high.thresholds = c(nGene, mito_cutoff),
                                     cells.use = sample_1_NO_LECs_BECs)
sample_1_seurat_minus <- FilterCells(object = sample_1_seurat_minus, subset.names = c("Ackr4","Pecam1"),
                                     low.thresholds = c(-Inf, -Inf), high.thresholds = c(1,1))
sampled_cells <- rownames(sample_1_seurat_minus@meta.data)

#extract count matrix
exprMat_SC <- sample_1_seurat_minus@raw.data[,sampled_cells]


##########################
#Analyse ONLY FSC
##########################
sample_1_1_seurat <- CreateSeuratObject(raw.data = sample_1_1.data, min.cells = min_cells, min.genes = min_gene_number)
sample_1_2_seurat <- CreateSeuratObject(raw.data = sample_1_2.data, min.cells = min_cells, min.genes = min_gene_number)
sample_1_seurat <- MergeSeurat(sample_1_1_seurat, sample_1_2_seurat, do.normalize = FALSE)
#Eliminate mitochondrial high cells and dublets
#Identify cells with mitochondrial readss
mito.genes <- grep(pattern = "^mt-", x = rownames(sample_1_seurat@data), value = TRUE)
percent.mito <- colSums(as.matrix(sample_1_seurat@raw.data[mito.genes, ]))/colSums(as.matrix(sample_1_seurat@raw.data))

# AddMetaData adds columns to object@data.info and good to store QC stats
sample_1_seurat <- AddMetaData(object = sample_1_seurat, metadata = percent.mito, col.name = "percent.mito")

#Identify cells with high rpl reads
rpl.genes <- grep(pattern = "^Rpl", x = rownames(sample_1_seurat@data), value = TRUE)
percent.rpl <- colSums(as.matrix(sample_1_seurat@raw.data[rpl.genes, ]))/colSums(as.matrix(sample_1_seurat@raw.data))

# AddMetaData for RPL percentage
sample_1_seurat <- AddMetaData(object = sample_1_seurat, metadata = percent.rpl, col.name = "percent.rpl")

#Filter for mitochondiral and duplets and LECs and BECs out
sample_1_seurat_minus <- FilterCells(object = sample_1_seurat, subset.names = c("nGene", "percent.mito"), 
                                     low.thresholds = c(min_gene_number, -Inf), high.thresholds = c(nGene, mito_cutoff),
                                     cells.use = sample_1_NO_LECs_BECs)
sample_1_seurat_minus <- FilterCells(object = sample_1_seurat_minus, subset.names = c("Ackr4","Pecam1"),
                                     low.thresholds = c(-Inf, -Inf), high.thresholds = c(1,1))
sample_1_seurat_minus <- NormalizeData(sample_1_seurat_minus)

#Find variable genes independent of expression level
sample_1_seurat_minus <- FindVariableGenes(object = sample_1_seurat_minus, mean.function = ExpMean, dispersion.function = LogVMR, 
                                           x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

#number of variable genes
length(x = sample_1_seurat_minus@var.genes)
#scale data
sample_1_seurat_minus <- ScaleData(object = sample_1_seurat_minus, vars.to.regress = c("nUMI", "percent.mito","percent.rpl"))

#Pick Number variable genes
var_genes <- rownames(head(sample_1_seurat_minus@hvg.info, n_variable_genes))

#perfrom PCA and print genes that define PCA
sample_1_seurat_minus <- RunPCA(object = sample_1_seurat_minus, pc.genes = var_genes, do.print = TRUE, pcs.print = 1:5, 
                                genes.print = 5, pcs.compute = 30)

#ProjectPCA scores each gene in the dataset 
sample_1_seurat_minus <- ProjectPCA(object = sample_1_seurat_minus, do.print = FALSE)
#PCHeatmap(object = sample_1_seurat_minus, pc.use = 1:pca_use, cells.use = 500, do.balanced = TRUE, 
#         label.columns = FALSE, use.full = FALSE)

#Group the cells into clusters
#perfrom Tsne
sample_1_seurat_minus <- RunTSNE(object = sample_1_seurat_minus, dims.use = 1:pca_use, do.fast = TRUE)
sample_1_seurat_minus <- FindClusters(object = sample_1_seurat_minus, reduction.type = "pca", dims.use = 1:pca_use, 
                                      resolution = resolution_use, print.output = 0, save.SNN = TRUE,
                                      force.recalc = TRUE)

PrintFindClustersParams(object = sample_1_seurat_minus)

#Change Cluster Names
current.cluster.ids <- c(0:12)

new.cluster.ids <- c("Ccl19highMadcam1+","Inmt+Cxcl12+","Inmt+","Cd34+Gdf10+",
                     "Cd34+Has1+","Cd34+Aldh1a2+","Il6+Cxcl1","Cd34+Ackr3+",
                     "Cd34+Cd248+","Ccl19+Il7+","pSC","PvC",
                     "SCx")

sample_1_seurat_minus@ident <- plyr::mapvalues(x = sample_1_seurat_minus@ident, from = current.cluster.ids, to = new.cluster.ids)

FeaturePlot(sample_1_seurat_minus, features.plot = c("Madcam1","Ackr3", "Aldh1a2",
                                                     "Ccl19","Il6","Vcam1",
                                                     "Icam1","Inmt","Nr4a1",
                                                     "Tnfsf11","Cxcl13"),
            min.cutoff = "q9",
            cols.use = c('gainsboro', 'darkred'),
            pt.size = 0.1)

TSNEPlot(object = sample_1_seurat_minus, pt.size = 0.3, do.label = TRUE)

#####
#Extract RCM for cell subsets
#####
#Grab CD34+ cluster IDs
Adventi_clusters <- c("Cd34+Gdf10+",
                     "Cd34+Has1+","Cd34+Aldh1a2+","Cd34+Ackr3+",
                     "Cd34+Cd248+")
nonAdventi_clusters <- c("Ccl19highMadcam1+","Inmt+Cxcl12+","Inmt+","Il6+Cxcl1","Cd34+Ackr3+",
                    "Ccl19+Il7+")
Adventi_cells <- names(sample_1_seurat_minus@ident[sample_1_seurat_minus@ident %in% Adventi_clusters])
nonAdventi_cells <- names(sample_1_seurat_minus@ident[sample_1_seurat_minus@ident %in% nonAdventi_clusters])



#extract count matrix
exprMat_Adventi <- as.matrix(sample_1_seurat_minus@raw.data[,Adventi_cells])
exprMat_nonAdventi <- as.matrix(sample_1_seurat_minus@raw.data[,nonAdventi_cells])
saveRDS(exprMat_Adventi, file=paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/mLN_SPF/Adventi/int/exprMat_Adventi.Rds", sep=""))
saveRDS(exprMat_nonAdventi, file=paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/mLN_SPF/nonAdventi/int/exprMat_nonAdventi.Rds", sep=""))


#Note: Use Seurat output for SCENIC
##################
#SCENIC
##################
#####
#Dataset to run
#####
exprMat <- exprMat_Adventi
#####
#Settings
#####
#set storage path
setwd(path_scenic_rds)
#Set Rcistarget location
org="mgi"
dbDir="/data/software/R-3.5.0/library/RcisTarget/db" # RcisTarget databases location
myDatasetTitle="mLN-SPF non-Endothelial SC"
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, datasetTitle=myDatasetTitle, nCores=12) 
#Set to mm10
scenicOptions@settings$dbs <- setNames("mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather", "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

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

#Delete matrix
rm(exprMat)
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
ex_network <- read.delim("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/mLN_SPF/Adventi/int/GRNBoost_linklist.tsv", header=FALSE)
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
colnames(ex_network) <- c("TF","Target","weight")
saveRDS(ex_network, file="int/1.4_GENIE3_linkList.Rds")

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

#####
#Reload Dataset of choice
#####
cell_type = "Adventi"
path_scenic_rds <- paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/mLN_SPF/", cell_type, "/int", sep="")
scenicOptions <- readRDS(file=paste(path_scenic_rds,"/scenicOptions.Rds",sep=""))
setwd(paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/mLN_SPF/", cell_type, "/", sep=""))

#chosen t-SNE can then be saved as default to use for plots
scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 10
scenicOptions@settings$defaultTsne$perpl <- 20
saveRDS(scenicOptions, file=paste(path_scenic_rds,"/scenicOptions.Rds",sep=""))

######
#Binarize and analyze
######

#Optional: Binarize AUC 
logMat <- exprMat # Better if it is logged/normalized
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
savedSelections <- shiny::runApp(aucellApp)

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
tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

# Show TF expression:
par(mfrow=c(1,1))
#Note: not functional
AUCell::AUCell_plotTSNE(tSNE_scenic$Y,
                        exprMat, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Nkx2-3", "Atf7", "Nfkb1", "Maf")],],
                        plots="Expression")

#Density plot to detect most likely stable states
library(KernSmooth)
library(RColorBrewer)
dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)

#Show several regulons simultaneously
par(mfrow=c(1,2))
#Note: only works with binarized data
regulonNames <- c( "Atf1","Atf2")
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)
text(-30,-25, attr(cellCol,"red"), col="red", cex=.7, pos=4)
text(-30,-25-4, attr(cellCol,"green"), col="green3", cex=.7, pos=4)

regulonNames <- list(red=c("Sox10", "Sox8"),
                     green=c("Maf"),
                     blue=c( "Neurod2"))
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC")
text(-30,-25, attr(cellCol,"red"), col="red", cex=.7, pos=4)
text(-30,-25-4, attr(cellCol,"green"), col="green3", cex=.7, pos=4)
text(-30,-25-8, attr(cellCol,"blue"), col="blue", cex=.7, pos=4)

#####
#Regulon targets
#####
#Only regulons with 10 genes or more are scored
regulons <- loadInt(scenicOptions, "regulons")
regulons[c("Atf3","Irf2","Rxra","Bach1","Erf","Stat3","Rarb")]

#Analyze Links of TF to genes and Motifs
regulons <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
regulonTargetsInfo <- RcisTarget::addLogo(regulonTargetsInfo, motifCol="bestMotif")
regulonTargetsInfo$Genie3Weight <- signif(regulonTargetsInfo$Genie3Weight, 2)
colsToShow <- c("TF", "gene", "nMotifs", "bestMotif", "logo", "NES", "highConfAnnot", "Genie3Weight")
DT::datatable(regulonTargetsInfo[TF=="Atf3" & highConfAnnot==TRUE, colsToShow, with=F], escape=FALSE, filter="top")

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
library(SingleCellExperiment)
#load("data/sceMouseBrain.RData")
#dgem <- counts(sceMouseBrain)
#head(colnames(dgem))  #should contain the Cell ID/name
# Export:
#scenicOptions@fileNames$output["loomFile",] <- "mouseBrain.loom"
#export2scope(scenicOptions, dgem)
