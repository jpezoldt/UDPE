#!/usr/bin/env Rscript
#args <- commandArgs(trailingOnly = T)
#blah = args[1]

# Autor: Joern Pezoldt
# 22.11.2018
# Function:
#1) Align Ontogeny samples using Seurat
#2) Identifiy DEGs per cluster


#Libraries
library("Seurat")
library("cowplot")
library("Matrix")
library("magrittr")
library("dplyr")
library("pryr")

#####
#PATHs & Global Variables
#####
#PATH
PATH_Output <- "/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/seurat/Ontogney_Multialign/D0_to_D300"
#Signatures from mLN
NC_Pezoldt_Pasztoi_2018 <- read.csv("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/2018_NatComm_Pezoldt_Pasztoi/NC_2018_mLN_Clusters.csv",sep = ";")
#Cell cycle genes
Cell_cycle_Regev_BLOOD <- read.delim("/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/Cell_Cycle_Regev/regev_lab_cell_cycle_genes.txt")
#Imprinted signature NC 2018 Pezoldt
NC_Pezoldt_Pasztoi_2018_repressed_mLN <- read.delim("/home/pezoldt/NAS2/pezoldt/Analysis/RNAseq/2017_HZI_FSC_Tx_DESeq2/NC_2018/NC_2018_Tx_FSC_imprinted_repressed_mLN.txt")
NC_Pezoldt_Pasztoi_2018_maintained_mLN <- read.delim("/home/pezoldt/NAS2/pezoldt/Analysis/RNAseq/2017_HZI_FSC_Tx_DESeq2/NC_2018/NC_2018_Tx_FSC_imprinted_maintained_mLN.txt")
#TFs from DARs
DAR_TF_existence_matrix <- read.table("/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/Motif/Homer_Output/allHomerResults/size_500bp/SPF/SPF_existence_matrix_known_homer_compile.txt")


#######
#Global variables
#######
#@User: Define sample names
sample_1 <- "d000"
sample_2 <- "d010"
sample_3 <- "d024"
sample_4 <- "d056"
sample_5 <- "d300"
v_sample <- c(sample_1, sample_2, sample_3, sample_4,sample_5)
#setwd("/Users/Pezoldt/PowerFolders/R/2018_Exp270_Ontogeny/Output/Multialign_SEURAT")

#Parameters Seurat
min_gene_number <- 1000
mito_cutoff <- 0.045
nGene <- 4700
min_cells <- 20
dims_use <- 12
resolution_single <- 1.2
resolution_merge <- 1.0
resolution_merge_MAIN <- 1.0

#######
#Load and label
#######
###
#mLN d0
###
#sample_1.data <- Read10X(data.dir = "C:/Users/jpe12/PowerFolders/R/2018_Exp270_Ontogeny/Data/neonatal/filtered_gene_bc_matrices/mm10")
sample_1.data <- Read10X(data.dir = "/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/270_2018-01-03_scRNA-Seq_mLN_ontogeny/Data/C11_neo_mLN/outs/filtered_gene_bc_matrices/mm10")
sample_1.data@Dimnames[[2]] <- paste("D000_", c(1:length(sample_1.data@Dimnames[[2]])), sep = "")
###
#mLN d10
###
#sample_2.data <- Read10X(data.dir = "C:/Users/jpe12/PowerFolders/R/2018_Exp270_Ontogeny/Data/d10/filtered_gene_bc_matrices/mm10")
sample_2.data <- Read10X(data.dir = "/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/270_2018-01-03_scRNA-Seq_mLN_ontogeny/Data/C12_d10_mLN/outs/filtered_gene_bc_matrices/mm10")
sample_2.data@Dimnames[[2]] <- paste("D010_", c(1:length(sample_2.data@Dimnames[[2]])), sep = "")

###
#mLN d24
###
#sample_3.data <- Read10X(data.dir = "C:/Users/jpe12/PowerFolders/R/2018_Exp270_Ontogeny/Data/d24/filtered_gene_bc_matrices/mm10")
sample_3.data <- Read10X(data.dir = "/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/270_2018-01-03_scRNA-Seq_mLN_ontogeny/Data/D1_d24_mLN/outs/filtered_gene_bc_matrices/mm10")
sample_3.data@Dimnames[[2]] <- paste("D024_", c(1:length(sample_3.data@Dimnames[[2]])), sep = "")

###
#mLN d60
###
#Experiment 1
# Cell IDs: 1 to nrow(Experiment1)
sample_4_1.data <- Read10X("/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/244_scRNA-Seq_mLN_pLN_SPF/Data/10X_results/L1700567_mLN_SPF/outs/filtered_gene_bc_matrices/mm10")
sample_4_1.data@Dimnames[[2]] <- paste("D056_", c(1:ncol(sample_4_1.data)), sep = "")
#Experiment 2
# Cell IDs: nrow(Experiment1) to nrow(Experiment2)
#sample_4.data <- Read10X(data.dir = "C:/Users/jpe12/PowerFolders/R/2017_252_scRNAseq/Data/mLN_SPF/B6/filtered_gene_bc_matrices/mm10")
sample_4_2.data <- Read10X(data.dir = "/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/252_2017-07-19_scRNA-Seq_mLN_pLN_SPF_GF/Data/10x/mLN_SPF/B6/outs/filtered_gene_bc_matrices/mm10")
sample_4_2.data@Dimnames[[2]] <- paste("D056_", c((ncol(sample_4_1.data)+1):(ncol(sample_4_1.data)+length(sample_4_2.data@Dimnames[[2]]))), sep = "")

###
#mLN d300
###
sample_5.data <- Read10X(data.dir = "/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/275_2018-04-23_scRNA-Seq_mLN_SPF_45wk/Data/E7_SPF_mLN_42wk/outs/filtered_gene_bc_matrices/mm10")
sample_5.data@Dimnames[[2]] <- paste("D300_", c(1:length(sample_5.data@Dimnames[[2]])), sep = "")


# Convert to sparse matrices for efficiency
sample_1.data <- as(as.matrix(sample_1.data), "dgCMatrix")
sample_2.data <- as(as.matrix(sample_2.data), "dgCMatrix")
sample_3.data <- as(as.matrix(sample_3.data), "dgCMatrix")
sample_4_1.data <- as(as.matrix(sample_4_1.data), "dgCMatrix")
sample_4_2.data <- as(as.matrix(sample_4_2.data), "dgCMatrix")
sample_5.data <- as(as.matrix(sample_5.data), "dgCMatrix")

#####
#Setup Samples
#####
# Create and setup Seurat objects for each dataset
#Sample_1
sample_1_seurat <- CreateSeuratObject(raw.data = sample_1.data, min.cells = min_cells, min.genes = min_gene_number)
sample_1_seurat <- NormalizeData(sample_1_seurat)
sample_1_seurat <- FindVariableGenes(sample_1_seurat, do.plot = F, display.progress = F)
sample_1_seurat <- ScaleData(sample_1_seurat)
sample_1_seurat@meta.data$tech <- sample_1

#sample_2
sample_2_seurat <- CreateSeuratObject(raw.data = sample_2.data, min.cells = min_cells, min.genes = min_gene_number)
sample_2_seurat <- NormalizeData(sample_2_seurat)
sample_2_seurat <- FindVariableGenes(sample_2_seurat, do.plot = F, display.progress = F)
sample_2_seurat <- ScaleData(sample_2_seurat)
sample_2_seurat@meta.data$tech <- sample_2

#sample_3
sample_3_seurat <- CreateSeuratObject(raw.data = sample_3.data, min.cells = min_cells, min.genes = min_gene_number)
sample_3_seurat <- NormalizeData(sample_3_seurat)
sample_3_seurat <- FindVariableGenes(sample_3_seurat, do.plot = F, display.progress = F)
sample_3_seurat <- ScaleData(sample_3_seurat)
sample_3_seurat@meta.data$tech <- sample_3

#sample_4
sample_4_1_seurat <- CreateSeuratObject(raw.data = sample_4_1.data, min.cells = min_cells, min.genes = min_gene_number)
sample_4_2_seurat <- CreateSeuratObject(raw.data = sample_4_2.data, min.cells = min_cells, min.genes = min_gene_number)
sample_4_seurat <- MergeSeurat(sample_4_1_seurat, sample_4_2_seurat, do.normalize = FALSE)
sample_4_seurat <- NormalizeData(sample_4_seurat)
sample_4_seurat <- FindVariableGenes(sample_4_seurat, do.plot = F, display.progress = F)
sample_4_seurat <- ScaleData(sample_4_seurat)
sample_4_seurat@meta.data$tech <- sample_4

#sample_5
sample_5_seurat <- CreateSeuratObject(raw.data = sample_5.data, min.cells = min_cells, min.genes = min_gene_number)
sample_5_seurat <- NormalizeData(sample_5_seurat)
sample_5_seurat <- FindVariableGenes(sample_5_seurat, do.plot = F, display.progress = F)
sample_5_seurat <- ScaleData(sample_5_seurat)
sample_5_seurat@meta.data$tech <- sample_5

######
#Determine the genes dominant in the PCAs
######
#Merge the seurat objects
mLN_Ontogeny.merged <- MergeSeurat(sample_1_seurat, sample_2_seurat, do.normalize = FALSE)
mLN_Ontogeny.merged <- MergeSeurat(mLN_Ontogeny.merged, sample_3_seurat, do.normalize = FALSE)
mLN_Ontogeny.merged <- MergeSeurat(mLN_Ontogeny.merged, sample_4_seurat, do.normalize = FALSE)
mLN_Ontogeny.merged <- MergeSeurat(mLN_Ontogeny.merged, sample_5_seurat, do.normalize = FALSE)
mLN_Ontogeny.merged <- NormalizeData(mLN_Ontogeny.merged)
mLN_Ontogeny.merged <- FindVariableGenes(mLN_Ontogeny.merged, do.plot = F, display.progress = F)
mLN_Ontogeny.merged <- ScaleData(mLN_Ontogeny.merged)
#Run PCA
mLN_Ontogeny.merged <- RunPCA(object = mLN_Ontogeny.merged, pc.genes = mLN_Ontogeny.merged@var.genes, do.print = TRUE, pcs.print = 1:10, 
                              genes.print = 10)

# Examine and visualize PCA results
PrintPCA(object = mLN_Ontogeny.merged, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = mLN_Ontogeny.merged, pcs.use = 1:6)
PCAPlot(object = mLN_Ontogeny.merged, dim.1 = 1, dim.2 = 2)

# Determine genes to use for CCA, must be highly variable in at least 2 datasets
ob.list <- list(sample_1_seurat, sample_2_seurat, sample_3_seurat, sample_4_seurat,sample_5_seurat)
genes.use <- c()
for (i in 1:length(ob.list)) {
  genes.use <- c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 1000))
}
genes.use <- names(which(table(genes.use) > 1))
for (i in 1:length(ob.list)) {
  genes.use <- genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
}
saveRDS(ob.list, paste(PATH_Output,"/obs.list_1000.Rds", sep=""))

# Run multi-set CCA
mLN_Ontogeny.integrated <- RunMultiCCA(ob.list, genes.use = genes.use, num.ccs = dims_use)

# CC Selection
MetageneBicorPlot(mLN_Ontogeny.integrated, grouping.var = "tech", dims.eval = 1:dims_use)

# Run rare non-overlapping filtering
mLN_Ontogeny.integrated <- CalcVarExpRatio(object = mLN_Ontogeny.integrated, reduction.type = "pca",
                                           grouping.var = "tech", dims.use = 1:dims_use)
mLN_Ontogeny.integrated <- SubsetData(mLN_Ontogeny.integrated, subset.name = "var.ratio.pca",
                                      accept.low = 0.5)

# Alignment
mLN_Ontogeny.integrated <- AlignSubspace(mLN_Ontogeny.integrated,
                                         reduction.type = "cca",
                                         grouping.var = "tech",
                                         dims.align = 1:dims_use)

#Plot the aligned CCA ACCA
#p1 <- VlnPlot(mLN_Ontogeny.integrated, features.plot = "ACC1", group.by = "tech", do.return = T)
#p2 <- VlnPlot(mLN_Ontogeny.integrated, features.plot = "ACC2", group.by = "tech", do.return = T)
#plot_grid(p1, p2)


# t-SNE and Clustering
mLN_Ontogeny.integrated <- FindClusters(mLN_Ontogeny.integrated, reduction.type = "cca",
                                        dims.use = 1:dims_use, save.SNN = T, resolution = resolution_merge)
mLN_Ontogeny.integrated <- RunTSNE(mLN_Ontogeny.integrated,
                                   reduction.use = "cca",
                                   dims.use = 1:dims_use, save.SNN = T,
                                   resolution = resolution_merge, force.recalc = TRUE)



#Plot the aligned CCA ACCA
#p1 <- VlnPlot(mLN_Ontogeny.integrated, features.plot = "ACC1", group.by = "tech", do.return = T)
#p2 <- VlnPlot(mLN_Ontogeny.integrated, features.plot = "ACC2", group.by = "tech", do.return = T)
#plot_grid(p1, p2)


#Check allocation of cells
#sample_1_cells <- grep(pattern = "^sample_1_", x = rownames(mLN_Ontogeny.integrated@meta.data), value = TRUE)
#sample_2_cells <- grep(pattern = "^sample_2_", x = rownames(mLN_Ontogeny.integrated@meta.data), value = TRUE)
#sample_3_cells <- grep(pattern = "^sample_3_", x = rownames(mLN_Ontogeny.integrated@meta.data), value = TRUE)
#sample_4_cells <- grep(pattern = "^sample_4_", x = rownames(mLN_Ontogeny.integrated@meta.data), value = TRUE)
#p1 <- TSNEPlot(mLN_Ontogeny.integrated, group.by = "tech", do.return = T, pt.size = 0.4, cells.use = sample_1_cells)
#p2 <- TSNEPlot(mLN_Ontogeny.integrated, group.by = "tech", do.return = T, pt.size = 0.4, cells.use = sample_2_cells)
#p3 <- TSNEPlot(mLN_Ontogeny.integrated, group.by = "tech", do.return = T, pt.size = 0.4, cells.use = sample_3_cells)
#p4 <- TSNEPlot(mLN_Ontogeny.integrated, group.by = "tech", do.return = T, pt.size = 0.4, cells.use = sample_4_cells)
#p5 <- TSNEPlot(mLN_Ontogeny.integrated, do.return = T, pt.size = 0.2, do.label = TRUE)

#p6 <- TSNEPlot(mLN_Ontogeny.integrated, group.by = "tech", do.return = T, pt.size = 0.2,
 #              colors.use = c("orange","green","blue","purple"))
#par(mfrow = c(3, 2))
#plot_grid(p1, p2, p3, p4, p5, p6)
#plot_grid(p5, p6)

#FeaturePlot(mLN_Ontogeny.integrated, features.plot = c("Pdpn","Pecam1","Lyve1","Ackr4", 
 #                                                      "Tnfsf11","Tnfsf13b","Acta2","Ccl19",
  #                                                     "Madcam1","Ltbr","Tnfsf13b","Il6", 
   #                                                    "Il7","Ackr3","Fn1","Col4a1"),
    #        min.cutoff = "q9",
     #       cols.use = c('lightgrey', 'brown'),
      #      pt.size = 0.2)

#####
#Eliminate LECs, BECs
#####
#For each sample 
#list for storing LECs and BECs per Sample
#VlnPlot(object = mLN_Ontogeny.integrated, features.plot = c("Pecam1"))
mLN_Ontogeny.integrated_aExp <- AverageExpression(mLN_Ontogeny.integrated)
aExp_Pecam1 <- subset(mLN_Ontogeny.integrated_aExp, rownames(mLN_Ontogeny.integrated_aExp) == c("Pecam1"))
#Grab colnames/clusters that have Pecam1 expression < 1 on average per cluster
cluster_keep <- as.numeric(rownames(subset(as.data.frame(t(aExp_Pecam1)), Pecam1 < 1)))

#Get cells that are not within the cluster BEC and LEC cluster
all_cells <- mLN_Ontogeny.integrated@meta.data
all_cells_keep <- subset(all_cells, res.1.2 == 1000)

#attain all cells per cluster ID and concatenate the tables
for(i in 1:length(cluster_keep)){
  cluster_id_i <- cluster_keep[i]
  print(cluster_id_i)
  cells_cluster_i <- subset(all_cells, res.1.2 == cluster_id_i)
  print(nrow(cells_cluster_i))
  all_cells_keep <- rbind(all_cells_keep, cells_cluster_i)
}
mLN_Ontogeny.integrated_NO_LECs_BECs <- as.character(rownames(all_cells_keep))

#l_NO_LECs_BECs <- list()
#for(i in 1:length(ob.list)){
#  #Identify LECs and BECs
#  sample_i <- ob.list[[i]]
#  sample_i.integrated_aExp <- AverageExpression(sample_i)
#  aExp_Pecam1 <- subset(sample_i.integrated_aExp, rownames(sample_i.integrated_aExp) == c("Pecam1"))
#  #Grab colnames/clusters that have Pecam1 expression < 1 on average per cluster
#  cluster_keep <- as.numeric(rownames(subset(as.data.frame(t(aExp_Pecam1)), Pecam1 < 1)))
#  
#  #Get cells that are not within the cluster BEC and LEC cluster
#  all_cells <- sample_i@meta.data
#  all_cells_keep <- subset(all_cells, res.1.2 == 1000)
#  
#  #attain all cells per cluster ID and concatenate the tables
#  for(i in 1:length(cluster_keep)){
#    cluster_id_i <- cluster_keep[i]
#    print(cluster_id_i)
#    cells_cluster_i <- subset(all_cells, res.1.2 == cluster_id_i)
#    print(nrow(cells_cluster_i))
#    all_cells_keep <- rbind(all_cells_keep, cells_cluster_i)
#  }
#  sample_i.integrated_NO_LECs_BECs <- as.character(rownames(all_cells_keep))
#  l_NO_LECs_BECs[[i]] <- sample_i.integrated_NO_LECs_BECs
#}

##########################
#Analyse ONLY FSC
##########################
#####
#Eliminate cells of poor quality
#####
mito.genes <- grep(pattern = "^mt-", x = rownames(mLN_Ontogeny.integrated@data), value = TRUE)
percent.mito <- colSums(as.matrix(mLN_Ontogeny.integrated@raw.data[mito.genes, ]))/colSums(as.matrix(mLN_Ontogeny.integrated@raw.data))

#Identify cells with low nGene/nUmi ratio
ratio.nUMITOnGene <- mLN_Ontogeny.integrated@meta.data$nUMI / mLN_Ontogeny.integrated@meta.data$nGene
names(ratio.nUMITOnGene) <- rownames(mLN_Ontogeny.integrated@meta.data)

# AddMetaData adds columns to object@data.info and good to store QC stats
mLN_Ontogeny.integrated <- AddMetaData(object = mLN_Ontogeny.integrated, metadata = percent.mito, col.name = "percent.mito")
mLN_Ontogeny.integrated <- AddMetaData(object = mLN_Ontogeny.integrated, metadata = ratio.nUMITOnGene, col.name = "ratio.nUMITOnGene")

#Filter for mitochondiral and duplets and LECs and BECs out
# We filter out cells according unique gene counts 
# and percentage mitochondrial reads
mLN_Ontogeny.integrated_minus <- FilterCells(object = mLN_Ontogeny.integrated, subset.names = c("nGene", "ratio.nUMITOnGene", "percent.mito"), 
                                     low.thresholds = c(min_gene_number, -Inf, -Inf), high.thresholds = c(nGene, 4, mito_cutoff),
                                     cells.use = mLN_Ontogeny.integrated_NO_LECs_BECs)

mLN_Ontogeny.integrated_minus <- NormalizeData(mLN_Ontogeny.integrated_minus)
mLN_Ontogeny.integrated_minus <- FindVariableGenes(mLN_Ontogeny.integrated_minus, do.plot = F, display.progress = F)
mLN_Ontogeny.integrated_minus <- ScaleData(mLN_Ontogeny.integrated_minus)
#Run PCA
mLN_Ontogeny.integrated_minus <- RunPCA(object = mLN_Ontogeny.integrated_minus, pc.genes = mLN_Ontogeny.integrated_minus@var.genes, do.print = TRUE, pcs.print = 1:10, 
                              genes.print = 10)

#####
#Align without LECs & BECs per sample
#####
obs.list_minus <- list()
for(i in 1:length(ob.list)){
  #Normalize Samples and cluster
  sample_i <- ob.list[[i]]

  mito.genes <- grep(pattern = "^mt-", x = rownames(sample_i@data), value = TRUE)
  percent.mito <- colSums(as.matrix(sample_i@raw.data[mito.genes, ]))/colSums(as.matrix(sample_i@raw.data))
  
  #Identify cells with low nGene/nUmi ratio
  ratio.nUMITOnGene <- sample_i@meta.data$nUMI / sample_i@meta.data$nGene
  names(ratio.nUMITOnGene) <- rownames(sample_i@meta.data)
  
  # AddMetaData adds columns to object@data.info and good to store QC stats
  sample_i <- AddMetaData(object = sample_i, metadata = percent.mito, col.name = "percent.mito")
  sample_i <- AddMetaData(object = sample_i, metadata = ratio.nUMITOnGene, col.name = "ratio.nUMITOnGene")
  #Normalize the expression and log-transform
  sample_i <- NormalizeData(object = sample_i, normalization.method = "LogNormalize", 
                                   scale.factor = 10000)
  
  #Find variable genes independent of expression level
  sample_i <- FindVariableGenes(object = sample_i, mean.function = ExpMean, dispersion.function = LogVMR, 
                                       x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot = FALSE)
  
  #scale data
  sample_i <- ScaleData(object = sample_i)
  
  #perfrom PCA and print genes that define PCA
  sample_i <- RunPCA(object = sample_i, pc.genes = sample_i@var.genes, do.print = TRUE, pcs.print = 1:5, 
                            genes.print = 5)
  
  #ProjectPCA scores each gene in the dataset 
  sample_i <- ProjectPCA(object = sample_i, do.print = FALSE)

  #Group the cells into clusters
  sample_i <- FindClusters(object = sample_i, reduction.type = "pca", dims.use = 1:dims_use, 
                                  resolution = resolution_single, print.output = 0, save.SNN = TRUE,
                                  force.recalc = TRUE)
  
  #perfrom Tsne
  sample_i <- RunTSNE(object = sample_i, dims.use = 1:dims_use, do.fast = TRUE)
  
  #Identify LECs and BECs
  sample_i.integrated_aExp <- AverageExpression(sample_i)
  aExp_Pecam1 <- subset(sample_i.integrated_aExp, rownames(sample_i.integrated_aExp) == c("Pecam1"))
  #Grab colnames/clusters that have Pecam1 expression < 1 on average per cluster
  cluster_keep <- as.numeric(rownames(subset(as.data.frame(t(aExp_Pecam1)), Pecam1 < 1)))
  
  #Get cells that are not within the cluster BEC and LEC cluster
  all_cells <- sample_i@meta.data
  all_cells_keep <- subset(all_cells, res.1.2 == 1000)
  
  #attain all cells per cluster ID and concatenate the tables
  for(j in 1:length(cluster_keep)){
    cluster_id_i <- cluster_keep[j]
    print(cluster_id_i)
    cells_cluster_i <- subset(all_cells, res.1.2 == cluster_id_i)
    print(nrow(cells_cluster_i))
    all_cells_keep <- rbind(all_cells_keep, cells_cluster_i)
  }
  sample_i.integrated_NO_LECs_BECs <- as.character(rownames(all_cells_keep))
  
  #Prep samples without LECs and BECs
  sample_minus_i <- FilterCells(object = sample_i, subset.names = c("nGene", "ratio.nUMITOnGene", "percent.mito"), 
                                       low.thresholds = c(min_gene_number, -Inf, -Inf), high.thresholds = c(nGene, 4, mito_cutoff),
                                       cells.use = sample_i.integrated_NO_LECs_BECs)
  sample_minus_i <- FilterCells(object = sample_minus_i, subset.names = c("Ackr4","Pecam1"),
                                       low.thresholds = c(-Inf, -Inf), high.thresholds = c(1,1))
  sample_minus_i <- NormalizeData(sample_minus_i)
  sample_minus_i <- FindVariableGenes(sample_minus_i, do.plot = F, display.progress = F)
  sample_minus_i <- ScaleData(sample_minus_i)
  sample_minus_i@meta.data$tech <- v_sample[i]
  obs.list_minus[[i]] <- sample_minus_i
}
saveRDS(obs.list_minus, paste(PATH_Output,"/obs.list_minus_1000.Rds", sep=""))
obs.list_minus <- readRDS(paste(PATH_Output,"/obs.list_minus_1000.Rds", sep=""))
#####
#Aligen Datasets  ONLY FSCs
#####
#Merge the seurat objects
# First two set the object
mLN_Ontogeny_minus.merged <- MergeSeurat(obs.list_minus[[1]], obs.list_minus[[2]], do.normalize = FALSE)
for(i in 1:(length(obs.list_minus)-2)){
  k = i + 2
  print(k)
  mLN_Ontogeny_minus.merged <- MergeSeurat(mLN_Ontogeny_minus.merged, obs.list_minus[[k]], do.normalize = FALSE)
}

#Normalize
mLN_Ontogeny_minus.merged <- NormalizeData(mLN_Ontogeny_minus.merged)
#Find variable genes
mLN_Ontogeny_minus.merged <- FindVariableGenes(mLN_Ontogeny_minus.merged, do.plot = F, display.progress = F)

#####
#Scale data and regress parameters
#####
# 1) Cell Cycle
#To lower case tranlating into mouse
Cell_cycle_Regev_BLOOD_minor <- tolower(Cell_cycle_Regev_BLOOD[,1])
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
Cell_cycle_genes_mouse <- firstup(Cell_cycle_Regev_BLOOD_minor)

#Phase specific
s.genes <- Cell_cycle_genes_mouse[1:43]
g2m.genes <- Cell_cycle_genes_mouse[44:97]

#Score each cell
mLN_Ontogeny_minus.merged <- CellCycleScoring(object = mLN_Ontogeny_minus.merged, s.genes = s.genes, g2m.genes = g2m.genes, 
                                              set.ident = TRUE)

# Visualize across states
#RidgePlot(object = mLN_Ontogeny_minus.merged, features.plot = c("Pcna", "Top2a", "Mcm6", "Mki67"), 
 #         nCol = 2)

# 2) Mitochondria
mito.genes <- grep(pattern = "^mt-", x = rownames(mLN_Ontogeny_minus.merged@data), value = TRUE)
percent.mito <- colSums(as.matrix(mLN_Ontogeny_minus.merged@raw.data[mito.genes, ]))/colSums(as.matrix(mLN_Ontogeny_minus.merged@raw.data))
# AddMetaData adds columns to object@data.info and good to store QC stats
mLN_Ontogeny_minus.merged <- AddMetaData(object = mLN_Ontogeny_minus.merged, metadata = percent.mito, col.name = "percent.mito")

# 3) Ribosomal
#Calculate % ribosomal genes
rpl.genes <- grep(pattern = "^Rpl", x = rownames(mLN_Ontogeny_minus.merged@data), value = TRUE)
rps.genes <- grep(pattern = "^Rps", x = rownames(mLN_Ontogeny_minus.merged@data), value = TRUE)
ribo.genes <- c(rpl.genes, rps.genes)
percent.ribo <- colSums(as.matrix(mLN_Ontogeny_minus.merged@raw.data[ribo.genes, ]))/colSums(as.matrix(mLN_Ontogeny_minus.merged@raw.data))
mLN_Ontogeny_minus.merged <- AddMetaData(object = mLN_Ontogeny_minus.merged, metadata = percent.ribo, col.name = "percent.ribo")

#Scale
mLN_Ontogeny_minus.merged <- ScaleData(object = mLN_Ontogeny_minus.merged, vars.to.regress = c("S.Score", "G2M.Score","nGene", "percent.mito", "percent.ribo"), 
                                       display.progress = FALSE)

#Run PCA
mLN_Ontogeny_minus.merged <- RunPCA(object = mLN_Ontogeny_minus.merged, pc.genes = mLN_Ontogeny_minus.merged@var.genes, do.print = TRUE, pcs.print = 1:10, 
                              genes.print = 10)


# Determine genes to use for CCA, must be highly variable in at least 2 datasets
genes.use <- c()
for (i in 1:length(obs.list_minus)) {
  genes.use <- c(genes.use, head(rownames(obs.list_minus[[i]]@hvg.info), 1500))
  print(names(obs.list_minus)[i])
  print(length(genes.use))
}
genes.use <- names(which(table(genes.use) > 1))
for (i in 1:length(obs.list_minus)) {
  genes.use <- genes.use[genes.use %in% rownames(obs.list_minus[[i]]@scale.data)]
  print(length(genes.use))
}

# Run multi-set CCA
mLN_Ontogeny_minus.merged <- RunMultiCCA(obs.list_minus, genes.use = genes.use, num.ccs = dims_use)

# CC Selection
MetageneBicorPlot(mLN_Ontogeny_minus.merged, grouping.var = "tech", dims.eval = 1:12)

# Run rare non-overlapping filtering
mLN_Ontogeny_minus.merged <- CalcVarExpRatio(object = mLN_Ontogeny_minus.merged, reduction.type = "pca",
                                           grouping.var = "tech", dims.use = 1:dims_use)
mLN_Ontogeny_minus.merged <- SubsetData(mLN_Ontogeny_minus.merged, subset.name = "var.ratio.pca",
                                      accept.low = 0.5)

# Alignment
mLN_Ontogeny_minus.merged <- AlignSubspace(mLN_Ontogeny_minus.merged,
                                         reduction.type = "cca",
                                         grouping.var = "tech",
                                         dims.align = 1:dims_use)

# t-SNE and Clustering
mLN_Ontogeny_minus.merged <- RunTSNE(mLN_Ontogeny_minus.merged,
                                   reduction.use = "cca",
                                   dims.use = 1:dims_use, save.SNN = T,
                                   resolution = 1.2, force.recalc = TRUE)
mLN_Ontogeny_minus.merged <- FindClusters(mLN_Ontogeny_minus.merged, reduction.type = "cca",
                                          dims.use = 1:dims_use, save.SNN = T, resolution = resolution_merge,
                                          force.recalc = TRUE)

saveRDS(mLN_Ontogeny_minus.merged, paste(PATH_Output,"/MultiAlign_12_1_1500_regressed_1000genes_D0_D10_D24_D56_D300.Rds", sep=""))

mLN_Ontogeny_minus.merged <- readRDS(paste(PATH_Output,"/MultiAlign_12_1_1500_regressed_1000genes_D0_D10_D24_D56_D300.Rds", sep=""))
#Change cluster ID
current.cluster.ids <- c(0:13)

new.cluster.ids <- c("feeder1","Inmt+","feeder2","CD34+CD248+Ackr3+Has1+",
                     "Ccl19+Il7+","Il6+Cxcl1+","feeder3","PvC",
                     "LTolike","CD34+Aldh1a2+","SCx","CD34+Cdk1+",
                     "pSC","CD34+Gdf10+")

mLN_Ontogeny_minus.merged@ident <- plyr::mapvalues(x = mLN_Ontogeny_minus.merged@ident, from = current.cluster.ids, to = new.cluster.ids)
mLN_Ontogeny_minus.merged@meta.data$Cluster <- mLN_Ontogeny_minus.merged@ident

#Check allocation of cells
sample_1_cells <- grep(pattern = "^D000_", x = rownames(mLN_Ontogeny_minus.merged@meta.data), value = TRUE)
sample_2_cells <- grep(pattern = "^D010_", x = rownames(mLN_Ontogeny_minus.merged@meta.data), value = TRUE)
sample_3_cells <- grep(pattern = "^D024_", x = rownames(mLN_Ontogeny_minus.merged@meta.data), value = TRUE)
sample_4_cells <- grep(pattern = "^D056_", x = rownames(mLN_Ontogeny_minus.merged@meta.data), value = TRUE)
sample_5_cells <- grep(pattern = "^D300_", x = rownames(mLN_Ontogeny_minus.merged@meta.data), value = TRUE)

p1 <- TSNEPlot(mLN_Ontogeny_minus.merged, group.by = "tech", do.return = T, pt.size = 0.4, cells.use = sample_1_cells)
p2 <- TSNEPlot(mLN_Ontogeny_minus.merged, group.by = "tech", do.return = T, pt.size = 0.4, cells.use = sample_2_cells)
p3 <- TSNEPlot(mLN_Ontogeny_minus.merged, group.by = "tech", do.return = T, pt.size = 0.4, cells.use = sample_3_cells)
p4 <- TSNEPlot(mLN_Ontogeny_minus.merged, group.by = "tech", do.return = T, pt.size = 0.4, cells.use = sample_4_cells)
p5 <- TSNEPlot(mLN_Ontogeny_minus.merged, group.by = "tech", do.return = T, pt.size = 0.4, cells.use = sample_5_cells)
p6 <- TSNEPlot(mLN_Ontogeny_minus.merged, do.return = T, pt.size = 0.2, do.label = TRUE)
#test <- LN@meta.data[order(-LN@meta.data$),]
p7 <- TSNEPlot(mLN_Ontogeny_minus.merged, group.by = "tech", do.return = T, pt.size = 0.2,
               colors.use = c("orange","green","blue","purple","black"))
par(mfrow = c(4,2))
plot_grid(p1, p2, p3, p4, p5, p6, p7)
#, p6, p7)
plot_grid(p7)

FeaturePlot(mLN_Ontogeny_minus.merged, features.plot = c("Tnfsf11","Tnfsf13b","Acta2","Ccl19",
                                                       "Madcam1","Ltbr","Tnfsf13b","Il6", 
                                                       "Il7","Ackr3","Fn1"), #"Col4a1", "Cxcl13",
                                                      # "Icam1","Cd34","Gdf10","Ptgis","Aldh1a2",
                                                       #"Tinagl1","Cdk1"),
            min.cutoff = "q9",
            cols.use = c('lightgrey', 'brown'),
            pt.size = 0.2)

FeaturePlot(mLN_Ontogeny_minus.merged, features.plot = c("Col4a1", "Cxcl13",
                                                         "Icam1","Cd34","Gdf10","Ptgis","Aldh1a2",
                                                         "Tinagl1","Cdk1"),
            min.cutoff = "q9",
            cols.use = c('lightgrey', 'brown'),
            pt.size = 0.2)

FeaturePlot(mLN_Ontogeny_minus.merged, features.plot = c("Igf1","Igf2",
                                                         "Igfbp2","Igfbp4","Igfbp5","Igfbp6","Igfbp7",
                                                         "Igf1r","Igf2r"),
            min.cutoff = "q9",
            cols.use = c('lightgrey', 'brown'),
            pt.size = 0.2)
FeaturePlot(mLN_Ontogeny_minus.merged, features.plot = c("Cdk1","Cxcl13","Cd34","Acta2"),
            min.cutoff = "q9",
            cols.use = c('lightgrey', 'brown'),
            pt.size = 0.2)

#saveRDS(mLN_Ontogeny_minus.merged, paste(PATH_Output,"/MultiAlign_D0_D10_D24_D56.Rds", sep=""))
#mLN_Ontogeny_minus.merged <- readRDS(paste(PATH_Output,"/MultiAlign_20_1_regressed_1000genes_D0_D10_D24_D56_D300.Rds", sep=""))
#saveRDS(mLN_Ontogeny_minus.merged, paste(PATH_Output,"/MultiAlign_D0_D10_D24_D56_D300.Rds", sep=""))

#####
#Obtain cell IDs minus
#####
#Define clusters to be excluded
unique_clusters <- c("PvC","feeder1","feeder3")
All_meta <- subset(mLN_Ontogeny_minus.merged@meta.data)
number_of_cells = 0
for(i in 1:length(v_sample)){
  v_sample_i <- v_sample[i]
  print(v_sample_i)
  day_i_meta <- subset(All_meta, tech == v_sample_i)
  print(nrow(day_i_meta))
  cells_include_i <- row.names(subset(day_i_meta, !(Cluster %in% unique_clusters)))
  print(length(cells_include_i))
  number_of_cells = number_of_cells + length(cells_include_i)
  saveRDS(cells_include_i, paste(PATH_Output,"/",v_sample_i,"_1500_1000_1_12_SC_minus_feeder1_3_PvC.Rds", sep=""))
}
print(number_of_cells)
#Write Metadata
saveRDS(All_meta, paste(PATH_Output,"/","MetaData_1500_1000_1_12_SC_minus_feeder1_3_PvC.Rds", sep=""))

#####
#Cluster identification via signature
#####
#Use list of DEGs from merged analysis of mLN and pLN multiple samples

DEGs_core_40 <- NC_Pezoldt_Pasztoi_2018

#Load DEG tables and make list of clusters
DEGs_list <- split(DEGs_core_40, DEGs_core_40$cluster)
#Average gene expression 
mLN_Ontogeny_minus.merged_minus_aExp <- AverageExpression(mLN_Ontogeny_minus.merged)

#calculate cZscore for each signature across all clusters
#number of rows for output
i = 1
out = NULL
for(i in i:length(DEGs_list)){
  cluster_DEG_i <- as.character(DEGs_list[[i]]$gene)
  cluster_DEG_n_i <- length(cluster_DEG_i)
  
  cluster_ID_i <- names(DEGs_list)[i]
  y_label_i <- paste("cZ: ", cluster_ID_i)
  
  print(cluster_ID_i)
  print(i)
  print(cluster_DEG_n_i)
  
  mLN_Ontogeny_minus.merged_minus_aExp_cluster_i <- subset(mLN_Ontogeny_minus.merged_minus_aExp, rownames(mLN_Ontogeny_minus.merged_minus_aExp) %in% cluster_DEG_i)
  
  Scale_mLN_Ontogeny_minus.merged_minus_aExp_cluster_i <- apply(mLN_Ontogeny_minus.merged_minus_aExp_cluster_i, 1, scale)
  Zscore_cluster_i <- rowSums(Scale_mLN_Ontogeny_minus.merged_minus_aExp_cluster_i)
  Zscore_cluster_i_norm <- Zscore_cluster_i / (nrow(mLN_Ontogeny_minus.merged_minus_aExp_cluster_i)/10)
  
  #mLN maintained Scores
  Zscore_table_cluster_i <- data.frame(module = c(rep(cluster_ID_i,length(Zscore_cluster_i_norm))),
                                       cluster = colnames(mLN_Ontogeny_minus.merged_minus_aExp),
                                       cZscore = Zscore_cluster_i_norm)
  
  
  #Generate cZscore bargraph matrix
  p <- ggplot(data = Zscore_table_cluster_i, aes(x = cluster, y = cZscore, fill = module)) +
    geom_bar(stat = "identity", colour = "black") +
    theme(axis.text.x=element_text(angle=270,hjust=1), legend.position="none") +
    scale_fill_manual(values = c("deepskyblue1")) +
    ylab(y_label_i)
  #print(p)
  out[[i]] <- p
}

title <- paste("Ontogeny_MultiAlign_DEGs_core_40",sep = "_")
pdf(paste(PATH_Output,"/MultiAlign_12_1_1000_regressed_all_NC_Top40_D0_D10_D24_D56_D300.pdf", sep=""))
marrangeGrob(out, nrow = round(length(DEGs_list) / 3 - 1), ncol = 3, top = title)
dev.off()



#####
#Find differentially expressed genes
#####
#Find markers for every cluster compared to all remaining cells
mLN_Ontogeny_minus.merged.markers <- FindAllMarkers(object = mLN_Ontogeny_minus.merged, only.pos = TRUE, min.pct = 0.25, 
                                                    thresh.use = 0.25)
mLN_Ontogeny_minus.merged.markers_20 <- mLN_Ontogeny_minus.merged.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
mLN_Ontogeny_minus.merged.markers_40 <- mLN_Ontogeny_minus.merged.markers %>% group_by(cluster) %>% top_n(40, avg_logFC)
mLN_Ontogeny_minus.merged.markers_20 <- mLN_Ontogeny_minus.merged.markers_20$gene
#Order LN_minus meta.data according to cluster size

#Save DEGs
write.table(mLN_Ontogeny_minus.merged.markers, paste(PATH_Output,"/MultiAlign_12_1_regressed_1000genes_DEGs_D0_D10_D24_D56_D300.csv", sep=""), dec=".", sep=",")
mLN_Ontogeny_minus.merged.markers <- read.table(paste(PATH_Output,"/MultiAlign_12_1_regressed_1000genes_DEGs_D0_D10_D24_D56_D300.csv", sep=""), dec=".", sep=",")


#mLN_Ontogeny_minus.merged.markers <- read.table(paste(PATH_Output,"/MultiAlign_20_1.2_regressed_500genes_DEGs_D0_D10_D24_D56_D300.csv", sep=""), dec=".", sep=",")
#Plot differential expression of Marker genes
DoHeatmap(object = mLN_Ontogeny_minus.merged, use.scaled = TRUE, group.by = "ident",
          genes.use = mLN_Ontogeny_minus.merged.markers_20, title = "Merged FSC Ontogeny",
          remove.key = FALSE,
          col.low = "white",
          col.mid =  "oldlace",
          col.high = "brown")

mLN_Ontogeny_minus.merged.markers <- read.table(paste(PATH_Output,"/DEGs_D0_D10_D24_D56_D300.csv", sep=""), dec=".", sep=",")

#####
#Dotplot
#####

mLN_Ontogeny_minus.merged.markers_nonDup <- mLN_Ontogeny_minus.merged.markers[!duplicated(mLN_Ontogeny_minus.merged.markers$gene),]
mLN_Ontogeny_minus.merged.markers_nonDup_10_pValue <- mLN_Ontogeny_minus.merged.markers_nonDup %>% group_by(cluster) %>% top_n(10, 1/p_val_adj)
mLN_Ontogeny_minus.merged.markers_nonDup_20_pValue <- mLN_Ontogeny_minus.merged.markers_nonDup %>% group_by(cluster) %>% top_n(20, 1/p_val_adj)
mLN_Ontogeny_minus.merged.markers_nonDup_10_FC <- mLN_Ontogeny_minus.merged.markers_nonDup %>% group_by(cluster) %>% top_n(8, avg_logFC)
mLN_Ontogeny_minus.merged.markers_nonDup_20_FC <- mLN_Ontogeny_minus.merged.markers_nonDup %>% group_by(cluster) %>% top_n(20, avg_logFC)

#Plot differential expression of Marker genes
title <- "Dim20_Res1_100_Regressed_Dotplot_Top10_D0_10_24_56_300"
#pdf(paste(path_output,"/",title, ".pdf", sep = ""), width = 25, height = 10)
DotPlot(mLN_Ontogeny_minus.merged, mLN_Ontogeny_minus.merged.markers_nonDup_10_FC$gene,
        x.lab.rot = TRUE, plot.legend = FALSE,
        col.min = 0, col.max = 5, dot.scale = 5,
        cols.use = c("lightgrey", "brown"))
#dev.off()
#dev.off()
#####

#####
#Hierarchical clustering according to Top DEGs per cluster
#####
#Average expression
mLN_Ontogeny_minus.merged_aExp <- AverageExpression(mLN_Ontogeny_minus.merged)
#Check for DEGs
#genes <- as.character(LN_minus.markers_20[1:260,]$gene)
mLN_Ontogeny_minus.merged_aExp_TopDEGs <- subset(mLN_Ontogeny_minus.merged_aExp, rownames(mLN_Ontogeny_minus.merged_aExp) %in% mLN_Ontogeny_minus.merged.markers$gene)
mLN_Ontogeny_minus.merged_aExp_TopDEGs <- as.matrix(mLN_Ontogeny_minus.merged_aExp_TopDEGs)

title <- paste("Heatmap_hierarchical",sep = "_")
#pdf(paste(path_output,"/",title, ".pdf", sep = ""), width = 6, height = 7)
pheatmap(mLN_Ontogeny_minus.merged_aExp_TopDEGs, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 30, show_rownames = FALSE, cluster_cols = TRUE,
         scale = "row", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("blue", "white","red"), space="rgb")(128),
         main = title)
#dev.off()


#####
#Calculate % per cluster across age
#####
meta_data <- mLN_Ontogeny_minus.merged@meta.data
l_condition <- split(meta_data, meta_data$tech)
l_cells_per_condition <- lapply(l_condition, function(x) {
  #x <- l_condition[[1]]
  x_cells <- unlist(lapply(split(x, x$Cluster), nrow))
  x_cells
})
names(l_cells_per_condition) <- names(l_condition)
t_cells_per_condition <- do.call(rbind.fill.matrix, lapply(l_cells_per_condition, function(x) {t(as.data.frame(x))}))
t_cells_per_condition[is.na(t_cells_per_condition)] <- 0
rownames(t_cells_per_condition) <- names(l_condition)
#Normalize to 1000 cells
norm_factor_vector <- rowSums(t_cells_per_condition) / 1000
#Divide cell number per cluster by normalization factor
t_cells_per_condition_norm <- t_cells_per_condition / norm_factor_vector
#Calculate Frequencies
total <- colSums(t_cells_per_condition_norm)
freq_cells_per_condition <- (t(t_cells_per_condition_norm) / total) * 100
freq_cells_per_condition <- as.matrix(freq_cells_per_condition)
pheatmap(freq_cells_per_condition, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = TRUE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("white","grey","grey23","black"), space="rgb")(128))



#####################
#Run analysis on pre-selected subsets (non-feeder and non-PvC)
#####################
#####
#Eliminate feeder cells MAIN
#####
unique_clusters <- c("PvC","feeder1","feeder2","feeder3")
All_meta <- subset(mLN_Ontogeny_minus.merged@meta.data)
number_of_cells = 0
l_cells_include_CCA <- list()
for(i in 1:length(v_sample)){
  v_sample_i <- v_sample[i]
  print(v_sample_i)
  day_i_meta <- subset(All_meta, tech == v_sample_i)
  print(nrow(day_i_meta))
  cells_include_i <- row.names(subset(day_i_meta, !(Cluster %in% unique_clusters)))
  print(length(cells_include_i))
  number_of_cells = number_of_cells + length(cells_include_i)
  l_cells_include_CCA[[i]] <- cells_include_i
  saveRDS(cells_include_i, paste(PATH_Output,"/",v_sample_i,"_1500_1000_1_12_SC_minus_feeder1_2_3_PvC.Rds", sep=""))
}
print(number_of_cells)
#####
#Align with MAIN cells
#####
obs.list_minus_main <- list()
for(i in 1:length(ob.list)){
  #Normalize Samples and cluster
  sample_i <- ob.list[[i]]
  
  mito.genes <- grep(pattern = "^mt-", x = rownames(sample_i@data), value = TRUE)
  percent.mito <- colSums(as.matrix(sample_i@raw.data[mito.genes, ]))/colSums(as.matrix(sample_i@raw.data))
  
  #Identify cells with low nGene/nUmi ratio
  ratio.nUMITOnGene <- sample_i@meta.data$nUMI / sample_i@meta.data$nGene
  names(ratio.nUMITOnGene) <- rownames(sample_i@meta.data)
  
  # AddMetaData adds columns to object@data.info and good to store QC stats
  sample_i <- AddMetaData(object = sample_i, metadata = percent.mito, col.name = "percent.mito")
  sample_i <- AddMetaData(object = sample_i, metadata = ratio.nUMITOnGene, col.name = "ratio.nUMITOnGene")
  #Normalize the expression and log-transform
  sample_i <- NormalizeData(object = sample_i, normalization.method = "LogNormalize", 
                            scale.factor = 10000)
  
  #Find variable genes independent of expression level
  sample_i <- FindVariableGenes(object = sample_i, mean.function = ExpMean, dispersion.function = LogVMR, 
                                x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot = FALSE)
  
  #scale data
  sample_i <- ScaleData(object = sample_i)
  
  #perfrom PCA and print genes that define PCA
  sample_i <- RunPCA(object = sample_i, pc.genes = sample_i@var.genes, do.print = TRUE, pcs.print = 1:5, 
                     genes.print = 5)
  
  #ProjectPCA scores each gene in the dataset 
  sample_i <- ProjectPCA(object = sample_i, do.print = FALSE)
  
  #Group the cells into clusters
  sample_i <- FindClusters(object = sample_i, reduction.type = "pca", dims.use = 1:dims_use, 
                           resolution = resolution_single, print.output = 0, save.SNN = TRUE,
                           force.recalc = TRUE)
  
  #perfrom Tsne
  sample_i <- RunTSNE(object = sample_i, dims.use = 1:dims_use, do.fast = TRUE)
  
  #Cells to use
  sample_i.integrated_cells_chosen <- l_cells_include_CCA[[i]]
  
  #Prep samples without LECs and BECs
  sample_minus_i <- FilterCells(object = sample_i, subset.names = c("nGene", "ratio.nUMITOnGene", "percent.mito"), 
                                low.thresholds = c(min_gene_number, -Inf, -Inf), high.thresholds = c(nGene, 4, mito_cutoff),
                                cells.use = sample_i.integrated_cells_chosen)
  sample_minus_i <- NormalizeData(sample_minus_i)
  sample_minus_i <- FindVariableGenes(sample_minus_i, do.plot = F, display.progress = F)
  sample_minus_i <- ScaleData(sample_minus_i)
  sample_minus_i@meta.data$tech <- v_sample[i]
  obs.list_minus_main[[i]] <- sample_minus_i
}
saveRDS(obs.list_minus_main, paste(PATH_Output,"/obs.list_minus_main_minus_feeder123_PvC_1000.Rds", sep=""))
obs.list_minus_main <- readRDS(paste(PATH_Output,"/obs.list_minus_main_minus_feeder123_PvC_1000.Rds", sep=""))

#Merge the seurat objects
# First two set the object
obs.mLN_Ontogeny_minus_main.merged <- MergeSeurat(obs.list_minus_main[[1]], obs.list_minus_main[[2]], do.normalize = FALSE)
for(i in 1:(length(obs.list_minus_main)-2)){
  k = i + 2
  print(k)
  obs.mLN_Ontogeny_minus_main.merged <- MergeSeurat(obs.mLN_Ontogeny_minus_main.merged, obs.list_minus_main[[k]], do.normalize = FALSE)
}

#Normalize
obs.mLN_Ontogeny_minus_main.merged <- NormalizeData(obs.mLN_Ontogeny_minus_main.merged)
#Find variable genes
obs.mLN_Ontogeny_minus_main.merged <- FindVariableGenes(obs.mLN_Ontogeny_minus_main.merged, do.plot = F, display.progress = F)

#####
#Scale data and regress parameters MAIN
#####
# 1) Cell Cycle
#To lower case tranlating into mouse
Cell_cycle_Regev_BLOOD_minor <- tolower(Cell_cycle_Regev_BLOOD[,1])
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
Cell_cycle_genes_mouse <- firstup(Cell_cycle_Regev_BLOOD_minor)

#Phase specific
s.genes <- Cell_cycle_genes_mouse[1:43]
g2m.genes <- Cell_cycle_genes_mouse[44:97]

#Score each cell
obs.mLN_Ontogeny_minus_main.merged <- CellCycleScoring(object = obs.mLN_Ontogeny_minus_main.merged, s.genes = s.genes, g2m.genes = g2m.genes, 
                                              set.ident = TRUE)

# 2) Mitochondria
mito.genes <- grep(pattern = "^mt-", x = rownames(obs.mLN_Ontogeny_minus_main.merged@data), value = TRUE)
percent.mito <- colSums(as.matrix(obs.mLN_Ontogeny_minus_main.merged@raw.data[mito.genes, ]))/colSums(as.matrix(obs.mLN_Ontogeny_minus_main.merged@raw.data))
# AddMetaData adds columns to object@data.info and good to store QC stats
obs.mLN_Ontogeny_minus_main.merged <- AddMetaData(object = obs.mLN_Ontogeny_minus_main.merged, metadata = percent.mito, col.name = "percent.mito")

# 3) Ribosomal
#Calculate % ribosomal genes
rpl.genes <- grep(pattern = "^Rpl", x = rownames(obs.mLN_Ontogeny_minus_main.merged@data), value = TRUE)
rps.genes <- grep(pattern = "^Rps", x = rownames(obs.mLN_Ontogeny_minus_main.merged@data), value = TRUE)
ribo.genes <- c(rpl.genes, rps.genes)
percent.ribo <- colSums(as.matrix(obs.mLN_Ontogeny_minus_main.merged@raw.data[ribo.genes, ]))/colSums(as.matrix(obs.mLN_Ontogeny_minus_main.merged@raw.data))
obs.mLN_Ontogeny_minus_main.merged <- AddMetaData(object = obs.mLN_Ontogeny_minus_main.merged, metadata = percent.ribo, col.name = "percent.ribo")

#Scale
obs.mLN_Ontogeny_minus_main.merged <- ScaleData(object = obs.mLN_Ontogeny_minus_main.merged, vars.to.regress = c("S.Score", "G2M.Score","nGene", "percent.mito", "percent.ribo"), 
                                       display.progress = FALSE)

#Run PCA
obs.mLN_Ontogeny_minus_main.merged <- RunPCA(object = obs.mLN_Ontogeny_minus_main.merged, pc.genes = obs.mLN_Ontogeny_minus_main.merged@var.genes, do.print = TRUE, pcs.print = 1:10, 
                                    genes.print = 10)


# Determine genes to use for CCA, must be highly variable in at least 2 datasets
genes.use <- c()
for (i in 1:length(obs.list_minus_main)) {
  genes.use <- c(genes.use, head(rownames(obs.list_minus_main[[i]]@hvg.info), 1500))
  print(names(obs.list_minus_main)[i])
  print(length(genes.use))
}
genes.use <- names(which(table(genes.use) > 1))
for (i in 1:length(obs.list_minus_main)) {
  genes.use <- genes.use[genes.use %in% rownames(obs.list_minus_main[[i]]@scale.data)]
  print(length(genes.use))
}

# Run multi-set CCA
obs.mLN_Ontogeny_minus_main.merged <- RunMultiCCA(obs.list_minus_main, genes.use = genes.use, num.ccs = dims_use)

# CC Selection
#MetageneBicorPlot(obs.mLN_Ontogeny_minus_main.merged, grouping.var = "tech", dims.eval = 1:12)

# Run rare non-overlapping filtering
obs.mLN_Ontogeny_minus_main.merged <- CalcVarExpRatio(object = obs.mLN_Ontogeny_minus_main.merged, reduction.type = "pca",
                                             grouping.var = "tech", dims.use = 1:dims_use)
obs.mLN_Ontogeny_minus_main.merged <- SubsetData(obs.mLN_Ontogeny_minus_main.merged, subset.name = "var.ratio.pca",
                                        accept.low = 0.5)

# Alignment
obs.mLN_Ontogeny_minus_main.merged <- AlignSubspace(obs.mLN_Ontogeny_minus_main.merged,
                                           reduction.type = "cca",
                                           grouping.var = "tech",
                                           dims.align = 1:dims_use)

# t-SNE and Clustering
obs.mLN_Ontogeny_minus_main.merged <- RunTSNE(obs.mLN_Ontogeny_minus_main.merged,
                                     reduction.use = "cca",
                                     dims.use = 1:dims_use, save.SNN = T,
                                     resolution = 1.2, force.recalc = TRUE)
obs.mLN_Ontogeny_minus_main.merged <- FindClusters(obs.mLN_Ontogeny_minus_main.merged, reduction.type = "cca",
                                          dims.use = 1:dims_use, save.SNN = T, resolution = resolution_merge_MAIN,
                                          force.recalc = TRUE)

saveRDS(obs.mLN_Ontogeny_minus_main.merged, paste(PATH_Output,"/MultiAlign_12_1_1500_regressed_1000genes_minus_feeder123_PvC_D0_D10_D24_D56_D300.Rds", sep=""))
obs.mLN_Ontogeny_minus_main.merged <- readRDS(paste(PATH_Output,"/MultiAlign_12_1_1500_regressed_1000genes_minus_feeder123_PvC_D0_D10_D24_D56_D300.Rds", sep=""))
#obs.mLN_Ontogeny_minus_main.merged <- readRDS(paste(PATH_Output,"/MultiAlign_12_1_1500_regressed_1000genes_MAIN_D0_D10_D24_D56_D300.Rds", sep=""))
#Change cluster ID
current.cluster.ids <- c(0:11)

new.cluster.ids <- c("Ccl19+Il7+","Il6+Cxcl1+","CD34+Ackr3+","Inmt+",
                     "Inmt+Cxcl12+","LTolike","CD34+Aldh1a2+","Cdk1+",
                     "CD34+CD248+","SCx","Cxcl9+","pSC")

obs.mLN_Ontogeny_minus_main.merged@ident <- plyr::mapvalues(x = obs.mLN_Ontogeny_minus_main.merged@ident, from = current.cluster.ids, to = new.cluster.ids)
obs.mLN_Ontogeny_minus_main.merged@meta.data$Cluster <- obs.mLN_Ontogeny_minus_main.merged@ident

#Add Age score to metadata
#age_score <- list()
#for(i in 1:length(v_sample)){
#  age_i <- v_sample[i]
#  binary_i <- apply(obs.mLN_Ontogeny_minus_main.merged@meta.data, 1, function(r) any(r %in% c(age_i)))
#  binary_numeric_i =ifelse(binary_i=="TRUE",10,0)
#  age_score[[i]] <- binary_numeric_i
#}
#age_score <- do.call("cbind", age_score)
#colnames(age_score) <- v_sample
#append score to metadata
#meta_data_temp <- cbind(obs.mLN_Ontogeny_minus_main.merged@meta.data, age_score)
#obs.mLN_Ontogeny_minus_main.merged@meta.data <- meta_data_temp
#Age plotting
#FeaturePlot(obs.mLN_Ontogeny_minus_main.merged, features.plot = c("d056"), 
 #           nCol = 1, cols.use = c("lightgrey", "black"), pt.size = 0.5,
  #          cells.use = c(head(sample_1_cells,760),head(sample_2_cells,760),head(sample_3_cells,760),head(sample_4_cells,760),head(sample_5_cells,760)))


#Check allocation of cells
sample_1_cells <- grep(pattern = "^D000_", x = rownames(obs.mLN_Ontogeny_minus_main.merged@meta.data), value = TRUE)
sample_2_cells <- grep(pattern = "^D010_", x = rownames(obs.mLN_Ontogeny_minus_main.merged@meta.data), value = TRUE)
sample_3_cells <- grep(pattern = "^D024_", x = rownames(obs.mLN_Ontogeny_minus_main.merged@meta.data), value = TRUE)
sample_4_cells <- grep(pattern = "^D056_", x = rownames(obs.mLN_Ontogeny_minus_main.merged@meta.data), value = TRUE)
sample_5_cells <- grep(pattern = "^D300_", x = rownames(obs.mLN_Ontogeny_minus_main.merged@meta.data), value = TRUE)

p1 <- TSNEPlot(obs.mLN_Ontogeny_minus_main.merged, group.by = "tech",
               do.return = T, pt.size = 0.1,
               cells.use = head(sample_1_cells,760),
               colors.use = "grey30")
p2 <- TSNEPlot(obs.mLN_Ontogeny_minus_main.merged, group.by = "tech",
               do.return = T, pt.size = 0.1,
               cells.use = head(sample_2_cells,760),
               colors.use = "grey30")
p3 <- TSNEPlot(obs.mLN_Ontogeny_minus_main.merged, group.by = "tech",
               do.return = T, pt.size = 0.1,
               cells.use = head(sample_3_cells,760),
               colors.use = "grey30")
p4 <- TSNEPlot(obs.mLN_Ontogeny_minus_main.merged, group.by = "tech",
               do.return = T, pt.size = 0.1,
               cells.use = head(sample_4_cells,760),
               colors.use = "grey30")
p5 <- TSNEPlot(obs.mLN_Ontogeny_minus_main.merged, group.by = "tech",
               do.return = T, pt.size = 0.1,
               cells.use = head(sample_5_cells,760),
               colors.use = "grey30")
p6 <- TSNEPlot(obs.mLN_Ontogeny_minus_main.merged, do.return = T, pt.size = 0.2, do.label = TRUE)
p7 <- TSNEPlot(obs.mLN_Ontogeny_minus_main.merged, group.by = "tech", do.return = T, pt.size = 0.2,
               colors.use = c("goldenrod2","darkorange3","deepskyblue3","grey60","black"),
               cells.use = c(head(sample_1_cells,760),head(sample_2_cells,760),head(sample_3_cells,760),head(sample_4_cells,760),head(sample_5_cells,760)))
plot_grid(p7)
par(mfrow = c(4,2))
plot_grid(p1, p2, p3, p4, p5, p6, p7)
#Grey thumbnail
par(mfrow = c(5,1))
plot_grid(p1, p2, p3, p4, p5)
#Cell cycle plot
#To lower case tranlating into mouse
Cell_cycle_Regev_BLOOD_minor <- tolower(Cell_cycle_Regev_BLOOD[,1])
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
Cell_cycle_genes_mouse <- firstup(Cell_cycle_Regev_BLOOD_minor)

#Phase specific
s.genes <- Cell_cycle_genes_mouse[1:43]
g2m.genes <- Cell_cycle_genes_mouse[44:97]

#Score each cell
obs.mLN_Ontogeny_minus_main.merged <- CellCycleScoring(object = obs.mLN_Ontogeny_minus_main.merged,
                                                       s.genes = s.genes, g2m.genes = g2m.genes, 
                                                       set.ident = TRUE)
FeaturePlot(obs.mLN_Ontogeny_minus_main.merged, features.plot = c("S.Score", "G2M.Score"), min.cutoff = "q05", max.cutoff = "q95", 
            nCol = 2, cols.use = c("lightgrey", "blue"), pt.size = 0.5)


FeaturePlot(obs.mLN_Ontogeny_minus_main.merged, features.plot = c("Tnfsf11","Tnfsf13b","Acta2","Ccl19",
                                                         "Madcam1","Ltbr","Tnfsf13b","Il6", 
                                                         "Il7","Ackr3","Fn1","Col4a1", "Cxcl13",
                                                         "Icam1","Cd34","Gdf10","Ptgis","Aldh1a2",
                                                          "Cxcl9","Cdk1"),
            min.cutoff = "q9",
            cols.use = c('lightgrey', 'brown'),
            pt.size = 0.2)

FeaturePlot(obs.mLN_Ontogeny_minus_main.merged, features.plot = c("Igf1","Igf2",
                                                         "Igfbp2","Igfbp4","Igfbp5","Igfbp6","Igfbp7",
                                                         "Igf1r","Igf2r"),
            min.cutoff = "q9",
            cols.use = c('lightgrey', 'brown'),
            pt.size = 0.2)

#####
#Imprinted and DAR-TF Signature across clusters and time points MAIN
#####
#Signature lists
DEGs_list <- list(NC_Pezoldt_Pasztoi_2018_maintained_mLN, NC_Pezoldt_Pasztoi_2018_repressed_mLN)
names(DEGs_list) <- c("maintained","repressed")
#Average gene expression 
obs.mLN_Ontogeny_minus_main.merged_minus_aExp <- AverageExpression(obs.mLN_Ontogeny_minus_main.merged)

#calculate cZscore for each signature across all clusters
#number of rows for output
index = 1
cluster_ID_i <- names(DEGs_list)[index]
y_label_i <- paste("cZ: ", cluster_ID_i)

print(cluster_ID_i)
print(i)

obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster <- subset(obs.mLN_Ontogeny_minus_main.merged_minus_aExp, rownames(obs.mLN_Ontogeny_minus_main.merged_minus_aExp) %in% DEGs_list[[index]][,1])

Scale_obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster <- apply(obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster, 1, scale)
#Kick columns wiht NaN
condense_table <- Scale_obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster
cols_to_take <- apply(condense_table,2, function(x) all(is.nan(x)))
condense_table <- condense_table[,!cols_to_take]
#Input for Zscore
Scale_obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster <- condense_table

Zscore_cluster_i <- rowSums(Scale_obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster)
Zscore_cluster_i_norm <- Zscore_cluster_i / (nrow(Scale_obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster)/10)

#mLN maintained Scores
Zscore_table_cluster_i <- data.frame(module = c(rep(cluster_ID_i,length(Zscore_cluster_i_norm))),
                                     cluster = colnames(obs.mLN_Ontogeny_minus_main.merged_minus_aExp),
                                     cZscore = Zscore_cluster_i_norm)
#Generate cZscore bargraph matrix
#p <- 
ggplot(data = Zscore_table_cluster_i, aes(x = cluster, y = cZscore, fill = module)) +
  geom_bar(stat = "identity", colour = "black") +
  theme(axis.text.x=element_text(angle=270,hjust=1), legend.position="none") +
  scale_fill_manual(values = c("pink")) +
  ylab(y_label_i)

#Make Heatmap
data_heatmap <- t(Scale_obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster)
colnames(data_heatmap) <- colnames(obs.mLN_Ontogeny_minus_main.merged_minus_aExp)
pheatmap(data_heatmap, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = TRUE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("blue", "white","red"), space="rgb")(128),
         main = title, fontsize = 8, fontsize_col = 12, fontsize_row = 12, fontsize_number = 12)

#####
#DAR-TF Signature across clusters and time points MAIN
#####
DAR_TF_existence_matrix <- DAR_TF_existence_matrix
DAR_TF_common <- rownames(DAR_TF_existence_matrix[rowSums(DAR_TF_existence_matrix) >= 3,])
DAR_TF_reg <- rownames(subset(DAR_TF_existence_matrix, 
                              pLN_Open_UP == 0 &
                                mLN_Open_UP == 0 &
                                pLN_peak_UP == 1 &
                                mLN_peak_UP == 0 &
                                pLN_Open_None == 0 &
                                mLN_Open_None == 0))
DAR_TF_act <- rownames(subset(DAR_TF_existence_matrix, 
                              pLN_Open_UP == 0 &
                                mLN_Open_UP == 0 &
                                pLN_peak_UP == 0 &
                                mLN_peak_UP == 0 &
                                pLN_Open_None == 0 &
                                mLN_Open_None == 1))
DAR_TF_act_reg <- rownames(subset(DAR_TF_existence_matrix, 
                                  pLN_Open_UP == 0 &
                                    mLN_Open_UP == 0 &
                                    pLN_peak_UP == 1 &
                                    mLN_peak_UP == 0 &
                                    pLN_Open_None == 0 &
                                    mLN_Open_None == 1))
ALL_TFs <- rownames(DAR_TF_existence_matrix)
ALL_mLN_TFs <- rownames(subset(DAR_TF_existence_matrix, 
                               pLN_Open_UP == 0 &
                                 pLN_peak_UP == 0 &
                                 pLN_Open_None == 0 &
                                 (mLN_Open_None == 1 |
                                 mLN_Open_UP == 1 |
                                 mLN_peak_UP == 1)))
DEGs_list <- list(DAR_TF_common,DAR_TF_reg, DAR_TF_act,DAR_TF_act_reg, ALL_TFs, ALL_mLN_TFs)


names(DEGs_list) <- c("common","TF_reg","TF_act","TF_act_reg","ALL","ALL_mLN_TFs")
#Average gene expression 
obs.mLN_Ontogeny_minus_main.merged_minus_aExp <- AverageExpression(obs.mLN_Ontogeny_minus_main.merged)

#calculate cZscore for each signature across all clusters
#number of rows for output
index = 6
cluster_ID_i <- names(DEGs_list)[index]
y_label_i <- paste("cZ: ", cluster_ID_i)

print(cluster_ID_i)

obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster <- subset(obs.mLN_Ontogeny_minus_main.merged_minus_aExp, rownames(obs.mLN_Ontogeny_minus_main.merged_minus_aExp) %in% DEGs_list[[index]])

Scale_obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster <- apply(obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster, 1, scale)
Zscore_cluster_i <- rowSums(Scale_obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster)
Zscore_cluster_i_norm <- Zscore_cluster_i / (nrow(Scale_obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster)/10)

#mLN maintained Scores
Zscore_table_cluster_i <- data.frame(module = c(rep(cluster_ID_i,length(Zscore_cluster_i_norm))),
                                     cluster = colnames(obs.mLN_Ontogeny_minus_main.merged_minus_aExp_j_aExp_cluster_i),
                                     cZscore = Zscore_cluster_i_norm)


#Generate cZscore bargraph matrix
#p <- 
ggplot(data = Zscore_table_cluster_i, aes(x = cluster, y = cZscore, fill = module)) +
  geom_bar(stat = "identity", colour = "black") +
  theme(axis.text.x=element_text(angle=270,hjust=1), legend.position="none") +
  scale_fill_manual(values = c("darkgreen")) +
  ylab(y_label_i)

#Make Heatmap
data_heatmap <- t(Scale_obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster)
colnames(data_heatmap) <- colnames(obs.mLN_Ontogeny_minus_main.merged_minus_aExp)
pheatmap(data_heatmap, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = TRUE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("blue", "white","red"), space="rgb")(128),
         main = title, fontsize = 8, fontsize_col = 12, fontsize_row = 12, fontsize_number = 12)

#####################
#Run analysis on pre-selected subsets (only Cdk1+)
#####################
#####
#Eliminate feeder cells Cdk1+
#####
unique_clusters <- c("Cdk1+")
All_meta <- subset(obs.mLN_Ontogeny_minus_main.merged@meta.data)
number_of_cells = 0
l_cells_include_CCA <- list()
for(i in 1:length(v_sample)){
  v_sample_i <- v_sample[i]
  print(v_sample_i)
  day_i_meta <- subset(All_meta, tech == v_sample_i)
  print(nrow(day_i_meta))
  cells_include_i <- row.names(subset(day_i_meta, Cluster %in% unique_clusters))
  print(length(cells_include_i))
  number_of_cells = number_of_cells + length(cells_include_i)
  l_cells_include_CCA[[i]] <- cells_include_i
  saveRDS(cells_include_i, paste(PATH_Output,"/",v_sample_i,"_1500_1000_1_12_only_Cdk1.Rds", sep=""))
}
print(number_of_cells)
#####
#Align with Cdk1+ cells
#####
obs.list_minus_main_Cdk1 <- list()
for(i in 1:length(obs.list_minus_main)){
  #Normalize Samples and cluster
  sample_i <- obs.list_minus_main[[i]]
  
  mito.genes <- grep(pattern = "^mt-", x = rownames(sample_i@data), value = TRUE)
  percent.mito <- colSums(as.matrix(sample_i@raw.data[mito.genes, ]))/colSums(as.matrix(sample_i@raw.data))
  
  #Identify cells with low nGene/nUmi ratio
  ratio.nUMITOnGene <- sample_i@meta.data$nUMI / sample_i@meta.data$nGene
  names(ratio.nUMITOnGene) <- rownames(sample_i@meta.data)
  
  # AddMetaData adds columns to object@data.info and good to store QC stats
  sample_i <- AddMetaData(object = sample_i, metadata = percent.mito, col.name = "percent.mito")
  sample_i <- AddMetaData(object = sample_i, metadata = ratio.nUMITOnGene, col.name = "ratio.nUMITOnGene")
  #Normalize the expression and log-transform
  sample_i <- NormalizeData(object = sample_i, normalization.method = "LogNormalize", 
                            scale.factor = 10000)
  
  #Find variable genes independent of expression level
  sample_i <- FindVariableGenes(object = sample_i, mean.function = ExpMean, dispersion.function = LogVMR, 
                                x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot = FALSE)
  
  #scale data
  sample_i <- ScaleData(object = sample_i)
  
  #perfrom PCA and print genes that define PCA
  sample_i <- RunPCA(object = sample_i, pc.genes = sample_i@var.genes, do.print = TRUE, pcs.print = 1:5, 
                     genes.print = 5)
  
  #ProjectPCA scores each gene in the dataset 
  sample_i <- ProjectPCA(object = sample_i, do.print = FALSE)
  
  #Group the cells into clusters
  sample_i <- FindClusters(object = sample_i, reduction.type = "pca", dims.use = 1:dims_use, 
                           resolution = resolution_single, print.output = 0, save.SNN = TRUE,
                           force.recalc = TRUE)
  
  #perfrom Tsne
  sample_i <- RunTSNE(object = sample_i, dims.use = 1:dims_use, do.fast = TRUE)
  
  #Cells to use
  sample_i.integrated_cells_chosen <- l_cells_include_CCA[[i]]
  if(length(sample_i.integrated_cells_chosen) > 10){
  
  #Prep samples without LECs and BECs
    sample_minus_i <- FilterCells(object = sample_i, subset.names = c("nGene", "ratio.nUMITOnGene", "percent.mito"), 
                                  low.thresholds = c(min_gene_number, -Inf, -Inf), high.thresholds = c(nGene, 4, mito_cutoff),
                                  cells.use = sample_i.integrated_cells_chosen)
    sample_minus_i <- NormalizeData(sample_minus_i)
    sample_minus_i <- FindVariableGenes(sample_minus_i, do.plot = F, display.progress = F)
    sample_minus_i <- ScaleData(sample_minus_i)
    sample_minus_i@meta.data$tech <- v_sample[i]
    obs.list_minus_main_Cdk1[[i]] <- sample_minus_i
  }
  rm(sample_minus_i)
}
#Eliminate NULL objects from llist
obs.list_minus_main_Cdk1 <- Filter(Negate(is.null), obs.list_minus_main_Cdk1)


saveRDS(obs.list_minus_main_Cdk1, paste(PATH_Output,"/obs.list_minus_main_only_Cdk1.Rds", sep=""))
obs.list_minus_main_Cdk1 <- readRDS(paste(PATH_Output,"/obs.list_minus_main_only_Cdk1.Rds", sep=""))

#Merge the seurat objects
# First two set the object
obs.mLN_Ontogeny_Cdk1.merged <- MergeSeurat(obs.list_minus_main_Cdk1[[1]], obs.list_minus_main_Cdk1[[2]], do.normalize = FALSE)
for(i in 1:(length(obs.list_minus_main_Cdk1)-2)){
  k = i + 2
  print(k)
  obs.mLN_Ontogeny_Cdk1.merged <- MergeSeurat(obs.mLN_Ontogeny_Cdk1.merged, obs.list_minus_main_Cdk1[[k]], do.normalize = FALSE)
}

#Normalize
obs.mLN_Ontogeny_Cdk1.merged <- NormalizeData(obs.mLN_Ontogeny_Cdk1.merged)
#Find variable genes
obs.mLN_Ontogeny_Cdk1.merged <- FindVariableGenes(obs.mLN_Ontogeny_Cdk1.merged, do.plot = F, display.progress = F)

#####
#Scale data and regress parameters Cdk1+
#####
# 1) Cell Cycle
#To lower case tranlating into mouse
Cell_cycle_Regev_BLOOD_minor <- tolower(Cell_cycle_Regev_BLOOD[,1])
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
Cell_cycle_genes_mouse <- firstup(Cell_cycle_Regev_BLOOD_minor)

#Phase specific
s.genes <- Cell_cycle_genes_mouse[1:43]
g2m.genes <- Cell_cycle_genes_mouse[44:97]

#Score each cell
obs.mLN_Ontogeny_Cdk1.merged <- CellCycleScoring(object = obs.mLN_Ontogeny_Cdk1.merged, s.genes = s.genes, g2m.genes = g2m.genes, 
                                                       set.ident = TRUE)

# 2) Mitochondria
mito.genes <- grep(pattern = "^mt-", x = rownames(obs.mLN_Ontogeny_Cdk1.merged@data), value = TRUE)
percent.mito <- colSums(as.matrix(obs.mLN_Ontogeny_Cdk1.merged@raw.data[mito.genes, ]))/colSums(as.matrix(obs.mLN_Ontogeny_Cdk1.merged@raw.data))
# AddMetaData adds columns to object@data.info and good to store QC stats
obs.mLN_Ontogeny_Cdk1.merged <- AddMetaData(object = obs.mLN_Ontogeny_Cdk1.merged, metadata = percent.mito, col.name = "percent.mito")

# 3) Ribosomal
#Calculate % ribosomal genes
rpl.genes <- grep(pattern = "^Rpl", x = rownames(obs.mLN_Ontogeny_Cdk1.merged@data), value = TRUE)
rps.genes <- grep(pattern = "^Rps", x = rownames(obs.mLN_Ontogeny_Cdk1.merged@data), value = TRUE)
ribo.genes <- c(rpl.genes, rps.genes)
percent.ribo <- colSums(as.matrix(obs.mLN_Ontogeny_Cdk1.merged@raw.data[ribo.genes, ]))/colSums(as.matrix(obs.mLN_Ontogeny_Cdk1.merged@raw.data))
obs.mLN_Ontogeny_Cdk1.merged <- AddMetaData(object = obs.mLN_Ontogeny_Cdk1.merged, metadata = percent.ribo, col.name = "percent.ribo")

#Scale
obs.mLN_Ontogeny_Cdk1.merged <- ScaleData(object = obs.mLN_Ontogeny_Cdk1.merged, vars.to.regress = c("S.Score", "G2M.Score","nGene", "percent.mito", "percent.ribo"), 
                                                display.progress = FALSE)

#Run PCA
obs.mLN_Ontogeny_Cdk1.merged <- RunPCA(object = obs.mLN_Ontogeny_Cdk1.merged, pc.genes = obs.mLN_Ontogeny_Cdk1.merged@var.genes, do.print = TRUE, pcs.print = 1:10, 
                                             genes.print = 10)


# Determine genes to use for CCA, must be highly variable in at least 2 datasets
genes.use <- c()
for (i in 1:length(obs.list_minus_main_Cdk1)) {
  genes.use <- c(genes.use, head(rownames(obs.list_minus_main_Cdk1[[i]]@hvg.info), 1500))
  print(names(obs.list_minus_main_Cdk1)[i])
  print(length(genes.use))
}
genes.use <- names(which(table(genes.use) > 1))
for (i in 1:length(obs.list_minus_main_Cdk1)) {
  genes.use <- genes.use[genes.use %in% rownames(obs.list_minus_main_Cdk1[[i]]@scale.data)]
  print(length(genes.use))
}

# Run multi-set CCA
obs.mLN_Ontogeny_Cdk1.merged <- RunCCA(obs.list_minus_main_Cdk1[[1]], obs.list_minus_main_Cdk1[[2]], genes.use = genes.use, num.cc = 20)

# CC Selection
MetageneBicorPlot(obs.mLN_Ontogeny_Cdk1.merged, grouping.var = "tech", dims.eval = 1:20)

# Run rare non-overlapping filtering
obs.mLN_Ontogeny_Cdk1.merged <- CalcVarExpRatio(object = obs.mLN_Ontogeny_Cdk1.merged, reduction.type = "pca",
                                                      grouping.var = "tech", dims.use = 1:9)
obs.mLN_Ontogeny_Cdk1.merged <- SubsetData(obs.mLN_Ontogeny_Cdk1.merged, subset.name = "var.ratio.pca",
                                                 accept.low = 0.5)

# Alignment
obs.mLN_Ontogeny_Cdk1.merged <- AlignSubspace(obs.mLN_Ontogeny_Cdk1.merged,
                                                    reduction.type = "cca",
                                                    grouping.var = "tech",
                                                    dims.align = 1:9)

# t-SNE and Clustering
obs.mLN_Ontogeny_Cdk1.merged <- RunTSNE(obs.mLN_Ontogeny_Cdk1.merged,
                                              reduction.use = "cca",
                                              dims.use = 1:9, save.SNN = T,
                                              resolution = 1.2, force.recalc = TRUE)
obs.mLN_Ontogeny_Cdk1.merged <- FindClusters(obs.mLN_Ontogeny_Cdk1.merged, reduction.type = "cca",
                                                   dims.use = 1:9, save.SNN = T, resolution = 0.6,
                                                   force.recalc = TRUE)

saveRDS(obs.mLN_Ontogeny_Cdk1.merged, paste(PATH_Output,"/MultiAlign_12_1_1500_regressed_1000genes_Cdk1_D0_D10_D24_D56_D300.Rds", sep=""))
obs.mLN_Ontogeny_Cdk1.merged <- readRDS(paste(PATH_Output,"/MultiAlign_12_1_1500_regressed_1000genes_Cdk1_D0_D10_D24_D56_D300.Rds", sep=""))
#obs.mLN_Ontogeny_minus_main.merged <- readRDS(paste(PATH_Output,"/MultiAlign_12_1_1500_regressed_1000genes_MAIN_D0_D10_D24_D56_D300.Rds", sep=""))
#Change cluster ID
current.cluster.ids <- c(0:2)

new.cluster.ids <- c("CD34+Crip1+","Cxcl13+Ccl19+","Cenpa+Ptgr1")

obs.mLN_Ontogeny_Cdk1.merged@ident <- plyr::mapvalues(x = obs.mLN_Ontogeny_Cdk1.merged@ident, from = current.cluster.ids, to = new.cluster.ids)
obs.mLN_Ontogeny_Cdk1.merged@meta.data$Cluster <- obs.mLN_Ontogeny_Cdk1.merged@ident

#Add Age score to metadata
#age_score <- list()
#for(i in 1:length(v_sample)){
#  age_i <- v_sample[i]
#  binary_i <- apply(obs.mLN_Ontogeny_minus_main.merged@meta.data, 1, function(r) any(r %in% c(age_i)))
#  binary_numeric_i =ifelse(binary_i=="TRUE",10,0)
#  age_score[[i]] <- binary_numeric_i
#}
#age_score <- do.call("cbind", age_score)
#colnames(age_score) <- v_sample
#append score to metadata
#meta_data_temp <- cbind(obs.mLN_Ontogeny_minus_main.merged@meta.data, age_score)
#obs.mLN_Ontogeny_minus_main.merged@meta.data <- meta_data_temp
#Age plotting
#FeaturePlot(obs.mLN_Ontogeny_minus_main.merged, features.plot = c("d056"), 
#           nCol = 1, cols.use = c("lightgrey", "black"), pt.size = 0.5,
#          cells.use = c(head(sample_1_cells,760),head(sample_2_cells,760),head(sample_3_cells,760),head(sample_4_cells,760),head(sample_5_cells,760)))


#Check allocation of cells
sample_1_cells <- grep(pattern = "^D000_", x = rownames(obs.mLN_Ontogeny_Cdk1.merged@meta.data), value = TRUE)
sample_2_cells <- grep(pattern = "^D010_", x = rownames(obs.mLN_Ontogeny_Cdk1.merged@meta.data), value = TRUE)
#sample_3_cells <- grep(pattern = "^D024_", x = rownames(obs.mLN_Ontogeny_Cdk1.merged@meta.data), value = TRUE)
#sample_4_cells <- grep(pattern = "^D056_", x = rownames(obs.mLN_Ontogeny_Cdk1.merged@meta.data), value = TRUE)
#sample_5_cells <- grep(pattern = "^D300_", x = rownames(obs.mLN_Ontogeny_Cdk1.merged@meta.data), value = TRUE)

p1 <- TSNEPlot(obs.mLN_Ontogeny_Cdk1.merged, group.by = "tech",
               do.return = T, pt.size = 0.1,
               cells.use = head(sample_1_cells,760),
               colors.use = "grey30")
p2 <- TSNEPlot(obs.mLN_Ontogeny_Cdk1.merged, group.by = "tech",
               do.return = T, pt.size = 0.1,
               cells.use = head(sample_2_cells,760),
               colors.use = "grey30")

p6 <- TSNEPlot(obs.mLN_Ontogeny_Cdk1.merged, do.return = T, pt.size = 0.6, do.label = TRUE,
               colors.use = c("darkorange3","blue","grey45"))
p7 <- TSNEPlot(obs.mLN_Ontogeny_Cdk1.merged, group.by = "tech", do.return = T, pt.size = 0.2,
               colors.use = c("goldenrod2","darkorange3"))

plot_grid(p6)
par(mfrow = c(4,2))
plot_grid(p1, p2, p3, p4, p5, p6, p7)
#Grey thumbnail
par(mfrow = c(5,1))
plot_grid(p1, p2, p3, p4, p5)
#Cell cycle plot
#To lower case tranlating into mouse
#Cell_cycle_Regev_BLOOD_minor <- tolower(Cell_cycle_Regev_BLOOD[,1])
#firstup <- function(x) {
#  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
#  x
#}
#Cell_cycle_genes_mouse <- firstup(Cell_cycle_Regev_BLOOD_minor)

#Phase specific
#s.genes <- Cell_cycle_genes_mouse[1:43]
#g2m.genes <- Cell_cycle_genes_mouse[44:97]

#Score each cell
#obs.mLN_Ontogeny_Cdk1.merged <- CellCycleScoring(object = obs.mLN_Ontogeny_Cdk1.merged,
 #                                                      s.genes = s.genes, g2m.genes = g2m.genes, 
  #                                                     set.ident = TRUE)
#FeaturePlot(obs.mLN_Ontogeny_Cdk1.merged, features.plot = c("S.Score", "G2M.Score"), min.cutoff = "q05", max.cutoff = "q95", 
 #           nCol = 2, cols.use = c("lightgrey", "blue"), pt.size = 0.5)


FeaturePlot(obs.mLN_Ontogeny_Cdk1.merged, features.plot = c("Tnfsf11","Nfkbia","Crip1","Ccl19",
                                                            "Madcam1","Cxcl13","Ptgr1",
                                                            "Icam1","Cd34","Ptgis","Sparcl1",
                                                            "Enpp2","Cdk1","Des"),
            min.cutoff = "q9",
            cols.use = c('lightgrey', 'brown'),
            pt.size = 0.5)

FeaturePlot(obs.mLN_Ontogeny_Cdk1.merged, features.plot = c("Ebf4","Ebf3","Mybl2","E2f1",
                                                            "Elf1","Atf5","Hoxb4","Cd34",
                                                            "Ebf2","Isl1","Mybl1","Elk1"),
            min.cutoff = "q9",
            cols.use = c('lightgrey', 'brown'),
            pt.size = 0.5)

#####
#DEGs Cdk1+
#####
#Write Metadata
saveRDS(obs.mLN_Ontogeny_Cdk1.merged@meta.data, paste(PATH_Output,"/","MetaData_1500_1000_1_12_SC_Cdk1.Rds", sep=""))

#####
#Find differentially expressed genes Cdk1+
#####
#Find markers for every cluster compared to all remaining cells
mLN_Ontogeny_minus.merged.markers <- FindAllMarkers(object = obs.mLN_Ontogeny_Cdk1.merged, only.pos = TRUE, min.pct = 0.1,
                                                    thresh.use = 0.1)
mLN_Ontogeny_minus.merged.markers_20 <- mLN_Ontogeny_minus.merged.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
mLN_Ontogeny_minus.merged.markers_40 <- mLN_Ontogeny_minus.merged.markers %>% group_by(cluster) %>% top_n(40, avg_logFC)
mLN_Ontogeny_minus.merged.markers_60 <- mLN_Ontogeny_minus.merged.markers %>% group_by(cluster) %>% top_n(60, avg_logFC)
mLN_Ontogeny_minus.merged.markers_20 <- mLN_Ontogeny_minus.merged.markers$gene
#Order LN_minus meta.data according to cluster size

#Save DEGs
write.table(mLN_Ontogeny_minus.merged.markers, paste(PATH_Output,"/MultiAlign_12_1_regressed_1000genes_DEGs_Cdk1_.csv", sep=""), dec=".", sep=",")
mLN_Ontogeny_minus.merged.markers <- read.table(paste(PATH_Output,"/MultiAlign_12_1_regressed_1000genes_DEGs_Cdk1.csv", sep=""), dec=".", sep=",")
#Plot differential expression of Marker genes
DoHeatmap(object = obs.mLN_Ontogeny_minus_main.merged, use.scaled = TRUE, group.by = "ident",
          genes.use = mLN_Ontogeny_minus.merged.markers_20, title = "Merged FSC Ontogeny",
          remove.key = FALSE,
          col.low = "white",
          col.mid =  "oldlace",
          col.high = "brown")

#####
#Imprinted Signature across clusters and time points Cdk1+
#####
#Signature lists
DEGs_list <- list(NC_Pezoldt_Pasztoi_2018_maintained_mLN, NC_Pezoldt_Pasztoi_2018_repressed_mLN)
names(DEGs_list) <- c("maintained","repressed")
#Average gene expression 
obs.mLN_Ontogeny_minus_main.merged_minus_aExp <- AverageExpression(obs.mLN_Ontogeny_Cdk1.merged)

#calculate cZscore for each signature across all clusters
#number of rows for output
index = 1
cluster_ID_i <- names(DEGs_list)[index]
y_label_i <- paste("cZ: ", cluster_ID_i)

print(cluster_ID_i)
print(i)

obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster <- subset(obs.mLN_Ontogeny_minus_main.merged_minus_aExp, rownames(obs.mLN_Ontogeny_minus_main.merged_minus_aExp) %in% DEGs_list[[index]][,1])

Scale_obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster <- apply(obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster, 1, scale)
#Kick columns wiht NaN
condense_table <- Scale_obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster
cols_to_take <- apply(condense_table,2, function(x) all(is.nan(x)))
condense_table <- condense_table[,!cols_to_take]
#Input for Zscore
Scale_obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster <- condense_table

Zscore_cluster_i <- rowSums(Scale_obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster)
Zscore_cluster_i_norm <- Zscore_cluster_i / (nrow(Scale_obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster)/10)

#mLN maintained Scores
Zscore_table_cluster_i <- data.frame(module = c(rep(cluster_ID_i,length(Zscore_cluster_i_norm))),
                                     cluster = colnames(obs.mLN_Ontogeny_minus_main.merged_minus_aExp),
                                     cZscore = Zscore_cluster_i_norm)
#Generate cZscore bargraph matrix
#p <- 
ggplot(data = Zscore_table_cluster_i, aes(x = cluster, y = cZscore, fill = module)) +
  geom_bar(stat = "identity", colour = "black") +
  theme(axis.text.x=element_text(angle=270,hjust=1), legend.position="none") +
  scale_fill_manual(values = c("pink")) +
  ylab(y_label_i)

#Make Heatmap
data_heatmap <- t(Scale_obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster)
colnames(data_heatmap) <- colnames(obs.mLN_Ontogeny_minus_main.merged_minus_aExp)
pheatmap(data_heatmap, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = TRUE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("blue", "white","red"), space="rgb")(128),
         main = title, fontsize = 8, fontsize_col = 12, fontsize_row = 12, fontsize_number = 12)



#####
#DAR-TF Signature across clusters and time points Cdk1+
#####
DAR_TF_existence_matrix <- DAR_TF_existence_matrix
DAR_TF_common <- rownames(DAR_TF_existence_matrix[rowSums(DAR_TF_existence_matrix) >= 3,])
DAR_TF_reg <- rownames(subset(DAR_TF_existence_matrix, 
                              pLN_Open_UP == 0 &
                                mLN_Open_UP == 0 &
                                pLN_peak_UP == 1 &
                                mLN_peak_UP == 0 &
                                pLN_Open_None == 0 &
                                mLN_Open_None == 0))
DAR_TF_act <- rownames(subset(DAR_TF_existence_matrix, 
                              pLN_Open_UP == 0 &
                                mLN_Open_UP == 0 &
                                pLN_peak_UP == 0 &
                                mLN_peak_UP == 0 &
                                pLN_Open_None == 0 &
                                mLN_Open_None == 1))
DAR_TF_act_reg <- rownames(subset(DAR_TF_existence_matrix, 
                                  pLN_Open_UP == 0 &
                                    mLN_Open_UP == 0 &
                                    pLN_peak_UP == 1 &
                                    mLN_peak_UP == 0 &
                                    pLN_Open_None == 0 &
                                    mLN_Open_None == 1))
ALL_TFs <- rownames(DAR_TF_existence_matrix)
ALL_mLN_TFs <- rownames(subset(DAR_TF_existence_matrix, 
                               pLN_Open_UP == 0 &
                                 pLN_peak_UP == 0 &
                                 pLN_Open_None == 0 &
                                 (mLN_Open_None == 1 |
                                    mLN_Open_UP == 1 |
                                    mLN_peak_UP == 1)))
DEGs_list <- list(DAR_TF_common,DAR_TF_reg, DAR_TF_act,DAR_TF_act_reg, ALL_TFs, ALL_mLN_TFs)


names(DEGs_list) <- c("common","TF_reg","TF_act","TF_act_reg","ALL","ALL_mLN_TFs")
#Average gene expression 
obs.mLN_Ontogeny_minus_main.merged_minus_aExp <- AverageExpression(obs.mLN_Ontogeny_Cdk1.merged)

#calculate cZscore for each signature across all clusters
#number of rows for output
index = 6
cluster_ID_i <- names(DEGs_list)[index]
y_label_i <- paste("cZ: ", cluster_ID_i)

print(cluster_ID_i)

obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster <- subset(obs.mLN_Ontogeny_minus_main.merged_minus_aExp, rownames(obs.mLN_Ontogeny_minus_main.merged_minus_aExp) %in% DEGs_list[[index]])

Scale_obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster <- apply(obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster, 1, scale)
#Kick columns wiht NaN
condense_table <- Scale_obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster
cols_to_take <- apply(condense_table,2, function(x) all(is.nan(x)))
condense_table <- condense_table[,!cols_to_take]
#Input for Zscore
Scale_obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster <- condense_table
Zscore_cluster_i <- rowSums(Scale_obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster)
Zscore_cluster_i_norm <- Zscore_cluster_i / (nrow(Scale_obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster)/10)

#mLN maintained Scores
Zscore_table_cluster_i <- data.frame(module = c(rep(cluster_ID_i,length(Zscore_cluster_i_norm))),
                                     cluster = colnames(obs.mLN_Ontogeny_minus_main.merged_minus_aExp),
                                     cZscore = Zscore_cluster_i_norm)


#Generate cZscore bargraph matrix
#p <- 
ggplot(data = Zscore_table_cluster_i, aes(x = cluster, y = cZscore, fill = module)) +
  geom_bar(stat = "identity", colour = "black") +
  theme(axis.text.x=element_text(angle=270,hjust=1), legend.position="none") +
  scale_fill_manual(values = c("darkgreen")) +
  ylab(y_label_i)

#Make Heatmap
data_heatmap <- t(Scale_obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster)
colnames(data_heatmap) <- colnames(obs.mLN_Ontogeny_minus_main.merged_minus_aExp)
pheatmap(data_heatmap, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = TRUE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("blue", "white","red"), space="rgb")(128),
         main = names(DEGs_list)[index], fontsize = 8, fontsize_col = 12, fontsize_row = 12, fontsize_number = 12)

#Overlap with DEGs
genes_in_heatmap <- subset(mLN_Ontogeny_minus.merged.markers, gene %in% DEGs_list[[index]])

#####
#Cluster identification via signature
#####
#Use list of DEGs from merged analysis of mLN and pLN multiple samples

DEGs_core_40 <- NC_Pezoldt_Pasztoi_2018

#Load DEG tables and make list of clusters
DEGs_list <- split(DEGs_core_40, DEGs_core_40$cluster)
#Average gene expression 
obs.mLN_Ontogeny_minus_main.merged_minus_aExp <- AverageExpression(obs.mLN_Ontogeny_minus_main.merged)

#calculate cZscore for each signature across all clusters
#number of rows for output
i = 1
out = NULL
for(i in i:length(DEGs_list)){
  cluster_DEG_i <- as.character(DEGs_list[[i]]$gene)
  cluster_DEG_n_i <- length(cluster_DEG_i)
  
  cluster_ID_i <- names(DEGs_list)[i]
  y_label_i <- paste("cZ: ", cluster_ID_i)
  
  print(cluster_ID_i)
  print(i)
  print(cluster_DEG_n_i)
  
  obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster_i <- subset(obs.mLN_Ontogeny_minus_main.merged_minus_aExp, rownames(obs.mLN_Ontogeny_minus_main.merged_minus_aExp) %in% cluster_DEG_i)
  
  Scale_obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster_i <- apply(obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster_i, 1, scale)
  Zscore_cluster_i <- rowSums(Scale_obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster_i)
  Zscore_cluster_i_norm <- Zscore_cluster_i / (nrow(obs.mLN_Ontogeny_minus_main.merged_minus_aExp_cluster_i)/10)
  
  #mLN maintained Scores
  Zscore_table_cluster_i <- data.frame(module = c(rep(cluster_ID_i,length(Zscore_cluster_i_norm))),
                                       cluster = colnames(obs.mLN_Ontogeny_minus_main.merged_minus_aExp),
                                       cZscore = Zscore_cluster_i_norm)
  
  
  #Generate cZscore bargraph matrix
  p <- ggplot(data = Zscore_table_cluster_i, aes(x = cluster, y = cZscore, fill = module)) +
    geom_bar(stat = "identity", colour = "black") +
    theme(axis.text.x=element_text(angle=270,hjust=1), legend.position="none") +
    scale_fill_manual(values = c("deepskyblue1")) +
    ylab(y_label_i)
  #print(p)
  out[[i]] <- p
}

title <- paste("Ontogeny_MultiAlign_DEGs_core_40",sep = "_")
pdf(paste(PATH_Output,"/MultiAlign_12_1_1000_regressed_all_NC_Top40_D0_D10_D24_D56_D300_MAIN_minus_feeder123_PvC_.pdf", sep=""))
marrangeGrob(out, nrow = round(length(DEGs_list) / 3 - 1), ncol = 3, top = title)
dev.off()


#####
#Find differentially expressed genes
#####
#Find markers for every cluster compared to all remaining cells
mLN_Ontogeny_minus.merged.markers <- FindAllMarkers(object = obs.mLN_Ontogeny_minus_main.merged, only.pos = TRUE, min.pct = 0.25, 
                                                    thresh.use = 0.25)
mLN_Ontogeny_minus.merged.markers_20 <- mLN_Ontogeny_minus.merged.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
mLN_Ontogeny_minus.merged.markers_40 <- mLN_Ontogeny_minus.merged.markers %>% group_by(cluster) %>% top_n(40, avg_logFC)
mLN_Ontogeny_minus.merged.markers_60 <- mLN_Ontogeny_minus.merged.markers %>% group_by(cluster) %>% top_n(60, avg_logFC)
mLN_Ontogeny_minus.merged.markers_20 <- mLN_Ontogeny_minus.merged.markers$gene
#Order LN_minus meta.data according to cluster size

#Save DEGs
write.table(mLN_Ontogeny_minus.merged.markers, paste(PATH_Output,"/MultiAlign_12_1_regressed_1000genes_DEGs_D0_D10_D24_D56_D300_MAIN_minus_feeder123_PvC_.csv", sep=""), dec=".", sep=",")
mLN_Ontogeny_minus.merged.markers <- read.table(paste(PATH_Output,"/MultiAlign_12_1_regressed_1000genes_DEGs_D0_D10_D24_D56_D300_MAIN_minus_feeder123_PvC_.csv", sep=""), dec=".", sep=",")
#Plot differential expression of Marker genes
DoHeatmap(object = obs.mLN_Ontogeny_minus_main.merged, use.scaled = TRUE, group.by = "ident",
         genes.use = mLN_Ontogeny_minus.merged.markers_20, title = "Merged FSC Ontogeny",
         remove.key = FALSE,
         col.low = "white",
      col.mid =  "oldlace",
     col.high = "brown")

mLN_Ontogeny_minus.merged.markers <- read.table(paste(PATH_Output,"/DEGs_D0_D10_D24_D56_D300.csv", sep=""), dec=".", sep=",")

#####
#Dotplot
#####

mLN_Ontogeny_minus.merged.markers_nonDup <- mLN_Ontogeny_minus.merged.markers[!duplicated(mLN_Ontogeny_minus.merged.markers$gene),]
mLN_Ontogeny_minus.merged.markers_nonDup_10_pValue <- mLN_Ontogeny_minus.merged.markers_nonDup %>% group_by(cluster) %>% top_n(10, 1/p_val_adj)
mLN_Ontogeny_minus.merged.markers_nonDup_20_pValue <- mLN_Ontogeny_minus.merged.markers_nonDup %>% group_by(cluster) %>% top_n(20, 1/p_val_adj)
mLN_Ontogeny_minus.merged.markers_nonDup_10_FC <- mLN_Ontogeny_minus.merged.markers_nonDup %>% group_by(cluster) %>% top_n(8, avg_logFC)
mLN_Ontogeny_minus.merged.markers_nonDup_20_FC <- mLN_Ontogeny_minus.merged.markers_nonDup %>% group_by(cluster) %>% top_n(20, avg_logFC)

#Plot differential expression of Marker genes
title <- "Dim20_Res1_100_Regressed_Dotplot_Top10_D0_10_24_56_300"
#pdf(paste(path_output,"/",title, ".pdf", sep = ""), width = 25, height = 10)
DotPlot(obs.mLN_Ontogeny_minus_main.merged, mLN_Ontogeny_minus.merged.markers_nonDup_10_FC$gene,
        x.lab.rot = TRUE, plot.legend = FALSE,
        col.min = 0, col.max = 5, dot.scale = 5,
        cols.use = c("lightgrey", "brown"))
#dev.off()
#dev.off()
#####

#####
#Hierarchical clustering according to Top DEGs per cluster
#####
#Average expression
mLN_Ontogeny_minus.merged_aExp <- AverageExpression(obs.mLN_Ontogeny_minus_main.merged)
#Check for DEGs
#genes <- as.character(LN_minus.markers_20[1:260,]$gene)
mLN_Ontogeny_minus.merged_aExp_TopDEGs <- subset(mLN_Ontogeny_minus.merged_aExp, rownames(mLN_Ontogeny_minus.merged_aExp) %in% mLN_Ontogeny_minus.merged.markers_60$gene)
mLN_Ontogeny_minus.merged_aExp_TopDEGs <- as.matrix(mLN_Ontogeny_minus.merged_aExp_TopDEGs)

title <- paste("Heatmap_hierarchical",sep = "_")
#pdf(paste(path_output,"/",title, ".pdf", sep = ""), width = 6, height = 7)
pheatmap(mLN_Ontogeny_minus.merged_aExp_TopDEGs, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 20, show_rownames = FALSE, cluster_cols = TRUE,
         scale = "row", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("blue", "white","red"), space="rgb")(128))
#dev.off()


#####
#Calculate % per cluster across age
#####
meta_data <- obs.mLN_Ontogeny_minus_main.merged@meta.data
l_condition <- split(meta_data, meta_data$tech)
l_cells_per_condition <- lapply(l_condition, function(x) {
  #x <- l_condition[[1]]
  x_cells <- unlist(lapply(split(x, x$Cluster), nrow))
  x_cells
})
names(l_cells_per_condition) <- names(l_condition)
t_cells_per_condition <- do.call(rbind.fill.matrix, lapply(l_cells_per_condition, function(x) {t(as.data.frame(x))}))
t_cells_per_condition[is.na(t_cells_per_condition)] <- 0
rownames(t_cells_per_condition) <- names(l_condition)
#Normalize to 1000 cells
norm_factor_vector <- rowSums(t_cells_per_condition) / 1000
#Divide cell number per cluster by normalization factor
t_cells_per_condition_norm <- t_cells_per_condition / norm_factor_vector
#Calculate Frequencies
total <- colSums(t_cells_per_condition_norm)
freq_cells_per_condition <- (t(t_cells_per_condition_norm) / total) * 100
freq_cells_per_condition <- as.matrix(freq_cells_per_condition)
pheatmap(freq_cells_per_condition, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = FALSE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("white","grey","grey23","black"), space="rgb")(128))

###########################################
#GO analysis
###########################################
library(stringr)
library(pheatmap)
library(ggplot2)
library(foreign)
library(GO.db)
library(topGO)
library(org.Mm.eg.db)
#####
#Make a gene2GO list
#####
x <- org.Mm.egGO
# Get the entrez gene identifiers that are mapped to a GO ID
mapped_genes <- mappedkeys(x)
# Build a list GeneID and GOs
GeneID2GO <- list()

xx <- as.list(x[mapped_genes])

for(i in 1:length(xx)){
  #Initiate vector to collect GOIDs for geneID i
  GO_vector <- c()
  #Get geneID
  Gene_ID <- ls(xx[i])
  #grab data for geneID_i
  temp <- xx[[i]]
  
  #check Ontology category
  for(i in 1:length(temp)){
    category <- as.character(temp[[i]]["Ontology"])
    #if ontology category matches collect GOIDs  (flex)
    if(category == "BP"){
      temp_GOID <- as.character(temp[[i]]["GOID"])
      GO_vector <- c(GO_vector, temp_GOID)
    }
    #print(GO_vector)
  }
  #Generate list name geneID, content
  GO_IDs <- GO_vector 
  GeneID2GO[[Gene_ID]] <- GO_IDs 
}

GO2GeneID <- inverseList(GeneID2GO)



#####
#Input data
#####
myfiles <- list(mLN_Ontogeny_minus.merged.markers)
#####
#GeneSymbol to GeneID
#####
for(i in 1:length(myfiles)){
  idfound <- myfiles[[i]]$gene %in% mappedRkeys(org.Mm.egSYMBOL)
  SYMBOL <- toTable(org.Mm.egSYMBOL)
  head(SYMBOL)
  m <- match(myfiles[[i]]$gene, SYMBOL$symbol)
  GENE_ID <- SYMBOL$gene_id[m]
  myfiles[[i]] <- cbind(GENE_ID, myfiles[[i]])
}
#####
#Perform GO analysis single
#####
#foreach resolution make
#determine number of clusters: split by cluster
#link gene IDs with p-values
#perform GO
l_gene_pval <- list()
l_all_gene_pval <- list()

for(i in 1:length(myfiles)){
  split_clusters <- split(myfiles[[i]], myfiles[[i]]$cluster)
  print("outer")
  print(i)
  for(i in 1:length(split_clusters)){
    genes_of_interest_i <- as.character(split_clusters[[i]]$GENE_ID)
    genes_pval_i <- split_clusters[[i]]$p_val
    names(genes_pval_i) <- genes_of_interest_i
    l_gene_pval[[i]] <- genes_pval_i 
  }
  print(class(l_gene_pval))
  print(length(l_gene_pval))
  l_all_gene_pval[[length(l_all_gene_pval)+1]] <- l_gene_pval
}


#get one gene list
l_gene_pval <- l_all_gene_pval[[1]]

#Gene universe
geneNames <- myfiles[[1]]$GENE_ID

#gene lists
l_gene_List <- list()
for(i in 1:length(l_gene_pval)){
  geneList_i <- factor(as.integer(geneNames %in% names(l_gene_pval[[i]])))
  names(geneList_i) <- geneNames
  l_gene_List[[i]] <- geneList_i 
}


List_allRes <- list()
#Do  GO statistics for all gene lists
for(i in 1:length(l_gene_List)){
  #Access the gene lists and p-values for the differentially expressed genes
  geneList <- l_gene_List[[i]]
  pvalue_of_interest <- l_gene_pval[[i]]
  
  #build GOdata object
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, geneSel = pvalue_of_interest,
                annot = annFUN.gene2GO , gene2GO = GeneID2GO, nodeSize = 5)
  
  #get number of differentially expressed genes in the GOdata object
  sg <- sigGenes(GOdata)
  numSigGenes(GOdata)
  #get the number of GO_IDs that are within the applied GeneUniverse
  graph(GOdata)
  number_GOIDs <- usedGO(GOdata)
  number_nodes <- length(number_GOIDs)
  
  
  #Run statistics
  
  #Fisher's works with weight
  test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
  resultWeight <- getSigGroups(GOdata, test.stat)
  #KS works with elim but not with weight
  test.stat <- new("elimScore", testStatistic = GOKSTest, name = "Fisher test", cutOff = 0.01)
  resultElim <- getSigGroups(GOdata, test.stat)
  #runTest a high level interface for testing Fisher
  resultFis <- runTest(GOdata, statistic = "fisher")
  #Kolmogorov-Smirnov includes p-values
  test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
  resultKS <- getSigGroups(GOdata,  test.stat)
  #runTest a high level interface for testing KS
  elim.ks <- runTest(GOdata, algorithm = "elim", statistic = "ks")
  
  #make table
  allRes <- GenTable(GOdata, classic = resultFis, KS = resultKS, weight = resultWeight,
                     orderBy = "weight", ranksOf = "KS", topNodes = number_nodes)
  #make list of result tables
  List_allRes[[i]] <- allRes
}

saveRDS(List_allRes, paste(PATH_Output,"/GO_D0_D10_D24_D60_D300_List_allRes.Rds", sep=""))
List_allRes <- readRDS(paste(PATH_Output,"/GO_D0_D10_D24_D60_D300_List_allRes.Rds", sep=""))

####
#Find the Top GOs from the lists
####
#Build table with Top GOs according to "weight"
List_TopGOs <- list()
for(i in 1:length(List_allRes)){
  table_i <- List_allRes[[i]]
  table_tophits <- subset(table_i, weight < 0.005)
  #print(i)
  #print(nrow(table_tophits))
  List_TopGOs[[i]] <- table_tophits
}

#collect all TopGos
TopGOs_vector <- c()
for(i in 1:length(List_TopGOs)){
  table_i <- List_TopGOs[[i]]
  TopGOs_i <- as.character(table_i$GO.ID)
  TopGOs_vector <- c(TopGOs_vector, TopGOs_i) 
}

#Condense tables and sort by GO.ID
l_topGO <- list()
for(i in 1:length(List_allRes)){
  table_i <- List_allRes[[i]]
  table_i_subset <- subset(table_i, table_i$GO.ID %in% TopGOs_vector)
  table_i_subset_GO_weight <- table_i_subset[,c("GO.ID","weight", "Term")]
  k = i - 1
  cluster_weight_name <- paste("c_", k, "_weight", sep="")
  colnames(table_i_subset_GO_weight)[2] <- cluster_weight_name
  table_i_subset_GO_weight[,2] <- as.numeric(table_i_subset_GO_weight[,2])
  table_i_subset_GO_weight <- table_i_subset_GO_weight[order(table_i_subset_GO_weight$"GO.ID", decreasing=TRUE), ]
  l_topGO[[i]] <- table_i_subset_GO_weight
}
#merger
data_TopGO_weight = Reduce(function(...) merge(..., all=T), l_topGO)



####
#Make GO comparison heatmap
####
data_heatmap <- data_TopGO_weight
row.names(data_heatmap) <- paste(data_heatmap$"GO.ID", data_heatmap$Term, sep = "_")
data_heatmap <- data_heatmap[,3:ncol(data_heatmap)]

min(data_heatmap)
data_heatmap <- -log10(data_heatmap)
#Common Minimum at 5
data_heatmap[data_heatmap > 5] <- 5
data_heatmap_matrix <- data.matrix(data_heatmap)
colnames(data_heatmap_matrix) <- new.cluster.ids
title <- c("Differential GO scRNASeq")

pheatmap(data_heatmap_matrix, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = FALSE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("white", "grey","deepskyblue2"), space="rgb")(128),
         main = title)

#######
#Identify gene in cluster contributing to GO
#######
#Associate GeneID to DEGs
Data_ready <- myfiles[[1]]




#########
#Heatmap for set of genes
#########
Data_ready <- myfiles[[1]]
List_GOs_interest <- c("GO:0006306")
#List_GOs_interest_overrep_D0_developmental <- c("GO:1905562","GO:0001657","GO:0043627","GO:0048706","GO:0030326","GO:0090263","GO:0045668","GO:0072074","GO:2000288")
Genes_for_GOs_interest <- subset(GO2GeneID, ls(GO2GeneID) %in% List_GOs_interest)
number_input <- length(Genes_for_GOs_interest)
AllGenes_interest <- unlist(l_gene_pval)

#Save Output Tables in List
l_Genes_GO_per_cluster_interest <- list()
#Store tables
l_genes_per_cluster_GOinterest <- list()

for(i in 1:number_input){
  
  #Get genes in GO group
  GO_interest_Genes <- AllGenes_interest[names(AllGenes_interest) %in% Genes_for_GOs_interest[[1]]]
  #print(GO_interest_Genes)
  number_genes_in_GO <- length(GO_interest_Genes)
  print(number_genes_in_GO)
  #Get gene names
  gene_names <- GO_interest_Genes
  print(gene_names)
  #Generate list of genes of interest contained in expression data
  genes_expression <- subset(Data_ready, Data_ready$GENE_ID %in% names(gene_names))
  number_genes_expressed_in_GO <- nrow(genes_expression)
  
  if(number_genes_expressed_in_GO >= 3){
    #format data for pheatmap
    l_genes_per_cluster_GOinterest[[i]] <- genes_expression
  }
}
l_genes_per_cluster_GOinterest[[1]]




