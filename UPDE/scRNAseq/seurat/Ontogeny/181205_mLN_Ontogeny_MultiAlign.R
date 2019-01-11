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

#######
#Global variables
#######
#@User: Define sample names
#M == merge
sample_1 <- "d000"
sample_2 <- "d010"
sample_3 <- "d024"
sample_4 <- "d056"
sample_5 <- "d300"
v_sample <- c(sample_1, sample_2, sample_3, sample_4,sample_5)
#setwd("/Users/Pezoldt/PowerFolders/R/2018_Exp270_Ontogeny/Output/Multialign_SEURAT")

#Parameters Seurat
min_gene_number <- 250
mito_cutoff <- 0.045
nGene <- 4700
min_cells <- 20
dims_use <- 16
resolution_single <- 1.2
resolution_merge <- 1.2

#######
#Load and label
#######
###
#mLN d0
###
#sample_1.data <- Read10X(data.dir = "C:/Users/jpe12/PowerFolders/R/2018_Exp270_Ontogeny/Data/neonatal/filtered_gene_bc_matrices/mm10")
sample_1.data <- Read10X(data.dir = "/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/270_2018-01-03_scRNA-Seq_mLN_ontogeny/Data/C11_neo_mLN/outs/filtered_gene_bc_matrices/mm10")
sample_1.data@Dimnames[[2]] <- paste("sample_1_", c(1:length(sample_1.data@Dimnames[[2]])), sep = "")
###
#mLN d10
###
#sample_2.data <- Read10X(data.dir = "C:/Users/jpe12/PowerFolders/R/2018_Exp270_Ontogeny/Data/d10/filtered_gene_bc_matrices/mm10")
sample_2.data <- Read10X(data.dir = "/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/270_2018-01-03_scRNA-Seq_mLN_ontogeny/Data/C12_d10_mLN/outs/filtered_gene_bc_matrices/mm10")
sample_2.data@Dimnames[[2]] <- paste("sample_2_", c(1:length(sample_2.data@Dimnames[[2]])), sep = "")

###
#mLN d24
###
#sample_3.data <- Read10X(data.dir = "C:/Users/jpe12/PowerFolders/R/2018_Exp270_Ontogeny/Data/d24/filtered_gene_bc_matrices/mm10")
sample_3.data <- Read10X(data.dir = "/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/270_2018-01-03_scRNA-Seq_mLN_ontogeny/Data/D1_d24_mLN/outs/filtered_gene_bc_matrices/mm10")
sample_3.data@Dimnames[[2]] <- paste("sample_3_", c(1:length(sample_3.data@Dimnames[[2]])), sep = "")

###
#mLN d60
###
#Experiment 1
sample_4_1.data <- read.csv("/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/2018_NC_stroma/Exp_244/mLNSPF.csv")
row.names(sample_4_1.data) <- sample_4_1.data$X
sample_4_1.data <- sample_4_1.data[,c(2:ncol(sample_4_1.data))]
colnames(sample_4_1.data) <- paste("sample_4_1_", c(1:ncol(sample_4_1.data)), sep = "")
#Experiment 2
#sample_4.data <- Read10X(data.dir = "C:/Users/jpe12/PowerFolders/R/2017_252_scRNAseq/Data/mLN_SPF/B6/filtered_gene_bc_matrices/mm10")
sample_4_2.data <- Read10X(data.dir = "/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/252_2017-07-19_scRNA-Seq_mLN_pLN_SPF_GF/Data/10x/mLN_SPF/B6/outs/filtered_gene_bc_matrices/mm10")
sample_4_2.data@Dimnames[[2]] <- paste("sample_4_2_", c(1:length(sample_4_2.data@Dimnames[[2]])), sep = "")

###
#mLN d300
###
sample_5.data <- Read10X(data.dir = "/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/275_2018-04-23_scRNA-Seq_mLN_SPF_45wk/Data/E7_SPF_mLN_42wk/outs/filtered_gene_bc_matrices/mm10")
sample_5.data@Dimnames[[2]] <- paste("sample_5_", c(1:length(sample_5.data@Dimnames[[2]])), sep = "")


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

# Run multi-set CCA
mLN_Ontogeny.integrated <- RunMultiCCA(ob.list, genes.use = genes.use, num.ccs = dims_use)

# CC Selection
MetageneBicorPlot(mLN_Ontogeny.integrated, grouping.var = "tech", dims.eval = 1:15)

# Run rare non-overlapping filtering
mLN_Ontogeny.integrated <- CalcVarExpRatio(object = mLN_Ontogeny.integrated, reduction.type = "pca",
                                           grouping.var = "tech", dims.use = 1:dims_use)
mLN_Ontogeny.integrated <- SubsetData(mLN_Ontogeny.integrated, subset.name = "var.ratio.pca",
                                      accept.low = 0.5)

# Alignment
mLN_Ontogeny.integrated <- AlignSubspace(mLN_Ontogeny.integrated,
                                         reduction.type = "cca",
                                         grouping.var = "tech",
                                         dims.align = 1:16)

#Plot the aligned CCA ACCA
#p1 <- VlnPlot(mLN_Ontogeny.integrated, features.plot = "ACC1", group.by = "tech", do.return = T)
#p2 <- VlnPlot(mLN_Ontogeny.integrated, features.plot = "ACC2", group.by = "tech", do.return = T)
#plot_grid(p1, p2)


# t-SNE and Clustering
mLN_Ontogeny.integrated <- FindClusters(mLN_Ontogeny.integrated, reduction.type = "cca",
                                        dims.use = 1:16, save.SNN = T, resolution = resolution_merge)
mLN_Ontogeny.integrated <- RunTSNE(mLN_Ontogeny.integrated,
                                   reduction.use = "cca",
                                   dims.use = 1:16, save.SNN = T,
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
               colors.use = c("orange","green","blue","purple"))
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
  ### Continue here:
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
  
  # Examine and visualize PCA results
  #PrintPCA(object = sample_1_seurat, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
  #VizPCA(object = sample_1_seurat, pcs.use = 1:6)
  #PCAPlot(object = sample_1_seurat, dim.1 = 1, dim.2 = 2)
  
  #ProjectPCA scores each gene in the dataset 
  sample_i <- ProjectPCA(object = sample_i, do.print = FALSE)
  #PCHeatmap(object = sample_i, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, 
   #         label.columns = FALSE, use.full = FALSE)
  
  
  #Group the cells into clusters
  sample_i <- FindClusters(object = sample_i, reduction.type = "pca", dims.use = 1:15, 
                                  resolution = resolution_single, print.output = 0, save.SNN = TRUE,
                                  force.recalc = TRUE)
  #PrintFindClustersParams(object = sample_1_seurat)
  
  #perfrom Tsne
  sample_i <- RunTSNE(object = sample_i, dims.use = 1:13, do.fast = TRUE)
  
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
  #l_NO_LECs_BECs[[i]] <- sample_i.integrated_NO_LECs_BECs
  
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
saveRDS(obs.list_minus, paste(PATH_Output,"/obs.list_minus.Rds", sep=""))
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
#Scale data
mLN_Ontogeny_minus.merged <- ScaleData(mLN_Ontogeny_minus.merged)
#Run PCA
mLN_Ontogeny_minus.merged <- RunPCA(object = mLN_Ontogeny_minus.merged, pc.genes = mLN_Ontogeny_minus.merged@var.genes, do.print = TRUE, pcs.print = 1:10, 
                              genes.print = 10)

# Determine genes to use for CCA, must be highly variable in at least 2 datasets
genes.use <- c()
for (i in 1:length(obs.list_minus)) {
  genes.use <- c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 1000))
}
genes.use <- names(which(table(genes.use) > 1))
for (i in 1:length(obs.list_minus)) {
  genes.use <- genes.use[genes.use %in% rownames(obs.list_minus[[i]]@scale.data)]
}

# Run multi-set CCA
mLN_Ontogeny_minus.merged <- RunMultiCCA(obs.list_minus, genes.use = genes.use, num.ccs = dims_use)

# CC Selection
MetageneBicorPlot(mLN_Ontogeny_minus.merged, grouping.var = "tech")

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

#Plot the aligned CCA ACCA
p1 <- VlnPlot(mLN_Ontogeny_minus.merged, features.plot = "ACC1", group.by = "tech", do.return = T)
p2 <- VlnPlot(mLN_Ontogeny_minus.merged, features.plot = "ACC2", group.by = "tech", do.return = T)
plot_grid(p1, p2)


# t-SNE and Clustering
mLN_Ontogeny_minus.merged <- RunTSNE(mLN_Ontogeny_minus.merged,
                                   reduction.use = "cca",
                                   dims.use = 1:dims_use, save.SNN = T,
                                   resolution = resolution_merge, force.recalc = TRUE)
mLN_Ontogeny_minus.merged <- FindClusters(mLN_Ontogeny_minus.merged, reduction.type = "cca",
                                          dims.use = 1:dims_use, save.SNN = T, resolution = resolution_merge,
                                          force.recalc = TRUE)

#saveRDS(mLN_Ontogeny_minus.merged, paste(PATH_Output,"/MultiAlign_D0_D10_D24_D56_D300.Rds", sep=""))
mLN_Ontogeny_minus.merged <- readRDS(paste(PATH_Output,"/MultiAlign_D0_D10_D24_D56_D300.Rds", sep=""))
#Change cluster ID
current.cluster.ids <- c(0:14)


new.cluster.ids <- c("neo1","Inmt+Cxcl12+","Ccl19+Il7+","neo2",
                     "Il6+Cxcl1+","Cd34+Ackr3+","PvC","LTO",
                     "neo3","Cd34+Gdf10+","Cd34+Aldh1a2+","Cd34+Cd248+",
                     "?pSC?","ProgProf","SCx")

mLN_Ontogeny_minus.merged@ident <- plyr::mapvalues(x = mLN_Ontogeny_minus.merged@ident, from = current.cluster.ids, to = new.cluster.ids)
mLN_Ontogeny_minus.merged@meta.data$Cluster <- mLN_Ontogeny_minus.merged@ident

#Check allocation of cells
sample_1_cells <- grep(pattern = "^sample_1_", x = rownames(mLN_Ontogeny_minus.merged@meta.data), value = TRUE)
sample_2_cells <- grep(pattern = "^sample_2_", x = rownames(mLN_Ontogeny_minus.merged@meta.data), value = TRUE)
sample_3_cells <- grep(pattern = "^sample_3_", x = rownames(mLN_Ontogeny_minus.merged@meta.data), value = TRUE)
sample_4_cells <- grep(pattern = "^sample_4_", x = rownames(mLN_Ontogeny_minus.merged@meta.data), value = TRUE)
sample_5_cells <- grep(pattern = "^sample_5_", x = rownames(mLN_Ontogeny_minus.merged@meta.data), value = TRUE)

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
plot_grid(p1, p2, p3, p4, p5)
#, p6, p7)
plot_grid(p5)

FeaturePlot(mLN_Ontogeny_minus.merged, features.plot = c("Tnfsf11","Tnfsf13b","Acta2","Ccl19",
                                                       "Madcam1","Ltbr","Tnfsf13b","Il6", 
                                                       "Il7","Ackr3","Fn1","Col4a1", "Cxcl13",
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

#saveRDS(mLN_Ontogeny_minus.merged, paste(PATH_Output,"/MultiAlign_D0_D10_D24_D56.Rds", sep=""))
mLN_Ontogeny_minus.merged <- readRDS(paste(PATH_Output,"/MultiAlign_D0_D10_D24_D56_D300.Rds", sep=""))
saveRDS(mLN_Ontogeny_minus.merged, paste(PATH_Output,"/MultiAlign_D0_D10_D24_D56_D300.Rds", sep=""))

#####
#Obtain cell IDs minus Day0 excluded
#####
#Define clusters to be excluded
Day0_unique_clusters <- c("neo1","neo2","neo3")
Day_0_meta <- subset(mLN_Ontogeny_minus.merged@meta.data, tech == "d000")
cells_include <- row.names(subset(Day_0_meta, !(Cluster %in% Day0_unique_clusters)))
saveRDS(cells_include, paste(PATH_Output,"/Day0_cells_include_withoutNEOs.Rds", sep=""))


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
pdf(paste(PATH_Output,"/NC_Top40_D0_D10_D24_D56.pdf", sep=""))
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
write.table(mLN_Ontogeny_minus.merged.markers, paste(PATH_Output,"/DEGs_D0_D10_D24_D56_D300.csv", sep=""), dec=".", sep=",")
mLN_Ontogeny_minus.merged.markers <- read.table(paste(PATH_Output,"/DEGs_D0_D10_D24_D56_D300.csv", sep=""), dec=".", sep=",")
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
mLN_Ontogeny_minus.merged.markers_nonDup_10_FC <- mLN_Ontogeny_minus.merged.markers_nonDup %>% group_by(cluster) %>% top_n(10, avg_logFC)
mLN_Ontogeny_minus.merged.markers_nonDup_20_FC <- mLN_Ontogeny_minus.merged.markers_nonDup %>% group_by(cluster) %>% top_n(20, avg_logFC)

#Plot differential expression of Marker genes
title <- "Dotplot_Top10_D0_10_24_56_300"
#pdf(paste(path_output,"/",title, ".pdf", sep = ""), width = 25, height = 10)
DotPlot(mLN_Ontogeny_minus.merged, mLN_Ontogeny_minus.merged.markers_nonDup_10_FC$gene,
        x.lab.rot = TRUE, plot.legend = TRUE,
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

title <- paste(organ,condition,"Heatmap_hierarchical",sep = "_")
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
t_cells_per_condition <- do.call(rbind, lapply(l_cells_per_condition, function(x) {t(as.data.frame(x))}))
rownames(t_cells_per_condition) <- names(l_condition)
#Calculate Frequencies
total <- colSums(t_cells_per_condition)
freq_cells_per_condition <- (t(t_cells_per_condition) / total) * 100
freq_cells_per_condition <- as.matrix(freq_cells_per_condition)
pheatmap(freq_cells_per_condition, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = FALSE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("white", "grey","grey56","grey23","black"), space="rgb")(128),
         main = title)

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
Data <- Data_ready
List_GOs_interest <- c("GO:0006306")
List_GOs_interest_overrep_D0_developmental <- c("GO:1905562","GO:0001657","GO:0043627","GO:0048706","GO:0030326","GO:0090263","GO:0045668","GO:0072074","GO:2000288")
Genes_for_GOs_interest <- subset(GO2GeneID, ls(GO2GeneID) %in% List_GOs_interest)
number_input <- length(Genes_for_GOs_interest)
AllGenes_interest <- unlist(l_gene_pval)

#Save Output Tables in List
l_Genes_GO_per_cluster_interest <- list()
cluster_of_interest <- 0

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
  number_genes_expressed_in_GO <- nrow(subset(genes_expression, cluster == 0))
  
  if(number_genes_expressed_in_GO >= 1){
    #format data for pheatmap
    l_Genes_GO_per_cluster_interest[[i]] <- genes_expression
    names(l_Genes_GO_per_cluster_interest)[i] <- names(Genes_for_GOs_interest)[i]
  }
}





