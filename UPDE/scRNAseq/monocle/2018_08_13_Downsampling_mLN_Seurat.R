#Seurat 
library("Seurat")
library("cowplot")
library("Matrix")
library("magrittr")
library("dplyr")
library("pryr")
library("ggplot2")
library("gridExtra")
library("pheatmap")
library("DropletUtils")

#######
#Load and label
#######
###
#mLN
###
#Experiment 2
#4_... IDs mLN cells single
sample_1_2.data <- Read10X(data.dir = "/Users/Pezoldt/PowerFolders/R/2017_252_scRNAseq/Data/mLN_SPF/B6/outs/filtered_gene_bc_matrices/mm10")
#sample_1_2.data <- Read10X(data.dir = "C:/Users/jpe12/PowerFolders/R/2017_252_scRNAseq/Data/mLN_SPF/B6/filtered_gene_bc_matrices/mm10")
sample_1_2.data@Dimnames[[2]] <- paste("sample_4_2_", c(1:length(sample_1_2.data@Dimnames[[2]])), sep = "")
#Experiment 1
sample_1_1.data <- read.csv("/Users/Pezoldt/PowerFolders/R/2017_244_scRNASeq/Data/mLNSPF.csv")
#sample_1_1.data <- read.csv("C:/Users/jpe12/PowerFolders/R/2017_244_scRNASeq/Data/mLNSPF.csv")
row.names(sample_1_1.data) <- sample_1_1.data$X
sample_1_1.data <- sample_1_1.data[,c(2:ncol(sample_1_1.data))]
colnames(sample_1_1.data) <- paste("sample_4_1_", c(1:ncol(sample_1_1.data)), sep = "")

#######
#Global variables
#######
sample_1 <- "mLN"
#setwd("C:/Users/jpe12/PowerFolders/R/Paper/2018_Pezoldt_Pasztoi_Revision_NC/scRNA-seq/Output/mLN_single")

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


#########
#Setup sample_1
#########
#sample_1_seurat <- CreateSeuratObject(raw.data = sample_1.data, min.cells = min_cells, min.genes = min_gene_number)

#Merger
sample_1_1_seurat <- CreateSeuratObject(raw.data = sample_1_1.data, min.cells = min_cells, min.genes = min_gene_number)
sample_1_2_seurat <- CreateSeuratObject(raw.data = sample_1_2.data, min.cells = min_cells, min.genes = min_gene_number)
sample_1_seurat <- MergeSeurat(sample_1_1_seurat, sample_1_2_seurat, do.normalize = FALSE)
rm(sample_1_1_seurat, sample_1_2_seurat)
#Generate downsampled files
mLN_5 <- CreateSeuratObject(raw.data = downsampleMatrix(sample_1_seurat@raw.data, prop = 0.05, bycol = TRUE), 
                                            min.cells = min_cells, min.genes = min_gene_number)
mLN_10 <- CreateSeuratObject(raw.data = downsampleMatrix(sample_1_seurat@raw.data, prop = 0.10, bycol = TRUE), 
                                             min.cells = min_cells, min.genes = min_gene_number)
mLN_20 <- CreateSeuratObject(raw.data = downsampleMatrix(sample_1_seurat@raw.data, prop = 0.20, bycol = TRUE), 
                                             min.cells = min_cells, min.genes = min_gene_number)
mLN_50 <- CreateSeuratObject(raw.data = downsampleMatrix(sample_1_seurat@raw.data, prop = 0.50, bycol = TRUE), 
                                             min.cells = min_cells, min.genes = min_gene_number)
mLN_100 <- sample_1_seurat

#Run one sample at a time
sample_1_seurat <- mLN_100

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



TSNEPlot(object = sample_1_seurat_minus, pt.size = 0.3, do.label = TRUE)
#LN_minus_mLN <- readRDS("C:/Users/jpe12/Desktop/Seurat_objects/LN_minus_mLN.rds")
#sample_1_seurat_minus <- LN_minus_mLN

FeaturePlot(sample_1_seurat_minus, features.plot = c("Cxcl9","Tnfsf11","Madcam1","Nkain4","Vcam1",
                                        "Icam1","Bst1","Cxcl13","Ccl19","Krt18",
                                        "Cd248","Cd34","Il6","Inmt","Nr4a1","Aldh1a2","Gdf10"),
            no.legend = TRUE,
            cols.use = c('gainsboro', 'darkred'),
            pt.size = 0.2)


#####
#Find differentially expressed genes
#####
#Find markers for every cluster compared to all remaining cells
sample_1_seurat_minus.markers <- FindAllMarkers(object = sample_1_seurat_minus, only.pos = TRUE, min.pct = 0.25, 
                                   thresh.use = 0.25)
sample_1_seurat_minus.markers_10 <- sample_1_seurat_minus.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
sample_1_seurat_minus.markers_20 <- sample_1_seurat_minus.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
sample_1_seurat_minus.markers_40 <- sample_1_seurat_minus.markers %>% group_by(cluster) %>% top_n(40, avg_logFC)

#Save DEGs
#write.table(sample_1_seurat_minus.markers, paste(getwd(),"/",sample_1,"_DEGs.csv", sep=""), dec=".", sep=",")
write.table(sample_1_seurat_minus.markers_20, paste("/Users/Pezoldt/Desktop/Downsampling_mLN_Seurat","/",sample_1,"100percent_DEGs_Top20.csv", sep=""), dec=".", sep=",")
#write.table(sample_1_seurat_minus.markers_40, paste(getwd(),"/",sample_1,"_DEGs_Top40.csv", sep=""), dec=".", sep=",")



#####
#Dotplots
#####
#Top10
sample_1_seurat_minus.markers_nonDup <- sample_1_seurat_minus.markers_40[!duplicated(sample_1_seurat_minus.markers_40$gene),]
sample_1_seurat_minus.markers_nonDup_10_pValue <- sample_1_seurat_minus.markers_nonDup %>% group_by(cluster) %>% top_n(10, 1/p_val_adj)
sample_1_seurat_minus.markers_nonDup_20_pValue <- sample_1_seurat_minus.markers_nonDup %>% group_by(cluster) %>% top_n(20, 1/p_val_adj)
sample_1_seurat_minus.markers_nonDup_10_FC <- sample_1_seurat_minus.markers_nonDup %>% group_by(cluster) %>% top_n(10, avg_logFC)
sample_1_seurat_minus.markers_nonDup_20_FC <- sample_1_seurat_minus.markers_nonDup %>% group_by(cluster) %>% top_n(20, avg_logFC)


DotPlot(sample_1_seurat_minus, sample_1_seurat_minus.markers_nonDup_10_FC$gene,
        x.lab.rot = TRUE, plot.legend = TRUE,
        col.min = 0, col.max = 5, dot.scale = 5,
        cols.use = c("lightgrey", "brown"))

#####
#Hierarchical clustering according to Top DEGs per cluster
#####
#Average expression
sample_1_seurat_minus_aExp <- AverageExpression(sample_1_seurat_minus)
#LN_minus_aExp <- LN_minus_aExp[,1:12]
#Check for DEGs
#genes <- as.character(LN_minus.markers_20[1:260,]$gene)
sample_1_seurat_minus_aExp_TopDEGs <- subset(sample_1_seurat_minus_aExp, rownames(sample_1_seurat_minus_aExp) %in% sample_1_seurat_minus.markers$gene)
sample_1_seurat_minus_aExp_TopDEGs <- as.matrix(sample_1_seurat_minus_aExp_TopDEGs)
title <- paste("Hierarchical_clustering",sample_1, sep = "_")
pheatmap(sample_1_seurat_minus_aExp_TopDEGs, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 30, show_rownames = FALSE, cluster_cols = TRUE,
         scale = "row", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("blue", "white","red"), space="rgb")(128),
         main = title)
