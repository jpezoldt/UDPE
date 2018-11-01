#Seurat 
library("Seurat")
library("cowplot")
#library("Matrix")
#library("magrittr")
library("dplyr")
#not installed on Fameux
#library("pryr")
#library("svglite")

#######
#Load and label
#######
###
#mLN
###
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

###
#pLN
###
#Label mLN cells with 2_1 & 2_2
#Experiment 1 pLN_SPF
sample_2_1.data <- read.csv("/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/2018_NC_stroma/Exp_244/pLNSPF.csv")
#sample_2_1.data <- read.csv("/Users/Pezoldt/PowerFolders/R/2017_244_scRNASeq/Data/pLNSPF.csv")
row.names(sample_2_1.data) <- sample_2_1.data$X
sample_2_1.data <- sample_2_1.data[,c(2:ncol(sample_2_1.data))]
colnames(sample_2_1.data) <- paste("sample_2_1_", c(1:ncol(sample_2_1.data)), sep = "")

#Experiment 2 pLN_SPF
sample_2_2.data <- Read10X(data.dir = "/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/2018_NC_stroma/Exp_252/pLN_SPF/B5/filtered_gene_bc_matrices/mm10")
#sample_2_2.data <- Read10X(data.dir = "/Users/Pezoldt/PowerFolders/R/2017_252_scRNAseq/Data/pLN_SPF/B5/filtered_gene_bc_matrices/mm10")
sample_2_2.data@Dimnames[[2]] <- paste("sample_2_2_", c(1:length(sample_2_2.data@Dimnames[[2]])), sep = "")


#######
#Global variables
#######
#@User: Define sample names
#M == merge
sample_1 <- "mLN"
sample_2 <- "pLN"
PATH_output <- "/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/seurat/LN_Commnesals/mLN_pLN_SPF"
PATH_output_integrated <- "/home/pezoldt/NAS2/pezoldt/Analysis/01_Integrated/scRNA_TFs"

#Important: If changing global variable Setup of pLN and mLN needs to be manually adjusted as the cluster IDs get changed
min_gene_number <- 250
mito_cutoff <- 0.045
nGene <- 3500
min_cells <- 20
resolution_single <- 1.4

#Parameters FSC only
n_variable_genes <- 1500
resolution_merge <- 1.1
dims_use <- 22

#########
#Setup sample_1
#########
#Merger
sample_1_1_seurat <- CreateSeuratObject(raw.data = sample_1_1.data, min.cells = min_cells, min.genes = min_gene_number)
sample_1_2_seurat <- CreateSeuratObject(raw.data = sample_1_2.data, min.cells = min_cells, min.genes = min_gene_number)
sample_1_seurat <- MergeSeurat(sample_1_1_seurat, sample_1_2_seurat, do.normalize = FALSE)
rm(sample_1_1_seurat, sample_1_2_seurat, 
   sample_1_1.data, sample_1_2.data)
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

# We filter out cells according unique gene counts 
# and percentage mitochondrial reads
sample_1_seurat <- FilterCells(object = sample_1_seurat, subset.names = c("nGene", "percent.mito"), 
                   low.thresholds = c(min_gene_number, -Inf), high.thresholds = c(nGene, mito_cutoff))

#Normalize the expression and log-transform
sample_1_seurat <- NormalizeData(object = sample_1_seurat, normalization.method = "LogNormalize", 
                     scale.factor = 10000)

#Scale Data
sample_1_seurat <- ScaleData(sample_1_seurat)

#Find variable genes independent of expression level
sample_1_seurat <- FindVariableGenes(object = sample_1_seurat, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

#number of variable genes
length(x = sample_1_seurat@var.genes)

#scale data
sample_1_seurat <- ScaleData(object = sample_1_seurat, vars.to.regress = c("nUMI", "percent.mito"))

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

FeaturePlot(sample_1_seurat, features.plot = c("Pdpn", "Ackr4", "Pecam1","Ackr3"),
            min.cutoff = "q9",
            cols.use = c('lightgrey', 'blue'),
            pt.size = 0.1)

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
all_cells_keep <- subset(all_cells, res.1.4 == 1000)

#attain all cells per cluster ID and concatenate the tables
for(i in 1:length(cluster_keep)){
  cluster_id_i <- cluster_keep[i]
  print(cluster_id_i)
  cells_cluster_i <- subset(all_cells, res.1.4 == cluster_id_i)
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
all_cells_keep <- subset(all_cells, res.1.4 == 1000)

#attain all cells per cluster ID and concatenate the tables
for(i in 1:length(cluster_keep)){
  cluster_id_i <- cluster_keep[i]
  print(cluster_id_i)
  cells_cluster_i <- subset(all_cells, res.1.4 == cluster_id_i)
  print(nrow(cells_cluster_i))
  all_cells_keep <- rbind(all_cells_keep, cells_cluster_i)
}
sample_1_ONLY_LECs_BECs <- as.character(rownames(all_cells_keep))

#########
#Setup sample_2
#########
#Merger
sample_2_1_seurat <- CreateSeuratObject(raw.data = sample_2_1.data, min.cells = min_cells, min.genes = min_gene_number)
sample_2_2_seurat <- CreateSeuratObject(raw.data = sample_2_2.data, min.cells = min_cells, min.genes = min_gene_number)
sample_2_seurat <- MergeSeurat(sample_2_1_seurat, sample_2_2_seurat, do.normalize = FALSE)
rm(sample_2_1_seurat, sample_2_2_seurat,
   sample_2_1.data, sample_2_2.data)
#Eliminate mitochondrial high cells and dublets
#Identify cells with mitochondrial readss
mito.genes <- grep(pattern = "^mt-", x = rownames(sample_2_seurat@data), value = TRUE)
percent.mito <- colSums(as.matrix(sample_2_seurat@raw.data[mito.genes, ]))/colSums(as.matrix(sample_2_seurat@raw.data))

# AddMetaData adds columns to object@data.info and good to store QC stats
sample_2_seurat <- AddMetaData(object = sample_2_seurat, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = sample_2_seurat, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(object = sample_2_seurat, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = sample_2_seurat, gene1 = "nUMI", gene2 = "nGene")

#Filter out cells according unique gene counts and percentage mitochondrial reads
sample_2_seurat <- FilterCells(object = sample_2_seurat, subset.names = c("nGene", "percent.mito"), 
                   low.thresholds = c(min_gene_number, -Inf), high.thresholds = c(nGene, mito_cutoff))

#Normalize the expression and log-transform
sample_2_seurat <- NormalizeData(object = sample_2_seurat, normalization.method = "LogNormalize", 
                     scale.factor = 10000)

#Scale Data
sample_2_seurat <- ScaleData(sample_2_seurat)

#Find variable genes independent of expression level
sample_2_seurat <- FindVariableGenes(object = sample_2_seurat, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
#number of variable genes
length(x = sample_2_seurat@var.genes)

#scale data
sample_2_seurat <- ScaleData(object = sample_2_seurat, vars.to.regress = c("nUMI", "percent.mito"))

#perfrom PCA and print genes that define PCA
sample_2_seurat <- RunPCA(object = sample_2_seurat, pc.genes = sample_2_seurat@var.genes, do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5)

# Examine and visualize PCA results a few different ways
PrintPCA(object = sample_2_seurat, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = sample_2_seurat, pcs.use = 1:6)
PCAPlot(object = sample_2_seurat, dim.1 = 1, dim.2 = 2)

# ProjectPCA scores each gene in the dataset 
sample_2_seurat <- ProjectPCA(object = sample_2_seurat, do.print = FALSE)
PCHeatmap(object = sample_2_seurat, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

#Group the cells into clusters
sample_2_seurat <- FindClusters(object = sample_2_seurat, reduction.type = "pca", dims.use = 1:17, 
                    resolution = resolution_single, print.output = 0, save.SNN = TRUE,
                    force.recalc = TRUE)
PrintFindClustersParams(object = sample_2_seurat)

#perfrom Tsne
sample_2_seurat <- RunTSNE(object = sample_2_seurat, dims.use = 1:17, do.fast = TRUE)
TSNEPlot(object = sample_2_seurat, pt.size = 0.3, do.label = TRUE)

FeaturePlot(sample_2_seurat, features.plot = c("Pdpn", "Ackr4", "Pecam1","Ackr3"),
            min.cutoff = "q9",
            cols.use = c('lightgrey', 'blue'),
            pt.size = 0.1)

#####
#Eliminate LECs, BECs sample_2
#####
#Identify LECs and BECs
VlnPlot(object = sample_2_seurat, features.plot = c("Pecam1"))
sample_2_aExp <- AverageExpression(sample_2_seurat)
aExp_Pecam1 <- subset(sample_2_aExp, rownames(sample_2_aExp) == c("Pecam1"))
#Grab colnames/clusters that have Pecam1 expression < 1 on average per cluster
cluster_keep <- as.numeric(rownames(subset(as.data.frame(t(aExp_Pecam1)), Pecam1 < 1)))

#Get cells that are not within the cluster BEC and LEC cluster
all_cells <- sample_2_seurat@meta.data
all_cells_keep <- subset(all_cells, res.1.4 == 1000)

#attain all cells per cluster ID and concatenate the tables
for(i in 1:length(cluster_keep)){
  cluster_id_i <- cluster_keep[i]
  print(cluster_id_i)
  cells_cluster_i <- subset(all_cells, res.1.4 == cluster_id_i)
  print(nrow(cells_cluster_i))
  all_cells_keep <- rbind(all_cells_keep, cells_cluster_i)
}
sample_2_NO_LECs_BECs <- as.character(rownames(all_cells_keep))

#####
#Select LECs and BECs sample_2
#####
#Identify LECs and BECs
VlnPlot(object = sample_2_seurat, features.plot = c("Pecam1"))
sample_2_aExp <- AverageExpression(sample_2_seurat)
aExp_Pecam1 <- subset(sample_2_aExp, rownames(sample_2_aExp) == c("Pecam1"))
#Grab colnames/clusters that have Pecam1 expression < 1 on average per cluster
cluster_keep <- as.numeric(rownames(subset(as.data.frame(t(aExp_Pecam1)), Pecam1 > 1)))

#Get cells that are not within the cluster BEC and LEC cluster
all_cells <- sample_2_seurat@meta.data
all_cells_keep <- subset(all_cells, res.1.4 == 1000)

#attain all cells per cluster ID and concatenate the tables
for(i in 1:length(cluster_keep)){
  cluster_id_i <- cluster_keep[i]
  print(cluster_id_i)
  cells_cluster_i <- subset(all_cells, res.1.4 == cluster_id_i)
  print(nrow(cells_cluster_i))
  all_cells_keep <- rbind(all_cells_keep, cells_cluster_i)
}
sample_2_ONLY_LECs_BECs <- as.character(rownames(all_cells_keep))

########################################
#Merged analysis without LECs and BECs
#########################################
#####
#sample_1
#####
#Filter for mitochondiral and duplets and LECs and BECs out
sample_1_seurat_minus <- FilterCells(object = sample_1_seurat, subset.names = c("nGene", "percent.mito"), 
                                     low.thresholds = c(min_gene_number, -Inf), high.thresholds = c(nGene, mito_cutoff),
                                     cells.use = sample_1_NO_LECs_BECs)

#Eliminate cells that have expression
sample_1_seurat_minus <- FilterCells(object = sample_1_seurat_minus, subset.names = c("Ackr4","Pecam1"),
                                     low.thresholds = c(-Inf, -Inf), high.thresholds = c(1,1))
sample_1_seurat_minus <- NormalizeData(sample_1_seurat_minus)


#Identify cells with high rpl reads
rpl.genes <- grep(pattern = "^Rpl", x = rownames(sample_1_seurat_minus@data), value = TRUE)
percent.rpl <- colSums(as.matrix(sample_1_seurat_minus@raw.data[rpl.genes, ]))/colSums(as.matrix(sample_1_seurat_minus@raw.data))

# AddMetaData for RPL percentage
sample_1_seurat_minus <- AddMetaData(object = sample_1_seurat_minus, metadata = percent.rpl, col.name = "percent.rpl")

#Identify cells with mitochondrial readss
mito.genes <- grep(pattern = "^mt-", x = rownames(sample_1_seurat_minus@data), value = TRUE)
percent.mito <- colSums(as.matrix(sample_1_seurat_minus@raw.data[mito.genes, ]))/colSums(as.matrix(sample_1_seurat_minus@raw.data))

# AddMetaData for percent mitochoncrial reads
sample_1_seurat_minus <- AddMetaData(object = sample_1_seurat_minus, metadata = percent.mito, col.name = "percent.mito")

#Scale data
sample_1_seurat_minus <- ScaleData(object = sample_1_seurat_minus, vars.to.regress = c("nUMI", "percent.mito","percent.rpl"))
cells.to.sample <- length(sample_1_seurat_minus@cell.names)

#####
#sample_2
#####
#Filter for mitochondiral and duplets and LECs and BECs out
sample_2_seurat_minus <- FilterCells(object = sample_2_seurat, subset.names = c("nGene", "percent.mito"), 
                   low.thresholds = c(min_gene_number, -Inf), high.thresholds = c(nGene, mito_cutoff),
                   cells.use = sample_2_NO_LECs_BECs)

#Eliminate cells that have expression
sample_2_seurat_minus <- FilterCells(object = sample_2_seurat_minus, subset.names = c("Ackr4","Pecam1"),
                   low.thresholds = c(-Inf, -Inf), high.thresholds = c(1,1))
sample_2_seurat_minus <- NormalizeData(sample_2_seurat_minus)


#Identify cells with high rpl reads
rpl.genes <- grep(pattern = "^Rpl", x = rownames(sample_2_seurat_minus@data), value = TRUE)
percent.rpl <- colSums(as.matrix(sample_2_seurat_minus@raw.data[rpl.genes, ]))/colSums(as.matrix(sample_2_seurat_minus@raw.data))

# AddMetaData for RPL percentage
sample_2_seurat_minus <- AddMetaData(object = sample_2_seurat_minus, metadata = percent.rpl, col.name = "percent.rpl")

#Identify cells with mitochondrial readss
mito.genes <- grep(pattern = "^mt-", x = rownames(sample_2_seurat_minus@data), value = TRUE)
percent.mito <- colSums(as.matrix(sample_2_seurat_minus@raw.data[mito.genes, ]))/colSums(as.matrix(sample_2_seurat_minus@raw.data))

# AddMetaData for percent mitochoncrial reads
sample_2_seurat_minus <- AddMetaData(object = sample_2_seurat_minus, metadata = percent.mito, col.name = "percent.mito")

#Scale data
sample_2_seurat_minus <- ScaleData(object = sample_2_seurat_minus, vars.to.regress = c("nUMI", "percent.mito","percent.rpl"))

#Find variable genes
sample_2_seurat_minus <- FindVariableGenes(sample_2_seurat_minus, do.plot = F)

#Align cell numbers across samples
set.seed(111)
sampled.cells <- sample(x = sample_2_seurat_minus@cell.names, size = cells.to.sample, replace = F)
sample_2_seurat_minus <- SubsetData(object = sample_2_seurat_minus, cells.use = sampled.cells)

#####
#Delete seurat Objects
#####
#rm(sample_1_seurat, sample_2_seurat)
#####
#Alignment
#####
###
#take top variable genes
hvg.sample_1_minus <- rownames(head(sample_1_seurat_minus@hvg.info, n_variable_genes))
hvg.sample_2_minus <- rownames(head(sample_2_seurat_minus@hvg.info, n_variable_genes))
####
hvg.union <- union(hvg.sample_1_minus, hvg.sample_2_minus)

#set labels
sample_1_seurat_minus@meta.data[,"protocol"] <- paste(sample_1, "_minus", sep="")
sample_2_seurat_minus@meta.data[,"protocol"] <- paste(sample_2, "_minus", sep="")

#Canonical Correlation Vectors
LN_minus <- RunCCA(sample_2_seurat_minus, sample_1_seurat_minus, genes.use = hvg.union, num.cc = 30)

#visualize results of CCA
p1 <- DimPlot(LN_minus, reduction.use = "cca", group.by = "protocol", pt.size = 0.5, do.return = T)
p2 <- VlnPlot(LN_minus, features.plot = "CC1", group.by = "protocol", do.return = T)
plot_grid(p1, p2)

#Search for cells whose expression profile cannot be well-explained by low-dimensional CCA
#LN_minus <- CalcVarExpRatio(LN_minus, reduction.type = "ica", grouping.var = "protocol", dims.use = 1:dims_use)

#align acording to CCA subspaces
LN_minus <- AlignSubspace(LN_minus, reduction.type = "cca", grouping.var = "protocol", dims.align = 1:dims_use)

#Plot the aligned CCA ACCA
p1 <- VlnPlot(LN_minus, features.plot = "ACC3", group.by = "protocol", do.return = T)
p2 <- VlnPlot(LN_minus, features.plot = "ACC4", group.by = "protocol", do.return = T)
plot_grid(p1, p2)

#Run single integrated analysis on all cells
LN_minus <- RunTSNE(LN_minus, reduction.use = "cca.aligned", dims.use = 1:22)
LN_minus <- FindClusters(LN_minus, reduction.type = "cca.aligned",
                         dims.use = 1:dims_use, save.SNN = T, resolution = resolution_merge,
                         force.recalc = TRUE)
#Change cluster ID
current.cluster.ids <- c(0:13)
  

new.cluster.ids <- c("Ccl19+Il7+","Nr4a1+","Inmt+","CD34+Gdf10+",
                      "Ccl19+","CD34+Aldh1a2+","CD34+Ackr3+","Il6+Cxcl1+",
                      "pSC","Ccl19high","CD34+CD248+","Cxcl9+",
                      "PvC","mSC")

LN_minus@ident <- plyr::mapvalues(x = LN_minus@ident, from = current.cluster.ids, to = new.cluster.ids)

#Save Object
saveRDS(LN_minus, file = filename)

sample_1_minus_cells <- grep(pattern = "^sample_1_", x = rownames(LN_minus@meta.data), value = TRUE)
sample_2_minus_cells <- grep(pattern = "^sample_2_", x = rownames(LN_minus@meta.data), value = TRUE)
p1 <- TSNEPlot(LN_minus, group.by = "protocol", do.return = T, pt.size = 0.4, cells.use = sample_2_minus_cells)
p2 <- TSNEPlot(LN_minus, group.by = "protocol", do.return = T, pt.size = 0.4, cells.use = sample_1_minus_cells)
p3 <- TSNEPlot(LN_minus, do.return = F, pt.size = 0.2, do.label = TRUE)
p4 <- TSNEPlot(LN_minus, group.by = "protocol", do.return = T, pt.size = 0.2, colors.use = c("green","blue"))
par(mfrow = c(2, 2))
plot_grid(p1, p2, p3, p4)

FeaturePlot(LN_minus, features.plot = c("Cxcl9","Tnfsf11","Madcam1","Il7","Vcam1",
                                        "Icam1","Bst1","Cxcl13","Ccl19","Fabp7",
                                        "Cd248","Cd34","Il6","Inmt","Nr4a1","Aldh1a2","Gdf10"),
            no.legend = TRUE,
            cols.use = c('gainsboro', 'darkred'),
            pt.size = 0.2)

#mLN Open not a DEG (TF motif)
FeaturePlot(LN_minus, features.plot = c("Atf5","E2f1","Ebf2","Ebf3","Ebf4",
                                        "Elf1","Elk1","Hoxb4","Isl1"),
            no.legend = TRUE,
            cols.use = c('gainsboro', 'darkred'),
            pt.size = 0.2)

#pLN Peak UP (TF motif)
FeaturePlot(LN_minus, features.plot = c("Bhlhe40","Egr1",
                                        "Gata4","Irf3",
                                        "Irf5","Irf7","Irf8","Nfkb1",
                                        "Nfkb2","Stat2",
                                        "Stat6"),
            no.legend = TRUE,
            cols.use = c('gainsboro', 'darkred'),
            pt.size = 0.2)



VlnPlot(LN_minus,c("Ccl19","Icam1"), x.lab.rot = TRUE)

saveRDS(LN_minus, paste(PATH_output,"/LN_minus_mLN_pLN_SPF.Rds",sep = ""))

#####
#Find differentially expressed genes
#####
#Find markers for every cluster compared to all remaining cells
LN_minus.markers <- FindAllMarkers(object = LN_minus, only.pos = TRUE, min.pct = 0.25, 
                                   thresh.use = 0.25)
LN_minus.markers_10 <- LN_minus.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
LN_minus.markers_20 <- LN_minus.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
LN_minus.markers_40 <- LN_minus.markers %>% group_by(cluster) %>% top_n(40, avg_logFC)
saveRDS(LN_minus.markers, paste(PATH_output,"/mLN_pLN_SPF_DEG_all.Rds",sep = ""))


#top 20 non-duplicated
LN_minus.markers_dup

#Top10
LN_minus.markers_nonDup <- LN_minus.markers_40[!duplicated(LN_minus.markers_40$gene),]
LN_minus.markers_nonDup_10 <- LN_minus.markers_nonDup %>% group_by(cluster) %>% top_n(10, 1/p_val_adj)
LN_minus.markers_nonDup_15 <- LN_minus.markers_nonDup %>% group_by(cluster) %>% top_n(15, 1/p_val_adj)
LN_minus.markers_nonDup_20 <- LN_minus.markers_nonDup %>% group_by(cluster) %>% top_n(20, 1/p_val_adj)


#Save DEGs
#write.table(LN_minus.markers, paste(getwd(),"/",sample_1,"_",sample_2,"_DEGs.csv", sep=""), dec=".", sep=",")
#write.table(LN_minus.markers_20, paste(getwd(),"/",sample_1,"_",sample_2,"_DEGs_Top20.csv", sep=""), dec=".", sep=",")
#write.table(LN_minus.markers_40, paste(getwd(),"/",sample_1,"_",sample_2,"_DEGs_Top40.csv", sep=""), dec=".", sep=",")


#Plot differential expression of Marker genes
subset_order <- c("PvC",
                  "CD34+CD248+","CD34+Aldh1a2+","CD34+Ackr3+","CD34+Gdf10+",
                  "Il6+Cxcl1+",
                  "Ccl19high","Ccl19+Il7+","Ccl19+","Nr4a1+",
                  "Cxcl9+",
                  "Inmt+","pSC","mSC")

DoHeatmap(object = LN_minus, use.scaled = TRUE, group.by = "ident",
          genes.use = LN_minus.markers_nonDup_20$gene, title = paste(sample_1, "&", sample_2, sep=" "),
          remove.key = FALSE,
          group.order = subset_order,
          col.low = "white",
          col.mid =  "oldlace",
          col.high = "brown")



####################
#Obtain DEGs per cluster
####################
PATH_input_scRNAseq <- "/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/seurat/LN_Commnesals/mLN_pLN_SPF"
LN_minus <- readRDS(paste(PATH_input_scRNAseq,"/LN_minus_mLN_pLN_SPF.Rds",sep = ""))
#Setup SEURAT Object where mLN and pLN Clusters are separate
LN_minus@meta.data$res.0.9 <- as.character(LN_minus@ident)
mLN_minus_cells <- grep(pattern = "^sample_1_", x = rownames(LN_minus@meta.data), value = TRUE)
pLN_minus_cells <- grep(pattern = "^sample_2_", x = rownames(LN_minus@meta.data), value = TRUE)
LN_minus_p <- SubsetData(object = LN_minus,
                         cells.use = pLN_minus_cells)
LN_minus_m <- SubsetData(object = LN_minus,
                         cells.use = mLN_minus_cells)

#Merge mLN and pLN Seurat objects
LN_minus_p_m <- MergeSeurat(LN_minus_p, LN_minus_m, do.normalize = FALSE)
#change identity
LN_minus_p_m <- SetIdent(LN_minus_p_m, ident.use = c(paste(sample_2,"_",LN_minus_p@meta.data$res.0.9, sep=""), paste(sample_1, "_",LN_minus_m@meta.data$res.0.9, sep="")))

#Normalize the expression and log-transform
LN_minus_p_m <- NormalizeData(object = LN_minus_p_m, normalization.method = "LogNormalize", 
                              scale.factor = 10000)

#Scale Data
LN_minus_p_m <- ScaleData(LN_minus_p_m)

#Find variable genes independent of expression level
LN_minus_p_m <- FindVariableGenes(object = LN_minus_p_m, mean.function = ExpMean, dispersion.function = LogVMR, 
                                  x.low.cutoff = 0.0125, do.plot = FALSE, x.high.cutoff = 3, y.cutoff = 0.5)
#number of variable genes
length(x = LN_minus_p_m@var.genes)

#scale data
LN_minus_p_m <- ScaleData(object = LN_minus_p_m, vars.to.regress = c("nUMI", "percent.mito"))

#perfrom PCA and print genes that define PCA
LN_minus_p_m <- RunPCA(object = LN_minus_p_m, pc.genes = LN_minus_p_m@var.genes, do.print = TRUE, pcs.print = 1:5, 
                       genes.print = 5)

# Examine and visualize PCA results a few different ways
PrintPCA(object = LN_minus_p_m, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = LN_minus_p_m, pcs.use = 1:6)

# ProjectPCA scores each gene in the dataset 
LN_minus_p_m <- ProjectPCA(object = LN_minus_p_m, do.print = FALSE)



#Compare Cluster by Cluster
#Store DEGs in List
n_cluster = length(levels(as.factor(LN_minus_p_m@meta.data$res.0.9)))
DEG_l <- list()
DEG_l_adjPval <- list()
i = 1

#test = "bimod"
# Note: no calls
#test = "roc"

# Note: relatively consistent results
#test = "MAST"
#test = "negbinom"
#test = "wilcox"
#test = "tobit"

# Note: more DEGs, but core similar to other clusters
test = "DESeq2"
for(i in i:n_cluster){
  #cluster ID
  #k = i - 1
  cluster_ID_m <- paste(sample_1,"_",levels(as.factor(LN_minus_m@meta.data$res.0.9))[i], sep="")
  cluster_ID_p <- paste(sample_2,"_",levels(as.factor(LN_minus_m@meta.data$res.0.9))[i], sep="")
  
  #Identify DEGs avg_diff>0 (high in mLN) avg_Diff<0 (high in pLN)
  DEGs <- FindMarkers(LN_minus_p_m, cluster_ID_m, cluster_ID_p, logfc.threshold = 0.5, 
                      test.use = test)
  GeneSymbol <- rownames(DEGs)
  DEGs <- cbind(DEGs, GeneSymbol)
  print(cluster_ID_m)
  DEGs_adjPval <- subset(DEGs, p_val_adj < 0.05) 
  print(nrow(DEGs))
  print(nrow(DEGs_adjPval))
  
  DEG_l[[i]] <- DEGs
  DEG_l_adjPval[[i]]<- DEGs_adjPval
}
names(DEG_l) <- levels(as.factor(LN_minus_p_m@meta.data$res.0.9))
names(DEG_l_adjPval) <- levels(as.factor(LN_minus_p_m@meta.data$res.0.9))

#write Table
t_DEG <- data.frame(matrix(ncol=7,nrow=0))
colnames(t_DEG) <- c("p_val","avg_logFC","pct.1","pct.2","p_val_adj","cluster_i","GeneSymbol")

for(i in 1:length(DEG_l_adjPval)){
  DEG_l_adjPval_i <- DEG_l_adjPval[[i]]
  print(i)
  cluster_i <- rep(names(DEG_l_adjPval)[i], nrow(DEG_l_adjPval_i))
  print(cluster_i[1])
  DEG_l_adjPval_i <- cbind(DEG_l_adjPval_i, cluster_i)
  head(DEG_l_adjPval_i)
  t_DEG <- rbind(t_DEG,DEG_l_adjPval_i)
  
}

write.table(t_DEG, paste(PATH_output,"/",sample_1,"_",sample_2,"_DEGs_per_cluster_",test,".csv", sep=""), dec=".", sep=",")



#####
#Plot TSNE Maps
#####
#Make standard TSNE maps
gene_list_map_on_TSNE <- c("Acta2","Ackr3","Aldh1a2","Apoe",
                           "Bst1",
                           "Ccl2","Ccl5","Ccl7","Ccl8","Ccl9","Ccl11","Ccl19","Ccl21a",
                           "Cxcl1","Cxcl2","Cxcl13","Cxcl9","Cxcl10",
                           "Cd34","Cd55","Cd248",
                           "Ch25h",
                           "Csf1",
                           "Col15a1","Col4a1","Col3a1",
                           "Des",
                           "Enpp2",
                           "Gdf10",
                           "Il7","Il6","Icam1","Inmt",
                           "Ly6a","Ly6c1",
                           "Ntrk1","Nr4a1",
                           "Madcam1","Malat1","Mfge8",
                           "Nos2",
                           "Ptgis",
                           "Sfrp4",
                           "Tnfsf13b","Tnfsf11", "Tinagl1",
                           "Timd4",
                           "Vcam1")
                           
for(i in 1:length(gene_list_map_on_TSNE)){
  setwd("C:/Users/jpe12/Desktop/scRNAseq_NewInput/TSNE_map")
  gene_i <- gene_list_map_on_TSNE[i]
  print(gene_i)
  filename_i <- paste(sample_1, "_", sample_2, "_", gene_i, ".eps", sep="")
  print(filename_i)
  FeaturePlot(LN_minus, features.plot = gene_i,
            no.legend = TRUE,
            cols.use = c('gainsboro', 'darkred'),
            pt.size = 0.2)
  ggsave(filename_i, width = 5.8, height = 6, units = "cm")
}
for(i in 1:length(gene_list_map_on_TSNE)){
  setwd("C:/Users/jpe12/Desktop/scRNAseq_NewInput/TSNE_map")
  gene_i <- gene_list_map_on_TSNE[i]
  print(gene_i)
  filename_i <- paste(sample_1, "_", sample_2, "_", gene_i,"_scale" , ".eps", sep="")
  print(filename_i)
  FeaturePlot(LN_minus, features.plot = gene_i,
              no.legend = FALSE,
              cols.use = c('gainsboro', 'darkred'),
              pt.size = 0.2)
  ggsave(filename_i, width = 5.8, height = 6, units = "cm")
}

#####
#Dotplots
#####
#Top10
LN_minus.markers_nonDup <- LN_minus.markers_40[!duplicated(LN_minus.markers_40$gene),]
LN_minus.markers_nonDup_10 <- LN_minus.markers_nonDup %>% group_by(cluster) %>% top_n(10, 1/p_val_adj)
LN_minus.markers_nonDup_20 <- LN_minus.markers_nonDup %>% group_by(cluster) %>% top_n(20, 1/p_val_adj)

DotPlot(LN_minus, LN_minus.markers_nonDup_10$gene, x.lab.rot = TRUE, plot.legend = TRUE,
        col.min = 0, col.max = 2.5, dot.scale = 4)

#####
#Hierarchical clustering according to Top DEGs per cluster
#####
#Average expression
LN_minus_aExp <- AverageExpression(LN_minus)
#LN_minus_aExp <- LN_minus_aExp[,1:12]
#Check for DEGs
#genes <- as.character(LN_minus.markers_20[1:260,]$gene)
LN_minus_aExp_TopDEGs <- subset(LN_minus_aExp, rownames(LN_minus_aExp) %in% LN_minus.markers$gene)
LN_minus_aExp_TopDEGs <- as.matrix(LN_minus_aExp_TopDEGs)
pheatmap(LN_minus_aExp_TopDEGs, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 30, show_rownames = FALSE, cluster_cols = TRUE,
         scale = "row", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("blue", "white","red"), space="rgb")(128))

#####
#FeaturePlot
#####
#Genes were chosen on their putative or described impact on environment
genes_interest <- c("Pdgfrb", "Aldh1a3","Ptgis","Sfrp4","Bmper","Cd34","Tcf21")
#genes_interest = c("Aldh1a2", "Aldh1a3", "Tcf21","Bmper","Sfrp4","Ptgis","Fez1")
FeatureHeatmap(object = LN_minus, features.plot = genes_interest, 
               group.by = "protocol", sep.scale = FALSE, pt.size = 0.1,
               cols.use = c("gray87","darkred"),
               min.exp = 0.5, max.exp = 5)





####################
#Expression per mLN and pLN Cluster
####################
#Setup SEURAT Object where mLN and pLN Clusters are separate
LN_minus@meta.data$res.0.9 <- as.character(LN_minus@ident)
mLN_minus_cells <- grep(pattern = "^sample_1_", x = rownames(LN_minus@meta.data), value = TRUE)
pLN_minus_cells <- grep(pattern = "^sample_2_", x = rownames(LN_minus@meta.data), value = TRUE)
LN_minus_p <- SubsetData(object = LN_minus,
                         cells.use = pLN_minus_cells)
LN_minus_m <- SubsetData(object = LN_minus,
                         cells.use = mLN_minus_cells)

#Merge mLN and pLN Seurat objects
LN_minus_p_m <- MergeSeurat(LN_minus_p, LN_minus_m, do.normalize = FALSE)
#change identity
LN_minus_p_m <- SetIdent(LN_minus_p_m, ident.use = c(paste("p_",LN_minus_p@meta.data$res.0.9, sep=""), paste("m_",LN_minus_m@meta.data$res.0.9, sep="")))

#Normalize the expression and log-transform
LN_minus_p_m <- NormalizeData(object = LN_minus_p_m, normalization.method = "LogNormalize", 
                              scale.factor = 10000)

#Scale Data
LN_minus_p_m <- ScaleData(LN_minus_p_m)

#Find variable genes independent of expression level
LN_minus_p_m <- FindVariableGenes(object = LN_minus_p_m, mean.function = ExpMean, dispersion.function = LogVMR, 
                                  x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
#number of variable genes
length(x = LN_minus_p_m@var.genes)

#scale data
LN_minus_p_m <- ScaleData(object = LN_minus_p_m, vars.to.regress = c("nUMI", "percent.mito"))

#perfrom PCA and print genes that define PCA
LN_minus_p_m <- RunPCA(object = LN_minus_p_m, pc.genes = LN_minus_p_m@var.genes, do.print = TRUE, pcs.print = 1:5, 
                       genes.print = 5)

# Examine and visualize PCA results a few different ways
PrintPCA(object = LN_minus_p_m, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = LN_minus_p_m, pcs.use = 1:6)

# ProjectPCA scores each gene in the dataset 
LN_minus_p_m <- ProjectPCA(object = LN_minus_p_m, do.print = FALSE)



#Compare Cluster by Cluster
#Store DEGs in List
n_cluster = length(levels(as.factor(LN_minus_p_m@meta.data$res.0.9)))
DEG_l <- list()
DEG_l_adjPval <- list()
i = 1
for(i in i:n_cluster){
  #cluster ID
  #k = i - 1
  cluster_ID_m <- paste("m_",levels(as.factor(LN_minus_m@meta.data$res.0.9))[i], sep="")
  cluster_ID_p <- paste("p_",levels(as.factor(LN_minus_m@meta.data$res.0.9))[i], sep="")
  
  #Identify DEGs avg_diff>0 (high in mLN) avg_Diff<0 (high in pLN)
  DEGs <- FindMarkers(LN_minus_p_m, cluster_ID_m, cluster_ID_p, logfc.threshold = 0.5, 
                      test.use = "tobit")
  print(cluster_ID_m)
  DEGs_adjPval <- subset(DEGs, p_val_adj < 0.05) 
  print(nrow(DEGs))
  print(nrow(DEGs_adjPval))
  
  DEG_l[[i]] <- DEGs
  DEG_l_adjPval[[i]]<- DEGs_adjPval
}
names(DEG_l) <- levels(as.factor(LN_minus_p_m@meta.data$res.0.9))
names(DEG_l_adjPval) <- levels(as.factor(LN_minus_p_m@meta.data$res.0.9))
#####
#Check overlap with imprinted genes
#####
#Imprinted
#mLN maintained
#mLN_main <- read.delim("/Users/Pezoldt/PowerFolders/R/2017_scRNASeq/Data/Input_Tx/mLN_maintained.txt")
mLN_main <- read.delim("C:/Users/jpe12/PowerFolders/R/Paper/2018_Pezoldt_Pasztoi_Revision_NC/Input/Input_Tx/mLN_maintained.txt")
mLN_main <- as.character(mLN_main$GeneSymbol)
#mLN repressed
#mLN_rep <- read.delim("/Users/Pezoldt/PowerFolders/R/2017_scRNASeq/Data/Input_Tx/mLN_repressed.txt")
mLN_rep <- read.delim("C:/Users/jpe12/PowerFolders/R/Paper/2018_Pezoldt_Pasztoi_Revision_NC/Input/Input_Tx/mLN_repressed.txt")
mLN_rep <- as.character(mLN_rep$GeneSymbol)

#Check overlap of DEGs with imprinted genes
#mLN maintained
DEG_l_mLN_main <- list()
for(i in 1:length(DEG_l)){
  DEG_l_i <- DEG_l[[i]]
  DEG_l_mLN_main_i <- subset(DEG_l_i, rownames(DEG_l_i) %in% mLN_main )
  DEG_l_mLN_main[[i]] <- DEG_l_mLN_main_i
  print(paste("cluster", i, sep=""))
  print(nrow(DEG_l_mLN_main_i))
}

DEG_l_mLN_rep <- list()
for(i in 1:length(DEG_l)){
  DEG_l_i <- DEG_l[[i]]
  DEG_l_mLN_rep_i <- subset(DEG_l_i, rownames(DEG_l_i) %in% mLN_rep )
  DEG_l_mLN_rep[[i]] <- DEG_l_mLN_rep_i
  print(paste("cluster", i, sep=""))
  print(nrow(DEG_l_mLN_rep_i))
}

####################
#Zscore over expression of all genes from genemodules
####################
library(pheatmap)

#Average Expression per cluster
LN_minus_p_m_aExp <- AverageExpression(LN_minus_p_m)

#Maintained genes mLN
LN_minus_p_m_aExp_mLN_main <- subset(LN_minus_p_m_aExp, rownames(LN_minus_p_m_aExp) %in% mLN_main)
pheatmap(LN_minus_p_m_aExp_mLN_main, scale = "row", cluster_cols = TRUE, border_color = "black", cellwidth = 10,
         treeheight_row = 0, treeheight_col = 20,
         cellheigth = 10, color = colorRampPalette(c("white", "oldlace","brown"), space="rgb")(128))
Scale_LN_minus_p_m_aExp_mLN_main <- apply(LN_minus_p_m_aExp_mLN_main, 1, scale)
Zscore_mLN_main <- rowSums(Scale_LN_minus_p_m_aExp_mLN_main)
Zscore_mLN_main_norm <- Zscore_mLN_main / (nrow(LN_minus_p_m_aExp_mLN_main)/10)
Zscore_mLN_main_norm_mLN <- Zscore_mLN_main_norm[1:((length(Zscore_mLN_main_norm))/2)]
Zscore_mLN_main_norm_pLN <- Zscore_mLN_main_norm[((length(Zscore_mLN_main_norm))/2+1):length(Zscore_mLN_main_norm)]

#Only mLN
LN_minus_p_m_aExp_mLN_main_mLNONly <- subset(LN_minus_p_m_aExp, rownames(LN_minus_p_m_aExp) %in% mLN_main)[,c(1:14)]
LN_minus_p_m_aExp_mLN_main_mLNONly <- as.matrix(LN_minus_p_m_aExp_mLN_main_mLNONly)

#Eliminate rows with only zeros
LN_minus_p_m_aExp_mLN_main_mLNONly <- LN_minus_p_m_aExp_mLN_main_mLNONly[apply(LN_minus_p_m_aExp_mLN_main_mLNONly[,-1], 1, function(x) !all(x==0)),]

LN_minus_p_m_aExp_mLN_main_mLNONly <- LN_minus_p_m_aExp_mLN_main_mLNONly[,paste("m_",subset_order, sep="")]
pheatmap(LN_minus_p_m_aExp_mLN_main_mLNONly, scale = "row", cluster_cols = FALSE, border_color = "black", cellwidth = 10,
         treeheight_row = 0, treeheight_col = 20,
         cellheigth = 10, color = colorRampPalette(c("white", "oldlace","brown"), space="rgb")(128))

#Repressed genes mLN
LN_minus_p_m_aExp_mLN_rep <- subset(LN_minus_p_m_aExp, rownames(LN_minus_p_m_aExp) %in% mLN_rep)
pheatmap(LN_minus_p_m_aExp_mLN_rep, scale = "row", cluster_cols = TRUE)
Scale_LN_minus_p_m_aExp_mLN_rep <- apply(LN_minus_p_m_aExp_mLN_rep, 1, scale)
Zscore_mLN_rep <- rowSums(Scale_LN_minus_p_m_aExp_mLN_rep)
Zscore_mLN_rep_norm <- Zscore_mLN_rep / (nrow(LN_minus_p_m_aExp_mLN_rep)/10)
Zscore_mLN_rep_norm_mLN <- Zscore_mLN_rep_norm[1:((length(Zscore_mLN_rep_norm))/2)]
Zscore_mLN_rep_norm_pLN <- Zscore_mLN_rep_norm[((length(Zscore_mLN_rep_norm))/2+1):length(Zscore_mLN_rep_norm)]

#All Scores
Zscore_table <- data.frame(module = c(rep("mLN_main",length(Zscore_mLN_main_norm)),
                                      rep("mLN_rep",length(Zscore_mLN_rep_norm))),
                           cluster = rep(colnames(LN_minus_p_m_aExp_mLN_main),2),
                           cZscore = c(Zscore_mLN_main_norm, Zscore_mLN_rep_norm))

#All clusters
ggplot(data = Zscore_table, aes(x = cluster, y = cZscore, fill = module)) +
  geom_bar(stat = "identity", colour = "black") +
  theme(axis.text.x=element_text(angle=270,hjust=1)) +
  scale_fill_manual(values = c("deepskyblue1", "deeppink"))

#mLN_scores Scores
Zscore_table_mLN <- data.frame(module = c(rep("mLN_main",length(Zscore_mLN_main_norm)/2),
                                          rep("mLN_rep",length(Zscore_mLN_rep_norm)/2)),
                               cluster = rep(colnames(LN_minus_p_m_aExp_mLN_main)[1:(length(colnames(LN_minus_p_m_aExp_mLN_main))/2)]),
                               cZscore = c(Zscore_mLN_main_norm_mLN, Zscore_mLN_rep_norm_mLN))

#Sort Zscore table according to DEG Top20 Heatmap
subset_names <- paste("m_",subset_order,
                      sep = "")
Zscore_table_mLN_upper <- Zscore_table_mLN[c(1:(nrow(Zscore_table_mLN)/2)),]
Zscore_table_mLN_upper <- Zscore_table_mLN_upper[match(subset_names, Zscore_table_mLN_upper$cluster),]
Zscore_table_mLN_upper$cluster <- paste(letters[1:nrow(Zscore_table_mLN_upper)],"_",Zscore_table_mLN_upper$cluster,sep="")

Zscore_table_mLN_lower <- Zscore_table_mLN[c((nrow(Zscore_table_mLN)/2+1):nrow(Zscore_table_mLN)),]
Zscore_table_mLN_lower <- Zscore_table_mLN_lower[match(subset_names, Zscore_table_mLN_lower$cluster),]
Zscore_table_mLN_lower$cluster <- paste(letters[1:nrow(Zscore_table_mLN_lower)],"_",Zscore_table_mLN_lower$cluster,sep="")

Zscore_table_mLN <- rbind(Zscore_table_mLN_upper, Zscore_table_mLN_lower)

#Only mLN
ggplot(data = Zscore_table_mLN, aes(x = cluster, y = cZscore, fill = module)) +
  geom_bar(stat = "identity", colour = "black") +
  theme(axis.text.x=element_text(angle=270,hjust=1)) +
  scale_fill_manual(values = c("deepskyblue1", "deeppink"))

###################################
#Merged analysis with LECs and BECs
###################################
#####
#Alignment
#####
#resolution set
resolution_merge = 1.4
dims_use = 26
#####
#take top variable genes
hvg.sample_1 <- rownames(head(sample_1_seurat@hvg.info, 1500))
hvg.sample_2 <- rownames(head(sample_2_seurat@hvg.info, 1500))
#####
hvg.union <- union(hvg.sample_1, hvg.sample_2)

#set labels
sample_1_seurat@meta.data[,"protocol"] <- sample_1
sample_2_seurat@meta.data[,"protocol"] <- sample_2

#Canonical Correlation Vectors
LN <- RunCCA(sample_2_seurat, sample_1_seurat, genes.use = hvg.union, num.cc = 40)

#visualize results of CCA
p1 <- DimPlot(LN, reduction.use = "cca", group.by = "protocol", pt.size = 0.5, do.return = T)
p2 <- VlnPlot(LN, features.plot = "CC5", group.by = "protocol", do.return = T)
plot_grid(p1, p2)

#Identify the CCs that give a decent resolution
DimHeatmap(LN, reduction.type = "cca", cells.use = 500, dim.use = 1:40, do.balanced = T)

#Search for cells whose expression profile cannot be well-explained by low-dimensional CCA
LN <- CalcVarExpRatio(LN, reduction.type = "ica", grouping.var = "protocol", dims.use = 1:dims_use)

#align acording to CCA subspaces
LN <- AlignSubspace(LN, reduction.type = "cca", grouping.var = "protocol", dims.align = 1:dims_use)

#Plot the aligned CCA ACCA
p1 <- VlnPlot(LN, features.plot = "ACC1", group.by = "protocol", do.return = T)
p2 <- VlnPlot(LN, features.plot = "ACC2", group.by = "protocol", do.return = T)
plot_grid(p1, p2)

#Run single integrated analysis on all cells
LN <- RunTSNE(LN, reduction.use = "cca.aligned", dims.use = 1:dims_use)
LN <- FindClusters(LN, reduction.type = "cca.aligned", dims.use = 1:dims_use, save.SNN = T, resolution = resolution_merge, force.recalc = TRUE)

sample_1_cells <- grep(pattern = "^sample_1_", x = rownames(LN@meta.data), value = TRUE)
sample_2_cells <- grep(pattern = "^sample_2_", x = rownames(LN@meta.data), value = TRUE)
p1 <- TSNEPlot(LN, group.by = "protocol", do.return = T, pt.size = 0.4, cells.use = sample_1_cells)
p2 <- TSNEPlot(LN, group.by = "protocol", do.return = T, pt.size = 0.4, cells.use = sample_2_cells)
p3 <- TSNEPlot(LN, do.return = T, pt.size = 0.2, do.label = FALSE)
p4 <- TSNEPlot(LN, group.by = "protocol", do.return = T, pt.size = 0.2, colors.use = c("green","blue"))
par(mfrow = c(2, 2))
plot_grid(p1, p2, p3, p4)
#450x370
plot_grid(p3)
#450x355
plot_grid(p4)

FeaturePlot(LN, features.plot = c("Pdpn","Pecam1","Lyve1","Ackr4","Madcam1"),
            min.cutoff = "q9",
            cols.use = c('lightgrey', 'brown'),
            pt.size = 0.2)

#####
#Plot TSNE Maps
#####
#Make standard TSNE maps
gene_list_map_on_TSNE <- c("Lyve1",
                           "Acta2","Ackr3","Aldh1a2","Apoe",
                           "Bst1",
                           "Ccl2","Ccl5","Ccl7","Ccl8","Ccl9","Ccl11","Ccl19","Ccl21a",
                           "Cxcl1","Cxcl2","Cxcl13","Cxcl9","Cxcl10",
                           "Cd34","Cd55","Cd248",
                           "Ch25h",
                           "Csf1",
                           "Col15a1","Col4a1","Col3a1",
                           "Des",
                           "Enpp2",
                           "Gdf10",
                           "Il7","Il6","Icam1","Inmt",
                           "Ly6a","Ly6c1",
                           "Ntrk1","Nr4a1",
                           "Madcam1","Malat1","Mfge8",
                           "Nos2",
                           "Ptgis","Pdpn","Pecam1",
                           "Sfrp4",
                           "Tnfsf13b","Tnfsf11", "Tinagl1",
                           "Timd4",
                           "Vcam1")

for(i in 1:length(gene_list_map_on_TSNE)){
  setwd("C:/Users/jpe12/Desktop/scRNAseq_NewInput/TSNE_map_LN")
  gene_i <- gene_list_map_on_TSNE[i]
  print(gene_i)
  filename_i <- paste(sample_1, "_", sample_2, "_", gene_i, ".eps", sep="")
  print(filename_i)
  FeaturePlot(LN, features.plot = gene_i,
              no.legend = TRUE,
              cols.use = c('gainsboro', 'darkred'),
              pt.size = 0.2)
  ggsave(filename_i, width = 5.8, height = 6, units = "cm")
}
for(i in 1:length(gene_list_map_on_TSNE)){
  setwd("C:/Users/jpe12/Desktop/scRNAseq_NewInput/TSNE_map_LN")
  gene_i <- gene_list_map_on_TSNE[i]
  print(gene_i)
  filename_i <- paste(sample_1, "_", sample_2, "_", gene_i,"_scale" , ".eps", sep="")
  print(filename_i)
  FeaturePlot(LN, features.plot = gene_i,
              no.legend = FALSE,
              cols.use = c('gainsboro', 'darkred'),
              pt.size = 0.2)
  ggsave(filename_i, width = 5.8, height = 6, units = "cm")
}


