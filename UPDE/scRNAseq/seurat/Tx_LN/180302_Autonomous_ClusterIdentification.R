#Seurat 
#mLN GF0
library("Seurat")
library("cowplot")
library("Matrix")
library("magrittr")
library("dplyr")


#######
#Load and label
#######
#Load 10x Data
sample_1.data <- Read10X(data.dir = "C:/Users/jpe12/PowerFolders/R/2018_Exp268_scRNAseq_Tx/Data/Tx_pLN_SPF/filtered_gene_bc_matrices/mm10")
#sample_1.data <- Read10X(data.dir = "/Users/Pezoldt/PowerFolders/R/2018_Exp268_scRNAseq_Tx/Data/Tx_mLN_SPF/filtered_gene_bc_matrices/mm10")
sample_1.data@Dimnames[[2]] <- paste("sample_1_", c(1:length(sample_1.data@Dimnames[[2]])), sep = "")

sample_2.data <- Read10X(data.dir = "C:/Users/jpe12/PowerFolders/R/2018_Exp268_scRNAseq_Tx/Data/Tx_mLN_SPF/filtered_gene_bc_matrices/mm10")
#sample_2.data <- Read10X(data.dir = "/Users/Pezoldt/PowerFolders/R/2017_252_scRNAseq/Data/mLN_SPF/B6/filtered_gene_bc_matrices/mm10")
sample_2.data@Dimnames[[2]] <- paste("sample_2_", c(1:length(sample_2.data@Dimnames[[2]])), sep = "")

#######
#Global variables
#######
#@User: Define sample names
sample_1 <- "Tx_pLN_SPF"
sample_2 <- "Tx_mLN_SPF"
setwd("C:/Users/jpe12/PowerFolders/R/2018_Exp268_scRNAseq_Tx/Output/pLN_SPF_Tx_vs_mLN_SPF_Tx")

min_gene_number <- 250
mito_cutoff <- 0.045
nGene <- 3500
min_cells <- 20
resolution_single <- 1.4


#########
#Setup sample_1
#########
sample_1_seurat <- CreateSeuratObject(raw.data = sample_1.data, min.cells = min_cells, min.genes = min_gene_number)

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
#@Ehsan: Here I take default parameters and am currently not adjusting the output. Would you also optimize here?
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

FeaturePlot(sample_1_seurat, features.plot = c("mt-Nd1","mt-Nd2","mt-Co1","mt-Co2",
                                   "mt-Atp8","mt-Atp6","mt-Co3","mt-Nd3",
                                   "mt-Nd4l","mt-Nd4","mt-Nd5"),
            min.cutoff = "q9",
            cols.use = c('lightgrey', 'blue'),
            pt.size = 0.1)

VlnPlot(object = sample_1_seurat, features.plot = c("mt-Nd1","mt-Nd2","mt-Co1","mt-Co2"))

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

#########
#Setup sample_2
#########
sample_2_seurat <- CreateSeuratObject(raw.data = sample_2.data, min.cells = min_cells, min.genes = min_gene_number)

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

FeaturePlot(sample_2_seurat, features.plot = c("mt-Nd1","mt-Nd2","mt-Co1","mt-Co2"),
            min.cutoff = "q9",
            cols.use = c('lightgrey', 'blue'),
            pt.size = 0.1)

VlnPlot(object = sample_2_seurat, features.plot = c("mt-Nd1","mt-Nd2","mt-Co1","mt-Co2"))

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

#LN_x <- RunPCA(object = LN, pc.genes = LN@var.genes, do.print = TRUE, pcs.print = 1:5, 
#             genes.print = 5, pcs.compute = 50)

#LN_x <- JackStraw(object = LN, num.replicate = 100, do.print = FALSE, num.pc = 50)
#!!! Obviously more PCs than are calculated are significant
#JackStrawPlot(object = LN_x, PCs = 1:40)

#visualize results of CCA
p1 <- DimPlot(LN, reduction.use = "cca", group.by = "protocol", pt.size = 0.5, do.return = T)
p2 <- VlnPlot(LN, features.plot = "CC5", group.by = "protocol", do.return = T)
plot_grid(p1, p2)

#Identify the CCs that give a decent resolution
#DimHeatmap(LN, reduction.type = "cca", cells.use = 500, dim.use = 1:40, do.balanced = T)

#Search for cells whose expression profile cannot be well-explained by low-dimensional CCA
LN <- CalcVarExpRatio(LN, reduction.type = "ica", grouping.var = "protocol", dims.use = 1:dims_use)

#Use until CC chosen
#Discard cells where variance explained by CCA is <2-fold (ratio < 0.5) compared to PCA
#LN.all.save <- LN
#Discards cells that are dataset specific
#This might not be useful
#Check how many cells get discarded
#LN <- SubsetData(LN, subset.name = "var.ratio.ica", accept.low = 0.5)
#LN.discard <- SubsetData(LN.all.save, subset.name = "var.ratio.ica", accept.high = 0.5)
#median(LN.discard@meta.data[,"nGene"])
#median(LN@meta.data[, "nGene"])
#VlnPlot(LN.discard, features.plot = "Pdpn", group.by = "protocol")
#rm(LN.all.save)

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
#test <- LN@meta.data[order(-LN@meta.data$),]
p4 <- TSNEPlot(LN, group.by = "protocol", do.return = T, pt.size = 0.2, colors.use = c("green","blue"))
par(mfrow = c(2, 2))
plot_grid(p1, p2, p3, p4)
#450x370
plot_grid(p3)
#450x355
plot_grid(p4)

FeaturePlot(LN, features.plot = c("Pdpn","Pecam1","Lyve1","Ackr4", 
                                  "Tnfsf11","Tnfsf13b","Icam1",
                                  "Vcam1","Ccl21a","Ccl19","Mfge8",
                                  "Des","Madcam1","Ltbr","Tnfsf13b","Il6", "Il7",
                                  "Ackr3","Fn1", "Col15a1","Kit","Col4a1","Col14a1"),
            min.cutoff = "q9",
            cols.use = c('lightgrey', 'brown'),
            pt.size = 0.2)

FeaturePlot(LN, features.plot = c("Kitl","Col14a1","Aldh1a2","Aldh1a3","Ptgis"),
            min.cutoff = "q9",
            cols.use = c('lightgrey', 'brown'),
            pt.size = 0.2)

#not to be detected:
#CD3e
#"Ncr1" Nkp46
#"Ptprc"CD45R
#CD127, Il7r

########################################
#Merged analysis without LECs and BECs
#########################################
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
sample_2_seurat_minus <- ScaleData(sample_2_seurat_minus)
sample_2_seurat_minus <- FindVariableGenes(sample_2_seurat_minus, do.plot = F)

#####
#sample_1
#####
#Filter for mitochondiral and duplets and LECs and BECs out
sample_1_seurat_minus <- FilterCells(object = sample_1_seurat, subset.names = c("nGene", "percent.mito"), 
                   low.thresholds = c(min_gene_number, -Inf), high.thresholds = c(nGene, mito_cutoff),
                   cells.use = sample_1_NO_LECs_BECs)
sample_1_seurat_minus <- FilterCells(object = sample_1_seurat_minus, subset.names = c("Ackr4","Pecam1"),
                   low.thresholds = c(-Inf, -Inf), high.thresholds = c(1,1))
sample_1_seurat_minus <- NormalizeData(sample_1_seurat_minus)
sample_1_seurat_minus <- ScaleData(sample_1_seurat_minus)
sample_1_seurat_minus <- FindVariableGenes(sample_1_seurat_minus, do.plot = F)

#####
#Alignment
#####
#resolution set
resolution_merge = 1.4
dims_use = 26
###
#take top variable genes
hvg.sample_1_minus <- rownames(head(sample_1_seurat_minus@hvg.info, 1500))
hvg.sample_2_minus <- rownames(head(sample_2_seurat_minus@hvg.info, 1500))
####
hvg.union <- union(hvg.sample_1_minus, hvg.sample_2_minus)

#set labels
sample_1_seurat_minus@meta.data[,"protocol"] <- paste(sample_1, "_minus", sep="")
sample_2_seurat_minus@meta.data[,"protocol"] <- paste(sample_2, "_minus", sep="")

#Canonical Correlation Vectors
LN_minus <- RunCCA(sample_2_seurat_minus, sample_1_seurat_minus, genes.use = hvg.union, num.cc = 40)

#visualize results of CCA
p1 <- DimPlot(LN_minus, reduction.use = "cca", group.by = "protocol", pt.size = 0.5, do.return = T)
p2 <- VlnPlot(LN_minus, features.plot = "CC1", group.by = "protocol", do.return = T)
plot_grid(p1, p2)

#Search for cells whose expression profile cannot be well-explained by low-dimensional CCA
LN_minus <- CalcVarExpRatio(LN_minus, reduction.type = "ica", grouping.var = "protocol", dims.use = 1:dims_use)

#Use until CC chosen
#Discard cells where variance explained by CCA is <2-fold (ratio < 0.5) compared to PCA
LN_minus.all.save <- LN_minus

#align acording to CCA subspaces
LN_minus <- AlignSubspace(LN_minus, reduction.type = "cca", grouping.var = "protocol", dims.align = 1:dims_use)

#Plot the aligned CCA ACCA
p1 <- VlnPlot(LN_minus, features.plot = "ACC1", group.by = "protocol", do.return = T)
p2 <- VlnPlot(LN_minus, features.plot = "ACC2", group.by = "protocol", do.return = T)
plot_grid(p1, p2)

#Run single integrated analysis on all cells
LN_minus <- RunTSNE(LN_minus, reduction.use = "cca.aligned", dims.use = 1:dims_use)
LN_minus <- FindClusters(LN_minus, reduction.type = "cca.aligned", dims.use = 1:dims_use, save.SNN = T, resolution = resolution_merge, force.recalc = TRUE)
#Change cluster ID
#current.cluster.ids <- c(10, 6, 0, 3, 1, 11, 2, 9, 8, 5, 7, 4, 12)

#new.cluster.ids <- c("Pery","MRC","TFSC1","TFSC2","TFSC3","Th1FSC",
       #              "FSCx1","FSCx2","FSCy1","FSCy2","FSCy3",
        #             "FSCy4","FSCx3")
#LN_minus@ident <- plyr::mapvalues(x = LN_minus@ident, from = current.cluster.ids, to = new.cluster.ids)


sample_1_minus_cells <- grep(pattern = "^sample_1_", x = rownames(LN_minus@meta.data), value = TRUE)
sample_2_minus_cells <- grep(pattern = "^sample_2_", x = rownames(LN_minus@meta.data), value = TRUE)
p1 <- TSNEPlot(LN_minus, group.by = "protocol", do.return = T, pt.size = 0.4, cells.use = sample_2_minus_cells)
p2 <- TSNEPlot(LN_minus, group.by = "protocol", do.return = T, pt.size = 0.4, cells.use = sample_1_minus_cells)
p3 <- TSNEPlot(LN_minus, do.return = T, pt.size = 0.2, do.label = TRUE)
#test <- LN@meta.data[order(-LN@meta.data$),]
p4 <- TSNEPlot(LN_minus, group.by = "protocol", do.return = T, pt.size = 0.2, colors.use = c("green","blue"))
list_tsne_plots <- list(p1, p3, p2, p4)

title <- paste("Overlay",sample_1, sample_2, sep = "_")
pdf(paste(title, ".pdf", sep = ""), height = 5, width = 9)
marrangeGrob(list_tsne_plots, nrow = 2, ncol = 2, top = title)
dev.off()


#####
#save(LN,file="C:/Users/jpe12/Desktop/LN_stromalCells.Robj")

FeaturePlot(LN_minus, features.plot = c("Mfge8","Madcam1","Icam1",
                                  "Vcam1","Ccl21a","Ccl19",
                                  "Tnfsf13b","Tnfsf11","Il6",
                                  "Ackr3","Bst1", "Il7"),
            min.cutoff = "q9",
            no.legend = TRUE,
            cols.use = c('lightgrey', 'brown'),
            pt.size = 0.2)


#Annotate clusters according to canonical markers
FeaturePlot(LN_minus, features.plot = c("Pdpn","Ptgis","Sfrp4",
                                  "Ackr3", "Csf1", "Ccl21a",
                                  "Il6", "Il7",
                                  "Notch2","mt-Co3","mt-Nd2"),
            min.cutoff = "q9",
            cols.use = c('lightgrey', 'blue'),
            pt.size = 0.1
)
FeaturePlot(LN_minus, features.plot = c("Tnfsf11","Tnfsf13b","Icam1",
                                  "Vcam1","Ccl21a","Ccl19","Mfge8",
                                  "Des","Madcam1","Ltbr","Tnfsf13b","Il6", "Il7",
                                  "Ptgs2","Ptgis"),
            min.cutoff = "q9",
            cols.use = c('lightgrey', 'brown'),
            pt.size = 0.2
)


FeaturePlot(LN_minus, features.plot = c("Tnfsf13b","Timd4","Ly6c1",
                                  "Cd34","Icam1","Il7",
                                  "Ackr3","Col15a1"),
            min.cutoff = "q9",
            cols.use = c('lightgrey', 'brown'),
            pt.size = 0.2
)

#Talk Jason Syster
FeaturePlot(LN_minus, features.plot = c("Cxcl9","Cxcl1","Ccl21a","Mfge8",
                                  "Tnfsf13b","Tnfsf11","Cd34",
                                  "Des","Acta2","Ch25h","Tmem119","Bst1",
                                  "Cd248","Enpp2","Ccl9"),
            min.cutoff = "q9",
            cols.use = c('lightgrey', 'brown'),
            pt.size = 0.2
)

#Imprinted genes
#FeatureHeatmap(object = LN_minus, features.plot = c("Aldh1a2", "Aldh1a3", "Ptgis",
 #                                                   "Rspo3","Sfrp4"), 
  #             group.by = "protocol", sep.scale = TRUE,
   #            pt.size = 0.5, cols.use = c('lightgrey', 'brown'))


#####
#Find differentially expressed genes
#####
#Find markers for every cluster compared to all remaining cells
LN_minus.markers <- FindAllMarkers(object = LN_minus, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)
LN_minus.markers_20 <- LN_minus.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)

LN_minus@ident <- sort(LN_minus@ident)

#Plot differential expression of Marker genes
p <- DoHeatmap(object = LN_minus, use.scaled = TRUE, group.by = "ident",
          genes.use = LN_minus.markers_20$gene, title = paste(sample_1, "&", sample_2, sep=" "),
          remove.key = FALSE,
          col.low = "white",
          col.mid =  "oldlace",
          col.high = "brown")

title <- paste("DEG_Marker",sample_1, sample_2, sep = "_")
pdf(paste(title, ".pdf", sep = ""), height = 5, width = 9)
plot(p)
dev.off()


#####
#Cluster identification via signature
#####
library(ggplot2)
library(gridExtra)
#Use list of DEGs per cluster extracted from previous analysis
DEGs_NC_submitted <- read.csv("C:/Users/jpe12/PowerFolders/R/2018_Exp268_scRNAseq_Tx/Input/DEG_Top20_NC_submission_2017.csv",
                              sep = ";")
#to identify known clusters in current analyis

#Load DEG tables and make list of clusters
DEGs_list <- split(DEGs_NC_submitted, DEGs_NC_submitted$cluster)
#Average gene expression 
LN_minus_aExp <- AverageExpression(LN_minus)

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
  
  LN_minus_aExp_cluster_i <- subset(LN_minus_aExp, rownames(LN_minus_aExp) %in% cluster_DEG_i)
  
  Scale_LN_minus_aExp_cluster_i <- apply(LN_minus_aExp_cluster_i, 1, scale)
  Zscore_cluster_i <- rowSums(Scale_LN_minus_aExp_cluster_i)
  Zscore_cluster_i_norm <- Zscore_cluster_i / (nrow(LN_minus_aExp_cluster_i)/10)
  
  #mLN maintained Scores
  Zscore_table_cluster_i <- data.frame(module = c(rep(cluster_ID_i,length(Zscore_cluster_i_norm))),
                                      cluster = colnames(LN_minus_aExp),
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

title <- paste("NC_Top20_DEGs",sample_1, "OVER", sample_2, sep = "_")
pdf(paste(title, ".pdf", sep = ""))
marrangeGrob(out, nrow = round(length(DEGs_list) / 3 - 1), ncol = 3, top = title)
dev.off()

######################################
#Calculate Z-score for imprinted genes
#######################################
#load genes lists to assess enrichment in clusters
#Imprinted
#pLN_SPF maintained
mLN_main <- read.delim("C:/Users/jpe12/PowerFolders/R/2017_244_scRNASeq/Data/Input_Tx/mLN_maintained.txt")
mLN_main <- as.character(mLN_main$GeneSymbol)
#mLN repressed
mLN_rep <- read.delim("C:/Users/jpe12/PowerFolders/R/2017_244_scRNASeq/Data/Input_Tx/mLN_repressed.txt")
mLN_rep <- as.character(mLN_rep$GeneSymbol)

#####
#Expression per sample_1 and sample_2 Cluster
#####
#Setup SEURAT Object where sample_1 and sample_2 Clusters are separate
LN_minus@meta.data$res.1.4 <- as.character(LN_minus@ident)
sample_1_minus_cells <- rownames(subset(LN_minus@meta.data, protocol == paste(sample_1, "_minus", sep="")))
sample_2_minus_cells <- rownames(subset(LN_minus@meta.data, protocol == paste(sample_2, "_minus", sep="")))

LN_minus_sample_1 <- SubsetData(object = LN_minus,
                   cells.use = sample_1_minus_cells)
LN_minus_sample_2 <- SubsetData(object = LN_minus,
                   cells.use = sample_2_minus_cells)

#Merge mLN and pLN_GF Seurat objects
LN_minus_1_2 <- MergeSeurat(LN_minus_sample_1, LN_minus_sample_2, do.normalize = FALSE)
#change identity
LN_minus_1_2 <- SetIdent(LN_minus_1_2, ident.use = c(paste(sample_1, as.numeric(LN_minus_sample_1@ident)-1, sep="_"),
                                                     paste(sample_2, as.numeric(LN_minus_sample_2@ident)-1, sep="_")))

#Normalize the expression and log-transform
LN_minus_1_2 <- NormalizeData(object = LN_minus_1_2, normalization.method = "LogNormalize", 
                     scale.factor = 10000)

#Scale Data
LN_minus_1_2 <- ScaleData(LN_minus_1_2)

#Find variable genes independent of expression level
LN_minus_1_2 <- FindVariableGenes(object = LN_minus_1_2, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
#number of variable genes
length(x = LN_minus_1_2@var.genes)

#scale data
LN_minus_1_2 <- ScaleData(object = LN_minus_1_2, vars.to.regress = c("nUMI", "percent.mito"))

#perfrom PCA and print genes that define PCA
LN_minus_1_2 <- RunPCA(object = LN_minus_1_2, pc.genes = LN_minus_1_2@var.genes, do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5)

# Examine and visualize PCA results a few different ways
PrintPCA(object = LN_minus_1_2, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = LN_minus_1_2, pcs.use = 1:6)
#PCAPlot(object = mLN, dim.1 = 1, dim.2 = 2)

# ProjectPCA scores each gene in the dataset 
LN_minus_1_2 <- ProjectPCA(object = LN_minus_1_2, do.print = FALSE)

#Compare Cluster by Cluster
#Store DEGs in List
n_cluster = length(levels(as.factor(LN_minus_1_2@meta.data$res.1.4)))
DEG_l <- list()
i = 1
for(i in i:n_cluster){
  #cluster ID
  k = i - 1
  cluster_ID_1 <- paste(sample_1,"_",k,sep="")
  cluster_ID_2 <- paste(sample_2,"_",k,sep="")
  
  #Identify DEGs avg_diff>0 (high in mLN) avg_Diff<0 (high in pLN_GF)
  DEGs <- FindMarkers(LN_minus_1_2, cluster_ID_1, cluster_ID_2, thresh.use = 0.25, 
                                     test.use = "roc")
  print(cluster_ID_1)
  print(nrow(DEGs))
  DEG_l[[i]] <- DEGs
}
names(DEG_l) <- paste("cluster_", 0:(length(levels(as.factor(LN_minus_1_2@meta.data$res.1.4)))-1), sep="")

###
#Check overlap with imprinted genes
###
#Imprinted
#mLN maintained
#load genes lists to assess enrichment in clusters
#Imprinted
#pLN_SPF maintained
mLN_main <- read.delim("C:/Users/jpe12/PowerFolders/R/2017_244_scRNASeq/Data/Input_Tx/mLN_maintained.txt")
mLN_main <- as.character(mLN_main$GeneSymbol)
#mLN repressed
mLN_rep <- read.delim("C:/Users/jpe12/PowerFolders/R/2017_244_scRNASeq/Data/Input_Tx/mLN_repressed.txt")
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

######
#Zscore over expression of all genes from genemodules
######
library(pheatmap)
#Average Expression per cluster
LN_minus_1_2_aExp <- AverageExpression(LN_minus_1_2)

#Maintained genes mLN
LN_minus_1_2_aExp_mLN_main <- subset(LN_minus_1_2_aExp, rownames(LN_minus_1_2_aExp) %in% mLN_main)
pheatmap(LN_minus_1_2_aExp_mLN_main, scale = "row", cluster_cols = TRUE, border_color = "black", cellwidth = 10,
         treeheight_row = 0, treeheight_col = 20,
         cellheigth = 10, color = colorRampPalette(c("white", "oldlace","brown"), space="rgb")(128))

Scale_LN_minus_1_2_aExp_mLN_main <- apply(LN_minus_1_2_aExp_mLN_main, 1, scale)
Zscore_mLN_main <- rowSums(Scale_LN_minus_1_2_aExp_mLN_main)
Zscore_mLN_main_norm <- Zscore_mLN_main / (nrow(LN_minus_1_2_aExp_mLN_main)/10)
Zscore_mLN_main_norm_1 <- Zscore_mLN_main_norm[1:((length(Zscore_mLN_main_norm))/2)]
Zscore_mLN_main_norm_2 <- Zscore_mLN_main_norm[((length(Zscore_mLN_main_norm))/2+1):length(Zscore_mLN_main_norm)]

#mLN maintained Scores
Zscore_table_mLN_main <- data.frame(module = c(rep("mLN_main",length(Zscore_mLN_main_norm))),
                               cluster = colnames(LN_minus_1_2_aExp),
                               cZscore = Zscore_mLN_main_norm)

#Only mLN
ggplot(data = Zscore_table_mLN_main, aes(x = cluster, y = cZscore, fill = module)) +
  geom_bar(stat = "identity", colour = "black") +
  theme(axis.text.x=element_text(angle=270,hjust=1)) +
  scale_fill_manual(values = c("deepskyblue1", "deeppink"))



#Only sample_1
LN_minus_1_2_aExp_mLN_main_1_only <- subset(LN_minus_1_2_aExp, rownames(LN_minus_1_2_aExp) %in% mLN_main)[,c(1:(ncol(LN_minus_1_2_aExp)/2))]
LN_minus_1_2_aExp_mLN_main_1_only <- subset(LN_minus_1_2_aExp_mLN_main_1_only, !(rowSums(LN_minus_1_2_aExp_mLN_main_1_only) == 0))
pheatmap(LN_minus_1_2_aExp_mLN_main_1_only, scale = "row", cluster_cols = FALSE, border_color = "black", cellwidth = 10,
         treeheight_row = 0, treeheight_col = 20,
         cellheigth = 10, color = colorRampPalette(c("white", "oldlace","brown"), space="rgb")(128))

#Repressed genes mLN
LN_minus_1_2_aExp_mLN_rep <- subset(LN_minus_1_2_aExp, rownames(LN_minus_1_2_aExp) %in% mLN_rep)
LN_minus_1_2_aExp_mLN_rep <- subset(LN_minus_1_2_aExp_mLN_rep, !(rowSums(LN_minus_1_2_aExp_mLN_rep) == 0))
pheatmap(LN_minus_1_2_aExp_mLN_rep, scale = "row", cluster_cols = FALSE, border_color = "black", cellwidth = 10,
         treeheight_row = 0, treeheight_col = 20,
         cellheigth = 10, color = colorRampPalette(c("white", "oldlace","brown"), space="rgb")(128))
Scale_LN_minus_1_2_aExp_mLN_rep <- apply(LN_minus_1_2_aExp_mLN_rep, 1, scale)
Zscore_mLN_rep <- rowSums(Scale_LN_minus_1_2_aExp_mLN_rep)
Zscore_mLN_rep_norm <- Zscore_mLN_rep / (nrow(LN_minus_1_2_aExp_mLN_rep)/10)
Zscore_mLN_rep_norm_1 <- Zscore_mLN_rep_norm[1:((length(Zscore_mLN_rep_norm))/2)]
Zscore_mLN_rep_norm_2 <- Zscore_mLN_rep_norm[((length(Zscore_mLN_rep_norm))/2+1):length(Zscore_mLN_rep_norm)]


#mLN repressed Scores
Zscore_table_mLN_rep <- data.frame(module = c(rep("mLN_rep",length(Zscore_mLN_rep_norm))),
                                    cluster = colnames(LN_minus_1_2_aExp),
                                    cZscore = Zscore_mLN_rep_norm)
ggplot(data = Zscore_table_mLN_rep, aes(x = cluster, y = cZscore, fill = module)) +
  geom_bar(stat = "identity", colour = "black") +
  theme(axis.text.x=element_text(angle=270,hjust=1)) +
  scale_fill_manual(values = c("deeppink"))

#mLN repressed and maintained scores
Zscore_table_mLN_rep_main <- data.frame(module = c(rep("mLN_main",length(Zscore_mLN_main_norm)),
                                              rep("mLN_rep",length(Zscore_mLN_rep_norm))),
                                   cluster = rep(colnames(LN_minus_1_2_aExp),2),
                                   cZscore = c(Zscore_mLN_main_norm,Zscore_mLN_rep_norm))
p <- ggplot(data = Zscore_table_mLN_rep_main, aes(x = cluster, y = cZscore, fill = module)) +
  geom_bar(stat = "identity", colour = "black") +
  theme(axis.text.x=element_text(angle=270,hjust=1)) +
  scale_fill_manual(values = c("deepskyblue1","deeppink"))


title <- paste("Tx_mLNFSC_rep_main_ON",sample_1, "OVER", sample_2, sep = "_")
pdf(paste(title, ".pdf", sep = ""))
plot(p, title = title)
dev.off()


###################################
#Dotplots FeaturePlot
###################################
#Imprinted
#mLN maintained
mLN_main <- read.delim("C:/Users/jpe12/PowerFolders/R/2017_scRNASeq/Data/Input_Tx/mLN_maintained.txt")
mLN_main <- as.character(mLN_main$GeneSymbol)

#Expressed
LN_minus_p_m_aExp <- AverageExpression(LN_minus_p_m)
#Maintained genes mLN
LN_minus_p_m_aExp_mLN_main <- subset(LN_minus_p_m_aExp, rownames(LN_minus_p_m_aExp) %in% mLN_main)



genes_interest <- rownames(LN_minus_p_m_aExp_mLN_main)
#Genes were chosen on their putative or described impact on environment
genes_interest <- c("Aldh1a2", "Aldh1a3","Ptgis","Rspo3","Sfrp4")
#genes_interest = c("Aldh1a2", "Aldh1a3", "Tcf21","Bmper","Sfrp4","Ptgis","Fez1")
p <- FeatureHeatmap(object = LN_minus, features.plot = genes_interest, 
               group.by = "protocol", sep.scale = FALSE, pt.size = 0.1,
               cols.use = c("gray87","darkred"),
               min.exp = 0.5, max.exp = 5)

title <- paste("NC_Feature_Genes",sample_1, "OVER", sample_2, sep = "_")
pdf(paste(title, ".pdf", sep = ""), height = 6.5, width = 5)
FeatureHeatmap(object = LN_minus, features.plot = genes_interest, 
               group.by = "protocol", sep.scale = FALSE, pt.size = 0.1,
               cols.use = c("gray87","darkred"),
               min.exp = 0.5, max.exp = 5)
dev.off()


###########################################
#Histogram plots JoyPlot
##########################################
#####
#Surface markers
#####
#Load Uniprot file
UP_CellMem <- read.delim("C:/Users/jpe12/PowerFolders/R/2017_scRNASeq/Data/Input_Uniprot/uniprot_2017-10-09_CellMembrane.tab")
UP_CellMem <- as.character(UP_CellMem$Gene.names...primary..)
#Overlap with LN data
LN_minus.markers_CellMem <- subset(LN_minus.markers, gene %in% UP_CellMem)
LN_minus.markers_CellMem_FC1 <- subset(LN_minus.markers_CellMem, avg_diff >= 1)
LN_minus.markers_CellMem_FC1_gene <- as.character(LN_minus.markers_CellMem_FC1$gene)
JoyPlot(object = LN_minus_m, features.plot = LN_minus.markers_CellMem_FC1_gene, nCol = 5)
#Overlap with LN_m data
LN_minus.markers_m <- FindAllMarkers(object = LN_minus_m, only.pos = TRUE, min.pct = 0.25, 
                                     thresh.use = 0.25)
LN_minus.markers_m_CellMem <- subset(LN_minus.markers_m, gene %in% UP_CellMem)
LN_minus.markers_m_CellMem_FC1 <- subset(LN_minus.markers_m_CellMem, avg_diff >= 1)
LN_minus.markers_m_CellMem_FC1_gene <- as.character(LN_minus.markers_m_CellMem_FC1$gene)
JoyPlot(object = LN_minus_m, features.plot = LN_minus.markers_m_CellMem_FC1_gene, nCol = 5)
#Overlap with LN_m data
LN_minus.markers_p <- FindAllMarkers(object = LN_minus_p, only.pos = TRUE, min.pct = 0.25, 
                                     thresh.use = 0.25)
LN_minus.markers_p_CellMem <- subset(LN_minus.markers_p, gene %in% UP_CellMem)
LN_minus.markers_p_CellMem_FC1 <- subset(LN_minus.markers_p_CellMem, avg_diff >= 1.2)
LN_minus.markers_p_CellMem_FC1_gene <- as.character(LN_minus.markers_p_CellMem_FC1$gene)
JoyPlot(object = LN_minus_p, features.plot = LN_minus.markers_p_CellMem_FC1_gene, nCol = 5)


#Load Uniprot file
UP_Sec <- read.delim("C:/Users/jpe12/PowerFolders/R/2017_scRNASeq/Data/Input_Uniprot/uniprot_2017-10-09_Secreted.tab")
UP_Sec <- as.character(UP_Sec$Gene.names...primary..)
#Overlap with LN data
LN_minus.markers_Sec <- subset(LN_minus.markers, gene %in% UP_Sec)
LN_minus.markers_Sec_FC1 <- subset(LN_minus.markers_Sec, avg_diff >= 2)
LN_minus.markers_Sec_FC1_gene <- as.character(LN_minus.markers_Sec_FC1$gene)
JoyPlot(object = LN_minus_m, features.plot = LN_minus.markers_Sec_FC1_gene, nCol = 5)
#Overlap with LN_m data
#LN.markers_m <- FindAllMarkers(object = LN_m, only.pos = TRUE, min.pct = 0.25, 
#                              thresh.use = 0.25)
LN_minus.markers_m_Sec <- subset(LN_minus.markers_m, gene %in% UP_Sec)
LN_minus.markers_m_Sec_FC1 <- subset(LN_minus.markers_m_Sec, avg_diff >= 1.5)
LN_minus.markers_m_Sec_FC1_gene <- as.character(LN_minus.markers_m_Sec_FC1$gene)
JoyPlot(object = LN_minus_m, features.plot = LN_minus.markers_m_Sec_FC1_gene, nCol = 5)
#Overlap with LN_m data
#LN.markers_p <- FindAllMarkers(object = LN_p, only.pos = TRUE, min.pct = 0.25, 
#                            thresh.use = 0.25)
LN_minus.markers_p_Sec <- subset(LN_minus.markers_p, gene %in% UP_Sec)
LN_minus.markers_p_Sec_FC1 <- subset(LN_minus.markers_p_Sec, avg_diff >= 2)
LN_minus.markers_p_Sec_FC1_gene <- as.character(LN_minus.markers_p_Sec_FC1$gene)
JoyPlot(object = LN_minus_p, features.plot = LN_minus.markers_p_Sec_FC1_gene, nCol = 5)




#Setup LN object to plot clusters of interest
mLN_cells <- grep(pattern = "^m_", x = rownames(LN_minus@meta.data), value = TRUE)
pLN_GF_cells <- grep(pattern = "^p_", x = rownames(LN_minus@meta.data), value = TRUE)
LN_minus_p <- SubsetData(object = LN_minus,
                         cells.use = pLN_GF_minus_cells)
LN_minus_m <- SubsetData(object = LN_minus,
                         cells.use = mLN_minus_cells)

genes_interest = c("Flt3l","Csf1")
#"Ltb","Jag1")
JoyPlot(object = LN_minus_p, features.plot = genes_interest, nCol = 2, do.sort = FALSE, y.log = FALSE)
#VlnPlot(object = LN_m, features.plot = genes_interest, nCol = 2, y.max = 3)
DotPlot(object = LN_minus_p,genes.plot = genes_interest, plot.legend = TRUE, cols.use = c("grey","darkred"),
        col.min = 0, col.max = 5, dot.scale = 6, dot.min = 0.2)


#####
#PieCharts
#####
sample_1_minus_cells <- row.names(subset(LN_minus@meta.data, protocol == paste(sample_1, "_minus", sep="")))
sample_2_minus_cells <- row.names(subset(LN_minus@meta.data, protocol == paste(sample_2, "_minus", sep="")))

LN_minus_sample_1 <- SubsetData(object = LN_minus,
                                cells.use = sample_1_minus_cells)
LN_minus_sample_2 <- SubsetData(object = LN_minus,
                                cells.use = sample_2_minus_cells)
S1 <- table(LN_minus_sample_1@meta.data$res.1.4)
S2 <- table(LN_minus_sample_2@meta.data$res.1.4)

count_1_2 <- rbind(S1, S2)
rownames(count_1_2) <- c(sample_1, sample_2)
par(mfcol = c(1, ncol(count_1_2)))

for(i in 1:ncol(count_1_2)){
  pie(count_1_2[,i], main = colnames(count_1_2)[i], col = c("blue","green"))
}




#####
#Differential expression between sample_1 and sample_2 per Cluster
#####
#Make Seurat object for only one subset
subset_0 <- SubsetData(LN_minus, ident.use = "0")
#make new identities
#count cells per cluster per sample
n_sample_1 <- nrow(subset(subset_0@meta.data, protocol == paste(sample_1, "_minus", sep="")))
n_sample_2 <- nrow(subset(subset_0@meta.data, protocol == paste(sample_2, "_minus", sep="")))
#generate factor to rename
cluster_sorting <- as.factor(c(rep("sample_1_cluster_0", n_sample_1),rep("sample_2_cluster_0", n_sample_2)))
cell_names <- names(subset_0@ident)
new_ids <- setNames(cluster_sorting,cell_names)
subset_0@ident <- new_ids
#Calculate DEGs for one subset
diff_0 <- FindMarkers(subset_0, ident.1 = "sample_1_cluster_0", ident.2 = "sample_2_cluster_0",
                      test.use = "tobit", thresh.use = 0.5)

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
myfiles <- list(LN_minus.markers)
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
title <- c("Differential GO scRNASeq")

pheatmap(data_heatmap, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = TRUE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("white", "grey","deepskyblue2"), space="rgb")(128),
         main = title)


