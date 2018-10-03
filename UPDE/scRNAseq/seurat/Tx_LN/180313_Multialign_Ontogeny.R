#Seurat 
library("Seurat")
library("cowplot")
library("Matrix")
library("magrittr")
library("dplyr")
library("pryr")


#######
#Load and label
#######
###
#mLN d0
###
#sample_1.data <- Read10X(data.dir = "C:/Users/jpe12/PowerFolders/R/2018_Exp270_Ontogeny/Data/neonatal/filtered_gene_bc_matrices/mm10")
sample_1.data <- Read10X(data.dir = "/Users/Pezoldt/PowerFolders/R/2018_Exp270_Ontogeny/Data/neonatal/filtered_gene_bc_matrices/mm10")
sample_1.data@Dimnames[[2]] <- paste("sample_1_", c(1:length(sample_1.data@Dimnames[[2]])), sep = "")
###
#mLN d10
###
#sample_2.data <- Read10X(data.dir = "C:/Users/jpe12/PowerFolders/R/2018_Exp270_Ontogeny/Data/d10/filtered_gene_bc_matrices/mm10")
sample_2.data <- Read10X(data.dir = "/Users/Pezoldt/PowerFolders/R/2018_Exp270_Ontogeny/Data/d10/filtered_gene_bc_matrices/mm10")
sample_2.data@Dimnames[[2]] <- paste("sample_2_", c(1:length(sample_2.data@Dimnames[[2]])), sep = "")

###
#mLN d24
###
#sample_3.data <- Read10X(data.dir = "C:/Users/jpe12/PowerFolders/R/2018_Exp270_Ontogeny/Data/d24/filtered_gene_bc_matrices/mm10")
sample_3.data <- Read10X(data.dir = "/Users/Pezoldt/PowerFolders/R/2018_Exp270_Ontogeny/Data/d24/filtered_gene_bc_matrices/mm10")
sample_3.data@Dimnames[[2]] <- paste("sample_3_", c(1:length(sample_3.data@Dimnames[[2]])), sep = "")

###
#mLN d60
###
#Experiment 2
#sample_4.data <- Read10X(data.dir = "C:/Users/jpe12/PowerFolders/R/2017_252_scRNAseq/Data/mLN_SPF/B6/filtered_gene_bc_matrices/mm10")
sample_4.data <- Read10X(data.dir = "/Users/Pezoldt/PowerFolders/R/2017_252_scRNAseq/Data/mLN_SPF/B6/filtered_gene_bc_matrices/mm10")
sample_4.data@Dimnames[[2]] <- paste("sample_4_", c(1:length(sample_4.data@Dimnames[[2]])), sep = "")

# Convert to sparse matrices for efficiency
sample_1.data <- as(as.matrix(sample_1.data), "dgCMatrix")
sample_2.data <- as(as.matrix(sample_2.data), "dgCMatrix")
sample_3.data <- as(as.matrix(sample_3.data), "dgCMatrix")
sample_4.data <- as(as.matrix(sample_4.data), "dgCMatrix")

#######
#Global variables
#######
#@User: Define sample names
#M == merge
sample_1 <- "d0"
sample_2 <- "d10"
sample_3 <- "d24"
sample_4 <- "d60"
setwd("/Users/Pezoldt/PowerFolders/R/2018_Exp270_Ontogeny/Output/Multialign_SEURAT")

#Important: If changing global variable Setup of pLN and mLN needs to be manually adjusted as the cluster IDs get changed
min_gene_number <- 250
mito_cutoff <- 0.045
nGene <- 3500
min_cells <- 20
resolution_single <- 1.4
dims_use <- 24

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
sample_4_seurat <- CreateSeuratObject(raw.data = sample_4.data, min.cells = min_cells, min.genes = min_gene_number)
sample_4_seurat <- NormalizeData(sample_4_seurat)
sample_4_seurat <- FindVariableGenes(sample_4_seurat, do.plot = F, display.progress = F)
sample_4_seurat <- ScaleData(sample_4_seurat)
sample_4_seurat@meta.data$tech <- sample_4

######
#Determine the genes dominant in the PCAs
######
#Merge the seurat objects
mLN_Ontogeny.merged <- MergeSeurat(sample_1_seurat, sample_2_seurat, do.normalize = FALSE)
mLN_Ontogeny.merged <- MergeSeurat(mLN_Ontogeny.merged, sample_3_seurat, do.normalize = FALSE)
mLN_Ontogeny.merged <- MergeSeurat(mLN_Ontogeny.merged, sample_4_seurat, do.normalize = FALSE)
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
ob.list <- list(sample_1_seurat, sample_2_seurat, sample_3_seurat, sample_4_seurat)
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
p1 <- VlnPlot(mLN_Ontogeny.integrated, features.plot = "ACC1", group.by = "tech", do.return = T)
p2 <- VlnPlot(mLN_Ontogeny.integrated, features.plot = "ACC2", group.by = "tech", do.return = T)
plot_grid(p1, p2)


# t-SNE and Clustering
mLN_Ontogeny.integrated <- FindClusters(mLN_Ontogeny.integrated, reduction.type = "cca",
                                    dims.use = 1:16, save.SNN = T, resolution = 0.8)
mLN_Ontogeny.integrated <- RunTSNE(mLN_Ontogeny.integrated,
                               reduction.use = "cca",
                               dims.use = 1:16, save.SNN = T,
                               resolution = resolution_merge, force.recalc = TRUE)



#Plot the aligned CCA ACCA
p1 <- VlnPlot(mLN_Ontogeny.integrated, features.plot = "ACC1", group.by = "tech", do.return = T)
p2 <- VlnPlot(mLN_Ontogeny.integrated, features.plot = "ACC2", group.by = "tech", do.return = T)
plot_grid(p1, p2)


#Check allocation of cells
sample_1_cells <- grep(pattern = "^sample_1_", x = rownames(mLN_Ontogeny.integrated@meta.data), value = TRUE)
sample_2_cells <- grep(pattern = "^sample_2_", x = rownames(mLN_Ontogeny.integrated@meta.data), value = TRUE)
sample_3_cells <- grep(pattern = "^sample_3_", x = rownames(mLN_Ontogeny.integrated@meta.data), value = TRUE)
sample_4_cells <- grep(pattern = "^sample_4_", x = rownames(mLN_Ontogeny.integrated@meta.data), value = TRUE)
p1 <- TSNEPlot(mLN_Ontogeny.integrated, group.by = "tech", do.return = T, pt.size = 0.4, cells.use = sample_1_cells)
p2 <- TSNEPlot(mLN_Ontogeny.integrated, group.by = "tech", do.return = T, pt.size = 0.4, cells.use = sample_2_cells)
p3 <- TSNEPlot(mLN_Ontogeny.integrated, group.by = "tech", do.return = T, pt.size = 0.4, cells.use = sample_3_cells)
p4 <- TSNEPlot(mLN_Ontogeny.integrated, group.by = "tech", do.return = T, pt.size = 0.4, cells.use = sample_4_cells)
p5 <- TSNEPlot(mLN_Ontogeny.integrated, do.return = T, pt.size = 0.2, do.label = FALSE)
#test <- LN@meta.data[order(-LN@meta.data$),]
p6 <- TSNEPlot(mLN_Ontogeny.integrated, group.by = "tech", do.return = T, pt.size = 0.2,
               colors.use = c("orange","green","blue","purple"))
par(mfrow = c(3, 2))
plot_grid(p1, p2, p3, p4, p5, p6)


FeaturePlot(mLN_Ontogeny.integrated, features.plot = c("Pdpn","Pecam1","Lyve1","Ackr4", 
                                  "Tnfsf11","Tnfsf13b","Acta2","Ccl19",
                                  "Madcam1","Ltbr","Tnfsf13b","Il6", 
                                  "Il7","Ackr3","Fn1","Col4a1"),
            min.cutoff = "q9",
            cols.use = c('lightgrey', 'brown'),
            pt.size = 0.2)




#####
#Find differentially expressed genes
#####
#Find markers for every cluster compared to all remaining cells
mLN_Ontogeny.integrated.markers <- FindAllMarkers(object = mLN_Ontogeny.integrated, only.pos = TRUE, min.pct = 0.25, 
                                   thresh.use = 0.25)
mLN_Ontogeny.integrated.markers_20 <- mLN_Ontogeny.integrated.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
mLN_Ontogeny.integrated.markers_20 <- mLN_Ontogeny.integrated.markers_20$gene
#Order LN_minus meta.data according to cluster size
#LN_minus_diff <- LN_minus
#LN_minus_diff@meta.data$res.1.4 <- reorder(as.factor(LN_minus_diff@meta.data$res.1.4), new.order = c("10","6","0","3","1",
#                                                                                        "11","2","9","8","5","7","4","12"))
#LN_minus_diff@meta.data <- LN_minus_diff@meta.data[order(LN_minus_diff@meta.data$res.1.4),]
#Save DEGs
#write.table(LN_minus.markers, paste("C:/Users/jpe12/PowerFolders/R/2017_252_scRNAseq/20180305_Denominator/Output/mLN/FSC/",sample_1,"_",sample_2,"_DEGs.csv", sep=""), dec=".", sep=",")


#LN_minus@ident <- sort(LN_minus@ident)

#Plot differential expression of Marker genes
DoHeatmap(object = mLN_Ontogeny.integrated, use.scaled = TRUE, group.by = "ident",
          genes.use = mLN_Ontogeny.integrated.markers_20, title = paste(sample_1, "&", sample_2, sep=" "),
          remove.key = FALSE,
          col.low = "white",
          col.mid =  "oldlace",
          col.high = "brown")



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
LN_x <- RunPCA(object = LN, pc.genes = LN@var.genes, do.print = TRUE, pcs.print = 1:5, 
             genes.print = 5, pcs.compute = 50)

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


FeaturePlot(LN, features.plot = c("Pdpn","Pecam1","Lyve1","Ackr4", 
                                  "Tnfsf11","Tnfsf13b","Icam1","Ccl19",
                                  "Madcam1","Ltbr","Tnfsf13b","Il6", 
                                  "Il7","Ackr3","Fn1","Col4a1"),
            min.cutoff = "q9",
            cols.use = c('lightgrey', 'brown'),
            pt.size = 0.2)

FeaturePlot(LN, features.plot = c("Pdpn","Pecam1","Lyve1","Ackr4", 
                                  "Tnfsf13b","Icam1","Ccl19",
                                  "Madcam1","Ltbr","Tnfsf13b","Il6", 
                                  "Il7","Ackr3","Fn1","Col4a1"),
            min.cutoff = "q9",
            cols.use = c('lightgrey', 'brown'),
            pt.size = 0.2)

FeatureHeatmap(object = LN, features.plot = c("Pdpn","Pecam1","Ackr4", "Ly6c1","Madcam1","Ltbr","Tnfsf13b","Il6", 
                                              "Ackr3","Fn1","Col4a1"),
               group.by = "protocol", sep.scale = TRUE,
               pt.size = 0.5, cols.use = c('lightgrey', 'brown'))

#not to be detected:
#CD3e
#"Ncr1" Nkp46
#"Ptprc"CD45R
#CD127, Il7r

########################################
#Merged analysis with ONLY LECs and BECs
#########################################
#####
#sample_2
#####
#Filter for mitochondiral and duplets and LECs and BECs out
sample_2_seurat_LEC_BEC <- FilterCells(object = sample_2_seurat, subset.names = c("nGene", "percent.mito"), 
                                     low.thresholds = c(min_gene_number, -Inf), high.thresholds = c(nGene, mito_cutoff),
                                     cells.use = sample_2_ONLY_LECs_BECs)
sample_2_seurat_LEC_BEC <- NormalizeData(sample_2_seurat_LEC_BEC)
sample_2_seurat_LEC_BEC <- ScaleData(sample_2_seurat_LEC_BEC)
sample_2_seurat_LEC_BEC <- FindVariableGenes(sample_2_seurat_LEC_BEC, do.plot = F)

#####
#sample_1
#####
#Filter for mitochondiral and duplets and LECs and BECs out
sample_1_seurat_LEC_BEC <- FilterCells(object = sample_1_seurat, subset.names = c("nGene", "percent.mito"), 
                                     low.thresholds = c(min_gene_number, -Inf), high.thresholds = c(nGene, mito_cutoff),
                                     cells.use = sample_1_ONLY_LECs_BECs)
sample_1_seurat_LEC_BEC <- NormalizeData(sample_1_seurat_LEC_BEC)
sample_1_seurat_LEC_BEC <- ScaleData(sample_1_seurat_LEC_BEC)
sample_1_seurat_LEC_BEC <- FindVariableGenes(sample_1_seurat_LEC_BEC, do.plot = F)

#####
#Alignment
#####
#resolution set
resolution_merge = 1.4
dims_use = 26
###
#take top variable genes
hvg.sample_1_LEC_BEC <- rownames(head(sample_1_seurat_LEC_BEC@hvg.info, 1500))
hvg.sample_2_LEC_BEC <- rownames(head(sample_2_seurat_LEC_BEC@hvg.info, 1500))
####
hvg.union <- union(hvg.sample_1_LEC_BEC, hvg.sample_2_LEC_BEC)

#set labels
sample_1_seurat_LEC_BEC@meta.data[,"protocol"] <- paste(sample_1, "_LEC_BEC", sep="")
sample_2_seurat_LEC_BEC@meta.data[,"protocol"] <- paste(sample_2, "_LEC_BEC", sep="")

#Canonical Correlation Vectors
LN_LEC_BEC <- RunCCA(sample_2_seurat_LEC_BEC, sample_1_seurat_LEC_BEC, genes.use = hvg.union, num.cc = 40)

#visualize results of CCA
p1 <- DimPlot(LN_LEC_BEC, reduction.use = "cca", group.by = "protocol", pt.size = 0.5, do.return = T)
p2 <- VlnPlot(LN_LEC_BEC, features.plot = "CC1", group.by = "protocol", do.return = T)
plot_grid(p1, p2)

#Search for cells whose expression profile cannot be well-explained by low-dimensional CCA
LN_LEC_BEC <- CalcVarExpRatio(LN_LEC_BEC, reduction.type = "ica", grouping.var = "protocol", dims.use = 1:dims_use)

#Use until CC chosen
#Discard cells where variance explained by CCA is <2-fold (ratio < 0.5) compared to PCA
LN_LEC_BEC.all.save <- LN_LEC_BEC

#align acording to CCA subspaces
LN_LEC_BEC <- AlignSubspace(LN_LEC_BEC, reduction.type = "cca", grouping.var = "protocol", dims.align = 1:dims_use)

#Plot the aligned CCA ACCA
p1 <- VlnPlot(LN_LEC_BEC, features.plot = "ACC1", group.by = "protocol", do.return = T)
p2 <- VlnPlot(LN_LEC_BEC, features.plot = "ACC2", group.by = "protocol", do.return = T)
plot_grid(p1, p2)

#Run single integrated analysis on all cells
LN_LEC_BEC <- RunTSNE(LN_LEC_BEC, reduction.use = "cca.aligned", dims.use = 1:dims_use)
LN_LEC_BEC <- FindClusters(LN_LEC_BEC, reduction.type = "cca.aligned", dims.use = 1:dims_use, save.SNN = T, resolution = 0.8, force.recalc = TRUE)
#Change cluster ID
#current.cluster.ids <- c(0:15)

#new.cluster.ids <- c("TFSC1","TFSC3","FSCx1/2","FSCx?",
#                    "FSCx1","FSCx2","TFSC2","FSCy4",
#                   "FSCy3","FSCy2","MRC","FSCy1",
#                  "CD8FSC","Pery", "?1","?2")
#LN_LEC_BEC@ident <- plyr::mapvalues(x = LN_LEC_BEC@ident, from = current.cluster.ids, to = new.cluster.ids)


sample_1_LEC_BEC_cells <- grep(pattern = "^sample_1_", x = rownames(LN_LEC_BEC@meta.data), value = TRUE)
sample_2_LEC_BEC_cells <- grep(pattern = "^sample_2_", x = rownames(LN_LEC_BEC@meta.data), value = TRUE)
p1 <- TSNEPlot(LN_LEC_BEC, group.by = "protocol", do.return = T, pt.size = 0.4, cells.use = sample_2_LEC_BEC_cells)
p2 <- TSNEPlot(LN_LEC_BEC, group.by = "protocol", do.return = T, pt.size = 0.4, cells.use = sample_1_LEC_BEC_cells)
p3 <- TSNEPlot(LN_LEC_BEC, do.return = T, pt.size = 0.2, do.label = TRUE)
#test <- LN@meta.data[order(-LN@meta.data$),]
p4 <- TSNEPlot(LN_LEC_BEC, group.by = "protocol", do.return = T, pt.size = 0.2, colors.use = c("green","blue"))
par(mfrow = c(2, 2))
plot_grid(p1, p2, p3, p4)

FeaturePlot(LN_LEC_BEC, features.plot = c("Pdpn","Pecam1","Lyve1","Ackr4","Icam1",
                                  "Ccl19","Ccl21a","Madcam1","Il7",
                                  "Ackr3","Col4a1","Cd74","H2-Ab1","Tap2"),
            min.cutoff = "q9",
            cols.use = c('lightgrey', 'brown'),
            pt.size = 0.2)

#Find markers for every cluster compared to all remaining cells
LN_LEC_BEC.markers <- FindAllMarkers(object = LN_LEC_BEC, only.pos = TRUE, min.pct = 0.25, 
                                   thresh.use = 0.25)
LN_LEC_BEC.markers_40 <- LN_LEC_BEC.markers %>% group_by(cluster) %>% top_n(40, avg_logFC)

#Save DEGs
write.table(LN_LEC_BEC.markers_40, paste(getwd(),"/",sample_1,"_",sample_2,"_DEGs.csv", sep=""), dec=".", sep=",")

DoHeatmap(object = LN_LEC_BEC, use.scaled = TRUE, group.by = "ident",
          genes.use = LN_LEC_BEC.markers_40$gene, title = paste("LEB_BEC",sample_1, "&", sample_2, sep=" "),
          remove.key = FALSE,
          col.low = "white",
          col.mid =  "oldlace",
          col.high = "brown")

VlnPlot(object = LN_LEC_BEC, features.plot = c("Pdpn","Pecam1"))

rm(LN_LEC_BEC)
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
sample_1_seurat_minus <- FilterCells(object = sample_1_seurat_minus, subset.names = c("Ackr4","Pecam1"),
                                     low.thresholds = c(-Inf, -Inf), high.thresholds = c(1,1))
sample_1_seurat_minus <- NormalizeData(sample_1_seurat_minus)
sample_1_seurat_minus <- ScaleData(sample_1_seurat_minus)
sample_1_seurat_minus <- FindVariableGenes(sample_1_seurat_minus, do.plot = F)
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
sample_2_seurat_minus <- ScaleData(sample_2_seurat_minus)
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
LN_minus <- RunCCA(sample_2_seurat_minus, sample_1_seurat_minus, genes.use = hvg.union, num.cc = 30)

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
#current.cluster.ids <- c(0:13)

#new.cluster.ids <- c("TFSC1","FSCx2/1","TFSC3","FSCx3",
 #                    "FSCy2","FSCx1/2","TFSC2",
  #                   "FSCy4","FSCy3","FSCy1/2","MRC",
   #                  "Pery","CD8FSC", "FSC?z")
#LN_minus@ident <- plyr::mapvalues(x = LN_minus@ident, from = current.cluster.ids, to = new.cluster.ids)


sample_1_minus_cells <- grep(pattern = "^sample_1_", x = rownames(LN_minus@meta.data), value = TRUE)
sample_2_minus_cells <- grep(pattern = "^sample_2_", x = rownames(LN_minus@meta.data), value = TRUE)
p1 <- TSNEPlot(LN_minus, group.by = "protocol", do.return = T, pt.size = 0.4, cells.use = sample_2_minus_cells)
p2 <- TSNEPlot(LN_minus, group.by = "protocol", do.return = T, pt.size = 0.4, cells.use = sample_1_minus_cells)
p3 <- TSNEPlot(LN_minus, do.return = T, pt.size = 0.2, do.label = TRUE)
#test <- LN@meta.data[order(-LN@meta.data$),]
p4 <- TSNEPlot(LN_minus, group.by = "protocol", do.return = T, pt.size = 0.2, colors.use = c("green","blue"))
par(mfrow = c(2, 2))
plot_grid(p1, p2, p3, p4)

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
FeaturePlot(LN_minus, features.plot = c("Cxcl9","Cxcl10","Nos2",
                                  "Vcam1","Ccl21a","Col15a1","Col4a1",
                                  "Il6","Il7","Timd4","Cd55","Icam1"),
            min.cutoff = "q9",
            cols.use = c('lightgrey', 'brown'),
            pt.size = 0.2
)

FeaturePlot(LN_minus, features.plot = c("Ccl2","Ccl8",
                                        "Ccl11","Ccl7",
                                        "Ccl5","Ly6c1","Bst1",
                                        "Ccl20","Ccl9","Ccl11"),
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

#Order LN_minus meta.data according to cluster size
#LN_minus_diff <- LN_minus
#LN_minus_diff@meta.data$res.1.4 <- reorder(as.factor(LN_minus_diff@meta.data$res.1.4), new.order = c("10","6","0","3","1",
 #                                                                                        "11","2","9","8","5","7","4","12"))
#LN_minus_diff@meta.data <- LN_minus_diff@meta.data[order(LN_minus_diff@meta.data$res.1.4),]
#Save DEGs
write.table(LN_minus.markers, paste("C:/Users/jpe12/PowerFolders/R/2017_252_scRNAseq/20180305_Denominator/Output/mLN/FSC/",sample_1,"_",sample_2,"_DEGs.csv", sep=""), dec=".", sep=",")


#LN_minus@ident <- sort(LN_minus@ident)

#Plot differential expression of Marker genes
DoHeatmap(object = LN_minus, use.scaled = TRUE, group.by = "ident",
          genes.use = LN_minus.markers_20$gene, title = paste(sample_1, "&", sample_2, sep=" "),
          remove.key = FALSE,
          col.low = "white",
          col.mid =  "oldlace",
          col.high = "brown")
#group.by = paste("res.1.4")

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


#####
#Cluster identification via signature
#####
library(ggplot2)
library(gridExtra)
#Use list of DEGs per cluster extracted from previous analysis
DEGs_NC_submitted <- read.csv("C:/Users/jpe12/PowerFolders/R/2018_Exp268_scRNAseq_Tx/Input/DEG_Top20_NC_submission_2017.csv",
                              sep = ";")
#Use list of DEGs from merged analysis of mLN and pLN multiple samples
DEGs_core <- read.csv("C:/Users/jpe12/PowerFolders/R/2017_252_scRNAseq/20180305_Denominator/Output/mLN_pLN/FSC/mLN_M_pLN_M_DEGs.csv",
                              sep = ",")
DEGs_core_40 <- DEGs_core %>% group_by(cluster) %>% top_n(40, avg_logFC)

#Load DEG tables and make list of clusters
DEGs_list <- split(DEGs_core_40, DEGs_core_40$cluster)
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

title <- paste("DEGs_core_40",sample_1, "OVER", sample_2, sep = "_")
pdf(paste(title, ".pdf", sep = ""))
marrangeGrob(out, nrow = round(length(DEGs_list) / 3 - 1), ncol = 3, top = title)
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



genes_interest <- rownames(LN_minus_1_2_aExp_mLN_main)
#Genes were chosen on their putative or described impact on environment
genes_interest <- c("Aldh1a2", "Aldh1a3","Ptgis","Rspo3","Sfrp4")
#genes_interest = c("Aldh1a2", "Aldh1a3", "Tcf21","Bmper","Sfrp4","Ptgis","Fez1")
FeatureHeatmap(object = LN_minus, features.plot = genes_interest, 
               group.by = "protocol", sep.scale = FALSE, pt.size = 0.1,
               cols.use = c("gray87","darkred"),
               min.exp = 0.5, max.exp = 5)


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


