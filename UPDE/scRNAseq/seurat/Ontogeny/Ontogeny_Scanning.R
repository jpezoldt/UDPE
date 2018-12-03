#Author: Joern Pezoldt
#Date: 13.11.2018
#Function: 
# 1) Select non-endothelial SCs
# 2) Compare impact of PCS, Resolution and number of variable genes on Clusters

#Seurat 
#mLN GF
library("Seurat")
library("cowplot")
library("Matrix")
library("magrittr")
library("dplyr")


#######
#Load and label
#######
#######
#Global variables
#######
# Note: Input required
organ <- "mLN"
condition <- "D0"
path_input <- "/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/270_2018-01-03_scRNA-Seq_mLN_ontogeny/Data/C11_neo_mLN/outs/filtered_gene_bc_matrices/mm10"

path_output <- paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/seurat/Subsets_to_scenic/", organ, "_", condition, sep = "")
dir.create(path_output, recursive = TRUE)

#Important: If changing global variable Setup of pLN and mLN needs to be manually adjusted as the cluster IDs get changed
min_gene_number <- 500
min_UMI <- 2000
mito_cutoff <- 0.045
nGene <- 4700
min_cells <- 20
resolution_single <- 1.2

#######
#Load and label
#######
sample_1.data <- Read10X(data.dir = path_input)
sample_1.data@Dimnames[[2]] <- paste(paste("sample_",organ,"_",condition,"_", sep=""), c(1:length(sample_1.data@Dimnames[[2]])), sep = "")

#########
#Setup sample_1
#########
sample_1_seurat <- CreateSeuratObject(raw.data = sample_1.data, min.cells = min_cells, min.genes = min_gene_number)

#Eliminate mitochondrial high cells and dublets
#Identify cells with mitochondrial reads
mito.genes <- grep(pattern = "^mt-", x = rownames(sample_1_seurat@data), value = TRUE)
percent.mito <- colSums(as.matrix(sample_1_seurat@raw.data[mito.genes, ]))/colSums(as.matrix(sample_1_seurat@raw.data))

#Identify cells with low nGene/nUmi ratio
ratio.nUMITOnGene <- sample_1_seurat@meta.data$nUMI / sample_1_seurat@meta.data$nGene
names(ratio.nUMITOnGene) <- rownames(sample_1_seurat@meta.data)

# AddMetaData adds columns to object@data.info and good to store QC stats
sample_1_seurat <- AddMetaData(object = sample_1_seurat, metadata = percent.mito, col.name = "percent.mito")
sample_1_seurat <- AddMetaData(object = sample_1_seurat, metadata = ratio.nUMITOnGene, col.name = "ratio.nUMITOnGene")
#VlnPlot(object = sample_1_seurat, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
#par(mfrow = c(1, 2))
#GenePlot(object = sample_1_seurat, gene1 = "nUMI", gene2 = "percent.mito")
#GenePlot(object = sample_1_seurat, gene1 = "nUMI", gene2 = "nGene")

# We filter out cells according unique gene counts 
# and percentage mitochondrial reads
sample_1_seurat <- FilterCells(object = sample_1_seurat, subset.names = c("nGene", "ratio.nUMITOnGene", "percent.mito"), 
                               low.thresholds = c(min_gene_number, -Inf, -Inf), high.thresholds = c(nGene, 4, mito_cutoff))

#Normalize the expression and log-transform
sample_1_seurat <- NormalizeData(object = sample_1_seurat, normalization.method = "LogNormalize", 
                                 scale.factor = 10000)

#Find variable genes independent of expression level
sample_1_seurat <- FindVariableGenes(object = sample_1_seurat, mean.function = ExpMean, dispersion.function = LogVMR, 
                                     x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5,
                                     do.plot = FALSE)

#scale data
sample_1_seurat <- ScaleData(object = sample_1_seurat, vars.to.regress = c("nUMI", "percent.mito"))

#perfrom PCA and print genes that define PCA
sample_1_seurat <- RunPCA(object = sample_1_seurat, pc.genes = sample_1_seurat@var.genes, do.print = TRUE, pcs.print = 1:5, 
                          genes.print = 5)

#ProjectPCA scores each gene in the dataset 
sample_1_seurat <- ProjectPCA(object = sample_1_seurat, do.print = FALSE)

#Group the cells into clusters
sample_1_seurat <- FindClusters(object = sample_1_seurat, reduction.type = "pca", dims.use = 1:15, 
                                resolution = resolution_single, print.output = 0, save.SNN = TRUE,
                                force.recalc = TRUE)

#perfrom Tsne
sample_1_seurat <- RunTSNE(object = sample_1_seurat, dims.use = 1:13, do.fast = TRUE)

#####
#Eliminate LECs, BECs sample_1
#####
#Identify LECs and BECs
#VlnPlot(object = sample_1_seurat, features.plot = c("Pecam1"))
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


##########################
#Analyse ONLY FSC
##########################
#resolution set
resolution_single = c(0.9, 1.1, 1.3, 1.5)
dims_use = c(12,18,24,30)
variable_genes = c(750,1000,1500)

sample_1_seurat <- CreateSeuratObject(raw.data = sample_1.data, min.cells = min_cells, min.genes = min_gene_number)
#Eliminate mitochondrial high cells and dublets
#Identify cells with mitochondrial readss
mito.genes <- grep(pattern = "^mt-", x = rownames(sample_1_seurat@data), value = TRUE)
percent.mito <- colSums(as.matrix(sample_1_seurat@raw.data[mito.genes, ]))/colSums(as.matrix(sample_1_seurat@raw.data))

#Identify cells with low nGene/nUmi ratio
ratio.nUMITOnGene <- sample_1_seurat@meta.data$nUMI / sample_1_seurat@meta.data$nGene
names(ratio.nUMITOnGene) <- rownames(sample_1_seurat@meta.data)

# AddMetaData adds columns to object@data.info and good to store QC stats
sample_1_seurat <- AddMetaData(object = sample_1_seurat, metadata = percent.mito, col.name = "percent.mito")
sample_1_seurat <- AddMetaData(object = sample_1_seurat, metadata = ratio.nUMITOnGene, col.name = "ratio.nUMITOnGene")

#Filter for mitochondiral and duplets and LECs and BECs out
# We filter out cells according unique gene counts 
# and percentage mitochondrial reads
sample_1_seurat <- FilterCells(object = sample_1_seurat, subset.names = c("nGene", "ratio.nUMITOnGene", "percent.mito"), 
                               low.thresholds = c(min_gene_number, -Inf, -Inf), high.thresholds = c(nGene, 4, mito_cutoff),
                               cells.use = sample_1_NO_LECs_BECs)

sample_1_seurat_minus <- FilterCells(object = sample_1_seurat, subset.names = c("Ackr4","Pecam1"),
                                     low.thresholds = c(-Inf, -Inf), high.thresholds = c(1,1))
sample_1_seurat_minus <- NormalizeData(sample_1_seurat_minus)
sample_1_seurat_minus <- ScaleData(object = sample_1_seurat_minus, vars.to.regress = c("nUMI", "percent.mito"))
sample_1_seurat_minus <- FindVariableGenes(object = sample_1_seurat_minus, mean.function = ExpMean, dispersion.function = LogVMR, 
                                           x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0., do.plot = FALSE)

for(a in 1:length(variable_genes)){
  variable_genes_i = variable_genes[a]
  print(paste("VarGenes=", variable_genes_i, sep=""))
  #take top variable genes
  hvg.sample_1_minus <- rownames(head(sample_1_seurat_minus@hvg.info, variable_genes_i))
  ####
  for(b in 1:length(dims_use)){
    dims_use_i = dims_use[b]
    print(paste("Dims=", dims_use_i, sep=""))
    #set labels for datasets in meta.data
    sample_1_seurat_minus@meta.data[,"protocol"] <- paste(organ,"_", condition, "_SC", sep="")
    
    #perfrom PCA and print genes that define PCA
    sample_1_seurat_minus <- RunPCA(object = sample_1_seurat_minus, pc.genes = hvg.sample_1_minus, do.print = TRUE, pcs.print = 1:5, 
                              genes.print = 5, pcs.compute = 30)
    
    #ProjectPCA scores each gene in the dataset 
    sample_1_seurat_minus <- ProjectPCA(object = sample_1_seurat_minus, do.print = FALSE)
    sample_1_seurat_minus <- RunTSNE(object = sample_1_seurat_minus, dims.use = 1:dims_use, do.fast = TRUE)
    #Group the cells into clusters
    
    for(i in 1:length(resolution_single)){
      resolution_merge_i = resolution_single[i]
      file_name_i = paste(organ, condition, "nVarGenes",variable_genes_i, "nDims", dims_use_i,  sep = "_", "Res", resolution_merge_i)
      print(file_name_i)
      sample_1_seurat_minus <- FindClusters(object = sample_1_seurat_minus, reduction.type = "pca", dims.use = 1:dims_use_i, 
                                      resolution = resolution_merge_i, print.output = 0, save.SNN = TRUE,
                                      force.recalc = TRUE)

      #Generate TSNE map
      p3 <- TSNEPlot(sample_1_seurat_minus, do.return = T, pt.size = 1.0, do.label = TRUE)
      
      #Calculate DEGs
      sample_1_seurat_minus.markers <- FindAllMarkers(object = sample_1_seurat_minus, only.pos = TRUE, min.pct = 0.25, 
                                                      thresh.use = 0.25)
      sample_1_seurat_minus.markers_20 <- sample_1_seurat_minus.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
      print(paste(file_name_i, ".pdf", sep=""))
      #######
      #save in PDF
      #######
      setwd(path_output)
      pdf(paste(file_name_i, ".pdf", sep=""), height = 10, width = 12)
      #each plot one slide
      
      #TSNE Plots and clustering
      par(mfrow = c(2, 2))
      p <- plot_grid(p3, labels = c("A"))
      title <- ggdraw() + draw_label(file_name_i, fontface='bold')
      print(plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)))
      
      #Dotplot of Key Genes on TSNE
      par(mfrow = c(1, 1))
      #Continue here: Does not run smoothly subsequent to addint the FeaturePlot()
      FP1 <- FeaturePlot(sample_1_seurat_minus, features.plot = c("Mfge8","Madcam1","Icam1",
                                                     "Vcam1","Ccl21a","mt-Co2",
                                                     "Tnfsf13b","Tnfsf11","Il6",
                                                     "Ackr3","Bst1", "Il7",
                                                     "Aldh1a2","Atf3","Maff",
                                                     "Cdk1"),
                         min.cutoff = "q9", no.legend = TRUE,
                         cols.use = c('lightgrey', 'brown'), pt.size = 0.2)
      print(plot_grid(FP1))
      
      #Heatmap Top DEGs per Cluster both samples
      par(mfrow = c(2, 1))
      FP2 <- DoHeatmap(object = sample_1_seurat_minus, use.scaled = TRUE, group.by = paste("res.",resolution_merge_i, sep=""),
                       genes.use = sample_1_seurat_minus.markers_20$gene, title = file_name_i,
                       remove.key = FALSE,
                       col.low = "white",col.mid =  "oldlace",col.high = "brown")
      print(plot_grid(FP2))
      
      dev.off()
      #write table for Top20 DEG
      write.table(sample_1_seurat_minus.markers_20, paste(file_name_i,"_ALL_Top20_DEGs.csv", sep=""), dec=".", sep=",")
      print(paste("+++++Done with", paste(file_name_i, ".pdf", sep="", "+++++")))
    }
  }
}
