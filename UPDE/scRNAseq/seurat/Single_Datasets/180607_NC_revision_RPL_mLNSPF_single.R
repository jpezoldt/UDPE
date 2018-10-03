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

#######
#Load and label
#######
###
#mLN
###
#Experiment 2
#4_... IDs mLN cells single
#sample_1_2.data <- Read10X(data.dir = "/Users/Pezoldt/PowerFolders/R/2017_252_scRNAseq/Data/mLN_SPF/B6/filtered_gene_bc_matrices/mm10")
sample_1_2.data <- Read10X(data.dir = "C:/Users/jpe12/PowerFolders/R/2017_252_scRNAseq/Data/mLN_SPF/B6/filtered_gene_bc_matrices/mm10")
sample_1_2.data@Dimnames[[2]] <- paste("sample_4_2_", c(1:length(sample_1_2.data@Dimnames[[2]])), sep = "")
#Experiment 1
#sample_1_1.data <- read.csv("/Users/Pezoldt/PowerFolders/R/2017_244_scRNASeq/Data/mLNSPF.csv")
sample_1_1.data <- read.csv("C:/Users/jpe12/PowerFolders/R/2017_244_scRNASeq/Data/mLNSPF.csv")
row.names(sample_1_1.data) <- sample_1_1.data$X
sample_1_1.data <- sample_1_1.data[,c(2:ncol(sample_1_1.data))]
colnames(sample_1_1.data) <- paste("sample_4_1_", c(1:ncol(sample_1_1.data)), sep = "")

#######
#Global variables
#######
sample_1 <- "mLN"
setwd("C:/Users/jpe12/PowerFolders/R/Paper/2018_Pezoldt_Pasztoi_Revision_NC/scRNA-seq/Output/mLN_single/02_scmap")
#setwd("/Users/Pezoldt/PowerFolders/R/Paper/2018_Pezoldt_Pasztoi_Revision_NC/scRNA-seq/Output/mLN_single/02_scmap")

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

#Change Cluster Names
current.cluster.ids <- c(0:12)

new.cluster.ids <- c("Ccl19highMadcam1+","Inmt+Cxcl12+","Inmt+","Cd34+Gdf10+",
                     "Cd34+Has1+","Cd34+Aldh1a2+","Il6+Cxcl1","Cd34+Ackr3+",
                     "Cd34+Cd248+","Ccl19+Il7+","pSC","PvC",
                     "SCx")

sample_1_seurat_minus@ident <- plyr::mapvalues(x = sample_1_seurat_minus@ident, from = current.cluster.ids, to = new.cluster.ids)

FeaturePlot(sample_1_seurat_minus, features.plot = c("Madcam1","Il7", "Aldh1a2",
                                                     "Ccl19","Il6","Vcam1",
                                                     "Icam1","Inmt","Nr4a1",
                                                     "Tnfsf11","Cxcl13"),
            min.cutoff = "q9",
            cols.use = c('gainsboro', 'darkred'),
            pt.size = 0.1)


TSNEPlot(object = sample_1_seurat_minus, pt.size = 0.3, do.label = TRUE)
filename <- paste("C:/Users/jpe12/Desktop/Seurat_objects/Final/","mLN_FSC_","PCA_",pca_use,"_nVarGene_",n_variable_genes,"_Res_",resolution_use,".rds",sep="")
saveRDS(sample_1_seurat_minus, file = filename)
#LN_minus_mLN <- readRDS("C:/Users/jpe12/Desktop/Seurat_objects/LN_minus_mLN.rds")
#sample_1_seurat_minus <- LN_minus_mLN

FeaturePlot(sample_1_seurat_minus, features.plot = c("Cxcl9","Tnfsf11","Madcam1","Nkain4","Vcam1",
                                        "Icam1","Bst1","Cxcl13","Ccl19","Krt18",
                                        "Cd248","Cd34","Il6","Inmt","Nr4a1","Aldh1a2","Gdf10"),
            no.legend = TRUE,
            cols.use = c('gainsboro', 'darkred'),
            pt.size = 0.2)

#Output TSNE maps
gene_list_map_on_TSNE <- c("Acta2","Ackr3","Bst1",
                           "Ccl2","Ccl5","Ccl7","Ccl8","Ccl9","Ccl11","Ccl19","Ccl21a",
                           "Cxcl1","Cxcl2","Cxcl13","Cxcl9","Cxcl10",
                           "Cd34","Cd55","Cd248",
                           "Ch25h",
                           "Csf1",
                           "Col15a1","Col4a1",
                           "Des",
                           "Enpp2",
                           "Il7","Il6","Icam1",
                           "Ly6a","Ly6c1",
                           "Madcam1","Mfge8",
                           "Nos2",
                           "Ptgis",
                           "Sfrp4",
                           "Tnfsf13b","Tnfsf11", "Tinagl1",
                           "Timd4",
                           "Vcam1")

for(i in 1:length(gene_list_map_on_TSNE)){
  setwd("C:/Users/jpe12/Desktop/scRNA-seq_CoreHunting")
  gene_i <- gene_list_map_on_TSNE[i]
  print(gene_i)
  filename_i <- paste(sample_1, "_", gene_i, ".eps", sep="")
  print(filename_i)
  FeaturePlot(sample_1_seurat_minus, features.plot = gene_i,
              no.legend = TRUE,
              cols.use = c('gainsboro', 'darkred'),
              pt.size = 0.2)
  ggsave(filename_i, width = 5.8, height = 6, units = "cm")
}
for(i in 1:length(gene_list_map_on_TSNE)){
  setwd("C:/Users/jpe12/Desktop/scRNA-seq_CoreHunting")
  gene_i <- gene_list_map_on_TSNE[i]
  print(gene_i)
  filename_i <- paste(sample_1, "_", gene_i,"_scale" , ".eps", sep="")
  print(filename_i)
  FeaturePlot(sample_1_seurat_minus, features.plot = gene_i,
              no.legend = FALSE,
              cols.use = c('gainsboro', 'darkred'),
              pt.size = 0.2)
  ggsave(filename_i, width = 5.8, height = 6, units = "cm")
}

FeaturePlot(sample_1_seurat_minus, features.plot = c("Bmp1", "Bmp2", "Bmp3","Bmp4",
                                               "Bmp5","Bmp7"),
            min.cutoff = "q9",
            cols.use = c('gainsboro', 'darkred'),
            pt.size = 0.1)

FeaturePlot(sample_1_seurat_minus, features.plot = c("Pparg","Fabp4","Cd34",
                                                     "Col15a1","Nes","Atxn1"),
            min.cutoff = "q9",
            cols.use = c('gainsboro', 'darkred'),
            pt.size = 0.1)
FeaturePlot(sample_1_seurat_minus, features.plot = c("Bmp1", "Bmp2", "Bmp3","Bmp4",
                                                     "Bmp5","Bmp7"),
            min.cutoff = "q9",
            cols.use = c('gainsboro', 'darkred'),
            pt.size = 0.1)

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
write.table(sample_1_seurat_minus.markers, paste(getwd(),"/",sample_1,"_DEGs.csv", sep=""), dec=".", sep=",")
write.table(sample_1_seurat_minus.markers_20, paste(getwd(),"/",sample_1,"_DEGs_Top20.csv", sep=""), dec=".", sep=",")
write.table(sample_1_seurat_minus.markers_40, paste(getwd(),"/",sample_1,"_DEGs_Top40.csv", sep=""), dec=".", sep=",")


#Plot differential expression of Marker genes
#DoHeatmap(object = sample_1_seurat_minus, use.scaled = TRUE, group.by = "ident",
 #         genes.use = sample_1_seurat_minus.markers_20$gene, title = paste(sample_1, sep=" "),
  #        remove.key = FALSE,
   #       col.low = "white",
    #      col.mid =  "oldlace",
     #     col.high = "brown")
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
myfiles <- list(sample_1_seurat_minus.markers)
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
  colnames(table_i_subset_GO_weight)[2] <- names(split_clusters)[i]
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
title <- c("Differential GO scRNA mLN")

#reorder according to Heatmap in main figure
data_heatmap_matrix <- data_heatmap_matrix[,new.cluster.ids]



#reorder according to Heatmap in main figure
data_heatmap_matrix <- data_heatmap_matrix[,c("FSCx","BFSC","MRC",
                                              "mTFSC1","mTFSC2",
                                              "TFSCx","TFSC3",
                                              "Pery","CD8FSC",
                                              "Neuro_FSC",
                                              "FSCy1_1","FSCy1_2_3",
                                              "FSCy3","FSCy4")]

pheatmap(data_heatmap_matrix, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = FALSE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("white", "grey","deepskyblue2"), space="rgb")(128),
         main = title)


