#Author: Joern Pezoldt
#Function: Compare DEGs within cluser for two different comparison





#Seurat 
library("Seurat")
library("cowplot")
#library("Matrix")
#library("magrittr")
library("dplyr")
#not installed on Fameux
#library("pryr")
#library("svglite")

#####
#Load
#####
seurat_pLN_Tx_pLN <- readRDS(file="/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/seurat/pLN_TxpLN/pLN_Tx_pLN.Rds")
seurat_pLN_Tx_mLN <- readRDS(file="/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/seurat/pLN_TxmLN/pLN_Tx_mLN.Rds")

####################
#Expression over seurat_pLN_Tx_pLN
####################
LN_minus <- seurat_pLN_Tx_pLN
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

DEG_within_cluster_pLN_Tx_pLN <- DEG_l


####################
#Expression over seurat_pLN_Tx_mLN
####################
LN_minus <- seurat_pLN_Tx_mLN
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

DEG_within_cluster_pLN_Tx_mLN <- DEG_l

#####
#Compare intersect across all clusters
#####
t_overlap_ratio_up <- matrix(nrow = 0, ncol = length(DEG_within_cluster_pLN_Tx_pLN))
t_overlap_ratio_down <- matrix(nrow = 0, ncol = length(DEG_within_cluster_pLN_Tx_pLN))


for(i in 1:length(DEG_within_cluster_pLN_Tx_mLN)){
  #take significant rownames=GeneSymbols
  DEG_1_cluster_up_i <- rownames(subset(DEG_within_cluster_pLN_Tx_mLN[[i]], p_val_adj < 0.05 & avg_logFC > 0))
  DEG_1_cluster_down_i <- rownames(subset(DEG_within_cluster_pLN_Tx_mLN[[i]], p_val_adj < 0.05 & avg_logFC < 0))
  print("List 1")
  print(i)
  print(DEG_1_cluster_up_i)
  print(DEG_1_cluster_down_i)
  print(length(DEG_1_cluster_up_i))
  print(length(DEG_1_cluster_down_i))
  
  #compare against significant rownames=GeneSymbols of second list
  overlap_up_row_i <- c()
  overlap_down_row_i <- c()
  for(j in 1:length(DEG_within_cluster_pLN_Tx_pLN)){
    DEG_2_cluster_up_i <- rownames(subset(DEG_within_cluster_pLN_Tx_pLN[[j]], p_val_adj < 0.05 & avg_logFC > 0))
    DEG_2_cluster_down_i <- rownames(subset(DEG_within_cluster_pLN_Tx_pLN[[j]], p_val_adj < 0.05 & avg_logFC < 0))
    print("List 2")
    print(j)
    print(DEG_2_cluster_up_i)
    print(DEG_2_cluster_down_i)
    print(length(DEG_2_cluster_up_i))
    print(length(DEG_2_cluster_down_i))
    #Number of toatal genes in Up comparison
    total_n_up_i <- length(unique(c(DEG_1_cluster_up_i,DEG_2_cluster_up_i)))
    #Number of toal gene in Down comparison
    total_n_down_i <- length(unique(c(DEG_1_cluster_down_i,DEG_2_cluster_down_i)))
    
    #calculate % of 1 in total of 1 & 2
    fraction_1_in_12_up_i <- length(DEG_1_cluster_up_i) / total_n_up_i * 100
    fraction_1_in_12_down_i <- length(DEG_1_cluster_down_i) / total_n_down_i * 100
    
    #store in vector
    overlap_up_row_i <- c(overlap_up_row_i, fraction_1_in_12_up_i)
    overlap_down_row_i <- c(overlap_down_row_i, fraction_1_in_12_down_i)
  }
  t_overlap_ratio_up <- rbind(t_overlap_ratio_up, overlap_up_row_i)
  t_overlap_ratio_down <- rbind(t_overlap_ratio_down, overlap_down_row_i)
}
#name columns and rows according to Tx LN used to compare to pLN
colnames(t_overlap_ratio_up) <- c(paste(0:(ncol(t_overlap_ratio_up)-1),"_","Tx_pLN",sep=""))
rownames(t_overlap_ratio_up) <- c(paste(0:(nrow(t_overlap_ratio_up)-1),"_","Tx_mLN",sep=""))

#name columns and rows according to Tx LN used to compare to pLN
colnames(t_overlap_ratio_down) <- c(paste(0:(ncol(t_overlap_ratio_down)-1),"_","Tx_pLN",sep=""))
rownames(t_overlap_ratio_down) <- c(paste(0:(nrow(t_overlap_ratio_down)-1),"_","Tx_mLN",sep=""))

#Heatmap
pheatmap(t_overlap_ratio_up, cluster_cols = FALSE, cluster_rows = FALSE)
pheatmap(t_overlap_ratio_down, cluster_cols = FALSE, cluster_rows = FALSE)

#Note: Cluster 09 is CD34 positive adventitial cluster, but has very low overlap of DEGs
#Identify DEGs solely maintained in Tx mLN
DEGs_TxpLN_Cd34_up <- rownames(subset(DEG_within_cluster_pLN_Tx_pLN[[10]], p_val_adj < 0.05 & avg_logFC > 0))
DEGs_TxpLN_Cd34_down <- rownames(subset(DEG_within_cluster_pLN_Tx_pLN[[10]], p_val_adj < 0.05 & avg_logFC < 0))

DEGs_TxmLN_Cd34_up <- rownames(subset(DEG_within_cluster_pLN_Tx_mLN[[10]], p_val_adj < 0.05 & avg_logFC > 0))
DEGs_TxmLN_Cd34_down <- rownames(subset(DEG_within_cluster_pLN_Tx_mLN[[10]], p_val_adj < 0.05 & avg_logFC < 0))


