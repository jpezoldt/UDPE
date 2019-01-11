#!/usr/bin/env Rscript
#args <- commandArgs(trailingOnly = T)
#blah = args[1]

# Autor: Joern Pezoldt
# 08.08.2018
# Function:
#1) Compile several datasets from kinetic into on cde object
#2) Classify cells using transcriptional signatures
#3) Check classification using scmap
#4) Perform Trajectory analysis
#5) Extract DEGs (e.g. TFs) per key cluster at branching point

#Libraries
library(monocle)
library(pheatmap)
library(reshape2)
library(cellrangerRkit)
library(biomaRt)
library(RcisTarget)
set.seed(123)

#####
#Notes - Important
#####
#Script only works if all datasets have the same phenotypcial data (gene id annotation files)

#####
#Global variable
#####
condition = "day0_10_24_56_300"
datasets <- c("d000","d010","d024","d056","d300")
#Number of cells in which a gene should be expressed to be included in analysis
num_cells_exp <- 20
#Export PATH
PATH_output <- "/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/monocle/day0_10_24_56_300/Cleared"
PATH_CDS_objects <- "/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/monocle/day0_10_24_56_300/Cleared/CDS_objects"
#Signatures from mLN
NC_Pezoldt_Pasztoi_2018 <- read.csv("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/2018_NatComm_Pezoldt_Pasztoi/NC_2018_mLN_Clusters.csv",sep = ";")


#####
#Load and compile data into cde object
#####
#Timepoint 1: e.g. neonatal/d0
dir_d0 <- "/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/270_2018-01-03_scRNA-Seq_mLN_ontogeny/Data/C11_neo_mLN"
d0 <- load_cellranger_matrix(dir_d0, genome = "mm10")

#Timepoint 2: e.g day10
dir_d10 <- "/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/270_2018-01-03_scRNA-Seq_mLN_ontogeny/Data/C12_d10_mLN"
d10 <- load_cellranger_matrix(dir_d10, genome = "mm10")

#Timepoint 3: e.g. day24
dir_d24 <- "/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/270_2018-01-03_scRNA-Seq_mLN_ontogeny/Data/D1_d24_mLN"
d24 <- load_cellranger_matrix(dir_d24, genome = "mm10")

#Timepoint 4: e.g. day56 (merge two datasets)
dir_d56 <- "/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/244_scRNA-Seq_mLN_pLN_SPF/Data/10X_results/L1700567_mLN_SPF"
d56 <- load_cellranger_matrix(dir_d56, genome = "mm10")
#dir_d56_2 <- "~/NAS2/pezoldt/Data/scRNAseq/252_2017-07-19_scRNA-Seq_mLN_pLN_SPF_GF/Data/10x/mLN_SPF/B6"
#d56_2 <- load_cellranger_matrix(dir_d56_2, genome = "mm10")

#Timepoint 6: e.g. day300
dir_d300 <- "/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/275_2018-04-23_scRNA-Seq_mLN_SPF_45wk/Data/E7_SPF_mLN_42wk"
d300 <- load_cellranger_matrix(dir_d300, genome = "mm10")

#list datasets
l_GBCMs <- list(d0,d10,d24,d56,d300)
#,d24,d56)
names(l_GBCMs) <- datasets

#quick glance at data (first entry number of genes, second entry number of cells)
lapply(l_GBCMs, function(x) dim(exprs(x)))

#Rename first columns
l_cds <- lapply(l_GBCMs, function(x) {
  feat <- fData(x)
  names(feat) <- c('id', 'gene_short_name')
  newCellDataSet(exprs(x),
                 phenoData = new("AnnotatedDataFrame", data = pData(x)),
                 featureData = new("AnnotatedDataFrame", data = feat),
                 lowerDetectionLimit = 1,
                 expressionFamily = negbinomial.size())
})

#####
#Keep only specific cells from D0 dataset
#####
# Note: Selection performed via Seurat
# At this stage of data processing both d0 Seurat and monocle object have the same number of cells
# Choose cells over index of row
PATH_Output <- "/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/seurat/Ontogney_Multialign/D0_to_D56"
cells_include <- readRDS(paste(PATH_Output,"/Day0_cells_include.Rds", sep=""))
#Row index
index <- as.numeric(unlist(lapply(strsplit(cells_include, "_"), function(x){x[[3]]})))
#make cds object for selection
expression <- exprs(l_GBCMs[[1]])[,index]
phenoData <- data.frame(barcode = pData(l_GBCMs[[1]])[index,])
rownames(phenoData) <- phenoData$barcode
feat <- fData(l_GBCMs[[1]])
names(feat) <- c('id', 'gene_short_name')
cds_cleared_d0 <- newCellDataSet(expression,
                                 phenoData = new("AnnotatedDataFrame", data = phenoData),
                                 featureData = new("AnnotatedDataFrame", data = feat),
                                 lowerDetectionLimit = 1,
                                 expressionFamily = negbinomial.size())

#replace element in list of cds
l_cds[[1]] <- cds_cleared_d0


#Add time_point and Cell_id column to pData
for(i in seq(length(datasets))){
  print(datasets[i])
  pData(l_cds[[i]])$Time_point <- rep(datasets[i], nrow(pData(l_cds[[i]])))
  pData(l_cds[[i]])$Cell_id <- paste(datasets[i],"_",1:nrow(pData(l_cds[[i]])),sep="")
}

#####Merge datasets to obtain one cde object
#build lists for pData
l_pData <- list()
#l_fData <- list()
l_expData <- list()
for(i in seq(length(datasets))){
  l_pData[[i]] <- pData(l_cds[[i]])
  l_expData[[i]] <- l_cds[[i]]@assayData$exprs
}
lapply(l_pData, nrow)
lapply(l_expData, nrow)
lapply(l_expData, ncol)

#rbind pData
pData_all_t <- do.call("rbind", l_pData)
rownames(pData_all_t) <- pData_all_t$Cell_id
#sparse Matrices
#Note: Only possible if the colnames (gene identifiers are identical across all conditions)
expData_all_t <- do.call("cbind", l_expData)
colnames(expData_all_t) <- pData_all_t$Cell_id
nrow(pData_all_t)
ncol(expData_all_t)
nrow(expData_all_t)
nrow(fData(l_cds[[1]]))

#Full data cde object
cde_all <- newCellDataSet(expData_all_t,
               phenoData = new("AnnotatedDataFrame", data = pData_all_t),
               featureData = new("AnnotatedDataFrame", data = fData(l_cds[[1]])),
               lowerDetectionLimit = 1,
               #for UMI
               expressionFamily = negbinomial.size())
dim(cde_all@assayData$exprs)
dim(pData(cde_all))
dim(fData(cde_all))

#####
#QC
#####
#Estimate factor matrices
cde_all <- estimateSizeFactors(cde_all)
cde_all <- estimateDispersions(cde_all)

#Filter low qualitiy cells
cde_all <- detectGenes(cde_all, min_expr = 0.1)
print(head(fData(cde_all)))
nrow(fData(cde_all))

#Thresh by number of cells gene is expressed in
expressed_genes <- row.names(subset(fData(cde_all), num_cells_expressed >= num_cells_exp))
saveRDS(expressed_genes, paste(PATH_output,"/expressed_genes.Rds",sep = ""))
print(head(pData(cde_all)))

#Check distribution of expression across condition
pData(cde_all)$Total_mRNAs <- Matrix::colSums(exprs(cde_all))
pData(cde_all)$Total_n_genes <- nrow(fData(cde_all)) - Matrix::colSums(exprs(cde_all)==0)
print(head(pData(cde_all)))

#Estimate cutoff
hist(log2(Matrix::colSums(exprs(cde_all))))

#Plot detected gene number range over conditions
upper_bound <- 250
lower_bound <- 5000
qplot(Total_n_genes, data = pData(cde_all), color = Time_point, geom = "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)
levels(as.factor(pData(cde_all)$Time_point))

#Ribosomal protein genes
rpl.genes <- grep(pattern = "^Rpl", x = fData(cde_all)$gene_short_name, value = TRUE)
rps.genes <- grep(pattern = "^Rps", x = fData(cde_all)$gene_short_name, value = TRUE)
ribo.genes <- c(rpl.genes, rps.genes)
#get ensemble Ids
ribo.genes.ensembl <- subset(fData(cde_all), gene_short_name %in% ribo.genes)$id
#subset expression matrix and perform colmeans per cell
ribo.genes_expr <- colSums(as.matrix(exprs(cde_all)[rownames(exprs(cde_all))%in% ribo.genes.ensembl,]))/colSums(as.matrix(exprs(cde_all)))
pData(cde_all)$ribo_expression <- ribo.genes_expr
#Add Gene signature averages to pData
#Cdk cell cycle genes
# Note: For more precies exclusion of cell-cycle effects check out https://satijalab.org/seurat/cell_cycle_vignette.html
cdk.genes <- grep(pattern = "^Cdk", x = fData(cde_all)$gene_short_name, value = TRUE)
#get ensemble Ids
cdk.genes.ensembl <- subset(fData(cde_all), gene_short_name %in% cdk.genes)$id
#subset expression matrix and perform colmeans per cell
cdk.genes_expr <- colSums(as.matrix(exprs(cde_all)[rownames(exprs(cde_all))%in% cdk.genes.ensembl,]))/colSums(as.matrix(exprs(cde_all)))
pData(cde_all)$cdk_expression <- cdk.genes_expr
#Mitochondrial genes
mito.genes <- grep(pattern = "^mt-", x = fData(cde_all)$gene_short_name, value = TRUE)
#get ensemble Ids
mito.genes.ensembl <- subset(fData(cde_all), gene_short_name %in% mito.genes)$id
#subset expression matrix and perform colmeans per cell
mito.genes_expr <- colSums(as.matrix(exprs(cde_all)[rownames(exprs(cde_all))%in% mito.genes.ensembl,]))/colSums(as.matrix(exprs(cde_all)))
pData(cde_all)$mito_expression <- mito.genes_expr
#Heat-shock proteins
hsp.genes <- grep(pattern = "^Hsp", x = fData(cde_all)$gene_short_name, value = TRUE)
#get ensemble Ids
hsp.genes.ensembl <- subset(fData(cde_all), gene_short_name %in% hsp.genes)$id
#subset expression matrix and perform colmeans per cell
hsp.genes_expr <- colSums(as.matrix(exprs(cde_all)[rownames(exprs(cde_all))%in% hsp.genes.ensembl,]))/colSums(as.matrix(exprs(cde_all)))
pData(cde_all)$hsp_expression <- hsp.genes_expr

#Get an overview of distributions
#Take only cells were number of mRNA and genes correlate linear on main axis
plot(pData(cde_all)$Total_n_genes, pData(cde_all)$Total_mRNAs)
# to do so: calculate regression per data point -> eliminate cells with too steep regression -> see below
#Regression of n_gene and total
pData(cde_all)$exp_regression <- pData(cde_all)$Total_mRNAs / pData(cde_all)$Total_n_genes


d0 <- subset(pData(cde_all), Time_point == "d000")
#d10 <- subset(pData(cde_all), Time_point == "d010")
#d24 <- subset(pData(cde_all), Time_point == "d024")
#Plot ribo counts -> cutoff ~0.28
plot(d0$Total_n_genes, d0$ribo_expression)
#plot(d10$Total_n_genes, d10$ribo_expression)
#plot(d24$Total_n_genes, d24$ribo_expression)
#Plot cdk counts -> no cutoff
plot(pData(cde_all)$Total_n_genes, pData(cde_all)$cdk_expression)
plot(d0$Total_n_genes, d0$cdk_expression)
#plot(d10$Total_n_genes, d10$cdk_expression)
#plot(d24$Total_n_genes, d24$cdk_expression)
#Plot Hsp counts -> cutoff ~0.09
plot(pData(cde_all)$Total_n_genes, pData(cde_all)$hsp_expression)
plot(d0$Total_n_genes, d0$hsp_expression)
#plot(d10$Total_n_genes, d10$hsp_expression)
#plot(d24$Total_n_genes, d24$hsp_expression)
#Plot mito counts -> cutoff ~0.06
plot(pData(cde_all)$Total_n_genes, pData(cde_all)$mito_expression)
plot(d0$Total_n_genes, d0$mito_expression)
#plot(d10$Total_n_genes, d10$mito_expression)
#plot(d24$Total_n_genes, d24$mito_expression)


#Filter the dataset by content of 'pData'
valid_cells <- row.names(subset(pData(cde_all),
                                (Total_n_genes >= 250 &
                                Total_n_genes <= 4700) &
                                #hsp_expression <= 0.06 &
                                mito_expression <= 0.045 &
                                #ribo_expression <= 0.28 &
                                exp_regression <= 4.0))
cde_all <- cde_all[,valid_cells]
dim(cde_all@assayData$exprs)
dim(pData(cde_all))
dim(fData(cde_all))

#Double check QC output with implemented cutoffs
#expression
plot(pData(cde_all)$Total_n_genes, pData(cde_all)$Total_mRNAs)
#Plot ribo counts -> cutoff ~0.28
plot(pData(cde_all)$Total_n_genes, pData(cde_all)$ribo_expression)
#Plot cdk counts -> no cutoff
plot(pData(cde_all)$Total_n_genes, pData(cde_all)$cdk_expression)
#Plot Hsp counts -> cutoff ~0.09
plot(pData(cde_all)$Total_n_genes, pData(cde_all)$hsp_expression)
#Plot mito counts -> cutoff ~0.06
plot(pData(cde_all)$Total_n_genes, pData(cde_all)$mito_expression)

#####Check for lognormal distribution
# Log-transform each value in the expression matrix.
#L <- log(exprs(cde_all[expressed_genes,]))
# Scale each gene
#melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
# Plot the distribution of scaled expression
#qplot(value, geom = "density", data = melted_dens_df) +
#  stat_function(fun = dnorm, size = 0.5, color = 'red') +
#  xlab("Standardized log(FPKM)") +
#  ylab("Density")

#####
#Define subsets
#####
#2) Via variable genes
#Identify variable genes
disp_table <- dispersionTable(cde_all)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cde_all <- setOrderingFilter(cde_all, unsup_clustering_genes$gene_id)
#black dots are used for ordering
plot_ordering_genes(cde_all)
nrow(unsup_clustering_genes)
#plot PCs explaining variance
plot_pc_variance_explained(cde_all, return_all = F) # norm_method = 'log',
#Note around 16 PCs describe the difference
#Count how many dimensions to use
cde_all <- reduceDimension(cde_all, max_components = 2, num_dim = 12,
                        reduction_method = 'tSNE', verbose = T)
cde_all <- clusterCells(cde_all, num_clusters = 13)
plot_cell_clusters(cde_all, 1, 2, color = "CellType", markers = c("Pdpn", "Pecam1","Ccl19"))
#plot cell clusters according to media/condition column of pData(HSMM/cds)
#plot_cell_clusters(cde_all, 1, 2, color = "Time_point")
plot_cell_clusters(cde_all, 1, 2, color = "Cluster")

#substract/regress factors
#Notes
#ribo and number of expressed genes make no big impact on clustering
cde_all <- reduceDimension(cde_all, max_components = 2, num_dim = 16,
                        reduction_method = 'tSNE',
                        residualModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression  + Time_point",
                        verbose = T)
cde_all <- clusterCells(cde_all, num_clusters = 14)
plot_cell_clusters(cde_all, 1, 2, color = "Cluster")
plot_cell_clusters(cde_all, 1, 2, color = "Cluster", markers = c("Il15")) #+ facet_wrap(~Time_point)
plot_cell_clusters(cde_all, 1, 2, color = "Cluster", markers = c("Pecam1","Cd34","Tnfsf11","Aldh1a2","Cxcl13",
                                                                "Nfkb1","Ltbr","Icam1","Vcam1","Acta2",
                                                                "Cd248","Ackr3","Cxcl13","Nes","Cdk1"),cell_size = 0.5)
plot_cell_clusters(cde_all, 1, 2, color = "Cluster") + facet_wrap(~Time_point)

#Save RDS
saveRDS(cde_all, file=paste(PATH_CDS_objects,"/day0_10_24_56_300_ALL_withoutNEO.Rds",sep=""))

#####
#Endothelial cell subset identification
#####
#Calculate mean expression of Pecam1 ove and add to pData
#Endothelial gene
Endothelial.genes <- c("Pecam1")
#get ensemble Ids
Endothelial.genes.ensembl <- subset(fData(cde_all), gene_short_name %in% Endothelial.genes)$id
#subset expression matrix and perform colmeans per cell
#Endothelial.genes_expr <- colSums(as.matrix(exprs(cde_all)[rownames(exprs(cde_all))%in% Endothelial.genes.ensembl,]))/colSums(as.matrix(exprs(cde_all)))
Endothelial.genes_expr <- as.matrix(exprs(cde_all)[rownames(exprs(cde_all)) %in% Endothelial.genes.ensembl,])
pData(cde_all)$Endothelial_expression <- Endothelial.genes_expr

#Calculate mean expression of Endothelial genes
# For each cluster calculate mean expression for Endothelial genes
# Store cluster IDs in vector
Endothelial.classificator <- c()
for(i in seq(length(levels(pData(cde_all)$Cluster)))){
  cluster_ID_i <- levels(pData(cde_all)$Cluster)[i]
  print(cluster_ID_i)
  mean_Endothelial.genes <- mean(subset(pData(cde_all), Cluster %in% cluster_ID_i)$Endothelial_expression)
  print(mean_Endothelial.genes)
  Endothelial.classificator[i] <- mean_Endothelial.genes
}
names(Endothelial.classificator) <- levels(pData(cde_all)$Cluster)

#Get all clusters with expression above threshold
# Note: Manually check average expressions
#       For Ontogeny mLN Dataset cutoff is 1.5
Endothelial_subsets <- Endothelial.classificator[Endothelial.classificator >= 1.5]
Endothelial_subsets
SC_subsets <- Endothelial.classificator[Endothelial.classificator < 1.5]
SC_subsets

#Separate data in to Endothelial and SC
Endothelial_cells <- row.names(subset(pData(cde_all),
                                Cluster %in% names(Endothelial_subsets)))
SC_cells <- row.names(subset(pData(cde_all),
                                      Cluster %in% names(SC_subsets)))

#Generate cds object for Endothelial and SC cells
cde_SC <- cde_all[,SC_cells]
cde_Endothelial <- cde_all[,Endothelial_cells]

#Save RDS
saveRDS(cde_Endothelial, file=paste(PATH_CDS_objects,"/day0_10_24_56_300_Endothelial_withoutNEO.Rds",sep=""))


##############
#Extract RCM
##############
#Generate ReadCountMatrix with GeneSymbols for export
#Endothelial cells
#RCM_endothelial <- as.matrix(exprs(cde_Endothelial))
#mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
#genes <- rownames(RCM_endothelial)
#RCM_endothelial <- cbind(genes, RCM_endothelial)
#rownames(RCM_endothelial) <- c()
#G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol"),values=genes,mart= mart)
#RCM_endothelial_export <- merge(G_list,RCM_endothelial,by.x="ensembl_gene_id",by.y="genes")
#RCM_endothelial_export <- RCM_endothelial_export[,2:ncol(RCM_endothelial_export)]
#write.table(RCM_endothelial_export, paste(PATH_output,"/RCM_endothelial_export.txt",sep=""), sep = "\t")

#non-Endothelial SCs
#RCM_SC <- as.matrix(exprs(cde_SC))
#mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
#genes <- rownames(RCM_SC)
#RCM_SC <- cbind(genes, RCM_SC)
#rownames(RCM_endothelial) <- c()
#G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol"),values=genes,mart= mart)
#RCM_SC_export <- merge(G_list,RCM_SC,by.x="ensembl_gene_id",by.y="genes")
#RCM_SC_export <- RCM_SC_export[,2:ncol(RCM_SC_export)]
#write.table(RCM_SC_export, paste(PATH_output,"/RCM_SC_export.txt",sep=""), sep = "\t")


###############
#Analyze SC cell
###############
#####
#Regress and cluster SCs
#####
#plot_pc_variance_explained(cde_SC, return_all = F)
cde_SC <- reduceDimension(cde_SC, max_components = 2, num_dim = 11,
                          reduction_method = 'tSNE',
                          residualModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression + Time_point",
                          verbose = T)
cde_SC <- clusterCells(cde_SC, num_clusters = 12)

plot_cell_clusters(cde_SC, 1, 2, color = "Cluster") + facet_wrap(~Time_point)
#plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Cd34")) + facet_wrap(~Time_point)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Cd34","Tnfsf11","Aldh1a2","Cxcl13",
                                                                "Nfkb1","Ltbr","Icam1","Vcam1","Acta2",
                                                                "Cd248","Ackr3","Cdk1","Nes"),cell_size = 0.5)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Tagln"),
                   cell_size = 0.5)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Cxcl13"),
                   cell_size = 0.5)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Tnfsf11"),
                   cell_size = 0.5)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Il21"),
                   cell_size = 0.5)

plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Acta2","Tagln","Lmod1","Tinagl1"),
                   cell_size = 0.5)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Vcam1","Icam1","Tnfsf11","Cxcl13"),
                   cell_size = 0.5)

#Save RDS
saveRDS(cde_SC, file=paste(PATH_CDS_objects,"/day0_10_24_56_300_SC_withoutNEO.Rds",sep=""))
cde_SC <- readRDS(file=paste(PATH_CDS_objects,"/day0_10_24_56_300_SC_withoutNEO.Rds",sep=""))
#####
#Separate SCs into subsets
#####
#Calculate mean expression of Pecam1 ove and add to pData
#PvC gene------------------------------
PvC.genes <- c("Acta2","Tagln","Lmod1","Tinagl1")
#get ensemble Ids
PvC.genes.ensembl <- subset(fData(cde_SC), gene_short_name %in% PvC.genes)$id
#subset expression matrix and perform colmeans per cell
PvC.genes_expr <- as.matrix(exprs(cde_SC)[rownames(exprs(cde_all)) %in% PvC.genes.ensembl,])
#Add_expression as columns to pData()
pData(cde_SC) <- cbind(pData(cde_SC),t(PvC.genes_expr))
colnames(pData(cde_SC))[(ncol(pData(cde_SC))-3):ncol(pData(cde_SC))] <- PvC.genes
#Calculate mean expression of PvC genes
# For each cluster calculate mean expression for Endothelial genes
# Store cluster IDs in vector
PvC.classificator <- list()
PvC.classificator_j <- c()
for(j in 1:length(PvC.genes)){
  PvC.genes_j <- PvC.genes[j]
  for(i in seq(length(levels(pData(cde_SC)$Cluster)))){
    cluster_ID_i <- levels(pData(cde_SC)$Cluster)[i]
    print(cluster_ID_i)
    mean_PvC.genes <- mean(subset(pData(cde_SC), Cluster %in% cluster_ID_i)[,PvC.genes_j])
    print(mean_PvC.genes)
    PvC.classificator_j[i] <- mean_PvC.genes
  }
  names(PvC.classificator_j) <- levels(pData(cde_SC)$Cluster)
  PvC.classificator[[j]] <- PvC.classificator_j
}
names(PvC.classificator) <- rev(PvC.genes)

#Make table from list
PvC_score <- colSums(do.call("rbind", PvC.classificator))

#Get all clusters with expression above threshold
# Note: Manually check average expressions
PvC_subsets <- PvC_score[PvC_score >= 12]
PvC_subsets
nonPvC_subsets <- PvC_score[PvC_score < 12]
nonPvC_subsets

#Separate data in to Endothelial and SC
PvC_cells <- row.names(subset(pData(cde_SC),
                                  Cluster %in% names(PvC_subsets)))
nonPvC_cells <- row.names(subset(pData(cde_SC),
                                     Cluster %in% names(nonPvC_subsets)))

#Generate cds object for Endothelial and SC cells
cde_nonPvC <- cde_SC[,nonPvC_cells]
cde_PvC <- cde_SC[,PvC_cells]

print(paste("Number of nonPvC SC:", nrow(pData(cde_nonPvC))))
print(paste("Number of PvC SC:", nrow(pData(cde_PvC))))

#########
#Building trajectories SC------------------------------
#########
cde_SC <- reduceDimension(cde_SC, max_components = 2, num_dim = 10,
                              reduction_method = 'tSNE',
                              residualModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression + Time_point",
                              verbose = T)
cde_SC <- clusterCells(cde_SC, num_clusters = 12)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", show_cell_names = TRUE, cell_name_size = 10) + facet_wrap(~Time_point)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Cd34")) + facet_wrap(~Time_point)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Cd34","Tnfsf11","Aldh1a2","Cxcl13",
                                                                    "Nfkb1","Ltbr","Icam1","Vcam1","Acta2",
                                                                    "Cd248","Gdf10","Cxcl9","Has1","Ackr3",
                                                                    "Cdk1","Il6"), cell_size = 0.5)
#Reduce dimensionality
cde_SC <- reduceDimension(cde_SC, max_components = 2, method = 'DDRTree')
#Build trajectory
cde_SC <- orderCells(cde_SC)

#Inferring the genes----------------------------------------------
plot_cell_trajectory(cde_SC, color_by = "Pseudotime")
plot_cell_trajectory(cde_SC, markers = c("Cdk1","Aldh1a2","Cxcl13","Cd34","Acta2"), use_color_gradient = TRUE, cell_size = 0.5)
plot_cell_trajectory(cde_SC, markers = c("Bst1","Pdpn","Pdgfrb","Cd34"), use_color_gradient = TRUE, cell_size = 0.5)

#Set root state
cde_SC <- orderCells(cde_SC, root_state = 1)

plot_cell_trajectory(cde_SC, color_by = "State") + facet_wrap(~Cluster)
plot_cell_trajectory(cde_SC, color_by = "State") + facet_wrap(~Time_point)
plot_cell_trajectory(cde_SC, color_by = "State", markers = c("Cxcl13"), use_color_gradient = TRUE)  + facet_wrap(~Time_point)
plot_cell_trajectory(cde_SC, markers = c("Cdk1"), use_color_gradient = TRUE) + facet_wrap(~Time_point)

plot_complex_cell_trajectory(cde_SC, color_by = 'State', show_branch_points = T,
                             cell_size = 0.5, cell_link_size = 0.3, root_states = c(1))  + facet_wrap(~Time_point)
plot_complex_cell_trajectory(cde_SC, color_by = 'State', show_branch_points = T,
                             cell_size = 0.5, cell_link_size = 0.3, root_states = c(2))  + facet_wrap(~Cluster)





#########
#Building trajectories nonPvC------------------------------
#########
cde_nonPvC <- reduceDimension(cde_nonPvC, max_components = 2, num_dim = 10,
                              reduction_method = 'tSNE',
                              residualModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression + Time_point",
                              verbose = T)
cde_nonPvC <- clusterCells(cde_nonPvC, num_clusters = 14)
plot_cell_clusters(cde_nonPvC, 1, 2, color = "Cluster", show_cell_names = TRUE, cell_name_size = 10) + facet_wrap(~Cluster)
plot_cell_clusters(cde_nonPvC, 1, 2, color = "Cluster", show_cell_names = TRUE, cell_name_size = 10) + facet_wrap(~Time_point)
plot_cell_clusters(cde_nonPvC, 1, 2, color = "Cluster", markers = c("Cd34")) + facet_wrap(~Cluster)
plot_cell_clusters(cde_nonPvC, 1, 2, color = "Cluster", markers = c("Cd34","Tnfsf11","Aldh1a2","Cxcl13",
                                                                    "Ccl19","Vcam1",
                                                                    "Cd248","Gdf10","Madcam1","Ackr3",
                                                                    "Cdk1","Il6"), cell_size = 0.5)

#"Clock","Hoxb3","Prrx1","Gata6",
#####
#Rename Clusters
#####
#current.cluster.ids <- c("Cd34+Aldh1a2+","Ccl19+Il7+","NEOx","Il6+Cxcl1+",
 #                        "SCx","Cd34+Aldh1a2+","?SC?","Ccl19+Madcam+","ProgProf",
  #                       "Inmt+Cxcl12+","Inmt+","prePvC","LTO")
#current.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13")

#new.cluster.ids <- c("Cd34+Aldh1a2+","Ccl19+Il7+","NEOx","Il6+Cxcl1+",
 #                    "SCx","Cd34+Aldh1a2+","TransSC","Ccl19+Madcam+","ProgProf",
  #                   "Inmt+Cxcl12+","Inmt+","preAdv","LTO")

#pData(cde_nonPvC)$Cluster <- plyr::mapvalues(x = pData(cde_nonPvC)$Cluster, from = current.cluster.ids, to = new.cluster.ids)


#Reduce dimensionality
cde_nonPvC <- reduceDimension(cde_nonPvC, max_components = 2, method = 'DDRTree')
#Build trajectory
cde_nonPvC <- orderCells(cde_nonPvC)

#Inferring the genes----------------------------------------------
diff_test_res <- differentialGeneTest(cde_nonPvC[expressed_genes,],
                                      "~Total_mRNAs + mito_expression + ribo_expression + Time_point",
                                      cores = 8)
#diff_test_res_test <- differentialGeneTest(cde_nonLTO[expressed_genes, pData(cde_nonLTO)$Cluster %in% c(1, 9)], fullModelFormulaStr="~Cluster", cores = 8)
#DEG_1_9 <- subset(diff_test_res, qval < 0.0000001)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
# Note: time trajectories are usefull to identify key genes
#       but not required
#Set trajectory genes in cde
#cde_nonPvC <- readRDS(file=paste(PATH_CDS_objects,condition,"/CDS_objects","/day0_10_24_56_300_nonPvC.Rds",sep=""))
plot_cell_trajectory(cde_nonPvC, color_by = "Pseudotime")
plot_cell_trajectory(cde_nonPvC, markers = c("Cdk1","Cxcl13","Cd34"), use_color_gradient = TRUE, cell_size = 0.5)


#Set root state
cde_nonPvC <- orderCells(cde_nonPvC, root_state = 2)

plot_cell_trajectory(cde_nonPvC, color_by = "Cluster") + facet_wrap(~Cluster)
plot_cell_trajectory(cde_nonPvC, color_by = "State") + facet_wrap(~Time_point)
plot_cell_trajectory(cde_nonPvC, color_by = "State", markers = c("Cd34"), use_color_gradient = TRUE)
plot_cell_trajectory(cde_nonPvC, markers = c("Cxcl13"), use_color_gradient = TRUE) + facet_wrap(~Time_point)

plot_complex_cell_trajectory(cde_nonPvC, color_by = 'State', show_branch_points = F,
                             cell_size = 0.5, cell_link_size = 0.3, root_states = c(2))
plot_complex_cell_trajectory(cde_nonPvC, color_by = 'State', show_branch_points = F,
                             cell_size = 0.5, cell_link_size = 0.3, root_states = c(2))  + facet_wrap(~Time_point)
plot_complex_cell_trajectory(cde_nonPvC, show_branch_points = F,
                             markers = c("Cdk1"), cell_size = 0.5,
                             cell_link_size = 0.3, root_states = c(2))

#Save RDS
#saveRDS(cde_nonPvC, file=paste(PATH_CDS_objects,"/day0_10_24_56_300_nonPvC_withoutNEO.Rds",sep=""))
cde_nonPvC <- readRDS(file=paste(PATH_CDS_objects,"/day0_10_24_56_300_nonPvC_withoutNEO.Rds",sep=""))

#####
#Extract RCM of CDE per Stage
#####
# Note: Input required
cde_RCM_export <- cde_nonPvC
name <- "RCM_cde_nonPvC"

#Prep Biomart intel
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol"),values=genes,mart= mart)

#Stages
States <- levels(pData(cde_RCM_export)$State)
for(i in 1:length(States)){
  #Get State
  State_i <- as.numeric(States[i])
  print(State_i)
  #Get cells in state and obtain RCM
  Cells_state_i <- row.names(subset(pData(cde_RCM_export), Time_point == "d024" & State == State_i))
  RCM_export_i <- as.matrix(exprs(cde_RCM_export)[,Cells_state_i])
  
  #Rename RCM
  genes <- rownames(RCM_export_i)
  RCM_export_i <- cbind(genes, RCM_export_i)
  print(ncol(RCM_export_i))
  rownames(RCM_export_i) <- c()
  RCM_export_i <- merge(G_list,RCM_export_i,by.x="ensembl_gene_id",by.y="genes")
  RCM_export_i <- RCM_export_i[,2:ncol(RCM_export_i)]
  #Export RCM for usage in progenitor profile
  write.table(RCM_export_i, paste(PATH_output,"/",name,"_d24_",State_i,".txt",sep=""), sep = "\t")
}

save_pData <- subset(pData(cde_RCM_export), Time_point == "d024")
write.table(save_pData, paste(PATH_output,"/",name,"_d24_phenotypicalData",".txt",sep=""), sep = "\t")

#Generate ReadCountMatrix with GeneSymbols for export
RCM_export <- as.matrix(exprs(cde_RCM_export))
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- rownames(RCM_export)
RCM_export <- cbind(genes, RCM_export)
rownames(RCM_export) <- c()
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol"),values=genes,mart= mart)
RCM_export <- merge(G_list,RCM_export,by.x="ensembl_gene_id",by.y="genes")
RCM_export <- RCM_export[,2:ncol(RCM_export)]
write.table(RCM_export, paste(PATH_output,"/",name,"txt",sep=""), sep = "\t")

#####
#Extract RCM of CDE per Cluster
#####
# Note: Input required
cde_RCM_export <- cde_selected
name <- "RCM_Progenitors"

#Prep Biomart intel
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol"),values=genes,mart= mart)

#Stages
Clusters <- levels(as.factor(pData(cde_RCM_export)$Cluster))
for(i in 1:length(Clusters)){
  #Get State
  Clusters_i <- as.numeric(Clusters[i])
  #Get cells in state and obtain RCM
  Cells_cluster_i <- row.names(subset(pData(cde_RCM_export), Cluster == Clusters_i))
  RCM_export_i <- as.matrix(exprs(cde_RCM_export)[,Cells_cluster_i])
  print(Clusters_i)
  print(ncol(RCM_export_i))
  #Rename RCM
  genes <- rownames(RCM_export_i)
  RCM_export_i <- cbind(genes, RCM_export_i)
  rownames(RCM_export_i) <- c()
  RCM_export_i <- merge(G_list,RCM_export_i,by.x="ensembl_gene_id",by.y="genes")
  RCM_export_i <- RCM_export_i[,2:ncol(RCM_export_i)]
  #Export RCM for usage in progenitor profile
  write.table(RCM_export_i, paste(PATH_output,"/",name,"_Cluster_","_d56_d300_State7",".txt",sep=""), sep = "\t")
}

#combine defined clusters
Cluster_3 <- read.delim(paste(PATH_output,"/",name,"_Cluster_","3",".txt",sep=""))
Cluster_4 <- read.delim(paste(PATH_output,"/",name,"_Cluster_","4",".txt",sep=""))
Cluster_3_4 <- cbind(Cluster_3,Cluster_4)
write.table(Cluster_3_4, paste(PATH_output,"/",name,"_Cluster_","3_4",".txt",sep=""), sep = "\t")


#Generate ReadCountMatrix with GeneSymbols for export
RCM_export <- as.matrix(exprs(cde_RCM_export))
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- rownames(RCM_export)
RCM_export <- cbind(genes, RCM_export)
rownames(RCM_export) <- c()
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol"),values=genes,mart= mart)
RCM_export <- merge(G_list,RCM_export,by.x="ensembl_gene_id",by.y="genes")
RCM_export <- RCM_export[,2:ncol(RCM_export)]
write.table(RCM_export, paste(PATH_output,"/",name,".txt",sep=""), sep = "\t")

#####
#select Stage / Cluster and rerun Trajectory
#####
#Select stage 2 (earliest developmental stage)
stage_2_cells <- row.names(subset(pData(cde_nonPvC),
                            State %in% c("2")))

Cluster_9_cells <- row.names(subset(pData(cde_nonPvC),
                                    Cluster %in% c("9")))

cde_selected <- cde_nonPvC[,Cluster_9_cells]

# Clustering
cde_selected <- reduceDimension(cde_selected, max_components = 2, num_dim = 12,
                              reduction_method = 'tSNE',
                              residualModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression + Time_point",
                              verbose = T)
cde_selected <- clusterCells(cde_selected, num_clusters = 5)
plot_cell_clusters(cde_selected, 1, 2, color = "Cluster", show_cell_names = TRUE, cell_name_size = 10)
plot_cell_clusters(cde_selected, 1, 2, color = "Cluster", show_cell_names = TRUE, cell_name_size = 10) + facet_wrap(~Time_point)
plot_cell_clusters(cde_selected, 1, 2, color = "Cluster", markers = c("Cdk1")) + facet_wrap(~Time_point)
plot_cell_clusters(cde_selected, 1, 2, color = "Cluster", markers = c("Cd34","Tnfsf11","Aldh1a2","Cxcl13",
                                                                    "Nfkb1","Ltbr","Icam1","Vcam1","Acta2",
                                                                    "Cd248","Gdf10","Cxcl9","Has1","Ackr3",
                                                                    "Cdk1","Il6"), cell_size = 0.5)
#Reduce dimensionality
cde_selected <- reduceDimension(cde_selected, max_components = 2, method = 'DDRTree')
#Build trajectory
cde_selected <- orderCells(cde_selected)

#Inferring the genes----------------------------------------------
expressed_genes <- row.names(subset(fData(cde_selected), num_cells_expressed >= num_cells_exp))
diff_test_res <- differentialGeneTest(cde_selected[expressed_genes,],
                                      "~Total_mRNAs + mito_expression + ribo_expression + Time_point",
                                      cores = 8)
#diff_test_res_test <- differentialGeneTest(cde_nonLTO[expressed_genes, pData(cde_nonLTO)$Cluster %in% c(1, 9)], fullModelFormulaStr="~Cluster", cores = 8)
#DEG_1_9 <- subset(diff_test_res, qval < 0.0000001)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
# Note: time trajectories are usefull to identify key genes
#       but not required
#Set trajectory genes in cde
#cde_nonPvC <- readRDS(file=paste(PATH_CDS_objects,condition,"/CDS_objects","/day0_10_24_56_300_nonPvC.Rds",sep=""))
plot_cell_trajectory(cde_selected, color_by = "Pseudotime")
plot_cell_trajectory(cde_selected, markers = c("Cdk1","Aldh1a2","Cxcl13","Cd34"), use_color_gradient = TRUE, cell_size = 0.5)

#Set root state
cde_nonPvC <- orderCells(cde_selected, root_state = 4)

plot_cell_trajectory(cde_selected, color_by = "State") + facet_wrap(~Cluster)
plot_cell_trajectory(cde_selected, color_by = "State") + facet_wrap(~Time_point)
plot_cell_trajectory(cde_selected, color_by = "State", markers = c("Cd34"), use_color_gradient = TRUE)  + facet_wrap(~Time_point)
plot_cell_trajectory(cde_selected, markers = c("Cxcl13"), use_color_gradient = TRUE) + facet_wrap(~Time_point)

plot_complex_cell_trajectory(cde_selected, color_by = 'State', show_branch_points = T,
                             cell_size = 0.5, cell_link_size = 0.3, root_states = c(4))  + facet_wrap(~Time_point)
plot_complex_cell_trajectory(cde_selected, color_by = 'State', show_branch_points = T,
                             cell_size = 0.5, cell_link_size = 0.3, root_states = c(5))  + facet_wrap(~Cluster)

#Save RDS
saveRDS(cde_nonPvC, file=paste(PATH_CDS_objects,"/day0_10_24_56_300_nonPvC_cleared_stage_2.Rds",sep=""))

#####
#Compile Tables for dynGENIE3
#####

# Select Stages of trajectory
# Adventitial
stages_adventi_cells <- subset(pData(cde_nonPvC),
                                  State %in% c("2","3","4","5","6","7"))
# Nonadventitial
stages_nonadventi_cells <- subset(pData(cde_nonPvC),
                                            State %in% c("1","2"))

# Dissect into clusters
time_points <- levels(as.factor(pData(cde_nonPvC)$Time_point))
l_adventi_cells_per_Age <- split(stages_adventi_cells, stages_adventi_cells$Time_point)
names(l_adventi_cells_per_Age) <- time_points
l_nonadventi_cells_per_Age <- split(stages_nonadventi_cells, stages_nonadventi_cells$Time_point)
names(l_nonadventi_cells_per_Age) <- time_points
# Average expression over time point and compile table
#' average_expr
#' average expression across each element of list
#' @param l_cell_substs named list of pData tables containing the cells per subset
#' @param cds_interest CDS object to be sampled from
#' @return t_average_expression table of averaged expression with columns being time-point and rows gene names

f_average_expr <- function(l_cell_subsets, cds_interest){
  #Get gene names
  gene_names <- fData(cde_nonPvC)$gene_short_name
  #initiate empty table
  t_averaged <- data.frame(matrix(ncol = length(names(l_cell_subsets)), nrow = length(gene_names)))
  
  #for each table
  for(i in 1:length(l_cell_subsets)){
    cell_subsets_i <- l_cell_subsets[[i]]
    nrow(cell_subsets_i)
    expr_i <- as.matrix(exprs(cds_interest)[,colnames(exprs(cds_interest)) %in% rownames(cell_subsets_i)])
    expr_mean_i <- rowMeans(expr_i)
    t_averaged[,i] <- expr_mean_i
  }
  colnames(t_averaged) <- names(l_cell_subsets)
  t_averaged <- cbind(gene_names,t_averaged)
  t_averaged <- t_averaged[!duplicated(t_averaged[c("gene_names")]),]
  gene_names_final <- t_averaged[,1]
  t_averaged <- t_averaged[,2:ncol(t_averaged)]
  rownames(t_averaged) <- gene_names_final
  #return
  t_averaged
}

# run function
averaged_adventi_dynGENIE3 <- f_average_expr(l_adventi_cells_per_Age,cde_nonPvC)
averaged_nonadventi_dynGENIE3 <- f_average_expr(l_nonadventi_cells_per_Age,cde_nonPvC)
subset(averaged_adventi_dynGENIE3, rownames(averaged_adventi_dynGENIE3) %in% c("Cd34"))
subset(averaged_nonadventi_dynGENIE3, rownames(averaged_nonadventi_dynGENIE3) %in% c("Cd34"))
#export
write.table(averaged_adventi_dynGENIE3, paste(PATH_output,"/","averaged_adventi_dynGENIE3.txt",sep=""), dec = ".", sep = "\t")
write.table(averaged_adventi_dynGENIE3, paste(PATH_output,"/","averaged_nonadventi_dynGENIE3.txt",sep=""), dec = ".", sep = "\t")

#####
#Signature identification across subsets
#####
cde_to_check <- cde_nonPvC
GeneSymbol <- fData(cde_to_check)$gene_short_name
#subset expression matrix and perform colmeans per cell
Expression <- as.matrix(exprs(cde_to_check))
#Add GeneSymbol as rownames
rownames(Expression) <- GeneSymbol
#For each cluster calculate the mean expression across all genes
l_meaned_expression <- list()
for(i in seq(length(levels(pData(cde_to_check)$Cluster)))){
  cluster_ID_i <- levels(pData(cde_to_check)$Cluster)[i]
  print(cluster_ID_i) 
  cells_cluster_i <- rownames(subset(pData(cde_to_check), Cluster == cluster_ID_i))
  Expression_mean_i <- rowMeans(Expression[,colnames(Expression) %in% cells_cluster_i])
  l_meaned_expression[[i]] <- Expression_mean_i
}
#Prep table for meaned expression per cluster
t_meaned_expression <- as.data.frame(do.call("cbind", l_meaned_expression))
t_meaned_expression <- cbind(GeneSymbol, t_meaned_expression)
colnames(t_meaned_expression) <- c("GeneSymbol", levels(pData(cde_to_check)$Cluster))

#Load DEG tables and make list of clusters
DEGs_core_40 <- NC_Pezoldt_Pasztoi_2018
DEGs_list <- split(DEGs_core_40, DEGs_core_40$cluster)


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
  
  aExp_cluster_i <- subset(t_meaned_expression, rownames(t_meaned_expression) %in% cluster_DEG_i)[,2:ncol(t_meaned_expression)]
  
  Scale_aExp_cluster_i <- apply(aExp_cluster_i, 1, scale)
  Zscore_cluster_i <- rowSums(Scale_aExp_cluster_i)
  Zscore_cluster_i_norm <- Zscore_cluster_i / (nrow(aExp_cluster_i)/10)
  
  #mLN maintained Scores
  Zscore_table_cluster_i <- data.frame(module = c(rep(cluster_ID_i,length(Zscore_cluster_i_norm))),
                                       cluster = colnames(t_meaned_expression)[2:ncol(t_meaned_expression)],
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

title <- paste(condition,"DEGs_core_40",sep = "_")
pdf(paste(PATH_output,"/day0_10_24_56_300_nonPvC_NC_clusters.pdf", sep=""))
marrangeGrob(out, nrow = round(length(DEGs_list) / 3 - 1), ncol = 3, top = title)
dev.off()

#####
#BEAM nonPvC
#####
#Perform BEAM for all branches
#List to store BEAMs
l_BEAM_res <- list()
for(i in 1:3){
  BEAM_res_i <- BEAM(cde_nonPvC, branch_point = i, cores =12)
  l_BEAM_res[[i]] <- BEAM_res_i
}
#names(l_BEAM_res) <- paste("Branch", c(1:2), sep = "_")

plot_genes_branched_heatmap(cde_nonPvC[row.names(subset(l_BEAM_res[[2]],
                                                   qval < 5*1e-2)),],
                            branch_point = 1,
                            num_clusters = length(levels(as.factor(pData(cde_nonPvC)$Cluster))),
                            cores = 8,
                            use_gene_short_name = T,
                            show_rownames = T,
                            branch_colors = c("blue", "white", "red"))

saveRDS(l_BEAM_res, paste(PATH_output,"/","BEAM_day0_10_24_56_300_nonPvC_Named.rds",sep=""))


#Identify TFs at branching points
# Use TF list from Genomatix
Genomatix_murineTFs <- read.table("/home/pezoldt/NAS2/pezoldt/Analysis/WGBS/FSC/Genomatix_TF/ANNa/Mouse/Genomatix_murineTFs_OnlyNames.txt",
                                  sep = "\t", skip = 1)
# DMR TFs
mLN_hypo <- read.table("/home/pezoldt/NAS2/pezoldt/Analysis/WGBS/FSC/Genomatix_TF/ANNa/Mouse/To_Use/Output/Compiled_TFs_Genomatix_hypo_all.txt")
mLN_hypo_TF <- as.character(mLN_hypo$V1)
mLN_hyper <- read.table("/home/pezoldt/NAS2/pezoldt/Analysis/WGBS/FSC/Genomatix_TF/ANNa/Mouse/To_Use/Output/Compiled_TFs_Genomatix_hyper_all.txt")
mLN_hyper_TF <- as.character(mLN_hyper$V1)

# DAR TFs existence matrix
DAR_existence_matrix <- read.table("/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/Motif/Homer_Output/allHomerResults/size_500bp/SPF/SPF_existence_matrix_known_homer_compile.txt")
Unique_pLN_peak_UP <- rownames(subset(DAR_existence_matrix, 
                                pLN_Open_UP == 0 &
                                  mLN_Open_UP == 0 &
                                  pLN_peak_UP == 1 &
                                  mLN_peak_UP == 0 &
                                  pLN_Open_None == 0 &
                                  mLN_Open_None == 0))
Unique_mLN_Open_None <- rownames(subset(DAR_existence_matrix, 
                             pLN_Open_UP == 0 &
                               mLN_Open_UP == 0 &
                               pLN_peak_UP == 0 &
                               mLN_peak_UP == 0 &
                               pLN_Open_None == 0 &
                               mLN_Open_None == 1))

Common_used <- rownames(subset(DAR_existence_matrix, 
                                        !(pLN_Open_UP == 0 &
                                          mLN_Open_UP == 0 &
                                          pLN_peak_UP == 0 &
                                          mLN_peak_UP == 0 &
                                          pLN_Open_None == 0 &
                                          mLN_Open_None == 1) |
                                          !(pLN_Open_UP == 0 &
                                          mLN_Open_UP == 0 &
                                          pLN_peak_UP == 1 &
                                          mLN_peak_UP == 0 &
                                          pLN_Open_None == 0 &
                                          mLN_Open_None == 0)))

#Eliminate ribosomal genes from branch defining list
TFs_genomatix <- Genomatix_murineTFs$V2
#get ensemble Ids
BEAM_TFs_genomatix_2 <- subset(l_BEAM_res[[2]], gene_short_name %in% TFs_genomatix)
BEAM_mLN_hypo_TF <- subset(l_BEAM_res[[2]], gene_short_name %in% mLN_hypo_TF)
BEAM_mLN_hyper_TF <- subset(l_BEAM_res[[2]], gene_short_name %in% mLN_hyper_TF)
BEAM_mLN_Unique_pLN_peak_UP <- subset(l_BEAM_res[[2]], gene_short_name %in% Unique_pLN_peak_UP)
BEAM_mLN_Unique_mLN_Open_None <- subset(l_BEAM_res[[2]], gene_short_name %in% Unique_mLN_Open_None)
BEAM_Common_pLN_peak_UP_mLN_Open_None <- subset(l_BEAM_res[[2]], gene_short_name %in% Common_pLN_peak_UP_mLN_Open_None)
BEAM_Common_used <- subset(l_BEAM_res[[2]], gene_short_name %in% Common_used)
plot_genes_branched_heatmap(cde_nonPvC[row.names(subset(BEAM_mLN_hypo_TF,
                                                        qval < 5*1e-2)),],
                            branch_point = 2,
                            #num_clusters = length(levels(as.factor(pData(cde_nonPvC)$Cluster))),
                            cores = 8,
                            use_gene_short_name = T,
                            show_rownames = T,
                            branch_labels = c("NonAdv", "Adv"))

plot_genes_branched_heatmap(cde_nonPvC[row.names(subset(BEAM_mLN_hyper_TF,
                                                        qval < 5*1e-2)),],
                            branch_point = 2,
                            #num_clusters = length(levels(as.factor(pData(cde_nonPvC)$Cluster))),
                            cores = 8,
                            use_gene_short_name = T,
                            show_rownames = T,
                            branch_labels = c("NonAdv", "Adv"))

plot_genes_branched_heatmap(cde_nonPvC[row.names(subset(BEAM_mLN_Unique_pLN_peak_UP,
                                                        qval < 5*1e-2)),],
                            branch_point = 2,
                            #num_clusters = length(levels(as.factor(pData(cde_nonPvC)$Cluster))),
                            cores = 8,
                            use_gene_short_name = T,
                            show_rownames = T,
                            branch_labels = c("NonAdv", "Adv"))

plot_genes_branched_heatmap(cde_nonPvC[row.names(subset(BEAM_mLN_Unique_mLN_Open_None,
                                                        qval < 5*1e-2)),],
                            branch_point = 2,
                            #num_clusters = length(levels(as.factor(pData(cde_nonPvC)$Cluster))),
                            cores = 8,
                            use_gene_short_name = T,
                            show_rownames = T,
                            branch_labels = c("NonAdv", "Adv"))

plot_genes_branched_heatmap(cde_nonPvC[row.names(subset(BEAM_Common_pLN_peak_UP_mLN_Open_None,
                                                        qval < 5*1e-2)),],
                            branch_point = 2,
                            #num_clusters = length(levels(as.factor(pData(cde_nonPvC)$Cluster))),
                            cores = 8,
                            use_gene_short_name = T,
                            show_rownames = T,
                            branch_labels = c("NonAdv", "Adv"))

plot_genes_branched_heatmap(cde_nonPvC[row.names(subset(BEAM_Common_used,
                                                        qval < 5*1e-2)),],
                            branch_point = 2,
                            #num_clusters = length(levels(as.factor(pData(cde_nonPvC)$Cluster))),
                            cores = 8,
                            use_gene_short_name = T,
                            show_rownames = T,
                            branch_labels = c("NonAdv", "Adv"))




test <- rbind(BEAM_mLN_Unique_mLN_Open_None,BEAM_mLN_Unique_pLN_peak_UP)

plot_genes_branched_heatmap(cde_nonPvC[row.names(subset(test,
                                                        qval < 5*1e-2)),],
                            branch_point = 2,
                            #num_clusters = length(levels(as.factor(pData(cde_nonPvC)$Cluster))),
                            cores = 8,
                            use_gene_short_name = T,
                            show_rownames = T,
                            branch_labels = c("NonAdv", "Adv"))

plot_genes_branched_heatmap(cde_nonPvC[row.names(BEAM_mLN_hyper_TF),],
                            branch_point = 2,
                            num_clusters = length(levels(as.factor(pData(cde_nonPvC)$Cluster))),
                            cores = 8,
                            use_gene_short_name = T,
                            show_rownames = T,
                            branch_labels = c("NonAdv", "Adv"))

plot_genes_branched_heatmap(cde_nonPvC[row.names(BEAM_Common_pLN_peak_UP_mLN_Open_None),],
                            branch_point = 2,
                            num_clusters = length(levels(as.factor(pData(cde_nonPvC)$Cluster))),
                            cores = 8,
                            use_gene_short_name = T,
                            show_rownames = T,
                            branch_labels = c("NonAdv", "Adv"))

plot_genes_branched_heatmap(cde_nonPvC[row.names(subset(l_BEAM_res[[3]],
                                                   qval < 1e-10)),],
                            branch_point = 1,
                            num_clusters = length(levels(as.factor(pData(cde_nonPvC)$Cluster))),
                            cores = 8,
                            use_gene_short_name = T,
                            show_rownames = F,
                            branch_labels = c("NonAdv", "Adv"))

#####
#Plot Genes expression
#####
#set cde of interest
cde_x <- readRDS(file=paste(PATH_CDS_objects,condition,"/CDS_objects","/day0_10_24_56_300_nonLTO.Rds",sep=""))

#Jitter-----------------------------
blast_genes <- row.names(subset(fData(cde_x),
                                gene_short_name %in% c("Nkx2-3", "Ccl9", "Cd34","Gdf10","Ackr3")))
blast_genes <- row.names(subset(fData(cde_x),
                                gene_short_name %in% c("Ccnb2", "Myod1", "Myog","Cdk1")))

plot_genes_jitter(cde_x[blast_genes,],
                  grouping = "State",
                  color_by = "Cluster",
                  min_expr = 0.05)

#Plot gene expression in pseudotime
cde_x_expressed_genes <-  row.names(subset(fData(cde_x),
                                          num_cells_expressed >= 10))
cde_x_filtered <- cde_x[cde_x_expressed_genes,]
my_genes <- row.names(subset(fData(cde_x_filtered),
                             gene_short_name %in% c("Aldh1a2","Cd34","Ccl19","Ackr3")))
cde_x_subset <- cde_x_filtered[my_genes,]
plot_genes_in_pseudotime(cde_x_subset, color_by = "State")


#####
#DEGs across clusters
#####
cde_x_cluster_DEG <- differentialGeneTest(cde_x[cde_x_expressed_genes,],
                                              fullModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression + Time_point",
                                              cores = 8)

#order
cde_x <- orderCells(cde_x, root_state = 4)

#
plot_cell_trajectory(cde_x)
plot_cell_trajectory(cde_x, color_by = "State") + facet_wrap(~Time_point)
plot_cell_trajectory(cde_x, color_by = "State") + facet_wrap(~Cluster)
plot_cell_trajectory(cde_x, markers = "Cdk1")  + facet_wrap(~Cluster)
plot_cell_trajectory(cde_x, markers = "Aldh1a2", use_color_gradient = TRUE)
plot_complex_cell_trajectory(cde_x, color_by = 'State', show_branch_points = T,
                             cell_size = 0.5, cell_link_size = 0.3, root_states = c(2)) + facet_wrap(~Cluster)

#####
#Multifactorial DEG analyis
#####
#Identify DEGs independent of facter (e.g. type of cluster)
# Note: Preliminary analsis indicates that clusters align substantially with states
#cde_x_CS_minusC_DEG <- differentialGeneTest(cde_x,
 #                                     fullModelFormulaStr = "~Cluster + State",
  #                                    reducedModelFormulaStr = "~Cluster",
   #                                    cores = 8)

#####
#Branch Analysis - BEAM
#####
#Analyse branch
#Input: Ordered according to pseudotime
cde_x <- orderCells(cde_x)

#Perform BEAM for all branches
#List to store BEAMs
l_BEAM_res <- list()
for(i in 1:2){
  BEAM_res_i <- BEAM(cde_x, branch_point = i, cores =12)
  l_BEAM_res[[i]] <- BEAM_res_i
}
names(l_BEAM_res) <- paste("Branch", c(1:2), sep = "_")



saveRDS(l_BEAM_res, paste(PATH_output,"/","BEAM_D0_res_all.rds",sep=""))
l_BEAM_res <- readRDS(paste(PATH_output,"/","BEMA_day0_10_24_56_300_nonPvC.rds",sep=""))
#Eliminate ribosomal genes from branch defining list
rpl.genes <- grep(pattern = "^Rpl", x = l_BEAM_res[[2]]$gene_short_name, value = TRUE)
rps.genes <- grep(pattern = "^Rps", x = l_BEAM_res[[2]]$gene_short_name, value = TRUE)
ribo.genes <- c(rpl.genes, rps.genes)
#get ensemble Ids
l_BEAM_res_minus_ribo_2 <- subset(l_BEAM_res[[2]], !(gene_short_name %in% ribo.genes))

#Save BEAM_tables
#for(i in 1:length(l_BEAM_res))
#Branch 1
plot_genes_branched_heatmap(cde_x[row.names(subset(l_BEAM_res[[1]],
                                                   qval < 1e-20)),],
                            branch_point = 1,
                            num_clusters = length(levels(as.factor(pData(cde_x)$Cluster))),
                            cores = 4,
                            use_gene_short_name = T,
                            show_rownames = T)
#Branch 2
plot_genes_branched_heatmap(cde_x[row.names(subset(l_BEAM_res[[2]],
                                                   qval < 1e-5)),],
                            branch_point = 2,
                            num_clusters = length(levels(as.factor(pData(cde_x)$Cluster))),
                            cores = 4,
                            use_gene_short_name = T,
                            show_rownames = T)

#Branch 3
plot_genes_branched_heatmap(cde_x[row.names(subset(l_BEAM_res[[2]],
                                                   qval < 1e-40)),],
                            branch_point = 1,
                            num_clusters = length(levels(as.factor(pData(cde_x)$Cluster))),
                            cores = 16,
                            use_gene_short_name = T,
                            show_rownames = T)

#Branch 3
plot_genes_branched_heatmap(cde_x[row.names(subset(l_BEAM_res[[2]],
                                                   qval < 1e-7)),],
                            branch_point = 4,
                            num_clusters = length(levels(as.factor(pData(cde_x)$Cluster))),
                            cores = 16,
                            use_gene_short_name = T,
                            show_rownames = T)
#####
#Check clusters
#####
expressed_genes <- readRDS(paste(PATH_output,"/expressed_genes.Rds",sep = ""))
diff_test_res <- differentialGeneTest(cde_x[expressed_genes,],
                                      #fullModelFormulaStr = "~Clusters",
                                      cores = 4)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
# Note: time trajectories are usefull to identify key genes
#       but not required
#Set trajectory genes in cde
cde_x <- setOrderingFilter(cde_x, ordering_genes)
plot_ordering_genes(cde_x)

#plot by condition/timepoint
plot_cell_trajectory(cde_x, color_by = "Time_point")
#plot by state
plot_cell_trajectory(cde_x, color_by = "State")
#plot by cluster
plot_cell_trajectory(cde_x, color_by = "State") + facet_wrap(~Cluster)



test_Atf <- subset(l_BEAM_res[[1]], gene_short_name %in% c("Atf1","Atf2","Atf3","Atf4","Atf5",
                                                           "Atf6","Nfkb1","Nkx2-3","Batf","Egr1",
                                                           "Egr2","Bach2",
                                                           "Klf1","Klf2","Klf3","Klf4","Klf5",
                                                           "Stat1","Stat2","Stat3","Stat4",
                                                           "Ctcf","Junb") &
                     qval < 1e-1)

#####
#Ordering based on genes that differ between clusters
#####
# Pick genes that are expressed in >5% of the cell
cde_x <- detectGenes(cde_x, min_expr = 0.1)
fData(cde_x)$use_for_ordering <-
  fData(cde_x)$num_cells_expressed > 0.05 * ncol(cde_x)
#Plot variance explainedby picked genes
plot_pc_variance_explained(cde_x, return_all = F)
#Reduce dimensions using genes
cde_x <- reduceDimension(cde_x,
                         max_components = 2,
                         norm_method = 'log',
                         num_dim = 5,
                         reduction_method = 'tSNE',
                         verbose = T)

#cluster cells
cde_x <- clusterCells(cde_x, num_clusters = 13)
plot_cell_clusters(cde_x, color_by = 'as.factor(Cluster)')

#Decision plot to define Rho and P
plot_rho_delta(cde_x, rho_threshold = 2, delta_threshold = 4 )

#Re-cluster cells implementing cut-off
cde_x <- clusterCells(cde_x,
                      rho_threshold = 70,
                      delta_threshold = 4,
                      skip_rho_sigma = T,
                      verbose = F)

#####
#Distinguish Cluster or State
#####
#Calculation DEGs
# Determine how good a gene explains the state of a cell by checking/building a model
# 1) that knows about the cluster annotation
# 2) that does not know about the cluster annotation
#Compare output of models and choose highest scoring models
#DEG according to Cluster
cde_x_Cluster_DEG <- differentialGeneTest(cde_x,
                                          fullModelFormulaStr = "~Cluster",
                                          cores = 8)
#DEG according to State
cde_x_State_DEG <- differentialGeneTest(cde_x, 
                                        fullModelFormulaStr = "~State",
                                        cores = 8)
#DEG according to pseudotime
cde_x_PseuTim_DEG <- differentialGeneTest(cde_x,
                                          fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                          cores = 8)


# Select genes that are significant at an FDR < 10%
sig_genes_State <- subset(cde_x_State_DEG, qval < 0.01)
sig_genes_Cluster <- subset(cde_x_Cluster_DEG, qval < 0.01)
sig_genes_PseuTim <- subset(cde_x_PseuTim_DEG, qval < 0.01)
saveRDS(sig_genes_State, paste(PATH_output,"/","sig_genes_State.rds",sep=""))
saveRDS(sig_genes_Cluster, paste(PATH_output,"/","sig_genes_Cluster.rds",sep=""))
saveRDS(sig_genes_PseuTim, paste(PATH_output,"/","sig_genes_PseuTim.rds",sep=""))
sig_genes_State <- readRDS(paste(PATH_output,"/","sig_genes_State.rds",sep=""))
sig_genes_Cluster <- readRDS(paste(PATH_output,"/","sig_genes_Cluster.rds",sep=""))
sig_genes_PseuTim <- readRDS(paste(PATH_output,"/","sig_genes_PseuTim.rds",sep=""))

sig_genes_PseuTim_TRUE <- subset(sig_genes_PseuTim, use_for_ordering == "TRUE")
sig_genes_Cluster_TRUE <- subset(sig_genes_Cluster, use_for_ordering == "TRUE")
sig_genes_Cluster_TRUE <- subset(sig_genes_Cluster, use_for_ordering == "TRUE")

#Check intesections
sig_genes_SCP <- Reduce(intersect, list(sig_genes_State$gene_short_name,sig_genes_Cluster$gene_short_name,sig_genes_PseuTim$gene_short_name))
sig_genes_SC <- Reduce(intersect, list(sig_genes_State$gene_short_name,sig_genes_Cluster$gene_short_name))
sig_genes_CP <- Reduce(intersect, list(sig_genes_Cluster$gene_short_name,sig_genes_PseuTim$gene_short_name))
sig_genes_SP <- Reduce(intersect, list(sig_genes_State$gene_short_name,sig_genes_PseuTim$gene_short_name))
print(paste("Number of common sig. SCP:", length(sig_genes_SCP), sep=" "))
print(paste("Number of common sig. SC:", length(sig_genes_SC), sep=" "))
print(paste("Number of common sig. CP:", length(sig_genes_CP), sep=" "))
print(paste("Number of common sig. SP:", length(sig_genes_SP), sep=" "))

# Note: Use genes commonly identified for state and pseudotime
sig_genes_common_SP <- subset(sig_genes_PseuTim, gene_short_name %in% sig_genes_SP)
head(sig_genes_common_SP[,c("gene_short_name", "pval", "qval")])

#Plot genes that describe Pseudotime
#Get Ensemble IDs for genes of interest
sig_genes_SP_ensem <- rownames(subset(sig_genes_PseuTim, gene_short_name %in% sig_genes_SP))

cde_x <- setOrderingFilter(cde_x, sig_genes_SP_ensem)
plot_ordering_genes(cde_x)
plot_pseudotime_heatmap(cde_x[sig_genes_SP_ensem,],
                        num_clusters = length(levels(as.factor(pData(cde_x)$Cluster))),
                        cores = 4,
                        show_rownames = T)
#trend_formula = "~Cluster")


