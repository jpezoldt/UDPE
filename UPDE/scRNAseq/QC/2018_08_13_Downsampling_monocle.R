#Downsample and co-analyze using Monocle

# Autor: Joern Pezoldt
# 13.08.2018
# Function:
#1) Downsample expression matrix
#2) Co-analyze using Monocle


#Libraries
library(monocle)
library(pheatmap)
library(reshape2)
library(cellrangerRkit)
library(Matrix)
library(DropletUtils)

#####
#Global variable
#####
datasets <- c("d056_100","d056_050","d056_020","d056_010","d056_005")
downsample_factor <- c(1.0,0.5,0.2,0.1,0.05)
#Number of cells in which a gene should be expressed to be included in analysis
num_cells_exp <- 50

#####
#Load and compile data into cde object
#####

#day56
#dir_d56 <- "~/NAS2/pezoldt/Data/scRNAseq/244_scRNA-Seq_mLN_pLN_SPF/Data/10X_results/L1700567_mLN_SPF"
dir_d56 <- "/Users/Pezoldt/PowerFolders/R/2017_252_scRNAseq/Data/mLN_SPF/B6"
d56_100 <- load_cellranger_matrix(dir_d56, genome = "mm10")

#Generate dataslots
d56_050 <- d56_100
d56_020 <- d56_100
d56_010 <- d56_100
d56_005 <- d56_100

#list datasets
l_GBCMs <- list(d56_100,d56_050,d56_020,d56_010,d56_005)
names(l_GBCMs) <- datasets

#quick glance at data (first entry number of genes, second entry number of cells)
lapply(l_GBCMs, function(x) dim(exprs(x)))

#Dowsample each countmatrix and generate vector with counts/UMIs per cell
set.seed(1)
#Empty list to store counts per cell
l_count_per_cell <- list()
l_GBCMs_reduced <- list()
for(i in seq(length(l_GBCMs))){
  prop_i <- downsample_factor[i]
  print(prop_i)
  #Downsample and convert to dgCMatrix
  #counts_dgC_i <- as.matrix(l_GBCMs[[i]])
  counts_dgC_i <- as(exprs(l_GBCMs[[i]]), "dgCMatrix")
  #counts_dgC_i@Dimnames[2] <- c()
  new.counts_i <- downsampleMatrix(counts_dgC_i, prop = prop_i, bycol = TRUE)
  #new.counts_i <- as(new.counts_i, "dgTMatrix")
  print(dim(new.counts_i))
  print(head(new.counts_i))
  #Count cells
  count_per_cell_i <- colSums(as.matrix(new.counts_i))
  print(length(count_per_cell_i))
  print(head(count_per_cell_i))
  #Store and convert to dgCMatrix
  l_GBCMs_reduced[[i]] <- new.counts_i
  l_count_per_cell[[i]] <- count_per_cell_i
}
names(l_count_per_cell) <- datasets
lapply(l_count_per_cell, head)
#test <- l_GBCMs[[1]]

#Replace count-matrizes in l_GBCMs
#for(i in seq(length(l_GBCMs))){
 # l_GBCMs[[i]]@assayData$exprs <- l_GBCMs_reduced[[i]]
#}

#Rename first columns
#l_cds <- lapply(l_GBCMs, function(x) {
 # feat <- fData(x)
#  names(feat) <- c('id', 'gene_short_name')
#  newCellDataSet(exprs(x),
 #                phenoData = new("AnnotatedDataFrame", data = pData(x)),
  #               featureData = new("AnnotatedDataFrame", data = feat),
   #              lowerDetectionLimit = 1,
    #             expressionFamily = negbinomial.size())
#})

#Add time_point and Cell_id column to pData
for(i in seq(length(datasets))){
  print(datasets[i])
  pData(l_cds[[i]])$Time_point <- rep(datasets[i], nrow(pData(l_cds[[i]])))
  pData(l_cds[[i]])$Cell_id <- paste(datasets[i],"_",1:nrow(pData(l_cds[[i]])),sep="")
}

#Check if appropriate
#head(pData(l_cds[[1]]))

#####Merge datasets to obtain one cde object
#build lists for pData
l_pData <- list()
#l_fData <- list()
l_expData <- list()
for(i in seq(length(datasets))){
  l_pData[[i]] <- pData(l_cds[[i]])
  l_expData[[i]] <- l_GBCMs_reduced[[i]]
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
print(head(pData(cde_all)))
tail(pData(cde_all))
#Check distribution of expression across condition
pData(cde_all)$Total_mRNAs <- Matrix::colSums(exprs(cde_all))
pData(cde_all)$Total_n_genes <- nrow(fData(cde_all)) - Matrix::colSums(exprs(cde_all)==0)
print(head(pData(cde_all)))
print(tail(pData(cde_all)))

#Estimate cutoff
hist(log2(Matrix::colSums(exprs(cde_all))))

#Check pData
hist(subset(pData(cde_all), Time_point == "d056_100")$Total_mRNAs)

#Plot detected gene number range over conditions
upper_bound <- 250
lower_bound <- 20000
qplot(Total_mRNAs, data = pData(cde_all), color = Time_point) +
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


#Plot ribo counts -> cutoff ~0.28
plot(pData(cde_all)$Total_n_genes, pData(cde_all)$ribo_expression)
#Plot cdk counts -> no cutoff
plot(pData(cde_all)$Total_n_genes, pData(cde_all)$cdk_expression)
#Plot Hsp counts -> cutoff ~0.09
plot(pData(cde_all)$Total_n_genes, pData(cde_all)$hsp_expression)
#Plot mito counts -> cutoff ~0.06
plot(pData(cde_all)$Total_n_genes, pData(cde_all)$mito_expression)

#Filter the dataset by content of 'pData'
valid_cells <- row.names(subset(pData(cde_all),
                                (Total_n_genes >= 50 &
                                   Total_n_genes <= 3750) &
                                  hsp_expression <= 0.06 &
                                  mito_expression <= 0.045 &
                                  ribo_expression <= 0.28 &
                                  exp_regression <= 5.0))
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
cde_all <- reduceDimension(cde_all, max_components = 2, num_dim = 16,
                           reduction_method = 'tSNE', verbose = T)
cde_all <- clusterCells(cde_all, num_clusters = 20)
plot_cell_clusters(cde_all, 1, 2, color = "CellType", markers = c("Pdpn", "Pecam1","Ccl19"))
#plot cell clusters according to media/condition column of pData(HSMM/cds)
plot_cell_clusters(cde_all, 1, 2, color = "Time_point")
plot_cell_clusters(cde_all, 1, 2, color = "Cluster")

#substract/regress factors
#Notes
#ribo and number of expressed genes mak no big impact on clustering
cde_all <- reduceDimension(cde_all, max_components = 2, num_dim = 10,
                           reduction_method = 'tSNE',
                           #residualModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression",
                           residualModelFormulaStr = "~Time_point",
                           verbose = T)
cde_all <- clusterCells(cde_all, num_clusters = 15)
plot_cell_clusters(cde_all, 1, 2, color = "Time_point")
plot_cell_clusters(cde_all, 1, 2, color = "Cluster", markers = c("Madcam1")) + facet_wrap(~Time_point)
plot_cell_clusters(cde_all, 1, 2, color = "Cluster") + facet_wrap(~Time_point)
plot_cell_clusters(cde_all, 1, 2, color = "Cluster", markers = c("Pecam1")) + facet_wrap(~Time_point)
plot_cell_clusters(cde_all, 1, 2, color = "Cluster")

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
Endothelial_subsets <- Endothelial.classificator[Endothelial.classificator >= 0.5]
Endothelial_subsets
SC_subsets <- Endothelial.classificator[Endothelial.classificator < 0.5]
SC_subsets

#Separate data in to Endothelial and SC
Endothelial_cells <- row.names(subset(pData(cde_all),
                                      Cluster %in% names(Endothelial_subsets)))
SC_cells <- row.names(subset(pData(cde_all),
                             Cluster %in% names(SC_subsets)))

#Generate cds object for Endothelial and SC cells
cde_SC <- cde_all[,SC_cells]
cde_Endothelial <- cde_all[,Endothelial_cells]

#####
#Building trajectories SC
#####
#Inferring the genes----------------------------------------------
diff_test_res <- differentialGeneTest(cde_SC[expressed_genes,],
                                      fullModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression",
                                      cores = 4)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
# Note: time trajectories are usefull to identify key genes
#       but not required
#Set trajectory genes in cde
cde_SC <- setOrderingFilter(cde_SC, ordering_genes)
plot_ordering_genes(cde_SC)

#Reduce dimensionality

cde_SC <- reduceDimension(cde_SC, max_components = 2, method = 'DDRTree')

#Build trajectory--------------------------------------------
cde_SC <- orderCells(cde_SC)
#plot by condition/timepoint
plot_cell_trajectory(cde_SC, color_by = "Time_point")
#plot by state
plot_cell_trajectory(cde_SC, color_by = "State")
#pass state to pData
GM_state <- function(cds){
  if (length(unique(pData(cde_SC)$State)) > 1){
    T0_counts <- table(pData(cde_SC)$State, pData(cde_SC)$Time_point)[,"d056_100"]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  }else {
    return (1)
  }
}
cde_SC <- orderCells(cde_SC, root_state = GM_state(cde_SC))
plot_cell_trajectory(cde_SC, color_by = "Pseudotime")
#facet the trajectory according to states
plot_cell_trajectory(cde_SC, color_by = "State") + facet_wrap(~Time_point, nrow = 1)

#####
#DEG analysis
#####
#!!! Not functional
cde_SC_DEGs <- differentialGeneTest(cde_SC, cores = 4)

cde_SC_cells_5 <- row.names(subset(pData(cde_SC),
                             Time_point == "d056_005"))
cde_SC_5 <- cde_all[,cde_SC_cells_5]
diff_test_res <- differentialGeneTest(cde_SC_5[expressed_genes,],
                                      fullModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression",
                                      cores = 4)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
cde_SC <- setOrderingFilter(cde_SC_5, ordering_genes)
cde_SC <- reduceDimension(cde_SC_5, max_components = 2, method = 'DDRTree')
cde_SC <- orderCells(cde_SC_5)
cde_SC_DEGs_5 <- differentialGeneTest(cde_SC_5, cores = 4)

cde_SC_cells_10 <- row.names(subset(pData(cde_SC),
                                   Time_point == "d056_010"))
cde_SC_10 <- cde_all[,cde_SC_cells_10]
