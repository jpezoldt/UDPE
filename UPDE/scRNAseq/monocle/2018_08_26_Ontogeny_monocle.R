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

#####
#Notes - Important
#####
#Script only works if all datasets have the same phenotypcial data (gene id annotation files)

#####
#Global variable
#####
datasets <- c("d000","d010","d024","d056","d300")
#Number of cells in which a gene should be expressed to be included in analysis
num_cells_exp <- 50
#Export PATH
PATH_output <- "/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/monocle/Output"

#####
#Load and compile data into cde object
#####
#Timepoint 1: e.g. neonatal/d0
dir_d0 <- "~/NAS2/pezoldt/Data/scRNAseq/270_2018-01-03_scRNA-Seq_mLN_ontogeny/Data/C11_neo_mLN"
d0 <- load_cellranger_matrix(dir_d0, genome = "mm10")

#Timepoint 2: e.g day10
dir_d10 <- "~/NAS2/pezoldt/Data/scRNAseq/270_2018-01-03_scRNA-Seq_mLN_ontogeny/Data/C12_d10_mLN"
d10 <- load_cellranger_matrix(dir_d10, genome = "mm10")

#Timepoint 3: e.g. day24
dir_d24 <- "~/NAS2/pezoldt/Data/scRNAseq/270_2018-01-03_scRNA-Seq_mLN_ontogeny/Data/D1_d24_mLN"
d24 <- load_cellranger_matrix(dir_d24, genome = "mm10")

#Timepoint 4: e.g. day56 (merge two datasets)
dir_d56 <- "~/NAS2/pezoldt/Data/scRNAseq/244_scRNA-Seq_mLN_pLN_SPF/Data/10X_results/L1700567_mLN_SPF"
d56 <- load_cellranger_matrix(dir_d56, genome = "mm10")
#dir_d56_2 <- "~/NAS2/pezoldt/Data/scRNAseq/252_2017-07-19_scRNA-Seq_mLN_pLN_SPF_GF/Data/10x/mLN_SPF/B6"
#d56_2 <- load_cellranger_matrix(dir_d56_2, genome = "mm10")

#Timepoint 5: e.g. day84
#Being processed

#Timepoint 6: e.g. day300
dir_d300 <- "~/NAS2/pezoldt/Data/scRNAseq/275_2018-04-23_scRNA-Seq_mLN_SPF_45wk/Data/E7_SPF_mLN_42wk"
d300 <- load_cellranger_matrix(dir_d300, genome = "mm10")

#list datasets
l_GBCMs <- list(d0,d10,d24,d56,d300)
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
d10 <- subset(pData(cde_all), Time_point == "d010")
d24 <- subset(pData(cde_all), Time_point == "d024")
#Plot ribo counts -> cutoff ~0.28
plot(pData(cde_all)$Total_n_genes, pData(cde_all)$ribo_expression, type = "p", pch = 4)
plot(d0$Total_n_genes, d0$ribo_expression)
plot(d10$Total_n_genes, d10$ribo_expression)
plot(d24$Total_n_genes, d24$ribo_expression)
#Plot cdk counts -> no cutoff
plot(pData(cde_all)$Total_n_genes, pData(cde_all)$cdk_expression)
plot(d0$Total_n_genes, d0$cdk_expression)
plot(d10$Total_n_genes, d10$cdk_expression)
plot(d24$Total_n_genes, d24$cdk_expression)
#Plot Hsp counts -> cutoff ~0.09
plot(pData(cde_all)$Total_n_genes, pData(cde_all)$hsp_expression)
plot(d0$Total_n_genes, d0$hsp_expression)
plot(d10$Total_n_genes, d10$hsp_expression)
plot(d24$Total_n_genes, d24$hsp_expression)
#Plot mito counts -> cutoff ~0.06
plot(pData(cde_all)$Total_n_genes, pData(cde_all)$mito_expression)
plot(d0$Total_n_genes, d0$mito_expression)
plot(d10$Total_n_genes, d10$mito_expression)
plot(d24$Total_n_genes, d24$mito_expression)


#Filter the dataset by content of 'pData'
valid_cells <- row.names(subset(pData(cde_all),
                                (Total_n_genes >= 250 &
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
cde_all <- reduceDimension(cde_all, max_components = 2, num_dim = 16,
                        reduction_method = 'tSNE',
                        #residualModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression",
                        residualModelFormulaStr = "~Time_point",
                        verbose = T)
cde_all <- clusterCells(cde_all, num_clusters = 23)
plot_cell_clusters(cde_all, 1, 2, color = "Time_point")
plot_cell_clusters(cde_all, 1, 2, color = "Cluster", markers = c("Ackr4")) + facet_wrap(~Time_point)
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
#Analyze Endothelial cell
###############
#print(paste("Number of non-endothelial SC:", nrow(pData(Endothelial_cells))))
#####
#Regress and cluster
#####
#plot_pc_variance_explained(cde_Endothelial, return_all = F)
#cde_Endothelial <- reduceDimension(cde_Endothelial, max_components = 2, num_dim = 11,
 #                         reduction_method = 'tSNE',
#                          residualModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression",
 #                         #residualModelFormulaStr = "~Time_point",
  #                        verbose = T)
#cde_Endothelial <- clusterCells(cde_Endothelial, num_clusters = 13)
#plot_cell_clusters(cde_Endothelial, 1, 2, color = "Time_point")
#plot_cell_clusters(cde_Endothelial, 1, 2, color = "Cluster") + facet_wrap(~Time_point)
#plot_cell_clusters(cde_Endothelial, 1, 2, color = "Cluster", markers = c("Pecam1")) + facet_wrap(~Time_point)
#plot_cell_clusters(cde_Endothelial, 1, 2, color = "Cluster")

###############
#Analyze SC cell
###############
print(paste("Number of non-endothelial SC:", nrow(pData(cde_SC))))
#####
#Regress and cluster
#####
plot_pc_variance_explained(cde_SC, return_all = F)
cde_SC <- reduceDimension(cde_SC, max_components = 2, num_dim = 11,
                           reduction_method = 'tSNE',
                           residualModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression",
                           #residualModelFormulaStr = "~Time_point",
                           verbose = T)
cde_SC <- clusterCells(cde_SC, num_clusters = 18)
plot_cell_clusters(cde_SC, 1, 2, color = "Time_point")
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster") + facet_wrap(~Time_point)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Cd34")) + facet_wrap(~Time_point)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster")

#####
#Building trajectories SC
#####
#Inferring the genes----------------------------------------------
diff_test_res <- differentialGeneTest(cde_SC[expressed_genes,],
                                      fullModelFormulaStr = "~Time_point",
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
#plot by cluster
plot_cell_trajectory(cde_SC, color_by = "Cluster")

#pass intel of earliest time point in trajectory to state to pData
GM_state <- function(cds){
  if (length(unique(pData(cde_SC)$State)) > 1){
    T0_counts <- table(pData(cde_SC)$State, pData(cde_SC)$Time_point)[,"d000"]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  }else {
    return (1)
  }
}
cde_SC <- orderCells(cde_SC, root_state = GM_state(cde_SC))
plot_cell_trajectory(cde_SC, color_by = "Pseudotime") + facet_wrap(~Time_point, nrow = 1)
#facet the trajectory according to states
plot_cell_trajectory(cde_SC, color_by = "State") + facet_wrap(~Time_point, nrow = 1)

#####
#Save/Load R object
#####
#saveRDS(cde_SC, file="~/NAS2/pezoldt/Analysis/scRNAseq/monocle/20180826_monocle_ontogeny.Rds") 
cde_SC <- readRDS(file="~/NAS2/pezoldt/Analysis/scRNAseq/monocle/20180815_monocle_ontogeny.Rds")

#####
#Plot Genes expression
#####
#Jitter-----------------------------
blast_genes <- row.names(subset(fData(cde_SC),
                                gene_short_name %in% c("Nkx2-3", "Ccl9", "Cd34")))
plot_genes_jitter(cde_SC[blast_genes,],
                  grouping = "Time_point",
                  min_expr = 0.05)

#Plot gene expression in pseudotime
cde_SC_expressed_genes <-  row.names(subset(fData(cde_SC),
                                          num_cells_expressed >= 10))
cde_SC_filtered <- cde_SC[cde_SC_expressed_genes,]
my_genes <- row.names(subset(fData(cde_SC_filtered),
                             gene_short_name %in% c("Cxcl1","Cd34","Ccl19")))
cde_SC_subset <- cde_SC_filtered[my_genes,]
plot_genes_in_pseudotime(cde_SC_subset, color_by = "Time_point")

#####
#Ordering based on genes taht difer between clusters
#####
# Pick genes that are expressed in >5% of the cell
cde_SC <- detectGenes(cde_SC, min_expr = 0.1)
fData(cde_SC)$use_for_ordering <-
  fData(cde_SC)$num_cells_expressed > 0.05 * ncol(cde_SC)
#Plot variance explainedby picked genes
plot_pc_variance_explained(cde_SC, return_all = F)
#Reduce dimensions using genes
cde_SC <- reduceDimension(cde_SC,
                            max_components = 2,
                            norm_method = 'log',
                            num_dim = 5,
                            reduction_method = 'tSNE',
                            verbose = T)

#cluster cells
cde_SC <- clusterCells(cde_SC, num_clusters = 13)
plot_cell_clusters(cde_SC, color_by = 'as.factor(Cluster)')
plot_cell_clusters(cde_SC, color_by = 'as.factor(Time_point)')

#Decision plot to define Rho and P
plot_rho_delta(cde_SC, rho_threshold = 2, delta_threshold = 4 )

#Re-cluster cells implementing cut-off
cde_SC <- clusterCells(cde_SC,
                         rho_threshold = 70,
                         delta_threshold = 4,
                         skip_rho_sigma = T,
                         verbose = F)

#Differential gene expression across clusters
clustering_DEG_genes <- differentialGeneTest(cde_SC[cde_SC_expressed_genes,],
                                              fullModelFormulaStr = '~Cluster',
                                              cores = 4)

#Print number of DEGs


#Take the top 1000 DEGs according to qval
cde_SC_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
#Use top DEGs as ordering for trajectories
cde_SC <- setOrderingFilter(cde_SC,
                    ordering_genes = cde_SC_ordering_genes)
#Perfomr RF
cde_SC <- reduceDimension(cde_SC, method = 'DDRTree')
#Order cells
#cde_SC <- orderCells(cde_SC)
cde_SC <- orderCells(cde_SC, root_state = GM_state(cde_SC))
#
plot_cell_trajectory(cde_SC, color_by = "State") + facet_wrap(~Time_point, nrow = 1)


#####
#DEG analysis
#####
#Computations expensive
diff_test_res <- differentialGeneTest(cde_SC, fullModelFormulaStr = "~Time_point",
                                      cores = 4)

# Select genes that are significant at an FDR < 10%
sig_genes <- subset(diff_test_res, qval < 0.1)

head(sig_genes[,c("gene_short_name", "pval", "qval")])

#####
#Distinguish Cluster or State
#####
#Calculation DEGs
# Determine how good a gene explains the state of a cell by checking/building a model
# 1) that knows about the cluster annotation
# 2) that does not know about the cluster annotation
#Compare output of models and choose highest scoring models
#DEG according to cluster
diff_test_res_cluster <- differentialGeneTest(cde_SC,
                                      fullModelFormulaStr = "~CellType")
head(diff_test_res_cluster[,c("gene_short_name", "pval", "qval")])

#DEG according to state
diff_test_res_state <- differentialGeneTest(cde_SC,
                                              fullModelFormulaStr = "~CellType")
sig_gene_names_pseudotime <- rownames(subset(diff_test_res_pseudotime, qval <= 0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001))

plot_pseudotime_heatmap(cde_SC[sig_gene_names_pseudotime,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T,
                        trend_formula = "~Time_point")
 
#####
#Multifactorial DEG analyis
#####







