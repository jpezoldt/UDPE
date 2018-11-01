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
set.seed(123)

#####
#Notes - Important
#####
#Script only works if all datasets have the same phenotypcial data (gene id annotation files)

#####
#Global variable
#####
datasets <- c("d000") #,"d010","d024","d056")
#Number of cells in which a gene should be expressed to be included in analysis
num_cells_exp <- 50
#Export PATH
PATH_output <- "/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/monocle/D0_Signature_Clusters"

#####
#Load and compile data into cde object
#####
#Timepoint 1: e.g. neonatal/d0
dir_d0 <- "~/NAS2/pezoldt/Data/scRNAseq/270_2018-01-03_scRNA-Seq_mLN_ontogeny/Data/C11_neo_mLN"
d0 <- load_cellranger_matrix(dir_d0, genome = "mm10")

#Timepoint 2: e.g day10
#dir_d10 <- "~/NAS2/pezoldt/Data/scRNAseq/270_2018-01-03_scRNA-Seq_mLN_ontogeny/Data/C12_d10_mLN"
#d10 <- load_cellranger_matrix(dir_d10, genome = "mm10")

#Timepoint 3: e.g. day24
#dir_d24 <- "~/NAS2/pezoldt/Data/scRNAseq/270_2018-01-03_scRNA-Seq_mLN_ontogeny/Data/D1_d24_mLN"
#d24 <- load_cellranger_matrix(dir_d24, genome = "mm10")

#Timepoint 4: e.g. day56 (merge two datasets)
#dir_d56 <- "~/NAS2/pezoldt/Data/scRNAseq/244_scRNA-Seq_mLN_pLN_SPF/Data/10X_results/L1700567_mLN_SPF"
#d56 <- load_cellranger_matrix(dir_d56, genome = "mm10")
#dir_d56_2 <- "~/NAS2/pezoldt/Data/scRNAseq/252_2017-07-19_scRNA-Seq_mLN_pLN_SPF_GF/Data/10x/mLN_SPF/B6"
#d56_2 <- load_cellranger_matrix(dir_d56_2, genome = "mm10")

#Timepoint 5: e.g. day84
#Being processed

#Timepoint 6: e.g. day300
#dir_d300 <- "~/NAS2/pezoldt/Data/scRNAseq/275_2018-04-23_scRNA-Seq_mLN_SPF_45wk/Data/E7_SPF_mLN_42wk"
#d300 <- load_cellranger_matrix(dir_d300, genome = "mm10")

#list datasets
l_GBCMs <- list(d0) #,d10,d24,d56)
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
cde_all <- clusterCells(cde_all, num_clusters = 13)
plot_cell_clusters(cde_all, 1, 2, color = "CellType", markers = c("Pdpn", "Pecam1","Ccl19"))
#plot cell clusters according to media/condition column of pData(HSMM/cds)
#plot_cell_clusters(cde_all, 1, 2, color = "Time_point")
plot_cell_clusters(cde_all, 1, 2, color = "Cluster")

#substract/regress factors
#Notes
#ribo and number of expressed genes mak no big impact on clustering
cde_all <- reduceDimension(cde_all, max_components = 2, num_dim = 16,
                        reduction_method = 'tSNE',
                        #residualModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression",
                        #residualModelFormulaStr = "~Time_point",
                        verbose = T)
cde_all <- clusterCells(cde_all, num_clusters = 23)
plot_cell_clusters(cde_all, 1, 2, color = "Time_point")
plot_cell_clusters(cde_all, 1, 2, color = "Cluster", markers = c("Vcam1")) #+ facet_wrap(~Time_point)
#plot_cell_clusters(cde_all, 1, 2, color = "Cluster") + facet_wrap(~Time_point)
#plot_cell_clusters(cde_all, 1, 2, color = "Cluster", markers = c("Pecam1")) + facet_wrap(~Time_point)
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
#####
#Regress and cluster SCs
#####
#plot_pc_variance_explained(cde_SC, return_all = F)
cde_SC <- reduceDimension(cde_SC, max_components = 2, num_dim = 11,
                          reduction_method = 'tSNE',
                          residualModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression",
                          #residualModelFormulaStr = "~Time_point",
                          verbose = T)
cde_SC <- clusterCells(cde_SC, num_clusters = 15)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster")
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Cxcl13")) #+ facet_wrap(~Time_point)
#plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Cd34")) + facet_wrap(~Time_point)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Cd34","Tnfsf11","Aldh1a2","Cxcl13",
                                                                "Nfkb1","Ltbr","Icam1","Vcam1","Acta2",
                                                                "Cd248","Ackr3"),
                   cell_size = 0.5)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Acta2"),
                   cell_size = 0.5)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Cxcl13"),
                   cell_size = 0.5)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Tnfsf11"),
                   cell_size = 0.5)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Cd34"),
                   cell_size = 0.5)

#Building trajectories SC-----------------------------------------
#Reduce dimensionality
cde_SC <- reduceDimension(cde_SC, max_components = 2, method = 'DDRTree')
#Build trajectory
cde_SC <- orderCells(cde_SC)

#Inferring the genes----------------------------------------------
diff_test_res <- differentialGeneTest(cde_SC[expressed_genes,],
                                      #fullModelFormulaStr = "~Clusters",
                                      cores = 8)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
# Note: time trajectories are usefull to identify key genes
#       but not required
#Set trajectory genes in cde
cde_SC <- setOrderingFilter(cde_SC, ordering_genes)
plot_ordering_genes(cde_SC)

#plot by condition/timepoint
#plot_cell_trajectory(cde_SC, color_by = "Time_point")
#plot by state
plot_cell_trajectory(cde_SC, color_by = "State")
plot_cell_trajectory(cde_SC, markers = c("Aldh1a2"), use_color_gradient = TRUE)

markers = c("Vcam1")
#plot by cluster
plot_cell_trajectory(cde_SC, color_by = "State") + facet_wrap(~Cluster)
plot_cell_trajectory(cde_SC, color_by = "Cluster") + facet_wrap(~State)
plot_cell_trajectory(cde_SC, markers = c("Cd34")) + facet_wrap(~State)

#####
#Separate SCs into subsets
#####
#Calculate mean expression of Pecam1 ove and add to pData
#PvC gene------------------------------
PvC.genes <- c("Acta2")
#get ensemble Ids
PvC.genes.ensembl <- subset(fData(cde_SC), gene_short_name %in% PvC.genes)$id
#subset expression matrix and perform colmeans per cell
PvC.genes_expr <- as.matrix(exprs(cde_SC)[rownames(exprs(cde_all)) %in% PvC.genes.ensembl,])
pData(cde_SC)$PvC_expression <- PvC.genes_expr
#Calculate mean expression of PvC genes
# For each cluster calculate mean expression for Endothelial genes
# Store cluster IDs in vector
PvC.classificator <- c()
for(i in seq(length(levels(pData(cde_SC)$Cluster)))){
  cluster_ID_i <- levels(pData(cde_SC)$Cluster)[i]
  print(cluster_ID_i)
  mean_PvC.genes <- mean(subset(pData(cde_SC), Cluster %in% cluster_ID_i)$PvC_expression)
  print(mean_PvC.genes)
  PvC.classificator[i] <- mean_PvC.genes
}
names(PvC.classificator) <- levels(pData(cde_SC)$Cluster)

#Get all clusters with expression above threshold
# Note: Manually check average expressions
PvC_subsets <- PvC.classificator[PvC.classificator >= 20]
PvC_subsets
nonPvC_subsets <- PvC.classificator[PvC.classificator < 20]
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
setwd(PATH_output)

#Adventitial gene------------------------------
Adventi.genes <- c("Cd34")
#get ensemble Ids
Adventi.genes.ensembl <- subset(fData(cde_nonPvC), gene_short_name %in% Adventi.genes)$id
#subset expression matrix and perform colmeans per cell
Adventi.genes_expr <- as.matrix(exprs(cde_nonPvC)[rownames(exprs(cde_all)) %in% Adventi.genes.ensembl,])
pData(cde_nonPvC)$Adventi_expression <- Adventi.genes_expr
#Calculate mean expression of Adventi genes
# For each cluster calculate mean expression for Endothelial genes
# Store cluster IDs in vector
Adventi.classificator <- c()
for(i in seq(length(levels(pData(cde_nonPvC)$Cluster)))){
  cluster_ID_i <- levels(pData(cde_nonPvC)$Cluster)[i]
  print(cluster_ID_i)
  mean_Adventi.genes <- mean(subset(pData(cde_nonPvC), Cluster %in% cluster_ID_i)$Adventi_expression)
  print(mean_Adventi.genes)
  Adventi.classificator[i] <- mean_Adventi.genes
}
names(Adventi.classificator) <- levels(pData(cde_nonPvC)$Cluster)

#Get all clusters with expression above threshold
# Note: Manually check average expressions
Adventi_subsets <- Adventi.classificator[Adventi.classificator >= 1.5]
Adventi_subsets
nonAdventi_subsets <- Adventi.classificator[Adventi.classificator < 1.5]
nonAdventi_subsets

#Separate data in to different SC subsets
Adventi_cells <- row.names(subset(pData(cde_nonPvC),
                                      Cluster %in% names(Adventi_subsets)))
nonAdventi_cells <- row.names(subset(pData(cde_nonPvC),
                             Cluster %in% names(nonAdventi_subsets)))

#Generate cds object for Endothelial and SC cells
cde_nonAdventi <- cde_nonPvC[,nonAdventi_cells]
cde_Adventi <- cde_nonPvC[,Adventi_cells]

print(paste("Number of nonAdventi SC:", nrow(pData(cde_nonAdventi))))
print(paste("Number of Adventi SC:", nrow(pData(cde_Adventi))))
setwd(PATH_output)

#Adventitial separation by state------------------------------
cells_Adventi_plus_stage <- rownames(pData(cde_Adventi)[pData(cde_Adventi)$State %in% c(3), ])
cells_Adventi_minus_stage <- rownames(pData(cde_Adventi)[ !pData(cde_Adventi)$State %in% c(3), ])


#Generate cds object for Endothelial and SC cells
cde_Adventi_minus_stage <- cde_Adventi[,cells_Adventi_minus_stage]
cde_Adventi_plus_stage <- cde_Adventi[,cells_Adventi_plus_stage]

print(paste("Number of Adventi minus:", nrow(pData(cde_Adventi_minus_stage))))
print(paste("Number of Adventi plus:", nrow(pData(cde_Adventi_plus_stage))))
setwd(PATH_output)

#LTO gene------------------------------
LTO.genes <- c("Cxcl13")
#get ensemble Ids
LTO.genes.ensembl <- subset(fData(cde_nonPvC), gene_short_name %in% LTO.genes)$id
#subset expression matrix and perform colmeans per cell
LTO.genes_expr <- as.matrix(exprs(cde_nonPvC)[rownames(exprs(cde_all)) %in% LTO.genes.ensembl,])
pData(cde_nonPvC)$LTO_expression <- LTO.genes_expr
#Calculate mean expression of LTO genes
# For each cluster calculate mean expression for Endothelial genes
# Store cluster IDs in vector
LTO.classificator <- c()
for(i in seq(length(levels(pData(cde_nonPvC)$Cluster)))){
  cluster_ID_i <- levels(pData(cde_nonPvC)$Cluster)[i]
  print(cluster_ID_i)
  mean_LTO.genes <- mean(subset(pData(cde_nonPvC), Cluster %in% cluster_ID_i)$LTO_expression)
  print(mean_LTO.genes)
  LTO.classificator[i] <- mean_LTO.genes
}
names(LTO.classificator) <- levels(pData(cde_nonPvC)$Cluster)

#Get all clusters with expression above threshold
# Note: Manually check average expressions
LTO_subsets <- LTO.classificator[LTO.classificator >= 3]
LTO_subsets
nonLTO_subsets <- LTO.classificator[LTO.classificator < 3]
nonLTO_subsets

#Separate data in to different SC subsets
LTO_cells <- row.names(subset(pData(cde_nonPvC),
                                  Cluster %in% names(LTO_subsets)))
nonLTO_cells <- row.names(subset(pData(cde_nonPvC),
                                     Cluster %in% names(nonLTO_subsets)))

#Generate cds object for Endothelial and SC cells
cde_nonLTO <- cde_nonPvC[,nonLTO_cells]
cde_LTO <- cde_nonPvC[,LTO_cells]

print(paste("Number of nonLTO SC:", nrow(pData(cde_nonLTO))))
print(paste("Number of LTO SC:", nrow(pData(cde_LTO))))
setwd(PATH_output)


#####
#Regress and cluster PvC
#####
#plot_pc_variance_explained(cde_PvC, return_all = F)
cde_PvC <- reduceDimension(cde_PvC, max_components = 2, num_dim = 11,
                          reduction_method = 'tSNE',
                          residualModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression",
                          #residualModelFormulaStr = "~Time_point",
                          verbose = T)
cde_PvC <- clusterCells(cde_PvC, num_clusters = 4)
plot_cell_clusters(cde_PvC, 1, 2, color = "Cluster", show_cell_names = TRUE, cell_name_size = 10)
plot_cell_clusters(cde_PvC, 1, 2, color = "Cluster") + facet_wrap(~Time_point)
#plot_cell_clusters(cde_PvC, 1, 2, color = "Cluster", markers = c("Cd34")) + facet_wrap(~Time_point)
plot_cell_clusters(cde_PvC, 1, 2, color = "Cluster", markers = c("Cd34","Tnfsf11","Aldh1a2","Cxcl13",
                                                                "Nfkb1","Ltbr","Icam1","Vcam1","Acta2",
                                                                "Cd248","Gdf10","Cxcl9","Has1","Ackr3"),
                   cell_size = 0.5)

#Building trajectories PvC------------------------------
#Reduce dimensionality
cde_PvC <- reduceDimension(cde_PvC, max_components = 2, method = 'DDRTree')
#Build trajectory
cde_PvC <- orderCells(cde_PvC)

#Inferring the genes----------------------------------------------
diff_test_res <- differentialGeneTest(cde_PvC[expressed_genes,],
                                      #fullModelFormulaStr = "~Clusters",
                                      cores = 4)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
# Note: time trajectories are usefull to identify key genes
#       but not required
#Set trajectory genes in cde
cde_PvC <- setOrderingFilter(cde_PvC, ordering_genes)
plot_ordering_genes(cde_PvC)

#plot by condition/timepoint
plot_cell_trajectory(cde_PvC, color_by = "Time_point")
#plot by state
plot_cell_trajectory(cde_PvC, color_by = "State")
#plot by cluster
plot_cell_trajectory(cde_PvC, color_by = "State") + facet_wrap(~Cluster)
plot_cell_trajectory(cde_PvC, color_by = "Cluster") + facet_wrap(~Time_point)

#####
#Regress and cluster Adventi
#####
plot_pc_variance_explained(cde_Adventi, return_all = F)
cde_Adventi <- reduceDimension(cde_Adventi, max_components = 2, num_dim = 11,
                               reduction_method = 'tSNE',
                               residualModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression",
                               #residualModelFormulaStr = "~Time_point",
                               verbose = T)
cde_Adventi <- clusterCells(cde_Adventi, num_clusters = 9)
plot_cell_clusters(cde_Adventi, 1, 2, color = "Cluster", show_cell_names = TRUE, cell_name_size = 10)
plot_cell_clusters(cde_Adventi, 1, 2, color = "Cluster") + facet_wrap(~Time_point)
plot_cell_clusters(cde_Adventi, markers = c("Ccl19"))
#plot_cell_clusters(cde_Adventi, 1, 2, color = "Cluster", markers = c("Cd34")) + facet_wrap(~Time_point)
plot_cell_clusters(cde_Adventi, 1, 2, color = "Cluster", markers = c("Cd34","Tnfsf11","Aldh1a2","Cxcl13",
                                                                     "Nfkb1","Ltbr","Icam1","Vcam1","Acta2",
                                                                     "Cd248","Gdf10","Cxcl9","Has1","Ackr3"),cell_size = 0.5)

#Building trajectories Adventi------------------------------
#Reduce dimensionality
cde_Adventi <- reduceDimension(cde_Adventi, max_components = 2, method = 'DDRTree')
#Build trajectory
cde_Adventi <- orderCells(cde_Adventi)

#Inferring the genes----------------------------------------------
diff_test_res <- differentialGeneTest(cde_Adventi[expressed_genes,],
                                      #fullModelFormulaStr = "~Clusters",
                                      cores = 8)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
# Note: time trajectories are usefull to identify key genes
#       but not required
#Set trajectory genes in cde
cde_Adventi <- setOrderingFilter(cde_Adventi, ordering_genes)
plot_ordering_genes(cde_Adventi)

#plot by condition/timepoint
plot_cell_trajectory(cde_Adventi, color_by = "Time_point")
#plot by state
plot_cell_trajectory(cde_Adventi, color_by = "Cluster")
plot_cell_trajectory(cde_Adventi, color_by = "Pseudotime")
#plot by cluster
plot_cell_trajectory(cde_Adventi, markers = c("Aldh1a2"), use_color_gradient = TRUE)
plot_cell_trajectory(cde_Adventi, color_by = "State") + facet_wrap(~Cluster)
plot_cell_trajectory(cde_Adventi, color_by = "Cluster") + facet_wrap(~Time_point)

#pass intel of earliest time point in trajectory to state to pData
#GM_state <- function(cds){
#  if (length(unique(pData(cde_Adventi)$State)) > 1){
#    T0_counts <- table(pData(cde_Adventi)$State, pData(cde_Adventi)$Time_point)[,"d000"]
#    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
#  }else {
#    return (1)
#  }
#}
#cde_Adventi <- orderCells(cde_Adventi, root_state = GM_state(cde_Adventi))
#plot_cell_trajectory(cde_Adventi, color_by = "Pseudotime") + facet_wrap(~Time_point, nrow = 1)
#facet the trajectory according to states
#plot_cell_trajectory(cde_Adventi, color_by = "State") + facet_wrap(~Time_point, nrow = 1)

#####
#Regress and cluster Adventi Plus
#####
plot_pc_variance_explained(cde_Adventi_minus_stage, return_all = F)
cde_Adventi_minus_stage <- reduceDimension(cde_Adventi_minus_stage, max_components = 2, num_dim = 11,
                               reduction_method = 'tSNE',
                               residualModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression",
                               #residualModelFormulaStr = "~Time_point",
                               verbose = T)
cde_Adventi_minus_stage <- clusterCells(cde_Adventi_minus_stage, num_clusters = 9)
plot_cell_clusters(cde_Adventi_minus_stage, 1, 2, color = "Cluster", show_cell_names = TRUE, cell_name_size = 10)
plot_cell_clusters(cde_Adventi_minus_stage, 1, 2, color = "Cluster") + facet_wrap(~Time_point)
plot_cell_clusters(cde_Adventi_minus_stage, markers = c("Ccl19"))
#plot_cell_clusters(cde_Adventi_minus_stage, 1, 2, color = "Cluster", markers = c("Cd34")) + facet_wrap(~Time_point)
plot_cell_clusters(cde_Adventi_minus_stage, 1, 2, color = "Cluster", markers = c("Cd34","Tnfsf11","Aldh1a2","Cxcl13",
                                                                     "Nfkb1","Ltbr","Icam1","Vcam1","Acta2",
                                                                     "Cd248","Gdf10","Cxcl9","Has1","Ackr3"),
                   cell_size = 0.5)

#Building trajectories Adventi------------------------------
#Reduce dimensionality
cde_Adventi_minus_stage <- reduceDimension(cde_Adventi_minus_stage, max_components = 2, method = 'DDRTree')
#Build trajectory
cde_Adventi_minus_stage <- orderCells(cde_Adventi_minus_stage)

#Inferring the genes----------------------------------------------
diff_test_res <- differentialGeneTest(cde_Adventi_minus_stage[expressed_genes,],
                                      #fullModelFormulaStr = "~Clusters",
                                      cores = 8)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
# Note: time trajectories are usefull to identify key genes
#       but not required
#Set trajectory genes in cde
cde_Adventi_minus_stage <- setOrderingFilter(cde_Adventi_minus_stage, ordering_genes)
plot_ordering_genes(cde_Adventi_minus_stage)

#plot by condition/timepoint
plot_cell_trajectory(cde_Adventi_minus_stage, color_by = "Time_point")
#plot by state
plot_cell_trajectory(cde_Adventi_minus_stage, color_by = "Cluster")
plot_cell_trajectory(cde_Adventi_minus_stage, color_by = "Pseudotime")
#plot by cluster
plot_cell_trajectory(cde_Adventi_minus_stage, markers = c("Aldh1a2"), use_color_gradient = TRUE)
plot_cell_trajectory(cde_Adventi_minus_stage, color_by = "State") + facet_wrap(~Cluster)
plot_cell_trajectory(cde_Adventi_minus_stage, color_by = "Cluster") + facet_wrap(~Time_point)

#####
#Regress and cluster Adventi Minus
#####
plot_pc_variance_explained(cde_Adventi_minus_stage, return_all = F)
cde_Adventi_minus_stage <- reduceDimension(cde_Adventi_minus_stage, max_components = 2, num_dim = 11,
                                      reduction_method = 'tSNE',
                                      residualModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression",
                                      #residualModelFormulaStr = "~Time_point",
                                      verbose = T)
cde_Adventi_minus_stage <- clusterCells(cde_Adventi_minus_stage, num_clusters = 9)
plot_cell_clusters(cde_Adventi_minus_stage, 1, 2, color = "Cluster", show_cell_names = TRUE, cell_name_size = 10)
plot_cell_clusters(cde_Adventi_minus_stage, 1, 2, color = "Cluster") + facet_wrap(~Time_point)
plot_cell_clusters(cde_Adventi_minus_stage, markers = c("Ccl19"))
#plot_cell_clusters(cde_Adventi_minus_stage, 1, 2, color = "Cluster", markers = c("Cd34")) + facet_wrap(~Time_point)
plot_cell_clusters(cde_Adventi_minus_stage, 1, 2, color = "Cluster", markers = c("Cd34","Tnfsf11","Aldh1a2","Cxcl13",
                                                                            "Nfkb1","Ltbr","Icam1","Vcam1","Acta2",
                                                                            "Cd248","Gdf10","Cxcl9","Has1","Ackr3"),
                   cell_size = 0.5)

#Building trajectories Adventi------------------------------
#Reduce dimensionality
cde_Adventi_minus_stage <- reduceDimension(cde_Adventi_minus_stage, max_components = 2, method = 'DDRTree')
#Build trajectory
cde_Adventi_minus_stage <- orderCells(cde_Adventi_minus_stage)

#Inferring the genes----------------------------------------------
diff_test_res <- differentialGeneTest(cde_Adventi_minus_stage[expressed_genes,],
                                      #fullModelFormulaStr = "~Clusters",
                                      cores = 8)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
# Note: time trajectories are usefull to identify key genes
#       but not required
#Set trajectory genes in cde
cde_Adventi_minus_stage <- setOrderingFilter(cde_Adventi_minus_stage, ordering_genes)
plot_ordering_genes(cde_Adventi_minus_stage)

#plot by condition/timepoint
plot_cell_trajectory(cde_Adventi_minus_stage, color_by = "Time_point")
#plot by state
plot_cell_trajectory(cde_Adventi_minus_stage, color_by = "Cluster")
plot_cell_trajectory(cde_Adventi_minus_stage, color_by = "Pseudotime")
#plot by cluster
plot_cell_trajectory(cde_Adventi_minus_stage, markers = c("Cd34"), use_color_gradient = TRUE)
plot_cell_trajectory(cde_Adventi_minus_stage, color_by = "State") + facet_wrap(~Cluster)
plot_cell_trajectory(cde_Adventi_minus_stage, color_by = "Cluster") + facet_wrap(~Time_point)




#####
#Regress and cluster nonAdventi
#####
#plot_pc_variance_explained(cde_nonAdventi, return_all = F)
cde_nonAdventi <- reduceDimension(cde_nonAdventi, max_components = 2, num_dim = 11,
                               reduction_method = 'tSNE',
                               residualModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression",
                               #residualModelFormulaStr = "~Time_point",
                               verbose = T)
cde_nonAdventi <- clusterCells(cde_nonAdventi, num_clusters = 8)
plot_cell_clusters(cde_nonAdventi, 1, 2, color = "Cluster", show_cell_names = TRUE, cell_name_size = 10)
plot_cell_clusters(cde_nonAdventi, 1, 2, color = "Cluster") + facet_wrap(~Time_point)
#plot_cell_clusters(cde_nonAdventi, 1, 2, color = "Cluster", markers = c("Cd34")) + facet_wrap(~Time_point)
plot_cell_clusters(cde_nonAdventi, 1, 2, color = "Cluster", markers = c("Cd34","Tnfsf11","Aldh1a2","Cxcl13",
                                                                     "Nfkb1","Ltbr","Icam1","Vcam1","Acta2",
                                                                     "Cd248","Gdf10","Cxcl9","Has1","Ackr3"),
                   cell_size = 0.5)

#Building trajectories nonAdventi------------------------------
#Reduce dimensionality
cde_nonAdventi <- reduceDimension(cde_nonAdventi, max_components = 2, method = 'DDRTree')
#Build trajectory
cde_nonAdventi <- orderCells(cde_nonAdventi)

#Inferring the genes----------------------------------------------
diff_test_res <- differentialGeneTest(cde_nonAdventi[expressed_genes,],
                                      #fullModelFormulaStr = "~Clusters",
                                      cores = 4)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
# Note: time trajectories are usefull to identify key genes
#       but not required
#Set trajectory genes in cde
cde_nonAdventi <- setOrderingFilter(cde_nonAdventi, ordering_genes)
plot_ordering_genes(cde_nonAdventi)

#plot by condition/timepoint
plot_cell_trajectory(cde_nonAdventi, color_by = "Time_point")
#plot by state
plot_cell_trajectory(cde_nonAdventi, color_by = "State")
plot_cell_trajectory(cde_nonAdventi, color_by = "Pseudotime")
#plot by cluster
plot_cell_trajectory(cde_nonAdventi, color_by = "State") + facet_wrap(~Cluster)
plot_cell_trajectory(cde_nonAdventi, color_by = "Cluster") + facet_wrap(~Time_point)

#####
#Regress and cluster Endothelial
#####
#plot_pc_variance_explained(cde_nonAdventi, return_all = F)
cde_Endothelial <- reduceDimension(cde_Endothelial, max_components = 2, num_dim = 11,
                                  reduction_method = 'tSNE',
                                  residualModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression",
                                  #residualModelFormulaStr = "~Time_point",
                                  verbose = T)
cde_Endothelial <- clusterCells(cde_Endothelial, num_clusters = 8)
plot_cell_clusters(cde_Endothelial, 1, 2, color = "Cluster", show_cell_names = TRUE, cell_name_size = 10)
plot_cell_clusters(cde_Endothelial, 1, 2, color = "Cluster") + facet_wrap(~Time_point)
#plot_cell_clusters(cde_Endothelial, 1, 2, color = "Cluster", markers = c("Cd34")) + facet_wrap(~Time_point)
plot_cell_clusters(cde_Endothelial, 1, 2, color = "Cluster", markers = c("Cd34","Tnfsf11","Aldh1a2","Cxcl13",
                                                                        "Nfkb1","Ltbr","Icam1","Vcam1","Acta2",
                                                                        "Cd248","Gdf10","Cxcl9","Has1","Ackr3"),
                   cell_size = 0.5)

#Building trajectories Endothelial------------------------------
#Reduce dimensionality
cde_Endothelial <- reduceDimension(cde_Endothelial, max_components = 2, method = 'DDRTree')
#Build trajectory
cde_Endothelial <- orderCells(cde_Endothelial)

#Inferring the genes----------------------------------------------
diff_test_res <- differentialGeneTest(cde_Endothelial[expressed_genes,],
                                      #fullModelFormulaStr = "~Clusters",
                                      cores = 4)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
# Note: time trajectories are usefull to identify key genes
#       but not required
#Set trajectory genes in cde
cde_Endothelial <- setOrderingFilter(cde_Endothelial, ordering_genes)
plot_ordering_genes(cde_Endothelial)

#plot by condition/timepoint
plot_cell_trajectory(cde_Endothelial, color_by = "Time_point")
#plot by state
plot_cell_trajectory(cde_Endothelial, color_by = "State")
#plot by cluster
plot_cell_trajectory(cde_Endothelial, color_by = "State") + facet_wrap(~Cluster)
plot_cell_trajectory(cde_Endothelial, color_by = "Cluster") + facet_wrap(~Time_point)



#####
#Regress and cluster nonLTO
#####
#plot_pc_variance_explained(cde_nonAdventi, return_all = F)
cde_nonLTO <- reduceDimension(cde_nonLTO, max_components = 2, num_dim = 11,
                                   reduction_method = 'tSNE',
                                   residualModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression",
                                   #residualModelFormulaStr = "~Time_point",
                                   verbose = T)
cde_nonLTO <- clusterCells(cde_nonLTO, num_clusters = 8)
plot_cell_clusters(cde_nonLTO, 1, 2, color = "Cluster", show_cell_names = TRUE, cell_name_size = 10)
plot_cell_clusters(cde_nonLTO, 1, 2, color = "Cluster") + facet_wrap(~Time_point)
#plot_cell_clusters(cde_nonLTO, 1, 2, color = "Cluster", markers = c("Cd34")) + facet_wrap(~Time_point)
plot_cell_clusters(cde_nonLTO, 1, 2, color = "Cluster", markers = c("Cd34","Tnfsf11","Aldh1a2","Cxcl13",
                                                                         "Nfkb1","Ltbr","Icam1","Vcam1","Acta2",
                                                                         "Cd248","Gdf10","Cxcl9","Has1","Ackr3"),cell_size = 0.5)
plot_cell_clusters(cde_nonLTO, 1, 2, color = "Cluster", markers = c("Ackr3"),cell_size = 1.0)
plot_cell_clusters(cde_nonLTO, 1, 2, color = "Cluster", markers = c("Aldh1a2"),cell_size = 1.0)
plot_cell_clusters(cde_nonLTO, 1, 2, color = "Cluster", markers = c("Cd34"),cell_size = 1.0)
plot_cell_clusters(cde_nonLTO, 1, 2, color = "Cluster", markers = c("Il6"),cell_size = 1.0)

#Building trajectories nonLTO------------------------------
#Reduce dimensionality
cde_nonLTO <- reduceDimension(cde_nonLTO, max_components = 2, method = 'DDRTree')
#Build trajectory
cde_nonLTO <- orderCells(cde_nonLTO)

#Inferring the genes----------------------------------------------
diff_test_res <- differentialGeneTest(cde_nonLTO[expressed_genes,],
                                      #fullModelFormulaStr = "~Clusters",
                                      cores = 4)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
# Note: time trajectories are usefull to identify key genes
#       but not required
#Set trajectory genes in cde
cde_nonLTO <- setOrderingFilter(cde_nonLTO, ordering_genes)
plot_ordering_genes(cde_nonLTO)

#plot by condition/timepoint
plot_cell_trajectory(cde_nonLTO, color_by = "Time_point")
#plot by state
plot_cell_trajectory(cde_nonLTO, color_by = "State")
#plot by cluster
plot_cell_trajectory(cde_nonLTO, color_by = "State") + facet_wrap(~Cluster)
plot_cell_trajectory(cde_nonLTO, color_by = "Cluster") + facet_wrap(~Time_point)
#####
#Regress and cluster nonLTO
#####
#plot_pc_variance_explained(cde_nonAdventi, return_all = F)
cde_nonLTO <- reduceDimension(cde_nonLTO, max_components = 2, num_dim = 11,
                           reduction_method = 'tSNE',
                           residualModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression",
                           #residualModelFormulaStr = "~Time_point",
                           verbose = T)
cde_nonLTO <- clusterCells(cde_nonLTO, num_clusters = 8)
plot_cell_clusters(cde_nonLTO, 1, 2, color = "Cluster", show_cell_names = TRUE, cell_name_size = 10)
plot_cell_clusters(cde_nonLTO, 1, 2, color = "Cluster")
#plot_cell_clusters(cde_nonLTO, 1, 2, color = "Cluster", markers = c("Cd34")) + facet_wrap(~Time_point)
plot_cell_clusters(cde_nonLTO, 1, 2, color = "Cluster", markers = c("Cd34","Tnfsf11","Aldh1a2","Cxcl13",
                                                                 "Nfkb1","Ltbr","Icam1","Vcam1","Acta2",
                                                                 "Cd248","Gdf10","Cxcl9","Has1","Ackr3"),
                   cell_size = 0.5)

#Building trajectories nonLTO------------------------------
#Reduce dimensionality
cde_nonLTO <- reduceDimension(cde_nonLTO, max_components = 2, method = 'DDRTree')
#Build trajectory
cde_nonLTO <- orderCells(cde_nonLTO)

#Inferring the genes----------------------------------------------
diff_test_res <- differentialGeneTest(cde_nonLTO[expressed_genes,],
                                      #fullModelFormulaStr = "~Clusters",
                                      cores = 4)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
# Note: time trajectories are usefull to identify key genes
#       but not required
#Set trajectory genes in cde
cde_nonLTO <- setOrderingFilter(cde_nonLTO, ordering_genes)
plot_ordering_genes(cde_nonLTO)

#plot by condition/timepoint
plot_cell_trajectory(cde_nonLTO, color_by = "Time_point")
#plot by state
plot_cell_trajectory(cde_nonLTO, color_by = "State")
#plot by cluster
plot_cell_trajectory(cde_nonLTO, color_by = "State") + facet_wrap(~Cluster)
plot_cell_trajectory(cde_nonLTO, color_by = "Cluster") + facet_wrap(~Time_point)




#####
#Save/Load R object
#####
#saveRDS(cde_SC, file="~/NAS2/pezoldt/Analysis/scRNAseq/monocle/20180826_monocle_ontogeny.Rds") 
cde_SC <- readRDS(file="~/NAS2/pezoldt/Analysis/scRNAseq/monocle/20180826_d0_SC.Rds")
saveRDS(cde_Adventi, file="~/NAS2/pezoldt/Analysis/scRNAseq/monocle/20180915_d0_Adventi.Rds")
saveRDS(cde_nonAdventi, file="~/NAS2/pezoldt/Analysis/scRNAseq/monocle/20180826_d0_nonAdventi.Rds")
saveRDS(cde_PvC, file="~/NAS2/pezoldt/Analysis/scRNAseq/monocle/20180826_d0_PvC.Rds")
saveRDS(cde_Endothelial, file="~/NAS2/pezoldt/Analysis/scRNAseq/monocle/20180826_d0_Endothelial.Rds")
saveRDS(cde_Adventi_plus_stage, file="~/NAS2/pezoldt/Analysis/scRNAseq/monocle/20180826_d0_Adventi_plus_stage.Rds")
saveRDS(cde_nonLTO, file="~/NAS2/pezoldt/Analysis/scRNAseq/monocle/20181027_d0_nonLTO.Rds")
saveRDS(cde_LTO, file="~/NAS2/pezoldt/Analysis/scRNAseq/monocle/20181027_d0_LTO.Rds")
cde_nonAdventi <- readRDS(file="~/NAS2/pezoldt/Analysis/scRNAseq/monocle/20180826_d0_d10_nonAdventi.Rds")
cde_Adventi <- readRDS(file="~/NAS2/pezoldt/Analysis/scRNAseq/monocle/20180915_d0_Adventi.Rds")
cde_nonLTO <- readRDS(file="~/NAS2/pezoldt/Analysis/scRNAseq/monocle/20181027_d0_nonLTO.Rds")
cde_SC
#####
#Plot Genes expression
#####
#set cde of interest
cde_x <- cde_nonLTO
#Jitter-----------------------------
blast_genes <- row.names(subset(fData(cde_x),
                                gene_short_name %in% c("Nkx2-3", "Ccl9", "Cd34","Gdf10","Ackr3")))
plot_genes_jitter(cde_SC[blast_genes,],
                  grouping = "Time_point",
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
                                              fullModelFormulaStr = '~Cluster',
                                              cores = 8)
head(cde_x_cluster_DEG)

#Take the top 1000 DEGs according to qval
cde_x_cluster_DEG_genes <- row.names(cde_x_cluster_DEG)[order(cde_x_cluster_DEG$qval)][1:1000]
cde_x_cluster_DEG_top <- cde_x_cluster_DEG[order(cde_x_cluster_DEG$qval),]
head(cde_x_cluster_DEG_top,50)
#Use top DEGs as ordering for trajectories
cde_x <- setOrderingFilter(cde_x,
                    ordering_genes = cde_x_cluster_DEG_genes)
#Perfomr RandomForest
cde_x <- reduceDimension(cde_x, method = 'DDRTree')
#Order cells and change the root state
cde_x <- orderCells(cde_x)
cde_x <- orderCells(cde_x, root_state = 5)
#cde_SC <- orderCells(cde_SC, root_state = GM_state(cde_SC))
#
plot_cell_trajectory(cde_x)
plot_cell_trajectory(cde_x, color_by = "State") + facet_wrap(~Cluster)
plot_cell_trajectory(cde_x, markers = "Aldh1a2")  + facet_wrap(~Cluster)
plot_cell_trajectory(cde_x, markers = "Cdk1", use_color_gradient = TRUE)

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
saveRDS(cde_x, paste(PATH_output,"/","cde_nonLTO_right_direction.rds",sep=""))
cde_x <- readRDS(paste(PATH_output,"/","cde_nonLTO_right_direction.rds",sep=""))
BEAM_res <- BEAM(cde_x, branch_point = 5, cores = 12)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]


#Perform BEAM for all branches
#List to store BEAMs
l_BEAM_res <- list()
for(i in 1:2){
  BEAM_res_i <- BEAM(cde_x, branch_point = i, cores =12)
  l_BEAM_res[[i]] <- BEAM_res_i
}
names(l_BEAM_res) <- paste("Branch", c(1:2), sep = "_")
saveRDS(l_BEAM_res, paste(PATH_output,"/","BEAM_D0_res_all.rds",sep=""))
l_BEAM_res <- readRDS(paste(PATH_output,"/","BEAM_D0_res_all.rds",sep=""))
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



