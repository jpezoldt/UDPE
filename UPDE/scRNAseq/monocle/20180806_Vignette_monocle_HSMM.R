#Implement vignette Monocle with HSMM dataset

#Libraries
library(monocle)
library(pheatmap)
library(reshape2)
library(HSMMSingleCell)


#Load data
HSMM <- load_HSMM()
HSMM_maker <- load_HSMM_markers()


HSMM <- newCellDataSet(as(HSMM@assayData$exprs, "sparseMatrix"),
                       phenoData = HSMM@phenoData,
                       featureData = HSMM@featureData,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())

#Estimate factor matrices
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

#Filter low qualitiy cells
HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))
nrow(fData(HSMM))

#Number of expressed genes 
# Note: Thresh by number of cells gene is expressed in
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
print(head(pData(HSMM)))

#Filter the dataset by content of 'pData'
#valid_cells <- row.names(subset(pData(HSMM),
 #                                 Cells.in.Well == 1 &
  #                                Control == FALSE &
   #                               Clump == FALSE &
    #                              Debris == FALSE &
     #                             Mapped.Fragments > 1000000))
#HSMM <- HSMM[,valid_cells]

#Check distribution of expression across condition
pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))

HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]

upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) +
                     2*sd(log10(pData(HSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) -
                     2*sd(log10(pData(HSMM)$Total_mRNAs)))

qplot(Total_mRNAs, data = pData(HSMM), color = Hours, geom = "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)
levels(as.factor(pData(HSMM)$Hour))

#Eliminate cells by pct thres
HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound &
               pData(HSMM)$Total_mRNAs < upper_bound]
HSMM <- detectGenes(HSMM, min_expr = 0.1)

#####Check for lognormal distribution
# Log-transform each value in the expression matrix.
L <- log(exprs(HSMM[expressed_genes,]))
# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
# Plot the distribution of the standardized gene expression values.
qplot(value, geom = "density", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")

######
#Classify cells by type
######
#1) by single gene
MYF5_id <- row.names(subset(fData(HSMM), gene_short_name == "MYF5"))
ANPEP_id <- row.names(subset(fData(HSMM), gene_short_name == "ANPEP"))
cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "Myoblast", classify_func = function(x) { x[MYF5_id,] >= 1 })
cth <- addCellType(cth, "Fibroblast", classify_func = function(x) { x[MYF5_id,] < 1 & x[ANPEP_id,] > 1 })
#Cell definition are classified over functions here stored at cth
HSMM <- classifyCells(HSMM, cth, 0.1)
table(pData(HSMM)$CellType)
pie <- ggplot(pData(HSMM), aes(x = factor(1), fill = factor(CellType))) +
  geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

#2) Implement accession of predefined cell subsets from Seurat via pData
#Identify variable genes
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
#black dots are used for ordering
plot_ordering_genes(HSMM)
nrow(unsup_clustering_genes)
#plot PCs explaining variance
plot_pc_variance_explained(HSMM, return_all = F) # norm_method = 'log',
#Count how many dimensions to use
HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 6,
                        reduction_method = 'tSNE', verbose = T)
HSMM <- clusterCells(HSMM, num_clusters = 2)
plot_cell_clusters(HSMM, 1, 2, color = "CellType", markers = c("MYF5", "ANPEP"))
#plot cell clusters according to media/condition column of pData(HSMM/cds)
plot_cell_clusters(HSMM, 1, 2, color = "Media")
#substract/regress factors         
HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 2,
                        reduction_method = 'tSNE',
                        residualModelFormulaStr = "~Media + num_genes_expressed",
                        verbose = T)
HSMM <- clusterCells(HSMM, num_clusters = 2)
plot_cell_clusters(HSMM, 1, 2, color = "CellType")
#redo unsupervised clustering
HSMM <- clusterCells(HSMM, num_clusters = 2)
plot_cell_clusters(HSMM, 1, 2, color = "Cluster") + facet_wrap(~CellType)

#3) Clustering cells using marker genes
#calculate DEGs
marker_diff <- markerDiffTable(HSMM[expressed_genes,],
                               cth,
                               residualModelFormulaStr = "~Media + num_genes_expressed",
                               cores = 1)
#Pick top ones
candidate_clustering_genes <- row.names(subset(marker_diff, qval < 0.01))
marker_spec <- calculateMarkerSpecificity(HSMM[candidate_clustering_genes,], cth)
#column 'specificity' infers correctnes of marker per cell type 1 is highest
head(selectTopMarkers(marker_spec, 3))
#Use top 500
semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 500)$gene_id)
HSMM <- setOrderingFilter(HSMM, semisup_clustering_genes)
#Show the top 500 genes
# Note: This analysis is honed to pinpoint the top genes that distungish fibroblasts and myocytes
plot_ordering_genes(HSMM)
plot_pc_variance_explained(HSMM, return_all = F)
HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 3, norm_method = 'log',
                        reduction_method = 'tSNE',
                        residualModelFormulaStr = "~Media + num_genes_expressed",
                        verbose = T)
HSMM <- clusterCells(HSMM, num_clusters = 2)

plot_cell_clusters(HSMM, 1, 2, color = "CellType")

#4) Imputing unknown cell types
HSMM <- clusterCells(HSMM,
                     num_clusters = 5,
                     frequency_thresh = 0.1,
                     cell_type_hierarchy = cth)
plot_cell_clusters(HSMM, 1, 2, color = "CellType", markers = c("MYF5", "ANPEP"))
pie <- ggplot(pData(HSMM), aes(x = factor(1), fill = factor(CellType))) +
  geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
HSMM_myo <- HSMM[,pData(HSMM)$CellType == "Myoblast"]
HSMM_myo <- estimateDispersions(HSMM_myo)

################
#Building trajectories
################
#####
#Inferring the genes
#####
diff_test_res <- differentialGeneTest(HSMM_myo[expressed_genes,],
                                      fullModelFormulaStr = "~Media")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
# Note: time trajectories are usefull to identify key genes
#       but not required
#Set trajectory genes in cde
HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)
plot_ordering_genes(HSMM_myo)
#####
#Reduce dimensionality
#####
HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2, method = 'DDRTree')

#####
#Build trajectory
#####
HSMM_myo <- orderCells(HSMM_myo)
#plot by condition/timepoint
plot_cell_trajectory(HSMM_myo, color_by = "Hours")
#plot by state
plot_cell_trajectory(HSMM_myo, color_by = "State")
#pass state to pData
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$Hours)[,"0"]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  }else {
    return (1)
  }
}
HSMM_myo <- orderCells(HSMM_myo, root_state = GM_state(HSMM_myo))
plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime")
#facet the trajectory according to states
plot_cell_trajectory(HSMM_myo, color_by = "State") + facet_wrap(~State, nrow = 1)

# Note: If no timeseries is available one can set the root of the tree
#        based on progenitor genes
blast_genes <- row.names(subset(fData(HSMM_myo),
                                gene_short_name %in% c("CCNB2", "MYOD1", "MYOG")))
plot_genes_jitter(HSMM_myo[blast_genes,], grouping = "State", min_expr = 0.1)

#Check with key genes whether states are correct
HSMM_expressed_genes <-  row.names(subset(fData(HSMM_myo),
                                          num_cells_expressed >= 10))
HSMM_filtered <- HSMM_myo[HSMM_expressed_genes,]
my_genes <- row.names(subset(fData(HSMM_filtered),
                             gene_short_name %in% c("CDK1", "MEF2C", "MYH3")))
cds_subset <- HSMM_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "Hours")

#####
#Alternative for inferring order genes
#####
# Note: requires clusters
#expressed in 5% of the cells
HSMM_myo <- detectGenes(HSMM_myo, min_expr = 0.1)
fData(HSMM_myo)$use_for_ordering <- fData(HSMM_myo)$num_cells_expressed > 0.05 * ncol(HSMM_myo)
plot_pc_variance_explained(HSMM_myo, return_all = F)
#choose PCs describing vatiation and reduce
HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2, norm_method = 'log',
                            num_dim = 3, reduction_method = 'tSNE', verbose = T)
#Calculate local densities -> distance of cells to each other
plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Cluster)')
plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Hours)')

#Decision plot to infer distance P
plot_rho_delta(HSMM_myo, rho_threshold = 2, delta_threshold = 4 )
#Rerun selection by defining P and Rho
HSMM_myo <- clusterCells(HSMM_myo,
                         rho_threshold = 2,
                         delta_threshold = 4,
                         skip_rho_sigma = T,
                         verbose = F)
#Check results
plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Cluster)')
plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Hours)')

#Perform DEG analysis
# >possible to set cores
clustering_DEG_genes <- differentialGeneTest(HSMM_myo[HSMM_expressed_genes,],
                                             fullModelFormulaStr = '~Cluster',
                                             cores = 1)
#Select top1000 for ordering
HSMM_ordering_genes <-row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes = HSMM_ordering_genes)
HSMM_myo <- reduceDimension(HSMM_myo, method = 'DDRTree')
HSMM_myo <- orderCells(HSMM_myo)
HSMM_myo <- orderCells(HSMM_myo, root_state = GM_state(HSMM_myo))
plot_cell_trajectory(HSMM_myo, color_by = "Hours")

# Note: Also possible to use the genes with high dispresion, wathc out for 'decent' mean expression
#disp_table <- dispersionTable(HSMM_myo)
#ordering_genes <- subset(disp_table,
 #                        mean_expression >= 0.5 &
  #                         dispersion_empirical >= 1 * dispersion_fit)$gene_id
#####
#Semi-supervised ordering 
#####
#allows exclusion of e.g. cycling genes
CCNB2_id <- row.names(subset(fData(HSMM_myo), gene_short_name == "CCNB2"))
MYH3_id <- row.names(subset(fData(HSMM_myo), gene_short_name == "MYH3"))

cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "Cycling myoblast",
                   classify_func = function(x) { x[CCNB2_id,] >= 1 })
cth <- addCellType(cth, "Myotube", classify_func = function(x) { x[MYH3_id,] >= 1 })
cth <- addCellType(cth, "Reserve cell",
                   classify_func = function(x) { x[MYH3_id,] == 0 & x[CCNB2_id,] == 0 })
HSMM_myo <- classifyCells(HSMM_myo, cth)
#Select genes that co-vary
marker_diff <- markerDiffTable(HSMM_myo[HSMM_expressed_genes,],
                               cth,
                               cores = 1)
#semisup_clustering_genes <- row.names(subset(marker_diff, qval < 0.05))
semisup_clustering_genes <- row.names(marker_diff)[order(marker_diff$qval)][1:1000]

#More defined trajectory based on top 1000 DEGs semisupervised
HSMM_myo <- setOrderingFilter(HSMM_myo, semisup_clustering_genes)
#plot_ordering_genes(HSMM_myo)
HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2,
                            method = 'DDRTree', norm_method = 'log')
HSMM_myo <- orderCells(HSMM_myo)
HSMM_myo <- orderCells(HSMM_myo, root_state = GM_state(HSMM_myo))
plot_cell_trajectory(HSMM_myo, color_by = "CellType") +
  theme(legend.position = "right")

#Exclude not fully differentiated cells
HSMM_filtered <- HSMM_myo[HSMM_expressed_genes,]
my_genes <- row.names(subset(fData(HSMM_filtered),
                             gene_short_name %in% c("CDK1", "MEF2C", "MYH3")))
cds_subset <- HSMM_filtered[my_genes,]
plot_genes_branched_pseudotime(cds_subset,
                               branch_point = 1,
                               color_by = "Hours",
                               ncol = 1)

#####
#DEGs
#####
#Computational intensive
marker_genes <- row.names(subset(fData(HSMM_myo),
                                 gene_short_name %in% c("MEF2C", "MEF2D", "MYF5",
                                                        "ANPEP", "PDGFRA","MYOG",
                                                        "TPM1",  "TPM2",  "MYH2",
                                                        "MYH3",  "NCAM1", "TNNT1",
                                                        "TNNT2", "TNNC1", "CDK1",
                                                        "CDK2",  "CCNB1", "CCNB2",
                                                        "CCND1", "CCNA1", "ID1",
                                                        "ACTB")))
#Identify DEGs based on protocol/media change
diff_test_res <- differentialGeneTest(HSMM_myo[marker_genes,],
                                      fullModelFormulaStr = "~Media")
# Select genes that are significant at an FDR < 10%
sig_genes <- subset(diff_test_res, qval < 0.1)
sig_genes[,c("gene_short_name", "pval", "qval")]
#Plot expression
MYOG_ID1 <- HSMM_myo[row.names(subset(fData(HSMM_myo),
                                      gene_short_name %in% c("MYOG", "CCNB2"))),]
plot_genes_jitter(MYOG_ID1, grouping = "Media", ncol= 2)


#1) DEGs that distinguish celltypes
to_be_tested <- row.names(subset(fData(HSMM),
                                 gene_short_name %in% c("UBC", "NCAM1", "ANPEP")))
cds_subset <- HSMM[to_be_tested,]
# Note: Two model comparison system ranking prediction quality of gene of interest
# Note: Any column can be used for compariosn in pData
diff_test_res <- differentialGeneTest(cds_subset, fullModelFormulaStr = "~CellType")
diff_test_res[,c("gene_short_name", "pval", "qval")]
plot_genes_jitter(cds_subset, grouping = "CellType", color_by = "CellType",
                  nrow= 1, ncol = NULL, plot_trend = TRUE)
#Identical to above
full_model_fits <- fitModel(cds_subset,  modelFormulaStr = "~CellType")
reduced_model_fits <- fitModel(cds_subset, modelFormulaStr = "~1")
diff_test_res <- compareModels(full_model_fits, reduced_model_fits)
diff_test_res

#2) DEGs over pseudotime
to_be_tested <- row.names(subset(fData(HSMM),
                                 gene_short_name %in% c("MYH3", "MEF2C", "CCNB2", "TNNT1")))
cds_subset <- HSMM_myo[to_be_tested,]
#Make model testing against changes to predicted and observed changes in pseudotime
#Non-linear expression modelling
diff_test_res <- differentialGeneTest(cds_subset,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
diff_test_res[,c("gene_short_name", "pval", "qval")]
plot_genes_in_pseudotime(cds_subset, color_by = "Hours")

#Clustering genes by pseudotemporal expression pattern
diff_test_res <- differentialGeneTest(HSMM_myo[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
plot_pseudotime_heatmap(HSMM_myo[sig_gene_names,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T)

#4) Multi-factorial DEGs
#Substract factors like time
to_be_tested <- row.names(subset(fData(HSMM),
                                 gene_short_name %in% c("TPM1", "MYH3", "CCNB2", "GAPDH")))
cds_subset <- HSMM[to_be_tested,]
diff_test_res <- differentialGeneTest(cds_subset,
                                      fullModelFormulaStr = "~CellType + Hours",
                                      reducedModelFormulaStr = "~Hours")
diff_test_res[,c("gene_short_name", "pval", "qval")]
plot_genes_jitter(cds_subset,
                  grouping = "Hours", color_by = "CellType", plot_trend = TRUE) +
                  facet_wrap( ~ feature_label, scales= "free_y")


#####
#Identify DEGs at branch decisions
#####
lung <- load_lung()
plot_cell_trajectory(lung, color_by = "Time")
#Choose Branch point and obtain DEGs per branch point
BEAM_res <- BEAM(lung, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
plot_genes_branched_heatmap(lung[row.names(subset(BEAM_res, qval < 1e-4)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
#genes over time
lung_genes <- row.names(subset(fData(lung),
                               gene_short_name %in% c("Ccnd2", "Sftpb", "Pdpn")))
plot_genes_branched_pseudotime(lung[lung_genes,],
                               branch_point = 1,
                               color_by = "Time",
                               ncol = 1)

#####
#Notes
#####
#Check for DESeq2 implementation of UMI counts

