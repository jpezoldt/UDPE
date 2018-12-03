# dissect compiled and curated Output from Homer
# Author: Joern Pezoldt
# Date: 12.10.2018
# Type of script: Exploratory

# Functions:
# 1) Identify TFs from curated files 
# 
# 3) Plot expression of TFs/regulons across cells


########### Libraries
library(pheatmap)
library(ggplot2)



############ PATH & Global Variables ----------------------------

sample_ID <- "SPF"
condition <- "pLN_SPF"
# path to binary matrix of TFs detected using homer via motifs identified
PATH_input_TF_homer_matrix <- paste("/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/Motif/Homer_Output/allHomerResults/size_500bp/", sample_ID, sep = "")
#PATH_input_TF_homer_matrix <- paste("/Volumes/Pezoldt_HD_4TB/UPDE_Transfer/20181101/Analysis/ATACseq/ATAC_FSC_all/Motif/Homer_Output/allHomerResults/size_500bp/", sample_ID, sep = "")



# path to regulons identified
# Note: Input required
organ <- c("pLN_SPF")
cell_subsets <- c("KeySubsets")
time_point <- c("D56")
#"D0","D10","D24",   ,"D300"
#PATH_analysis_scenic <- "/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic"
PATH_analysis_scenic <- "/Volumes/Pezoldt_HD_4TB/UPDE_Transfer/20181101/Analysis/scRNAseq/scenic"

#PATH_input_regulons_scenic <- "/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/mLN_SPF"
PATH_input_regulons_scenic <- "/Volumes/Pezoldt_HD_4TB/UPDE_Transfer/20181101/Analysis/scRNAseq/scenic/mLN_SPF"

cell_subsets_1 <- "KeySubsets_D56"
cell_subsets_2 <- "nonAdventi"



############ Functions ----------------------------------------
#' loadCSVs
#' loads CSV tables containing the regulon score per cell
#' @param PATH_analysis_scenic path to the general folder
#' @param organ character vector with organs of choice
#' @param cell_subsets character vector with cell_subsets of choice
#' @param time_point character vector with time points of choice
#' @param type_regulon either "core" or "extended" regulons
#' @return l_Regulon_tables object contains only the cells that align with the predefined annotation of in cluster_of_interest 

loadCSVs <- function(PATH_analysis_scenic, organ,cell_subsets,time_point,type_regulon){
  #go along the nameing vectors and concatenate to obtain names and store in vector
  v_names <- c()
  v_PATHs <- c()
  for(k in 1:length(time_point)){
    time_point_k <- time_point[k]
    v_names_ijk <- paste(organ,cell_subsets,time_point_k,sep = "_")
    v_names <- c(v_names, v_names_ijk)
    v_PATHs_ijk <- paste(PATH_analysis_scenic,"/",organ,"/",cell_subsets,"_",time_point_k,"/",sep="")
    v_PATHs <- c(v_PATHs,v_PATHs_ijk)
  }
  
  #returned names list of tables
  # core regulons
  # extended regulons
  l_regulons <- list()
  for(i in 1:length(v_PATHs)){
    t_regulon_i <- read.table(paste(v_PATHs[i],"AUCell_scores_",type_regulon,"_",v_names[i],".csv", sep=""),dec = ".", sep = ",")
    l_regulons[[i]] <- t_regulon_i
    
  }
  names(l_regulons) <- v_names
  
  #return
  l_regulons
}

#' regulon_seuratNames
#' append cell names defined by Seurat via the scenic CellInfo file
#' @param PATH_analysis_scenic path to the general folder
#' @param l_regulon list of regulon AUC score tables
#' @param organ character vector with organs of choice
#' @param cell_subsets character vector with cell_subsets of choice
#' @param time_point character vector with time points of choice
#' @return l_Regulon_tables object contains only the cells that align with the predefined annotation of in cluster_of_interest 
regulon_seuratNames <- function(l_regulon,PATH_analysis_scenic, organ,cell_subsets,time_point){
  l_regulon_seuratNamed <- list()
  for(i in 1:length(l_regulon_ext)){
    #load CellInfo file
    CellInfo_i <- readRDS(paste(PATH_analysis_scenic,"/",organ,"/",cell_subsets,"_",time_point[i],"/int/cellInfo.Rds",sep=""))
    #get cell names
    cell_names_i <- paste(CellInfo_i$cell_type1,"_",rownames(CellInfo_i), sep = "")
    #replace column names from regulon table with cell names derived from Seurat
    l_regulon_ext_i <- l_regulon_ext[[i]]
    colnames(l_regulon_ext_i)[(ncol(l_regulon_ext_i)-length(cell_names_i)+1):length(l_regulon_ext_i)] <- cell_names_i
    l_regulon_seuratNamed[[i]] <- l_regulon_ext_i
  }
  names(l_regulon_seuratNamed) <- paste(organ, cell_subsets, time_point, sep = "_")
  l_regulon_seuratNamed
}


#' annotation_pheatmap
#' load and store annotation from CellInfo in a "pheatmap" compatible format
#' @param PATH_analysis_scenic path to the general folder
#' @param l_regulon list of regulon AUC score tables
#' @param organ character vector with organs of choice
#' @param cell_subsets character vector with cell_subsets of choice
#' @param time_point character vector with time points of choice
#' @return l_annotation_pheatmap object contains list of annotation files for pheatmap
annotation_pheatmap <- function(l_regulon,PATH_analysis_scenic, organ,cell_subsets,time_point){
  l_annotation_pheatmap <- list()
  for(i in 1:length(time_point)){
    #load CellInfo file
    CellInfo_i <- readRDS(paste(PATH_analysis_scenic,"/",organ,"/",cell_subsets,"_",time_point[i],"/int/cellInfo.Rds",sep=""))
    #get cell names
    annotation_heatmap_i <- data.frame(cell_type1 = CellInfo_i$cell_type1)
    rownames(annotation_heatmap_i) <- paste(CellInfo_i$cell_type1,row.names(CellInfo_i),sep="_")
    l_annotation_pheatmap[[i]] <- annotation_heatmap_i
  }
  l_annotation_pheatmap
}


#' regulons_significant
#' Identify significantly different regulons across subsets per conditions
#' @param PATH_analysis_scenic path to the general folder
#' @param l_regulon list of regulon AUC score tables
#' @param organ character vector with organs of choice
#' @param cell_subsets character vector with cell_subsets of choice
#' @param time_point character vector with time points of choice
#' @return l_regulons_significant object that contains a list of vectors containing the significant TFs

regulons_significant <- function(l_regulon_seuratNames_core,l_annotation_pheatmap){
  #Identify significant TFs
  #loop over regulon matrices
  l_regulons_significant <- list()
  for(j in 1:length(l_regulon_seuratNames_core)){
    #format matrix for extraction of regulons per TF
    regulon_Matrix_j <- l_regulon_seuratNames_core[[j]]
    rownames(regulon_Matrix_j) <- regulon_Matrix_j$GeneSymbol
    regulon_Matrix_j <- regulon_Matrix_j[,4:ncol(regulon_Matrix_j)]
    #load annotation table
    l_annotation_pheatmap_j <- l_annotation_pheatmap[[j]]
    #Generate list of tables per motif that can be processed with the aov()
    l_regulon_to_anova <- list()
    for(i in 1:nrow(regulon_Matrix_j)){
      TFmotif <- rep(as.character(rownames(regulon_Matrix_j)[i]),ncol(regulon_Matrix_j))
      AUCell_scores <- as.numeric(regulon_Matrix_j[i,])
      names(AUCell_scores) <- colnames(regulon_Matrix_j)
      Cluster <- as.character(l_annotation_pheatmap_j$cell_type1)
      regulon_to_anova_i <- data.frame(TFmotif = TFmotif, AUCell_scores = AUCell_scores, Cluster = Cluster)
      l_regulon_to_anova[[i]] <- regulon_to_anova_i
    }
    names(l_regulon_to_anova) <- rownames(regulon_Matrix_j)
    
    #Run anova per TF/motif
    l_res.aov <- lapply(l_regulon_to_anova, function(x){
      aov(AUCell_scores ~ Cluster, data = x)
    })
    names(l_res.aov) <- rownames(regulon_Matrix_j)
    
    test_aov <- l_res.aov[[1]]
    #Extract p-values
    p_values_aov <- lapply(l_res.aov, function(x){
      summary(x)[[1]][["Pr(>F)"]][1]
    })
    
    #FDR correct the aov() p-values and thresh
    p_values_aov_fdr <- p.adjust(p_values_aov, method = "fdr")
    
    #Select regulons with with padj < 0.05
    p_values_aov_fdr_thres <- p_values_aov_fdr[p_values_aov_fdr <= 0.05]
    
    #Select regulons that are FDR corrected
    l_res.aov_fdr <- l_res.aov[names(p_values_aov_fdr_thres)]
    
    #Run Tukey FDR pairwise-comparison
    l_res.tuk <- lapply(l_res.aov_fdr, function(x){
      TukeyHSD(x)
    })
    names(l_res.tuk) <- names(l_res.aov_fdr)
    #l_res.tuk[["Nfkb2"]]
    #l_res.aov[["Nfkb2"]]
    
    #Identify significant regulons with biggest difference among comparison
    diff_tuk <- unlist(lapply(l_res.tuk, function(x){
      max(abs(x$Cluster[,1]))
    }))
    #select regulons with a maximal diff
    diff_tuk_names <- diff_tuk[diff_tuk >= 0.04]
    print(paste(j, ":", names(diff_tuk_names), sep = " "))
    l_regulons_significant[[j]] <- diff_tuk_names
  }
  l_regulons_significant
}

#### Process data with function--------------------------------------------------

#Load tables
type_regulon <- "ext"
l_regulon_ext <- loadCSVs(PATH_analysis_scenic, organ,cell_subsets,time_point,type_regulon)
type_regulon <- "core"
l_regulon_core <- loadCSVs(PATH_analysis_scenic, organ,cell_subsets,time_point,type_regulon)

# Rename regulon tables
l_regulon_seuratNames_ext <- regulon_seuratNames(l_regulon_ext,PATH_analysis_scenic, organ,cell_subsets,time_point)
l_regulon_seuratNames_core <- regulon_seuratNames(l_regulon_core,PATH_analysis_scenic, organ,cell_subsets,time_point)

#Load annotation for pheatmap
l_annotation_pheatmap <- annotation_pheatmap(l_regulon,PATH_analysis_scenic, organ,cell_subsets,time_point)

#Calculate significant regulons
l_regulons_significant <- regulons_significant(l_regulon_seuratNames_core,l_annotation_pheatmap)

########################
#### Process exploratory
########################
#Make heatmap for all regulons
indexed <- 1
regulon_seuratNames_ext_x <- l_regulon_seuratNames_ext[[indexed]]
rownames(regulon_seuratNames_ext_x) <- regulon_seuratNames_ext_x$GeneSymbol
regulon_seuratNames_ext_heat_x <- as.matrix(regulon_seuratNames_ext_x[,(ncol(regulon_seuratNames_ext_x)-nrow(l_annotation_pheatmap[[indexed]])+1):ncol(regulon_seuratNames_ext_x)])
#Heatmap across all regulons and all cells
pheatmap(regulon_seuratNames_ext_heat_x, annotation = l_annotation_pheatmap[[indexed]],cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 15, show_rownames = FALSE, cluster_cols = FALSE,
         scale = "none", show_colnames = FALSE,
         color = colorRampPalette(c("black", "white","green"), space="rgb")(128),
         main = paste("All regulons ", time_point[[indexed]], sep = ""))

#Make heatmap for significant regulons
indexed <- 1
regulon_seuratNames_ext_x <- l_regulon_seuratNames_ext[[indexed]]
l_regulons_significant_x <- names(l_regulons_significant[[indexed]])
regulon_seuratNames_ext_sig_x <- subset(regulon_seuratNames_ext_x, GeneSymbol %in% l_regulons_significant_x)
regulon_seuratNames_ext_sig_heat_x <- as.matrix(regulon_seuratNames_ext_sig_x[,(ncol(regulon_seuratNames_ext_sig_x)-nrow(l_annotation_pheatmap[[indexed]])+1):ncol(regulon_seuratNames_ext_sig_x)])
pheatmap(regulon_seuratNames_ext_sig_heat_x, annotation = l_annotation_pheatmap[[indexed]],cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 15, cluster_cols = FALSE,
         scale = "none", show_colnames = FALSE,show_rownames = TRUE,
         color = colorRampPalette(c("black", "white","green"), space="rgb")(128),
         main = paste("Significant_TopDiff_regulons ", time_point[[indexed]], sep = ""))

#Boxplot of regulon
name_x <- names(diff_tuk_names)[1]
test_boxplot <- l_regulon_to_anova[[name_x]]
ggplot(test_boxplot,aes(factor(Cluster), AUCell_scores)) +
  geom_violin(aes(fill = Cluster)) +
  geom_jitter(size = 0.5) +
  ggtitle(name_x)

#motif TFs ATAC
TFs_ATAC <- read.table(paste(PATH_input_TF_homer_matrix,"/",sample_ID,"_existence_matrix_known_homer_compile.txt", sep = ""), dec = ".", sep = "\t")

############ Extract TFs of interest ----------------------------
TFs_interest <- subset(TFs_ATAC, 
                                  Open_DEG_SPF_lowRNAseq == 0 &  
                                  Closed_DEG_SPF_lowRNAseq == 0 &
                                  No_DAR_UP_SPF_lowRNAseq == 0 &
                                  No_DAR_DOWN_SPF_lowRNAseq == 0 &
                                  NoDEG_Closed_SPF_lowRNAseq == 0 &
                                  NoDEG_Open_SPF_lowRNAseq == 1)
upset(TFs_ATAC, sets = colnames(TFs_ATAC),
      mainbar.y.label = "TFs ATACseq",
      number.angles = 30, point.size = 5,
      text.scale = 2, keep.order = TRUE)


############ Overlap and heatmap -------------------------------
indexed <- 1
regulon_seuratNames_ext_x <- l_regulon_seuratNames_ext[[indexed]]
regulon_seuratNames_ext_x_TFs <- subset(regulon_seuratNames_ext_x, GeneSymbol %in% rownames(TFs_interest))
rownames(regulon_seuratNames_ext_x_TFs) <- regulon_seuratNames_ext_x_TFs$GeneSymbol
#Heatmap across all regulons and all cells
pheatmap(as.matrix(regulon_seuratNames_ext_x_TFs[,4:ncol(regulon_seuratNames_ext_x_TFs)]),
         annotation = l_annotation_pheatmap[[indexed]],cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 15, cluster_cols = FALSE,
         scale = "none", show_colnames = FALSE,show_rownames = TRUE,
         color = colorRampPalette(c("black", "white","green"), space="rgb")(128),
         main = paste("ATACseq_TF_regulons_mLN_ ", time_point[[indexed]], sep = ""))


########### Include content of regulon---------------------------
l_PATH_scenicOptions <- paste(PATH_analysis_scenic,"/",organ,"/",cell_subsets,"_",time_point,"/int/scenicOptions.Rds",sep="")
l_PATH_scenicOptions_NO_int <- paste(PATH_analysis_scenic,"/",organ,"/",cell_subsets,"_",time_point,"/",sep="")

indexed <- 2
scenicOptions <- readRDS(l_PATH_scenicOptions[[indexed]])
loadInt(scenicOptions)
setwd(paste(l_PATH_scenicOptions_NO_int[[indexed]]))
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
#Get regulon dat
regulonTargetsInfo <- RcisTarget::addLogo(regulonTargetsInfo, motifCol="bestMotif")
regulonTargetsInfo$Genie3Weight <- signif(regulonTargetsInfo$Genie3Weight, 2)
regulonTargetsInfo_highConfAnnot <- subset(regulonTargetsInfo, highConfAnnot == "TRUE")
#split by Regulon
l_regulonTargets <- split(regulonTargetsInfo, regulonTargetsInfo$TF)
#Get all the genes behind the regulons
all_genes_in_regulons <- unlist(lapply(l_regulonTargets, function(x){
  genes_i <- x$gene
  genes_i
}))

test <- subset(TFs_interest, rownames(TFs_interest) %in% all_genes_in_regulons)


#1) Identify the genes that are targeted by TFs present in regulons

#2) Plot expression of genes at the different stages and cell subsets





# Generate some data
test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
colnames(test) = paste("Test", 1:10, sep = "")
rownames(test) = paste("Gene", 1:20, sep = "")

# Draw heatmaps
pheatmap(test)
pheatmap(test, kmeans_k = 2)
pheatmap(test, scale = "row", clustering_distance_rows = "correlation")
pheatmap(test, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
pheatmap(test, cluster_row = FALSE)
pheatmap(test, legend = FALSE)
pheatmap(test, display_numbers = TRUE)
pheatmap(test, display_numbers = TRUE, number_format = "%.1e")
pheatmap(test, cluster_row = FALSE, legend_breaks = -1:4, legend_labels = c("0",
                                                                            "1e-4", "1e-3", "1e-2", "1e-1", "1"))
pheatmap(test, cellwidth = 12, cellheight = 12, main = "Example heatmap")
#pheatmap(test, cellwidth = 15, cellheight = 12, fontsize = 8, filename = "test.pdf")


# Generate column annotations
annotation = data.frame(Var1 = factor(1:10 %% 2 == 0,
                                      labels = c("Class1", "Class2")), Var2 = 1:10)
annotation$Var1 = factor(annotation$Var1, levels = c("Class1", "Class2", "Class3"))
rownames(annotation) = paste("Test", 1:10, sep = "")

pheatmap(test, annotation = annotation)
pheatmap(test, annotation = annotation, annotation_legend = FALSE)
pheatmap(test, annotation = annotation, annotation_legend = FALSE, drop_levels = FALSE)

# Specify colors
Var1 = c("navy", "darkgreen")
names(Var1) = c("Class1", "Class2")
Var2 = c("lightgreen", "navy")

ann_colors = list(Var1 = Var1, Var2 = Var2)

#Specify row annotations
row_ann <- data.frame(foo=gl(2,nrow(test)/2),`Bar`=relevel(gl(2,nrow(test)/2),"2"))
rownames(row_ann)<-rownames(test)
pheatmap(test, annotation = annotation, annotation_legend = FALSE,
         drop_levels = FALSE,row_annotation = row_ann)

#Using cytokine annotations
M<-matrix(rnorm(8*20),ncol=8)
row_annotation<-data.frame(A=gl(4,nrow(M)/4),B=gl(4,nrow(M)/4))
eg<-expand.grid(factor(c(0,1)),factor(c(0,1)),factor(c(0,1)))
colnames(eg)<-c("IFNg","TNFa","IL2")
rownames(eg)<-apply(eg,1,function(x)paste0(x,collapse=""))
rownames(M)<-1:nrow(M)
colnames(M)<-rownames(eg)
cytokine_annotation=eg
pheatmap(M,annotation=annotation,row_annotation=row_annotation,
         annotation_legend=TRUE,row_annotation_legend=TRUE,
         cluster_rows=FALSE,cytokine_annotation=cytokine_annotation,cluster_cols=FALSE)

# Specifying clustering from distance matrix
drows = dist(test, method = "minkowski")
dcols = dist(t(test), method = "minkowski")
pheatmap(test, clustering_distance_rows = drows, clustering_distance_cols = dcols)

