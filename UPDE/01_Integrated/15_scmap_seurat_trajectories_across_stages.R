# Author: Joern Pezoldt
# Functions:
# 1) Make annotated SCE object from monocle SEURAT object
# 2) scmap subset specification across other SEURAT/SCE objects

# Libraries
library(monocle)
library(scmap)
library(SingleCellExperiment)
library(plotly)

#####
#PATHs and global variables
#####
# Input required
time_point <- "mLN_D0"

PATH_Key_Subsets <- paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/seurat/Subsets_to_scenic/",time_point,"/",time_point,".Rds",sep="")
time_points <- c("mLN_D10","mLN_D24","mLN_D56","mLN_D300","pLN_D56")

v_PATH_mapped_to_Subsets <- paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/seurat/Subsets_to_scenic/",time_points,"/",time_points,".Rds",sep="")
PATH_store_output <- "/home/pezoldt/NAS2/pezoldt/Analysis/scmap/seurat_annotations"
# ToDo: add folder names for different stages containing the .Rds of the monocle files
# v_Ages <- ("")
  
#scmap Params
#Genes to include
n_feature <- 1000
#Threshold for percentage of match from cell-to-cell cluster projection
threshold <- 0.7


#####
#Functions
#####
#' transform_SEURAT_to_SCE
#' loads CSV tables containing the regulon score per cell
#' @param PATH_Key_Subsets path to the general folder
#' @param v_PATH_mapped_to_Subsets character vector with organs of choice
#' @return returns list of sce objects with the first one being the reference composition

transform_SEURAT_to_SCE <- function(PATH_Key_Subsets, v_PATH_mapped_to_Subsets){
  #make string vector containing all the directories with the SEURAT to be transformed
  l_transform_SEURAT_to_SCE <- list()
  #store all PATHs in one vector of which the first one 
  v_PATH_Key_Subsets <- c(PATH_Key_Subsets,v_PATH_mapped_to_Subsets)
  for(i in 1:length(v_PATH_Key_Subsets)){
    v_PATH_Key_Subsets_i <- v_PATH_Key_Subsets[i]
    seurat_i <- readRDS(v_PATH_Key_Subsets_i)
    # Obtain expression matrix
    seurat_count_matrix_i <- as.matrix(seurat_i@raw.data)
    #Annotate by cluster == cell_type1 or State of trajectory == State
    seurat_annotation_i <- data.frame(cell_type1=seurat_i@ident, row.names = rownames(seurat_i@meta.data))
    seurat_count_matrix_i <- seurat_count_matrix_i[,rownames(seurat_annotation_i)]
    #Generate SCE object
    SCE_i <- SingleCellExperiment(assays = list(counts = as.matrix(seurat_count_matrix_i)), colData = seurat_annotation_i)
    logcounts(SCE_i) <- log2(counts(SCE_i) + 1)
    #Genesymbols as feature symbols
    rowData(SCE_i)$feature_symbol <- rownames(SCE_i)
    #eliminate duplicated rows
    SCE_i <- SCE_i[!duplicated(rownames(SCE_i)), ]
    l_transform_SEURAT_to_SCE[[i]] <- SCE_i
  }
  #name list elements
  names(l_transform_SEURAT_to_SCE) <- v_PATH_Key_Subsets
  #return
  l_transform_SEURAT_to_SCE
}

#####Run Functions----------------------------------------------------------------------------
#Generate SCEs from SEURAT
l_transform_SEURAT_to_SCE <- transform_SEURAT_to_SCE(PATH_Key_Subsets, v_PATH_mapped_to_Subsets)

#Plot Sankey Plots ---------------------------------------------------------------------------
#Select features
for(i in 1:length(l_transform_SEURAT_to_SCE)){
  l_transform_SEURAT_to_SCE[[i]] <- selectFeatures(l_transform_SEURAT_to_SCE[[i]], suppress_plot = FALSE, n_features = n_feature)
  table(rowData(l_transform_SEURAT_to_SCE[[1]])$scmap_features)
  #Calculate mean expression per cluster based on ID stored in "cell_type1" column
  l_transform_SEURAT_to_SCE[[i]] <- indexCluster(l_transform_SEURAT_to_SCE[[i]])
  
}

#Map reference onto projection
# First element of SCE list is the reference
l_scmapCluster_mapping <- list()
for(i in 1:(length(l_transform_SEURAT_to_SCE)-1)){
  reference <- l_transform_SEURAT_to_SCE[[1]]
  mapping_i <- scmapCluster(projection = reference,
                            threshold = threshold,
                            index_list = list(
                            older = metadata(l_transform_SEURAT_to_SCE[[(i+1)]])$scmap_cluster_index))
  l_scmapCluster_mapping[[i]] <- mapping_i
}
names(l_scmapCluster_mapping) <- c(time_points)

#plot Sankeys opening a browser window each
for(i in 1:length(l_scmapCluster_mapping)){
  print(i)
  plot(
    getSankey(
      colData(l_transform_SEURAT_to_SCE[[1]])$cell_type1, 
      l_scmapCluster_mapping[[i]]$scmap_cluster_labs[,'older'],
      plot_height = 500, plot_width = 300)
    #main = names(l_scmapCluster_mapping)[i]
  )
}





