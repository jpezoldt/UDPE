# Author: Joern Pezoldt
# Date: 31.11.2018
# Function:
# 1) Plot DMRs


#libraries 
library(pheatmap)

#####
#PATH & Global
#####
PATH_output <- "/home/pezoldt/NAS2/pezoldt/Analysis/WGBS/FSC/Plots"

####################
#Load data
####################
#Load CpGs
CpGs <- read.csv("/home/pezoldt/NAS2/pezoldt/Data/WGBS/WGBS_FSC/CpGs_in_DMRs.csv")


#Load Bsmooth data
mLN_SPF_GF <- read.delim2("/home/pezoldt/NAS2/pezoldt/Data/WGBS/WGBS_FSC/BSmooth_DMRs_mLN_SPF_vs_mLN_GF.txt")
pLN_SPF_GF <- read.delim2("/home/pezoldt/NAS2/pezoldt/Data/WGBS/WGBS_FSC/BSmooth_DMRs_pLN_SPF_vs_pLN_GF.txt")
SPF_mLN_pLN <- read.delim2("/home/pezoldt/NAS2/pezoldt/Data/WGBS/WGBS_FSC/BSmooth_DMRs_mLN_SPF_vs_pLN_SPF.txt")
GF_mLN_pLN <- read.delim2("/home/pezoldt/NAS2/pezoldt/Data/WGBS/WGBS_FSC/BSmooth_DMRs_pLN_GF_vs_mLN_GF.txt")


##############
#Assemble essential DMR info
##############
Bsmooth_tables <- list(SPF_mLN_pLN, mLN_SPF_GF,
                       GF_mLN_pLN, pLN_SPF_GF)
Bsmooth_names_analysis <- c("SPF_mLN_pLN", "mLN_SPF_GF",
                            "GF_mLN_pLN","pLN_SPF_GF")

#DMR_identifier <- paste(mLN_SPF_pLN_SPF$chr, mLN_SPF_pLN_SPF$start, mLN_SPF_pLN_SPF$end, sep = "_")
Bsmooth_intel <- data.frame(DMR_identifier = character(),
                            closest_gene = character(),
                            distance_closest_TSS = numeric(),
                            DMR_comp_origin = character(),
                            meanDiff = numeric(),
                            areaStat = numeric(),
                            width = numeric())

#Assemble all necessary identfiers and calculations from Bsmooth pair-wisecomparisons into
# one list
for(i in 1:4){
  x <- Bsmooth_tables[[i]]
  DMR_identifier <- paste(x$chr, x$start, x$end, sep = "_")
  closest_gene <- as.character(x$nearest_gene_name)
  distance_closest_TSS <- as.numeric(x$distance_to_nearest_tss)
  name_table <- Bsmooth_names_analysis[1]
  DMR_comp_origin <- rep(name_table, len = length(closest_gene)) 
  meanDiff <- as.numeric(x$meanDiff)
  areaStat <- as.numeric(x$areaStat)
  width <- as.numeric(x$width)
  table_i <- cbind(DMR_identifier, closest_gene,
                   distance_closest_TSS, DMR_comp_origin,
                   meanDiff, areaStat, width)
  Bsmooth_intel <- rbind(Bsmooth_intel, table_i)
}  
Bsmooth_intel$distance_closest_TSS <- as.numeric(Bsmooth_intel$distance_closest_TSS)
Bsmooth_intel$meanDiff <- as.numeric(as.character(Bsmooth_intel$meanDiff))
Bsmooth_intel$areaStat <- as.numeric(as.character(Bsmooth_intel$areaStat))
Bsmooth_intel$width <- as.numeric(as.character(Bsmooth_intel$width))


######
#Pick or threshhold DMRs
######
Bsmooth_intel_subset <- subset(Bsmooth_intel, (Bsmooth_intel$meanDiff <= -0.5 | Bsmooth_intel$meanDiff >= 0.5) |
                                 (Bsmooth_intel$areaStat >= 300 | Bsmooth_intel$areaStat <= -300))
#Pick
# Genes of interest
Genes_DMRs_interest <- c("Gdf6","Ptger2","Rsrc1","Agr2","Ptgis","Aldh1a2","Cxcl13","Wnt4")
Bsmooth_intel_subset <- subset(Bsmooth_intel, closest_gene %in% Genes_DMRs_interest)

#############
#Grab Threshed DMRs from CpG list
#############
#split according to DMR_identifier
CpGs_DMRs <- split(CpGs, CpGs$DMR_identifier)

#get DMR_identifiers from theshed list
threshed_DMRs <- as.character(Bsmooth_intel_subset$DMR_identifier)
CpGs_DMRs_threshed <- CpGs_DMRs[threshed_DMRs]

#########
#heatmap
#########

#get data from Bsmooth analysis
path_variable = paste(PATH_output, "/DMR_per_Gene.pdf", sep="")
pdf(file = path_variable, height = 2.3, width = 20)

for(i in 1:length(CpGs_DMRs_threshed)){
  par(mfcol = c(4,1))
  
  #Get info for main/title
  CpGs_DMRs_i <- CpGs_DMRs_threshed[[5]]
  data_heatmap <- CpGs_DMRs_i
  data_heatmap <- unique(data_heatmap)
  DMR_identifier_i <- as.character(CpGs_DMRs_i$DMR_identifier)[1]
  DMR_meanDiff <- as.character(Bsmooth_intel[Bsmooth_intel$DMR_identifier == DMR_identifier_i,5])
  DMR_genename <- as.character(Bsmooth_intel[Bsmooth_intel$DMR_identifier == DMR_identifier_i,2])
  DMR_dist_TSS <- as.character(Bsmooth_intel[Bsmooth_intel$DMR_identifier == DMR_identifier_i,3])
  
  #get CpG positions
  print(i)
  rownames(data_heatmap) <- as.character(data_heatmap$CpG_start)
  
  #make matrix
  data_heatmap <- data_heatmap[,8:17]
  min(data_heatmap)
  max(data_heatmap)
  data_heatmap_matrix <- as.matrix(data_heatmap)
  data_heatmap_matrix_t <- t(data_heatmap_matrix)
  
  #replace NA with 100, results in white squares in heatmap
  data_heatmap_matrix[is.na(data_heatmap_matrix_t)] <- 100 
  #add column for standard scaling
  scaling <- c(1,1,1,1,1,0,0,0,0,0)
  data_heatmap_matrix <- cbind(data_heatmap_matrix_t,scaling)
  
  title <- paste(DMR_genename,"_", "meanDiff=", DMR_meanDiff,"_distTSS=",DMR_dist_TSS, sep = "")
  
  pheatmap(data_heatmap_matrix_t, cluster_rows = FALSE, legend = TRUE,
           treeheight_row = 0, treeheight_col = 30, show_rownames = TRUE, cluster_cols = FALSE, 
           scale = "none", border_color = "black", cellwidth = 10, main = title,
           fontsize = 6, fontsize_row = 8, fontsize_col = 8,
           cellheigth = 10, color = colorRampPalette(c("yellow", "green","blue"), space="rgb")(128))
  
  
}
dev.off()

#####
#Plot DMR meaned tracks per sample group
#####

#Get info for main/title
CpGs_DMRs_i <- CpGs_DMRs_threshed[[4]]
data_heatmap <- CpGs_DMRs_i
data_heatmap <- unique(data_heatmap)
DMR_identifier_i <- as.character(CpGs_DMRs_i$DMR_identifier)[1]
DMR_meanDiff <- as.character(Bsmooth_intel[Bsmooth_intel$DMR_identifier == DMR_identifier_i,5])
DMR_genename <- as.character(Bsmooth_intel[Bsmooth_intel$DMR_identifier == DMR_identifier_i,2])
DMR_dist_TSS <- as.character(Bsmooth_intel[Bsmooth_intel$DMR_identifier == DMR_identifier_i,3])

#get CpG positions
print(i)
rownames(data_heatmap) <- as.character(data_heatmap$CpG_start)

#make matrix
data_heatmap <- data_heatmap[,8:17]
#mean groups
pLN_GF <- rowMeans(data_heatmap[,1:2])
mLN_GF <- rowMeans(data_heatmap[,3:4])
pLN_SPF <- rowMeans(data_heatmap[,5:7])
mLN_SPF <- rowMeans(data_heatmap[,8:10])

# bind and name
data_heatmap_mean <- cbind(pLN_GF, mLN_GF, pLN_SPF, mLN_SPF)
rownames(data_heatmap_mean) <- as.character(data_heatmap$CpG_start)
  
min(data_heatmap_mean)
max(data_heatmap_mean)
data_heatmap_matrix <- as.matrix(data_heatmap_mean)
data_heatmap_matrix_t <- t(data_heatmap_matrix)

#replace NA with 100, results in white squares in heatmap
data_heatmap_matrix[is.na(data_heatmap_matrix_t)] <- 100 
#add column for standard scaling
scaling <- c(1,1,0,0)
data_heatmap_matrix <- cbind(data_heatmap_matrix_t,scaling)

title <- paste(DMR_genename,"_", "meanDiff=", DMR_meanDiff,"_distTSS=",DMR_dist_TSS, sep = "")

pheatmap(data_heatmap_matrix_t, cluster_rows = FALSE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 30, show_rownames = TRUE, cluster_cols = FALSE, 
         scale = "none", border_color = "black", cellwidth = 10, main = title,
         fontsize = 12, fontsize_row = 8, fontsize_col = 8,
         cellheigth = 10, color = colorRampPalette(c("yellow", "green","blue"), space="rgb")(128))
