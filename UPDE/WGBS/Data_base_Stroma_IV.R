################
#Database: Stromal cells
################
#1) Methylome data
#a) Bsmooth analysis
#b) CpGs per DMR analysis
#2) Transcriptome Analysis

########
#Libraries
########
library(GenomicRanges)
require(data.table)
require(ggplot2)
#require(ggfortify)
require(limma)
require(DESeq2)
library(pheatmap)
library(dplyr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ChIPpeakAnno)
library(org.Mm.eg.db)
library(biomaRt)
library(foreign)
library(GO.db)
library(topGO)
library(GenomicFeatures)

#####
#PATHs
#####
output_beds_for_motifs <-"/home/pezoldt/NAS2/pezoldt/Analysis/WGBS/FSC/Motif/Regions_of_Interest"
path_RNAseq_DESeq2 <- "/home/pezoldt/NAS2/pezoldt/Data/RNAseq/2017_FSC_LEC_BEC_mLN_pLN_GF_SPF/FSC"
Output_DARs_DMRs <- "/home/pezoldt/NAS2/pezoldt/Analysis/01_Integrated/DMRs_DARs"
Output_DARs <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/DARs"
Output_DMRs_DEGs <- "/home/pezoldt/NAS2/pezoldt/Analysis/01_Integrated/DMRs_DEGs"

#Param
log2FC_RNA = 1.0
padj = 0.05

#Global variables
#Param
padj = 0.05
log2FC_ATAC = 1.0

######################
#Bsmooth
######################
#load
DMR_mLN_SPF_GF <- read.delim2("/home/pezoldt/NAS2/pezoldt/Data/WGBS/WGBS_FSC/BSmooth_DMRs_mLN_SPF_vs_mLN_GF.txt")
DMR_pLN_SPF_GF <- read.delim2("/home/pezoldt/NAS2/pezoldt/Data/WGBS/WGBS_FSC/BSmooth_DMRs_pLN_SPF_vs_pLN_GF.txt")
DMR_SPF_mLN_pLN <- read.delim2("/home/pezoldt/NAS2/pezoldt/Data/WGBS/WGBS_FSC/BSmooth_DMRs_mLN_SPF_vs_pLN_SPF.txt")
DMR_GF_mLN_pLN <- read.delim2("/home/pezoldt/NAS2/pezoldt/Data/WGBS/WGBS_FSC/BSmooth_DMRs_pLN_GF_vs_mLN_GF.txt")
#CpGs per DMR
CpGs_per_DMR <- read.csv("/home/pezoldt/NAS2/pezoldt/Data/WGBS/WGBS_FSC/CpGs_in_DMRs.csv")


#####
#Calculate
#####
#split according to DMR_identifier
#for every DMR you get a table
CpGs_DMRs <- split(CpGs_per_DMR, CpGs_per_DMR$DMR_identifier)
Calc_CpGs_DMRs <- lapply(CpGs_DMRs, function(x){
  DMR_mean_replicate <- colMeans(x[,8:17])
})
#make dataframe
df_DMRs <- data.frame(t(sapply(Calc_CpGs_DMRs,c)))
#kick NAs, are a result as soon as one CpG has no/insufficient coverage
df_DMRs_NA_replicates <- na.omit(df_DMRs)

#SPF
CpGs_per_DMR_SPF <- subset(CpGs_per_DMR, DMR_comp_origin == "mLN_SPF_vs_pLN_SPF")
CpGs_DMRs_SPF <- split(CpGs_per_DMR_SPF, CpGs_per_DMR$DMR_identifier)
Calc_CpGs_DMRs_SPF <- lapply(CpGs_DMRs_SPF, function(x){
  DMR_mean_replicate <- colMeans(x[,8:17])
})
#make dataframe
df_DMRs_SPF <- data.frame(t(sapply(Calc_CpGs_DMRs_SPF,c)))
#kick NAs, are a result as soon as one CpG has no/insufficient coverage
df_DMRs_SPF_NA_replicates <- na.omit(df_DMRs_SPF)

#GF
CpGs_per_DMR_GF <- subset(CpGs_per_DMR, DMR_comp_origin == "pLN_GF_vs_mLN_GF") 
CpGs_DMRs_GF <- split(CpGs_per_DMR_GF, CpGs_per_DMR$DMR_identifier)
Calc_CpGs_DMRs_GF <- lapply(CpGs_DMRs_GF, function(x){
  DMR_mean_replicate <- colMeans(x[,8:17])
})
#make dataframe
df_DMRs_GF <- data.frame(t(sapply(Calc_CpGs_DMRs_GF,c)))
#kick NAs, are a result as soon as one CpG has no/insufficient coverage
df_DMRs_GF_NA_replicates <- na.omit(df_DMRs_GF)

#####
#Calculate Methylation values for each DMR for sample
#####
Calc_CpGs_DMRs <- lapply(CpGs_DMRs, function(x){
  DMR_mean_replicate <- colMeans(x[,8:17])
  DMR_mean_pLN_GF <- mean(DMR_mean_replicate[1:2])
  DMR_mean_mLN_GF <- mean(DMR_mean_replicate[3:4])
  DMR_mean_pLN_SPF <- mean(DMR_mean_replicate[5:7])
  DMR_mean_mLN_SPF <- mean(DMR_mean_replicate[8:10])
  DMR_mean_samples <- c(DMR_mean_pLN_GF,DMR_mean_mLN_GF,DMR_mean_pLN_SPF,DMR_mean_mLN_SPF)
})

Calc_CpGs_DMRs_SPF <- lapply(CpGs_DMRs, function(x){
  DMR_mean_replicate <- colMeans(x[,8:17])
  DMR_mean_pLN_SPF <- mean(DMR_mean_replicate[5:7])
  DMR_mean_mLN_SPF <- mean(DMR_mean_replicate[8:10])
  DMR_mean_samples <- c(DMR_mean_pLN_SPF,DMR_mean_mLN_SPF)
})

Calc_CpGs_DMRs_GF <- lapply(CpGs_DMRs, function(x){
  DMR_mean_replicate <- colMeans(x[,8:17])
  DMR_mean_pLN_GF <- mean(DMR_mean_replicate[1:2])
  DMR_mean_mLN_GF <- mean(DMR_mean_replicate[3:4])
  DMR_mean_samples <- c(DMR_mean_pLN_GF,DMR_mean_mLN_GF)
})

#make dataframe
df_DMRs <- data.frame(t(sapply(Calc_CpGs_DMRs,c)))
colnames(df_DMRs) <- c("pLN_GF", "mLN_GF", "pLN_SPF", "mLN_SPF")

#kick NAs, are a result as soon as one CpG has no/insufficient coverage
df_DMRs_NA_samples <- na.omit(df_DMRs)

#####
#Calculate Methylation values for each DMR for replicates
#####
Calc_CpGs_DMRs <- lapply(CpGs_DMRs, function(x){
  DMR_mean_replicate <- colMeans(x[,8:17])
})
#make dataframe
df_DMRs <- data.frame(t(sapply(Calc_CpGs_DMRs,c)))
#kick NAs, are a result as soon as one CpG has no/insufficient coverage
df_DMRs_NA_replicates <- na.omit(df_DMRs)

# Make heatmap
data_heatmap <- df_DMRs_NA_replicates
#SPF
data_heatmap <- df_DMRs_SPF_NA_replicates[,5:10]
#GF
data_heatmap <- df_DMRs_GF_NA_replicates[,1:4]
#GF_SPF only
data_heatmap <- rbind(df_DMRs_SPF_NA_replicates,df_DMRs_GF_NA_replicates)
#data_heatmap <- data_heatmap[2:nrow(data_heatmap),])
min(data_heatmap)
max(data_heatmap)
data_heatmap_matrix <- as.matrix(data_heatmap)

pheatmap(data_heatmap_matrix, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 30, show_rownames = FALSE, cluster_cols = TRUE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("yellow", "green","blue"), space="rgb")(128))

#####
#Process Bsmooth results for Granges
#####
#rename group.1 group.2
colnames(DMR_mLN_SPF_GF) <- c("X","chr","start","end","n","width","areaStat","meanDiff","meth_mLN_GF",       
                        "meth_mLN_SPF","direction","nearest_tss","distance_to_nearest_tss","nearest_gene_id",
                        "nearest_gene_name","nearest_gene_desc","preceding_tss","following_tss")
colnames(DMR_pLN_SPF_GF) <- c("X","chr","start","end","n","width","areaStat","meanDiff","meth_pLN_GF",       
                          "meth_pLN_SPF","direction","nearest_tss","distance_to_nearest_tss","nearest_gene_id",
                          "nearest_gene_name","nearest_gene_desc","preceding_tss","following_tss")
colnames(DMR_SPF_mLN_pLN) <- c("X","chr","start","end","n","width","areaStat","meanDiff","meth_pLN_SPF",       
                          "meth_mLN_SPF","direction","nearest_tss","distance_to_nearest_tss","nearest_gene_id",
                          "nearest_gene_name","nearest_gene_desc","preceding_tss","following_tss")
colnames(DMR_GF_mLN_pLN) <- c("X","chr","start","end","n","width","areaStat","meanDiff","meth_mLN_GF",       
                          "meth_pLN_GF","direction","nearest_tss","distance_to_nearest_tss","nearest_gene_id",
                          "nearest_gene_name","nearest_gene_desc","preceding_tss","following_tss")

#add columns indicating the comparison results
DMR_mLN_SPF_GF[,"Bsmooth_comp"] <- "mLN_SPF_GF"
DMR_pLN_SPF_GF[,"Bsmooth_comp"] <- "pLN_SPF_GF"
DMR_SPF_mLN_pLN[,"Bsmooth_comp"] <- "SPF_mLN_pLN"
DMR_GF_mLN_pLN[,"Bsmooth_comp"] <- "GF_mLN_pLN"
Bsmooth_list <- list(DMR_mLN_SPF_GF,DMR_pLN_SPF_GF, DMR_SPF_mLN_pLN, DMR_GF_mLN_pLN)
#reduce data content
for(i in 1:length(Bsmooth_list)){
  data_i <- Bsmooth_list[[i]]
  reduced_i <- data_i[,c(1:7,13,15,19)]
  Bsmooth_list[[i]] <- reduced_i
}
rm(data_i, reduced_i)
#add identifier for eadch DMR
Bsmooth_list <- lapply(Bsmooth_list, function(x) {
  k <- cbind(x, paste(x$chr, x$start, x$end, sep = "_"))
  colnames(k)[ncol(k)] <- "DMR_identifier"
  return(k)
})

#Get identifiers of unique DMRs
gr_DMR_mLN_SPF_GF <- with(DMR_mLN_SPF_GF, GRanges(chr, IRanges(start = start, end = end), 
                                          Gene_Symbol = nearest_gene_name, areastat = areaStat, width_dmr = width, mean_Diff = meanDiff,
                                          mLN_GF = meth_mLN_GF, mLN_SPF = meth_mLN_SPF,
                                          distance_TSS = distance_to_nearest_tss))
gr_DMR_pLN_SPF_GF <- with(DMR_pLN_SPF_GF, GRanges(chr, IRanges(start = start, end = end), 
                                          Gene_Symbol = nearest_gene_name, areastat = areaStat, width_dmr = width, mean_Diff = meanDiff,
                                          pLN_GF = meth_pLN_GF, pLN_SPF = meth_pLN_SPF,
                                          distance_TSS = distance_to_nearest_tss))
gr_DMR_SPF_mLN_pLN <- with(DMR_SPF_mLN_pLN, GRanges(chr, IRanges(start = start, end = end), 
                                            Gene_Symbol = nearest_gene_name, areastat = areaStat, width_dmr = width, mean_Diff = meanDiff,
                                            pLN_SPF = meth_pLN_SPF, mLN_SPF = meth_mLN_SPF,
                                            distance_TSS = distance_to_nearest_tss))
gr_DMR_GF_mLN_pLN <- with(DMR_GF_mLN_pLN, GRanges(chr, IRanges(start = start, end = end), 
                                          Gene_Symbol = nearest_gene_name, areastat = areaStat, width_dmr = width, mean_Diff = meanDiff,
                                          mLN_GF = meth_mLN_GF, pLN_GF = meth_pLN_GF,
                                          distance_TSS = distance_to_nearest_tss))
###Identify overlap location
hits_location <- findOverlaps(gr_DMR_SPF_mLN_pLN, gr_DMR_GF_mLN_pLN, minoverlap = 50)
#Get index
query_Hits <- queryHits(hits_location)
subject_Hits <- subjectHits(hits_location)
#Get metadata
overlap_location_SPF <- gr_DMR_SPF_mLN_pLN[query_Hits]
overlap_location_GF <- gr_DMR_GF_mLN_pLN[subject_Hits]
#nonoverlap
#index vector for DMRs not overlapping with location SPF
index_DMR_SPF_mLN_pLN <- c(1:nrow(DMR_SPF_mLN_pLN))
remove_DMR_SPF_mLN_pLN <- query_Hits
non_overlap_SPF <- index_DMR_SPF_mLN_pLN[! index_DMR_SPF_mLN_pLN %in% remove_DMR_SPF_mLN_pLN]
#index vector for DMRs not overlapping with location GF
index_DMR_GF_mLN_pLN <- c(1:nrow(DMR_GF_mLN_pLN))
remove_DMR_GF_mLN_pLN <- query_Hits
non_overlap_GF <- index_DMR_GF_mLN_pLN[! index_DMR_GF_mLN_pLN %in% remove_DMR_GF_mLN_pLN]
#Granges object for non_overlap
nonoverlap_location_SPF <- gr_DMR_SPF_mLN_pLN[non_overlap_SPF]
nonoverlap_location_GF <- gr_DMR_GF_mLN_pLN[non_overlap_GF]
#make data fram
overlap_location_DMR_SPF_mLN_pLN <- as.data.frame(overlap_location_SPF)
nonoverlap_location_DMR_SPF_mLN_pLN <- as.data.frame(nonoverlap_location_SPF)
nonoverlap_location_DMR_GF_mLN_pLN <- as.data.frame(nonoverlap_location_GF)

###Identify overlap commensal

#Note: no overlap
#make data fram
nonoverlap_commensal_DMR_mLN_SPF_GF <- as.data.frame(gr_DMR_mLN_SPF_GF)
nonoverlap_commensal_DMR_pLN_SPF_GF <- as.data.frame(gr_DMR_pLN_SPF_GF)

#merge column chr stop end of each dataframe to generate identifier for DMR CpG average heatmap
o_seq <- with(overlap_location_DMR_SPF_mLN_pLN, paste(seqnames, "_", start, "_", end, sep = ""))
non_SPF <- with(nonoverlap_location_DMR_SPF_mLN_pLN, paste(seqnames, "_", start, "_", end, sep = ""))
non_GF <- with(nonoverlap_location_DMR_GF_mLN_pLN, paste(seqnames, "_", start, "_", end, sep = ""))
non_mLN <- with(nonoverlap_commensal_DMR_mLN_SPF_GF, paste(seqnames, "_", start, "_", end, sep = ""))
non_pLN <- with(nonoverlap_commensal_DMR_pLN_SPF_GF, paste(seqnames, "_", start, "_", end, sep = ""))
#all hits
all_DMRs <- c(o_seq, non_SPF, non_GF, non_mLN, non_pLN)

#####
#CpGs per DMRs
#####
#Gene assciated with DMR
split_CpGs_per_DMR <- split(CpGs_per_DMR, CpGs_per_DMR$DMR_identifier)
Calc_CpGs_DMRs <- lapply(split_CpGs_per_DMR, function(x){
  DMR_mean_replicate <- colMeans(x[,8:17])
  DMR_mean_pLN_GF <- mean(DMR_mean_replicate[1:2])
  DMR_mean_mLN_GF <- mean(DMR_mean_replicate[3:4])
  DMR_mean_pLN_SPF <- mean(DMR_mean_replicate[5:7])
  DMR_mean_mLN_SPF <- mean(DMR_mean_replicate[8:10])
  DMR_mean_samples <- c(DMR_mean_pLN_GF,DMR_mean_mLN_GF,DMR_mean_pLN_SPF,DMR_mean_mLN_SPF)
})
df_DMRs <- data.frame(t(sapply(Calc_CpGs_DMRs,c)))
DMR_identifiers <- rownames(df_DMRs)
rownames(df_DMRs) <- NULL
df_DMRs <- cbind(DMR_identifiers,df_DMRs)
colnames(df_DMRs) <- c("DMR_identifier","pLN_GF", "mLN_GF", "pLN_SPF", "mLN_SPF")

#####
#CpG methylation and Bsmooth Merger
#####
#use only data where DMRs are unique
df_DMRs_unique <- subset(df_DMRs, DMR_identifiers %in% all_DMRs)
Bsmooth_list_unique <- lapply(Bsmooth_list, function(x) {
  k <- subset(x, DMR_identifier %in% all_DMRs)
  return(k)
})
#merge methylation over all samples with Bsmooth data
Bsmooth_unique_CpGmeaned <- lapply(Bsmooth_list_unique, function(x) {
  k <- merge(x, df_DMRs, by = "DMR_identifier", all = FALSE)
  return(k)
})

#make dataframe
Bsmooth_CpGmeaned <- rbind(Bsmooth_unique_CpGmeaned[[1]], Bsmooth_unique_CpGmeaned[[2]])
Bsmooth_CpGmeaned <- rbind(Bsmooth_CpGmeaned, Bsmooth_unique_CpGmeaned[[3]])
Bsmooth_CpGmeaned <- rbind(Bsmooth_CpGmeaned, Bsmooth_unique_CpGmeaned[[4]])

#Eliminate DMRs that 
#lack gene name
#lack methylation value
Bsmooth_CpGmeaned_NA <- na.omit(Bsmooth_CpGmeaned)
Bsmooth_CpGmeaned_NA_complete <- Bsmooth_CpGmeaned_NA[!(Bsmooth_CpGmeaned_NA$nearest_gene_name == ""), ]
#Bsmooth_CpGmeaned_NA_complete <- Bsmooth_CpGmeaned_NA[!(Bsmooth_CpGmeaned_NA$nearest_gene_name == "NaN"), ]

#Mean genes with more than one DMR
#for all numeric parameters
split_Bsmooth_CpGmeaned_NA <- split(Bsmooth_CpGmeaned_NA_complete, as.character(Bsmooth_CpGmeaned_NA_complete$nearest_gene_name))
data_frame <- data.frame(nearest_gene_name = character(),n_CpGs = numeric(),width = numeric(),
                        areaStat = numeric(), distance_to_nearest_tss = numeric(),
                        pLN_GF = numeric(),mLN_GF = numeric(),pLN_SPF = numeric(),
                        mLN_SPF = numeric())

DMR_mean_Bsmooth_CpGmeaned_NA <- lapply(split_Bsmooth_CpGmeaned_NA, function(x){
  DMR_mean_replicate <- colMeans(x[,c(6:9,12:15)])
  #print(DMR_mean_replicate)
  nearest_gene_name <- as.character(unique(x$nearest_gene_name))
  data_frame_row <- data.frame(nearest_gene_name = as.character(nearest_gene_name),
                               n_CpGs = as.numeric(DMR_mean_replicate[1]),
                               width = as.numeric(DMR_mean_replicate[2]),
                               areaStat = as.numeric(DMR_mean_replicate[3]),
                               distance_to_nearest_tss = as.numeric(DMR_mean_replicate[4]),
                               pLN_GF = as.numeric(DMR_mean_replicate[5]),
                               mLN_GF = as.numeric(DMR_mean_replicate[6]),
                               pLN_SPF = as.numeric(DMR_mean_replicate[7]),
                               mLN_SPF = as.numeric(DMR_mean_replicate[8]))
})
df_DMR_mean_Bsmooth_CpGmeaned_NA <- do.call("rbind", DMR_mean_Bsmooth_CpGmeaned_NA)

#Calculate delta Methylation
dMeth_DMR_mLN_SPF_GF <- df_DMR_mean_Bsmooth_CpGmeaned_NA$mLN_SPF - df_DMR_mean_Bsmooth_CpGmeaned_NA$mLN_GF
dMeth_DMR_pLN_SPF_GF <- df_DMR_mean_Bsmooth_CpGmeaned_NA$pLN_SPF - df_DMR_mean_Bsmooth_CpGmeaned_NA$pLN_GF
dMeth_DMR_SPF_mLN_pLN <- df_DMR_mean_Bsmooth_CpGmeaned_NA$mLN_SPF - df_DMR_mean_Bsmooth_CpGmeaned_NA$pLN_SPF
dMeth_DMR_GF_mLN_pLN <- df_DMR_mean_Bsmooth_CpGmeaned_NA$mLN_GF - df_DMR_mean_Bsmooth_CpGmeaned_NA$pLN_GF
df_DMR_mean_Bsmooth_CpGmeaned_dMeth <- cbind(df_DMR_mean_Bsmooth_CpGmeaned_NA, dMeth_DMR_mLN_SPF_GF)
df_DMR_mean_Bsmooth_CpGmeaned_dMeth <- cbind(df_DMR_mean_Bsmooth_CpGmeaned_dMeth, dMeth_DMR_pLN_SPF_GF)
df_DMR_mean_Bsmooth_CpGmeaned_dMeth <- cbind(df_DMR_mean_Bsmooth_CpGmeaned_dMeth, dMeth_DMR_SPF_mLN_pLN)
df_DMR_mean_Bsmooth_CpGmeaned_dMeth <- cbind(df_DMR_mean_Bsmooth_CpGmeaned_dMeth, dMeth_DMR_GF_mLN_pLN)
df_DMR_mean_Bsmooth_CpGmeaned_dMeth_TSS <- subset(df_DMR_mean_Bsmooth_CpGmeaned_dMeth, distance_to_nearest_tss < 10000)

#Demeth mLN
Demeth_mLN <- subset(df_DMR_mean_Bsmooth_CpGmeaned_dMeth_TSS, dMeth_SPF_mLN_pLN < -0.15)
Demeth_pLN <- subset(df_DMR_mean_Bsmooth_CpGmeaned_dMeth_TSS, dMeth_SPF_mLN_pLN > 0.15)
Demeth_pLN_mLN <- subset(df_DMR_mean_Bsmooth_CpGmeaned_dMeth_TSS, abs(dMeth_SPF_mLN_pLN) > 0.15)
#Heatmap
data <- Demeth_pLN_mLN
data_heatmap_matrix <- as.matrix(data[,c("pLN_GF","mLN_GF","pLN_SPF","mLN_SPF")])


pheatmap(data_heatmap_matrix, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 30, show_rownames = FALSE, cluster_cols = TRUE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("yellow", "green","blue"), space="rgb")(128))


####################
#RNAseq analysis low Input
####################
#Load data
setwd(path_RNAseq_DESeq2)
names = list.files(pattern="*.csv")
myfiles = lapply(names, read.csv)

#merrge Pair-wise comparisons
names <- unlist(strsplit(names, split='DESeq2_genes_diffexp_by_fc_', fixed=TRUE))
names <- unlist(strsplit(names, split='.csv', fixed=TRUE))
#vector of length 2 x
names_group <- unlist(strsplit(names, split='_vs_', fixed=TRUE))
#take even and uneven elements
#uneven
names_group_1 <- names_group[seq(1, length(names_group), 2)]
#even
names_group_2 <- names_group[seq(2, length(names_group), 2)]
#new names with propper divivsion order
names <- paste(names_group_2, "_VS_", names_group_1, sep = "")
#get relevant columns from list of tables (RPKM, FC, padj)
#assign list names
names(myfiles) <- names

#get list name and attach to column names
list_all <- lapply(seq_along(myfiles),function(i){
  name_i <- names(myfiles[i])
  #colnames(myfiles[[i]])
  table_i <- cbind(myfiles[[i]][,c("GeneSymbol","log2FoldChange","padj")],myfiles[[i]][,grepl("RPKMcounts", colnames(myfiles[[i]]))])
  colnames(table_i)[2] <- paste("log2FC_", name_i, sep = "")
  colnames(table_i)[3] <- paste("padj_", name_i, sep = "")
  table_i
}
)

#make table from list
DF <- list_all[[1]]
for (.df in list_all) {
  DF <-merge(DF,.df,by = "GeneSymbol", all=T)
  DF <- DF[!duplicated(DF$GeneSymbol),]
}
data_all <- cbind(DF$GeneSymbol, 
                  #DF[,grepl(c(".x"), colnames(DF))],
                  DF[,grepl(c("padj"), colnames(DF))],
                  DF[,grepl(c("log2FC"), colnames(DF))])
colnames(data_all)[1] <- "GeneSymbol"
#Eliminate columns
data_all <- data_all[,c(1,2,4:7,9:ncol(data_all))]
data_all <- data_all[!(is.na(data_all$GeneSymbol)),]
data_all <- data_all[!duplicated(data_all$GeneSymbol),]
#reorder
data_all <- data_all[,c(1,6:9,2:5)]
colnames(data_all) <- c("GeneSymbol",
                        "GF_mLN_pLN_log2FoldChange","mLN_SPF_GF_log2FoldChange",
                        "SPF_mLN_pLN_log2FoldChange","pLN_SPF_GF_log2FoldChange",
                        "GF_mLN_pLN_padj","mLN_SPF_GF_padj",
                        "SPF_mLN_pLN_padj","pLN_SPF_GF_padj")
#Invert FC to align to renaming of columns
data_all_test <- data_all[,2:6] * (-1)
#Reorder dataframe
data_all <- data_all[,c("GeneSymbol",
                        "SPF_mLN_pLN_log2FoldChange","GF_mLN_pLN_log2FoldChange",
                        "mLN_SPF_GF_log2FoldChange","pLN_SPF_GF_log2FoldChange",
                        "SPF_mLN_pLN_padj","GF_mLN_pLN_padj",
                        "mLN_SPF_GF_padj","pLN_SPF_GF_padj")]
#store
#Reverse Foldchange to fit to name
data_all$SPF_mLN_pLN_log2FoldChange <- data_all$SPF_mLN_pLN_log2FoldChange * (-1)
data_all$GF_mLN_pLN_log2FoldChange <- data_all$GF_mLN_pLN_log2FoldChange * (-1)


diff_Genes <- data_all


####################
#ATACseq
####################
# Input required: Set directory
# 1) Counts per peak per replicate
path_input <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/homer/Overlap_Group_Merged/Run_3_Outlier_excluded"
#path_input <- "/Volumes/Pezoldt_HD_4TB/UPDE_Transfer/ATAC_FSC_all/homer/Overlap_Group_Merged/Run_3_Outlier_excluded"
# 2) Common Peak regions
path_common_peak_regions <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/peaks/broad/MACS2_broad/Run_3/ATAC_FSC_all_broad_merged_peaks.bed"
#path_common_peak_regions <- "/Volumes/Pezoldt_HD_4TB/UPDE_Transfer/ATAC_FSC_all/peaks/broad/MACS2_broad/Run_3/ATAC_FSC_all_broad_merged_peaks.bed"
# 3) Output path
#path_output <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/DESeq2/Run_3_Outlier_excluded"
#path_output_Motif <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/Motif/Regions_of_Interest"
# Input required: 
name = "ATAC_FSC_all"
paste(path_input, "/",name,".txt",sep="")
#read count table
counts <- fread(paste(path_input, "/",name,".txt",sep=""))

#####
#QC
#####
# peaks never found
which(rowSums(counts)==0) ## none (expected)
# peaks found in < 50 % of samples
# Input required: Number of samples
number_of_samples <- 13

# peaks with zero count per sample 
zero <- colSums(counts[,-1]==0)
zero.df <- data.frame(colnames(counts[,-1]), zero)
colnames(zero.df) <- c("cell", "zero")

ggplot(zero.df) + geom_col(aes(x=(reorder(zero.df$cell, zero.df$zero)), y= zero.df$zero)) +
  theme(axis.text.x = element_text(angle = 50, hjust=1))

#Store rownames
row_counts <- rownames(counts)

# transform into count matrix
counts.mat <-(counts[,-1])
counts.mat <- as.matrix(counts.mat)
rownames(counts.mat) <- row_counts

#####
#DEseq2 ATAC-seq
#####
#Define conditons 
condition <- factor(substr(colnames(counts.mat),1,nchar(colnames(counts.mat))-1))
dds <- DESeqDataSetFromMatrix(counts.mat, DataFrame(condition),~condition)
dds.factors <- estimateSizeFactors(dds)

#Perform DESeq2
dds <- DESeq(dds)
res <- results( dds )

#Perform RPKM calculation
peak_regions <- read.delim(path_common_peak_regions, header = FALSE)
colnames(peak_regions) <- c("chr","start","end")
#make Granges object all samples are computed on the same Granges object thus all have the same 
gr_peak_regions <- toGRanges(peak_regions, seqnames.field=c("chr"), start.field="start", end.field="end")
rowRanges(dds) <- gr_peak_regions
fpkm_dds <- fpkm(dds)
fpm_dds <- fpm(dds)

#generate dds-resuls objects for relevant comparisons
mLN_SPF_GF <- results(dds, contrast=c("condition","ATAC_mLNSPF","ATAC_mLNGF"))
pLN_SPF_GF <- results(dds, contrast=c("condition","ATAC_pLNSPF","ATAC_pLNGF"))
SPF_mLN_pLN <- results(dds, contrast=c("condition","ATAC_mLNSPF","ATAC_pLNSPF"))
GF_mLN_pLN <- results(dds, contrast=c("condition","ATAC_mLNGF","ATAC_pLNGF"))

#QC
#MA
par(mfrow = c(2,2))
plotMA(mLN_SPF_GF, ylim = c(-6, 6), main = "mLN_SPF_GF")
plotMA(pLN_SPF_GF, ylim = c(-6, 6), main = "pLN_SPF_GF" )
plotMA(SPF_mLN_pLN, ylim = c(-6, 6), main = "SPF_mLN_pLN" )
plotMA(GF_mLN_pLN, ylim = c(-6, 6), main = "GF_mLN_pLN" )

#DispEsts
plotDispEsts( dds, ylim = c(1e-6, 1e1), main = "All" )


#Hits pValue
par(mfrow = c(2,2))
hist( mLN_SPF_GF$pvalue, breaks=20, col="grey", main = "mLN_SPF_GF"  )
hist( pLN_SPF_GF$pvalue, breaks=20, col="grey", main = "pLN_SPF_GF"  )
hist( SPF_mLN_pLN$pvalue, breaks=20, col="grey", main = "SPF_mLN_pLN"  )
hist( GF_mLN_pLN$pvalue, breaks=20, col="grey", main = "GF_mLN_pLN"  )

par(mfrow = c(1,1))
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

#Get log2FC and pValue for all peaks across all chosen comparisons
DAR_mLN_SPF_GF <- as.data.frame(results(dds, contrast=c("condition","ATAC_mLNSPF","ATAC_mLNGF"), format = c("DataFrame")))
DAR_mLN_SPF_GF <- DAR_mLN_SPF_GF[,c("log2FoldChange","padj")]
colnames(DAR_mLN_SPF_GF) <- c(paste("log2FC_","mLN_SPF_GF",sep=""),paste("padj_","mLN_SPF_GF",sep=""))

DAR_pLN_SPF_GF <- as.data.frame(results(dds, contrast=c("condition","ATAC_pLNSPF","ATAC_pLNGF"), format = c("DataFrame")))
DAR_pLN_SPF_GF <- DAR_pLN_SPF_GF[,c("log2FoldChange","padj")]
colnames(DAR_pLN_SPF_GF) <- c(paste("log2FC_","pLN_SPF_GF",sep=""),paste("padj_","pLN_SPF_GF",sep=""))

DAR_SPF_mLN_pLN <- as.data.frame(results(dds, contrast=c("condition","ATAC_mLNSPF","ATAC_pLNSPF"), format = c("DataFrame")))
DAR_SPF_mLN_pLN <- DAR_SPF_mLN_pLN[,c("log2FoldChange","padj")]
colnames(DAR_SPF_mLN_pLN) <- c(paste("log2FC_","SPF_mLN_pLN",sep=""),paste("padj_","SPF_mLN_pLN",sep=""))


DAR_GF_mLN_pLN <- as.data.frame(results(dds, contrast=c("condition","ATAC_mLNGF","ATAC_pLNGF"), format = c("DataFrame")))
DAR_GF_mLN_pLN <- DAR_GF_mLN_pLN[,c("log2FoldChange","padj")]
colnames(DAR_GF_mLN_pLN) <- c(paste("log2FC_","GF_mLN_pLN",sep=""),paste("padj_","GF_mLN_pLN",sep=""))

#Table of all conditions
DARs <- cbind(id = rownames(DAR_mLN_SPF_GF), DAR_mLN_SPF_GF, DAR_pLN_SPF_GF, DAR_SPF_mLN_pLN, DAR_GF_mLN_pLN)
rownames(DARs) <- c()

#####
#Associate genes with DESeq results
#####
#Load the genomic regions associated with ids
# Common peak track for Experiment is annotated using homer
regions <- read.delim(path_common_peak_regions, header = FALSE)
colnames(regions) <- c("chr","start","end")
DARs_regions <- cbind(regions, DARs)
#replace chromosome numbers with 1 -> chr1
DARs_regions$chr <- paste("chr",DARs_regions$chr, sep="")
rownames(DARs_regions) <- c()
#DARs with FPKM
DARs_regions <- cbind(DARs_regions,fpkm_dds)

#####
#Peak/gene location
#####
library("GenomicFeatures")
#All features
mm10_TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#TSS sites BioMart
mm10 = useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl") 
TSS.mouse.mm10 = getAnnotation(mart=mm10, featureType="TSS")
Exon.mouse.mm10 = getAnnotation(mart=mm10, featureType="Exon")
#mm10 whole database
#annoData_exon <- toGRanges(mm10_TxDb, feature=c("exon"))
#Grange object from peak file
gr_DARs_regions <- toGRanges(DARs_regions, names = DARs_regions$id)
#Annotate Peaks
gr_DARs_regions_anno_TSS <- annotatePeakInBatch(gr_DARs_regions, AnnotationData=TSS.mouse.mm10)
#add gene name
gr_DARs_regions_anno_TSS <- addGeneIDs(annotatedPeak=gr_DARs_regions_anno_TSS, 
                                       feature_id_type="ensembl_gene_id",
                                       orgAnn="org.Mm.eg.db", 
                                       IDs2Add="symbol")
#Generate dataframe that contains DAR location, FCs and padj.
DARs_features <- as.data.frame(gr_DARs_regions_anno_TSS)
#Number of features
length(DARs_features$insideFeature[DARs_features$insideFeature == "upstream"])

#Commensal specific DARs
DARs_features_mLN_Commensal <- subset(DARs_features, abs(log2FC_mLN_SPF_GF) >= log2FC_ATAC & padj_mLN_SPF_GF <= 0.05)
write.table(DARs_features_mLN_Commensal, paste(Output_DARs,"/DARs_features_mLN_Commensal.txt", sep = ""), dec = ".", sep = "\t")
unique(DARs_features_mLN_Commensal$symbol)
DARs_features_pLN_Commensal <- subset(DARs_features, abs(log2FC_pLN_SPF_GF) >= log2FC_ATAC & padj_pLN_SPF_GF <= 0.05)
write.table(DARs_features_pLN_Commensal, paste(Output_DARs,"/DARs_features_pLN_Commensal.txt", sep = ""), dec = ".", sep = "\t")
unique(DARs_features_pLN_Commensal$symbol)

####################
#Heatmaps for DARs
####################
#####
#SPF
#####
DARs_features_SPF <- subset(DARs_features, abs(log2FC_SPF_mLN_pLN) >= log2FC_ATAC & padj_SPF_mLN_pLN <= padj)
DARs_features_SPF_NA <- DARs_features_SPF[!is.na(DARs_features_SPF$symbol),]
#Heatmap
data_heatmap <- DARs_features_SPF_NA[,c("id","ATAC_mLNSPF1","ATAC_mLNSPF2","ATAC_mLNSPF3","ATAC_mLNSPF4",
                                        "ATAC_pLNSPF1","ATAC_pLNSPF2","ATAC_pLNSPF3")]
rownames(data_heatmap) <- data_heatmap$id
data_heatmap <- as.matrix(data_heatmap[,2:ncol(data_heatmap)])
data_heatmap[data_heatmap == 0] <- 0.4
data_heatmap <- log2(data_heatmap)

min(data_heatmap)
pheatmap(data_heatmap, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = TRUE, cellwidth = 10,
         scale = "row", border_color = "black", color = colorRampPalette(c("white","grey", "darkgreen"), space="rgb")(128),
         main = "DARs TSS associated SPF")

#####
#GF
#####
DARs_features_GF <- subset(DARs_features, abs(log2FC_GF_mLN_pLN) >= log2FC_ATAC & padj_GF_mLN_pLN <= padj)
DARs_features_GF_NA <- DARs_features_GF[!is.na(DARs_features_GF$symbol),]
#Heatmap
data_heatmap <- DARs_features_GF_NA[,c("id",
                                       "ATAC_mLNGF1","ATAC_mLNGF2","ATAC_mLNGF3",
                                       "ATAC_pLNGF1","ATAC_pLNGF3")]
rownames(data_heatmap) <- data_heatmap$id
data_heatmap <- as.matrix(data_heatmap[,2:ncol(data_heatmap)])
data_heatmap[data_heatmap == 0] <- 0.4
min(data_heatmap)
data_heatmap <- log2(data_heatmap)

min(data_heatmap)
pheatmap(data_heatmap, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = FALSE, cluster_cols = TRUE, cellwidth = 10,
         scale = "row", border_color = "black", color = colorRampPalette(c("white","grey", "darkgreen"), space="rgb")(128),
         main = "DARs TSS associated SPF")

#####
#GF & SPF
#####
DARs_features_GF_SPF <- subset(DARs_features, (abs(log2FC_GF_mLN_pLN) >= log2FC_ATAC & padj_GF_mLN_pLN <= padj) |
                                 (abs(log2FC_GF_mLN_pLN) >= log2FC_ATAC & padj_GF_mLN_pLN <= padj) )
DARs_features_GF_SPF_NA <- DARs_features_GF_SPF[!is.na(DARs_features_GF_SPF$symbol),]
#Heatmap
data_heatmap <- DARs_features_GF_SPF_NA[,c("id",
                                           "ATAC_mLNSPF1","ATAC_mLNSPF2","ATAC_mLNSPF3","ATAC_mLNSPF4",
                                           "ATAC_pLNSPF1","ATAC_pLNSPF2","ATAC_pLNSPF3",
                                       "ATAC_mLNGF1","ATAC_mLNGF2","ATAC_mLNGF3",
                                       "ATAC_pLNGF1","ATAC_pLNGF3")]
rownames(data_heatmap) <- data_heatmap$id
data_heatmap <- as.matrix(data_heatmap[,2:ncol(data_heatmap)])
data_heatmap[data_heatmap == 0] <- 0.4
min(data_heatmap)
data_heatmap <- log2(data_heatmap)

min(data_heatmap)
pheatmap(data_heatmap, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = FALSE, cluster_cols = TRUE, cellwidth = 10,
         scale = "row", border_color = "black", color = colorRampPalette(c("white","grey", "darkgreen"), space="rgb")(128),
         main = "DARs TSS associated SPF")




#####
#Meaned GF & SPF
#####
DARs_features_GF_SPF <- subset(DARs_features, (abs(log2FC_GF_mLN_pLN) >= log2FC_ATAC & padj_GF_mLN_pLN <= padj) |
                                 (abs(log2FC_SPF_mLN_pLN) >= log2FC_ATAC & padj_SPF_mLN_pLN <= padj) )
DARs_features_GF_SPF_NA <- DARs_features_GF_SPF[!is.na(DARs_features_GF_SPF$symbol),]
DARs_features_GF_SPF_NA_Prom <- subset(DARs_features_GF_SPF_NA, distancetoFeature >= -400 & distancetoFeature <= 10000)
DARs_features_GF_SPF_NA_Prom_unique <- DARs_features_GF_SPF_NA_Prom[!duplicated(DARs_features_GF_SPF_NA_Prom[,c("symbol")]),]


#Heatmap
data_heatmap <- DARs_features_GF_SPF_NA_Prom_unique[,c("symbol",
                                           "ATAC_mLNSPF1","ATAC_mLNSPF2","ATAC_mLNSPF3","ATAC_mLNSPF4",
                                           "ATAC_pLNSPF1","ATAC_pLNSPF2","ATAC_pLNSPF3",
                                           "ATAC_mLNGF1","ATAC_mLNGF2","ATAC_mLNGF3",
                                           "ATAC_pLNGF1","ATAC_pLNGF3")]

rownames(data_heatmap) <- data_heatmap$symbol
data_heatmap <- as.matrix(data_heatmap[,2:ncol(data_heatmap)])
mLN_SPF <- rowMeans(data_heatmap[,1:4])
pLN_SPF <- rowMeans(data_heatmap[,5:7])
mLN_GF <- rowMeans(data_heatmap[,8:10])
pLN_GF <- rowMeans(data_heatmap[,11:12])
data_heatmap <- cbind(mLN_SPF,pLN_SPF,mLN_GF,pLN_GF)
#data_heatmap[data_heatmap == 0] <- 0.4
min(data_heatmap)
data_heatmap <- log2(data_heatmap)

min(data_heatmap)
pheatmap(data_heatmap, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 20, show_rownames = FALSE, cluster_cols = TRUE, cellwidth = 10,
         scale = "none", border_color = "black", color = colorRampPalette(c("black","dimgrey","gray48","grey","white"), space="rgb")(128),
         main = "DARs TSS associated SPF")
outs <- pheatmap(data_heatmap, cluster_rows = TRUE, legend = TRUE,
                 treeheight_row = 0, treeheight_col = 20, show_rownames = TRUE, cluster_cols = TRUE, cellwidth = 10,
                 scale = "none", border_color = "black", color = colorRampPalette(c("black","dimgrey","gray48","grey","white"), space="rgb")(128),
                 main = "DARs TSS associated SPF")

####################
#Compare DEG vs. DMR
####################
colnames(DMR_SPF_mLN_pLN)
#significant DMRs
DMR_sig_SPF_mLN_pLN <- subset(DMR_SPF_mLN_pLN, distance_to_nearest_tss <= 10000 & n >= 3 & abs(meanDiff) >= 0.2)

DMR_sig_SPF_mLN_pLN$nearest_gene_name <- as.character(DMR_sig_SPF_mLN_pLN$nearest_gene_name) 
split_DMR_sig_SPF_mLN_pLN <- split(DMR_sig_SPF_mLN_pLN, DMR_sig_SPF_mLN_pLN$nearest_gene_name)
#Quick check if typ of methylation is unidirectional
#Note: Within TSS for mLN vs. pLN that is exclusivly the case

#Take first DMR per Gene Locus
Meaned_split_DMR_sig_SPF_mLN_pLN <- lapply(split_DMR_sig_SPF_mLN_pLN, function(x){
  mean_per_gene <- mean(x$meanDiff)
  #For adjacent data use first element
  data.frame(GeneSymbol = x$nearest_gene_name[1],
             mean_dMeth_gene = mean_per_gene,
             distance_tss = x$distance_to_nearest_tss[1])
})

Meaned_DMR_sig_SPF_mLN_pLN <- do.call(rbind.data.frame, Meaned_split_DMR_sig_SPF_mLN_pLN)

#DEGs 
diff_Genes_DEG <- subset(diff_Genes, abs(SPF_mLN_pLN_log2FoldChange) >= 1.0 & SPF_mLN_pLN_padj <= 0.05)

overlap_DMR_DEG <- subset(diff_Genes_DEG, GeneSymbol %in% Meaned_DMR_sig_SPF_mLN_pLN$GeneSymbol)
#Merge DMR and DEG table
merge_DMR_DEG <- merge(diff_Genes_DEG, Meaned_DMR_sig_SPF_mLN_pLN,
                       by = "GeneSymbol")
#Note: Not completely consistent with the dogma of Demeth makes expression go high,
#      but OK
#check quadrant
DeMeth_but_UP <- subset(merge_DMR_DEG, SPF_mLN_pLN_log2FoldChange >= 1.0 & mean_dMeth_gene <= -0.2)

plot(merge_DMR_DEG$SPF_mLN_pLN_log2FoldChange, merge_DMR_DEG$mean_dMeth_gene)
mLN_Demeth_UP <- subset(merge_DMR_DEG, mean_dMeth_gene <= -0.2 & SPF_mLN_pLN_log2FoldChange >= 1.0)
unique(mLN_Demeth_UP$GeneSymbol)
mLN_Meth_UP <- subset(merge_DMR_DEG, mean_dMeth_gene >= 0.2 & SPF_mLN_pLN_log2FoldChange >= 1.0)
pLN_Demeth_UP <- subset(merge_DMR_DEG, mean_dMeth_gene >= 0.2 & SPF_mLN_pLN_log2FoldChange <= -1.0)
unique(pLN_Demeth_UP$GeneSymbol)
pLN_Meth_UP <- subset(merge_DMR_DEG, mean_dMeth_gene <= -0.2 & SPF_mLN_pLN_log2FoldChange <= -1.0)

#Output
write.table(mLN_Demeth_UP, paste(Output_DMRs_DEGs,"/mLN_Demeth_UP.txt", sep = ""), dec = ".", sep = "\t")
write.table(mLN_Meth_UP, paste(Output_DMRs_DEGs,"/mLN_Meth_UP.txt", sep = ""), dec = ".", sep = "\t")
write.table(pLN_Demeth_UP, paste(Output_DMRs_DEGs,"/pLN_Demeth_UP.txt", sep = ""), dec = ".", sep = "\t")
write.table(pLN_Meth_UP, paste(Output_DMRs_DEGs,"/pLN_Meth_UP.txt", sep = ""), dec = ".", sep = "\t")

####################
#DMR / DAR
####################
#Keep only DARs with annotated GeneSymbol and close to the TSS -10000 +200
DARs_features_Genes <- DARs_features[!(is.na(DARs_features$symbol)),]
DARs_features_TSS_sig <- subset(DARs_features_Genes, padj_SPF_mLN_pLN <=  0.05 & abs(log2FC_SPF_mLN_pLN) >= 1.0)
#Merge DMR and DAR Table
#DAR_genes_in_DMRs <- subset(Meaned_DMR_sig_SPF_mLN_pLN, GeneSymbol %in% DARs_features_TSS_sig$symbol)
#DMRs_genes_in_DARs <- subset(DARs_features_TSS_sig, symbol %in% DMR_SPF_mLN_pLN$nearest_gene_name)
DARs_DMRs_TSS <- merge(DARs_features_TSS_sig, DAR_genes_in_DMRs,
                       by.x = "symbol", by.y = "GeneSymbol")

mLN_Demeth_Open <- subset(DARs_DMRs_TSS, mean_dMeth_gene <= -0.2 & log2FC_SPF_mLN_pLN >= 1.0)
unique(mLN_Demeth_Open$symbol)
mLN_Demeth_Closed <- subset(DARs_DMRs_TSS, mean_dMeth_gene <= -0.2 & log2FC_SPF_mLN_pLN <= -1.0)
pLN_Demeth_Closed <- subset(DARs_DMRs_TSS, mean_dMeth_gene >= 0.2 & log2FC_SPF_mLN_pLN >= 1.0)
pLN_Demeth_Open <- subset(DARs_DMRs_TSS, mean_dMeth_gene >= 0.2 & log2FC_SPF_mLN_pLN <= -1.0)
unique(pLN_Demeth_Closed$symbol)

#write tables
write.table(mLN_Demeth_Open, paste(Output_DARs_DMRs,"/mLN_Demeth_Open.txt", sep = ""), dec = ".", sep = "\t")
write.table(mLN_Demeth_Closed, paste(Output_DARs_DMRs,"/mLN_Demeth_Closed.txt", sep = ""), dec = ".", sep = "\t")
write.table(pLN_Demeth_Open, paste(Output_DARs_DMRs,"/pLN_Demeth_Open.txt", sep = ""), dec = ".", sep = "\t")
write.table(pLN_Demeth_Closed, paste(Output_DARs_DMRs,"/pLN_Demeth_Closed.txt", sep = ""), dec = ".", sep = "\t")



#FC plot
# Threshed
plot(DARs_DMRs_TSS$log2FC_SPF_mLN_pLN, DARs_DMRs_TSS$mean_dMeth_gene)

# All DARs for DMRs
DARs_DMRs_Genes <- merge(DARs_features_Genes, DMR_SPF_mLN_pLN,
                       by.x = "symbol", by.y = "nearest_gene_name")
plot(DARs_DMRs_Genes$log2FC_SPF_mLN_pLN, DARs_DMRs_Genes$meanDiff)

#####
#List of all DMRs that are in accessible regions
#####
###Identify overlap location
hits_location <- findOverlaps(gr_SPF_mLN_pLN, gr_GF_mLN_pLN, minoverlap = 50)
#Get index
query_Hits <- queryHits(hits_location)
subject_Hits <- subjectHits(hits_location)
#Get metadata
overlap_location_SPF <- gr_SPF_mLN_pLN[query_Hits]
overlap_location_GF <- gr_GF_mLN_pLN[subject_Hits]
test <- as.data.frame(overlap_location_GF)

#####
#Overlap between DMRs and Peaks/DARs
#####
#Make Granges object of DARs TSS
gr_DARs_TSS <- makeGRangesFromDataFrame(DARs_features_TSS,
                                         seqnames.field=c("seqnames", "seqname",
                                                          "chromosome", "chrom",
                                                          "chr", "chromosome_name",
                                                          "seqid"),
                                         start.field="start",
                                         end.field=c("end"),
                                         strand.field="strand",
                                         starts.in.df.are.0based=FALSE,
                                         keep.extra.columns = TRUE)
#Make Granges object of DARs all
gr_DARs_all <- makeGRangesFromDataFrame(DARs_features,
                                        seqnames.field=c("seqnames", "seqname",
                                                         "chromosome", "chrom",
                                                         "chr", "chromosome_name",
                                                         "seqid"),
                                        start.field="start",
                                        end.field=c("end"),
                                        strand.field="strand",
                                        starts.in.df.are.0based=FALSE,
                                        keep.extra.columns = TRUE)
#Make Granges object of DMRs
gr_DMRs_TSS <- makeGRangesFromDataFrame(DMR_SPF_mLN_pLN,
                                        seqnames.field=c("seqnames", "seqname",
                                                         "chromosome", "chrom",
                                                         "chr", "chromosome_name",
                                                         "seqid"),
                                        start.field="start",
                                        end.field=c("end"),
                                        strand.field="strand",
                                        starts.in.df.are.0based=FALSE,
                                        keep.extra.columns = TRUE)
#Using DARs TSS as reference
hits <- findOverlaps(gr_DARs_TSS, gr_DMRs_TSS, minoverlap = 10)
#check hits
gr_DARs_TSS_overlap <- gr_DARs_TSS[queryHits(hits)]
gr_DMRs_TSS_overlap <- gr_DMRs_TSS[subjectHits(hits)]
DARs_TSS_overlap <- as.data.frame(gr_DARs_TSS_overlap)
DMRs_TSS_overlap <- as.data.frame(gr_DMRs_TSS_overlap)

#Using DARs TSS as reference
hits <- findOverlaps(gr_DARs_all, gr_DMRs_TSS, minoverlap = 10)
#check hits
gr_DARs_all_overlap <- gr_DARs_all[queryHits(hits)]
gr_DMRs_TSS_overlap <- gr_DMRs_TSS[subjectHits(hits)]
DARs_All_overlap <- as.data.frame(gr_DARs_all_overlap)
DMRs_All_overlap <- as.data.frame(gr_DMRs_TSS_overlap)

#####
#List of all DMRs that are not in accessible regions
#####
NonOverlap_DMRs_with_DARs <- as.data.frame(gr_DMRs_TSS[!gr_DMRs_TSS %over% gr_DARs_TSS,])
length(unique(NonOverlap_DMRs_with_DARs$nearest_gene_name))
NonOverlap_DMRs_with_DARs_TSS <- subset(NonOverlap_DMRs_with_DARs, distance_to_nearest_tss <= 10000)
length(unique(NonOverlap_DMRs_with_DARs_TSS$nearest_gene_name))

write.table(NonOverlap_DMRs_with_DARs, file=paste(Output_DARs_DMRs, "/", "NonOverlap_DMRs_with_DARs.txt", sep = ""),
            quote=F, sep="\t", row.names=T, col.names=T)
write.table(NonOverlap_DMRs_with_DARs_TSS, file=paste(Output_DARs_DMRs, "/", "NonOverlap_DMRs_with_DARs_TSS.txt", sep = ""),
            quote=F, sep="\t", row.names=T, col.names=T)


####################
#Perform GO for DMRs
####################
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
#GeneSymbol to GeneID
#####
#DMRs
idfound <- DMR_SPF_mLN_pLN_TSS$symbol %in% mappedRkeys(org.Mm.egSYMBOL)
SYMBOL <- toTable(org.Mm.egSYMBOL)
head(SYMBOL)
m <- match(DMR_SPF_mLN_pLN_TSS$symbol, SYMBOL$symbol)
GENE_ID <- SYMBOL$gene_id[m]
DMR_SPF_mLN_pLN_TSS <- cbind(GENE_ID, DMR_SPF_mLN_pLN_TSS)

#DARs and Peaks
idfound <- DARs_features_TSS$symbol %in% mappedRkeys(org.Mm.egSYMBOL)
SYMBOL <- toTable(org.Mm.egSYMBOL)
head(SYMBOL)
m <- match(DARs_features_TSS$symbol, SYMBOL$symbol)
GENE_ID <- SYMBOL$gene_id[m]
DARs_features_TSS <- cbind(GENE_ID, DARs_features_TSS)

#####
#Get gene groups
#####
DMRs_hypo_mLN <- subset(DMR_SPF_mLN_pLN_TSS, mean_Diff <= 0)
DMRs_hypo_pLN <- subset(DMR_SPF_mLN_pLN_TSS, mean_Diff >= 0)



#Demethylated in mLN
genes_of_interest_1 <- as.character(DMRs_hypo_mLN$GENE_ID)
# Note: pvalue is not used for fisher's test and filled with 0.05 for lack of p-value from Bsmooth
pvalue_of_interest_1 <- rep(0.05, length(genes_of_interest_1))
names(pvalue_of_interest_1) <- genes_of_interest_1

#Demethylated in pLN
genes_of_interest_2 <- as.character(DMRs_hypo_pLN$GENE_ID)
# Note: pvalue is not used for fisher's test and filled with 0.05 for lack of p-value from Bsmooth
pvalue_of_interest_2 <- rep(0.05, length(genes_of_interest_2))
names(pvalue_of_interest_2) <- genes_of_interest_2

#####
#GO analysis
#####
#Gene universe is the entitiy of transcribed, accessible (at TSS) and demeth genes
#all_groups <- c(DARs_features_TSS$symbol,
 #               genes_of_interest_2,genes_of_interest_1,
  #              diff_Genes$GeneSymbol)
geneNames <- unique(c(genes_of_interest_2,genes_of_interest_1,DARs_features_TSS$GENE_ID))

#gene lists
geneList_1 <- factor(as.integer(geneNames %in% genes_of_interest_1))
geneList_2 <- factor(as.integer(geneNames %in% genes_of_interest_2))

#named gene lists
names(geneList_1) <- geneNames
names(geneList_2) <- geneNames


#list the lists for looping
list_gene_List <- list(geneList_1, geneList_2)
list_pvalue_of_interest <- list(pvalue_of_interest_1, pvalue_of_interest_2)
List_allRes <- list()

#Do  GO statistics for all gene lists
for(i in 1:2){
  #Access the gene lists and p-values for the differentially expressed genes
  geneList <- list_gene_List[[i]]
  pvalue_of_interest <- list_pvalue_of_interest[[i]]
  
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
  #Fishe test
  test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
  resultWeight <- getSigGroups(GOdata, test.stat)
  #KS 
  test.stat <- new("elimScore", testStatistic = GOKSTest, name = "Fisher test", cutOff = 0.01)
  resultElim <- getSigGroups(GOdata, test.stat)
  #runTest
  resultFis <- runTest(GOdata, statistic = "fisher")
  #Kolmogorov-Smirnov -> used for further downstream analysis
  test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
  resultKS <- getSigGroups(GOdata,  test.stat)
  #runTest
  elim.ks <- runTest(GOdata, algorithm = "elim", statistic = "ks")
  
  #make table
  allRes <- GenTable(GOdata, classic = resultFis, KS = resultKS, weight = resultWeight,
                     orderBy = "weight", ranksOf = "KS", topNodes = number_nodes)
  #make list of result tables
  List_allRes[[i]] <- allRes
}


#####
#Find the Top GOs from the lists
#####
#Build table with Top GOs according to "weight"
List_TopGOs <- list()
for(i in 1:2){
  table_i <- List_allRes[[i]]
  table_tophits <- subset(table_i, weight < 0.05)
  #print(i)
  #print(nrow(table_tophits))
  List_TopGOs[[i]] <- table_tophits
}

#collect all TopGos
TopGOs_vector <- c()
for(i in 1:2){
  table_i <- List_TopGOs[[i]]
  TopGOs_i <- as.character(table_i$GO.ID)
  TopGOs_vector <- c(TopGOs_vector, TopGOs_i) 
}

table_1 <- List_allRes[[1]]
table_1_subset <- subset(table_1, table_1$GO.ID %in% TopGOs_vector)
table_1_subset_GO_weight <- table_1_subset[,c("GO.ID","weight", "Term")]
table_1_subset_GO_weight[,2] <- as.numeric(table_1_subset_GO_weight$weight)

table_2 <- List_allRes[[2]]
table_2_subset <- subset(table_2, table_2$GO.ID %in% TopGOs_vector)
table_2_subset_GO_weight <- table_2_subset[,c("GO.ID","weight", "Term")]
table_2_subset_GO_weight[,2] <- as.numeric(table_2_subset_GO_weight$weight)

#merger
a <- merge(table_1_subset_GO_weight, table_2_subset_GO_weight, by = "GO.ID")
data_TopGOs_weight <- a
#sort table
data_TopGOs_weight <- data_TopGOs_weight[,c(1,3,2,4)]
colnames(data_TopGOs_weight) <- c("GOID","GO_Term",
                                  "DMR_hypo_mLN", "DMR_hypo_pLN")
#setwd("C:/Users/jpe12/PowerFolders/R/Paper/2017_Pezoldt_Pasztoi/GO_Tx_FRCs/Output")
#write.table(data_TopGOs_weight, "TopGOs_weight.txt", sep = "\t", dec = ",")

#####
#Make GO comparison heatmap
#####
GO_DMRs <- data_TopGOs_weight
row_names <- paste(GO_DMRs$"GOID", GO_DMRs$GO_Term, sep = "_")
row.names(GO_DMRs) <- row_names
GO_DMRs <- GO_DMRs[,3:4]

min(GO_DMRs)
GO_DMRs <- -log10(GO_DMRs)
#Common Minimum at 5
GO_DMRs[GO_DMRs > 5] <- 5
GOs_maintained_curated_matrix <- data.matrix(GO_DMRs)
title <- c("GOs DMRs")

pheatmap(GOs_maintained_curated_matrix, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = TRUE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("white", "grey","deepskyblue2"), space="rgb")(128),
         main = title)

####################
#Perform GO for DARs
####################
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
#GeneSymbol to GeneID
#####
DARs_features_NA <- DARs_features[!is.na(DARs_features$symbol),]
#DARs
idfound <- DARs_features_NA$symbol %in% mappedRkeys(org.Mm.egSYMBOL)
SYMBOL <- toTable(org.Mm.egSYMBOL)
head(SYMBOL)
m <- match(DARs_features_NA$symbol, SYMBOL$symbol)
GENE_ID <- SYMBOL$gene_id[m]
DARs_features_NA <- cbind(GENE_ID, DARs_features_NA)



#####
#Get gene groups
#####
DARs_Open_mLN_SPF <- subset(DARs_features_NA, log2FC_SPF_mLN_pLN >= 1.0 & padj_SPF_mLN_pLN <= 0.05)
DARs_Closed_mLN_SPF <- subset(DARs_features_NA, log2FC_SPF_mLN_pLN <= -1.0 & padj_SPF_mLN_pLN <= 0.05)



#Demethylated in mLN
genes_of_interest_1 <- unique(as.character(DARs_Open_mLN_SPF$GENE_ID))
# Note: pvalue is not used for fisher's test and filled with 0.05 for lack of p-value from Bsmooth
pvalue_of_interest_1 <- rep(0.05, length(genes_of_interest_1))
names(pvalue_of_interest_1) <- genes_of_interest_1

#Demethylated in pLN
genes_of_interest_2 <- unique(as.character(DARs_Closed_mLN_SPF$GENE_ID))
# Note: pvalue is not used for fisher's test and filled with 0.05 for lack of p-value from Bsmooth
pvalue_of_interest_2 <- rep(0.05, length(genes_of_interest_2))
names(pvalue_of_interest_2) <- genes_of_interest_2

#####
#GO analysis
#####
#Gene universe is the entitiy of transcribed, accessible (at TSS) and demeth genes
#all_groups <- c(DARs_features_TSS$symbol,
#               genes_of_interest_2,genes_of_interest_1,
#              diff_Genes$GeneSymbol)
geneNames <- unique(c(genes_of_interest_2,genes_of_interest_1,DARs_features_TSS$GENE_ID))

#gene lists
geneList_1 <- factor(as.integer(geneNames %in% genes_of_interest_1))
geneList_2 <- factor(as.integer(geneNames %in% genes_of_interest_2))

#named gene lists
names(geneList_1) <- geneNames
names(geneList_2) <- geneNames


#list the lists for looping
list_gene_List <- list(geneList_1, geneList_2)
list_pvalue_of_interest <- list(pvalue_of_interest_1, pvalue_of_interest_2)
List_allRes <- list()

#Do  GO statistics for all gene lists
for(i in 1:2){
  #Access the gene lists and p-values for the differentially expressed genes
  geneList <- list_gene_List[[i]]
  pvalue_of_interest <- list_pvalue_of_interest[[i]]
  
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
  #Fishe test
  test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
  resultWeight <- getSigGroups(GOdata, test.stat)
  #KS 
  test.stat <- new("elimScore", testStatistic = GOKSTest, name = "Fisher test", cutOff = 0.01)
  resultElim <- getSigGroups(GOdata, test.stat)
  #runTest
  resultFis <- runTest(GOdata, statistic = "fisher")
  #Kolmogorov-Smirnov -> used for further downstream analysis
  test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
  resultKS <- getSigGroups(GOdata,  test.stat)
  #runTest
  elim.ks <- runTest(GOdata, algorithm = "elim", statistic = "ks")
  
  #make table
  allRes <- GenTable(GOdata, classic = resultFis, KS = resultKS, weight = resultWeight,
                     orderBy = "weight", ranksOf = "KS", topNodes = number_nodes)
  #make list of result tables
  List_allRes[[i]] <- allRes
}


#####
#Find the Top GOs from the lists
#####
#Build table with Top GOs according to "weight"
List_TopGOs <- list()
for(i in 1:2){
  table_i <- List_allRes[[i]]
  table_tophits <- subset(table_i, weight < 0.05)
  #print(i)
  #print(nrow(table_tophits))
  List_TopGOs[[i]] <- table_tophits
}

#collect all TopGos
TopGOs_vector <- c()
for(i in 1:2){
  table_i <- List_TopGOs[[i]]
  TopGOs_i <- as.character(table_i$GO.ID)
  TopGOs_vector <- c(TopGOs_vector, TopGOs_i) 
}

table_1 <- List_allRes[[1]]
table_1_subset <- subset(table_1, table_1$GO.ID %in% TopGOs_vector)
table_1_subset_GO_weight <- table_1_subset[,c("GO.ID","weight", "Term")]
table_1_subset_GO_weight[,2] <- as.numeric(table_1_subset_GO_weight$weight)

table_2 <- List_allRes[[2]]
table_2_subset <- subset(table_2, table_2$GO.ID %in% TopGOs_vector)
table_2_subset_GO_weight <- table_2_subset[,c("GO.ID","weight", "Term")]
table_2_subset_GO_weight[,2] <- as.numeric(table_2_subset_GO_weight$weight)

#merger
a <- merge(table_1_subset_GO_weight, table_2_subset_GO_weight, by = "GO.ID")
data_TopGOs_weight <- a
#sort table
data_TopGOs_weight <- data_TopGOs_weight[,c(1,3,2,4)]
colnames(data_TopGOs_weight) <- c("GOID","GO_Term",
                                  "DAR_Open_mLN_SPF", "DMR_Closed_mLN_SPF")
#setwd("C:/Users/jpe12/PowerFolders/R/Paper/2017_Pezoldt_Pasztoi/GO_Tx_FRCs/Output")
#write.table(data_TopGOs_weight, "TopGOs_weight.txt", sep = "\t", dec = ",")

#####
#Make GO comparison heatmap
#####
GO_DARs <- data_TopGOs_weight
row_names <- paste(GO_DARs$"GOID", GO_DARs$GO_Term, sep = "_")
row.names(GO_DARs) <- row_names
GO_DARs <- GO_DARs[,3:4]

min(GO_DARs)
GO_DARs <- -log10(GO_DARs)
#Common Minimum at 5
GO_DARs[GO_DARs > 5] <- 5
GOs_maintained_curated_matrix <- data.matrix(GO_DARs)
title <- c("GOs DARs")

pheatmap(GOs_maintained_curated_matrix, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = TRUE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("white", "grey","deepskyblue2"), space="rgb")(128),
         main = title)


#####################
#Which genes are in identified GOs
#####################
AllGenes_interest <- c(genes_of_interest_1, genes_of_interest_2)


List_GOs_interest <- as.character(data_TopGOs_weight$GOID)
List_GOs_interest <- c("GO:0060348","GO:2000112","GO:0032872","GO:0014066")
Genes_for_GOs_interest <- subset(GO2GeneID, ls(GO2GeneID) %in% List_GOs_interest)
number_input <- length(Genes_for_GOs_interest)

for(i in 1:number_input){
  
  #Get genes in GO group
  GO_interest_Genes <- subset(AllGenes_interest, AllGenes_interest %in% Genes_for_GOs_interest[[4]])
  print(GO_interest_Genes)
  number_genes_in_GO <- length(Genes_for_GOs_interest[[1]])
  #Get gene names
  gene_names <- GO_interest_Genes
  #print(gene_names)
  #Generate list of genes of interest contained in expression data
  genes_expression <- subset(Data_ready, Data_ready$GENE_ID %in% gene_names)
  number_genes_expressed_in_GO <- nrow(genes_expression)
  
  if(number_genes_expressed_in_GO >= 2){
    #format data for pheatmap
    row.names(genes_expression) <- paste(round(log2(rowMeans(genes_expression[,3:10]))), "_",
                                         genes_expression$GeneSymbol, sep = "")
    data_heatmap <- genes_expression[, c(2,3:6)]
    data_heatmap <- data_heatmap[,2:5]
    data_heatmap <- log2(data_heatmap)
    min(data_heatmap)
    #Common Minimum
    data_heatmap[data_heatmap == -Inf] <- -9
    data_heatmap_matrix <- data.matrix(data_heatmap)
    print(i)
    #content of GO expressed
    ratio_expressed_to_GO <- paste(number_genes_expressed_in_GO, "/", number_genes_in_GO, sep = "" )
    
    #title
    list_GOs <- ls(Genes_for_GOs_interest)
    title <- list_GOs[i]
    Go_term <- as.character(GOs_interest[GOs_interest$"GO.ID" == title, c("GO_Term")])
    Go_term <- unlist(strsplit(Go_term, split=':', fixed=TRUE))[1]
    column_header <- gsub(" ", "_", Go_term, fixed = TRUE)
    title_figure <- paste(ratio_expressed_to_GO,title, Go_term, sep = "_")
    #filename
    file_name <- unlist(strsplit(title, split=':', fixed=TRUE))[2]
    file_name <- paste("GO_", file_name, "__", column_header, ".pdf", sep = "")
    
    #Determine heigth of file
    height_pdf = 2 + nrow(data_heatmap_matrix) * 10 / 72
    
    setwd("C:/Users/jpe12/Desktop/R_Output")
    #draw heatmap and save to folder
    pdf(file_name, height = height_pdf, width = 8)
    pheatmap(data_heatmap_matrix, cluster_rows = TRUE, legend = TRUE,
             treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = TRUE,
             scale = "row", border_color = "black", cellwidth = 10,
             cellheigth = 10, color = colorRampPalette(c("blue", "yellow","red"), space="rgb")(128),
             main = title_figure)
    dev.off()
  }
}































####################
#Generate a Heatmap for each DMR !!!!!Under construction
####################

#Windows-Powerfolder
#Load CpGs
CpGs <- read.csv("C:/Users/jpe12/PowerFolders/Stroma_Paper/Methylome/DMR_heatmap/Data/CpGs_in_DMRs.csv")

#Load Bsmooth data
mLN_SPF_vs_pLN_SPF <-  read.delim("C:/Users/jpe12/PowerFolders/Stroma_Paper/Methylome/DMR_heatmap/Data/DMRs_bsmooth/BSmooth_DMRs_mLN_SPF_vs_pLN_SPF.txt", dec=",")
mLN_SPF_vs_mLN_GF <- read.delim("C:/Users/jpe12/PowerFolders/Stroma_Paper/Methylome/DMR_heatmap/Data/DMRs_bsmooth/BSmooth_DMRs_mLN_SPF_vs_mLN_GF.txt", dec=",")
pLN_GF_vs_mLN_GF <-read.delim("C:/Users/jpe12/PowerFolders/Stroma_Paper/Methylome/DMR_heatmap/Data/DMRs_bsmooth/BSmooth_DMRs_pLN_GF_vs_mLN_GF.txt", dec=",")
pLN_SPF_vs_pLN_GF <- read.delim("C:/Users/jpe12/PowerFolders/Stroma_Paper/Methylome/DMR_heatmap/Data/DMRs_bsmooth/BSmooth_DMRs_pLN_SPF_vs_pLN_GF.txt", dec=",")

#Assemble essential DMR info -------------------------------------
Bsmooth_tables <- list(mLN_SPF_vs_pLN_SPF, mLN_SPF_vs_mLN_GF,
                       pLN_GF_vs_mLN_GF, pLN_SPF_vs_pLN_GF)
Bsmooth_names_analysis <- c("mLN_SPF_vs_pLN_SPF", "mLN_SPF_vs_mLN_GF",
                            "pLN_GF_vs_mLN_GF","pLN_SPF_vs_pLN_GF")

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
  x <- Bsmooth_tables[[1]]
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

#threshhold DMRs------------------------------
Bsmooth_intel_subset <- subset(Bsmooth_intel, (Bsmooth_intel$meanDiff <= -0.5 | Bsmooth_intel$meanDiff >= 0.5) |
                                 (Bsmooth_intel$areaStat >= 300 | Bsmooth_intel$areaStat <= -300))

#Grab Threshed DMRs from CpG list -----------------------------------------------------
#split according to DMR_identifier
CpGs_DMRs <- split(CpGs, CpGs$DMR_identifier)

#get DMR_identifiers from theshed list
threshed_DMRs <- as.character(Bsmooth_intel_subset$DMR_identifier)
CpGs_DMRs_threshed <- CpGs_DMRs[threshed_DMRs]


#heatmap -----------------------------------------------
#get data from Bsmooth analysis
path_variable = "C:/Users/jpe12/Desktop/R_Output/test_pdf.pdf"
pdf(file = path_variable, height = 2.3, width = 20)

for(i in 1:length(CpGs_DMRs_threshed)){
  par(mfcol = c(4,1))
  
  #Get info for main/title
  CpGs_DMRs_i <- CpGs_DMRs_threshed[[i]]
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




################
#Generate Output for home Motif analysis
################
#####
#Peak/gene location
#####
library("GenomicFeatures")
#All features
mm10_TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#TSS sites BioMart
mm10 = useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl") 
TSS.mouse.mm10 = getAnnotation(mart=mm10, featureType="TSS")

#Annotate Peaks SPF-------------------------------
gr_DMR_SPF_mLN_pLN_anno_TSS <- annotatePeakInBatch(gr_DMR_SPF_mLN_pLN, AnnotationData=TSS.mouse.mm10)

#add gene name
gr_DMR_SPF_mLN_pLN_anno_TSS <- addGeneIDs(annotatedPeak=gr_DMR_SPF_mLN_pLN_anno_TSS, 
                                          feature_id_type="ensembl_gene_id",
                                          orgAnn="org.Mm.eg.db", 
                                          IDs2Add="symbol")
#Generate dataframe
DMR_SPF_mLN_pLN_anno_TSS <- as.data.frame(gr_DMR_SPF_mLN_pLN_anno_TSS)

#Annotate Peaks GF-------------------------------
gr_DMR_GF_mLN_pLN_anno_TSS <- annotatePeakInBatch(gr_DMR_GF_mLN_pLN, AnnotationData=TSS.mouse.mm10)

#add gene name
gr_DMR_GF_mLN_pLN_anno_TSS <- addGeneIDs(annotatedPeak=gr_DMR_GF_mLN_pLN_anno_TSS, 
                                         feature_id_type="ensembl_gene_id",
                                         orgAnn="org.Mm.eg.db", 
                                         IDs2Add="symbol")
#Generate dataframe
DMR_GF_mLN_pLN_anno_TSS <- as.data.frame(gr_DMR_GF_mLN_pLN_anno_TSS)

#####
#Export for Motif analysis
#####
### Subset Tables
# SPF
DMR_SPF_mLN_pLN_TSS <- DMR_SPF_mLN_pLN_anno_TSS[!(is.na(DMR_SPF_mLN_pLN_anno_TSS$symbol)),]
DMR_SPF_mLN_pLN_TSS_hypo_thresh <- subset(DMR_SPF_mLN_pLN_TSS, mean_Diff <= -0.3)
DMR_SPF_mLN_pLN_TSS_hypo_all <- subset(DMR_SPF_mLN_pLN_TSS, mean_Diff < 0)
DMR_SPF_mLN_pLN_TSS_hyper_thresh <- subset(DMR_SPF_mLN_pLN_TSS, mean_Diff >= 0.3)
DMR_SPF_mLN_pLN_TSS_hyper_all <- subset(DMR_SPF_mLN_pLN_TSS, mean_Diff > 0)
# GF
DMR_GF_mLN_pLN_TSS <- DMR_GF_mLN_pLN_anno_TSS[!(is.na(DMR_GF_mLN_pLN_anno_TSS$symbol)),]
DMR_GF_mLN_pLN_TSS_hypo_thresh <- subset(DMR_GF_mLN_pLN_TSS, mean_Diff <= -0.3)
DMR_GF_mLN_pLN_TSS_hypo_all <- subset(DMR_GF_mLN_pLN_TSS, mean_Diff < 0)
DMR_GF_mLN_pLN_TSS_hyper_thresh <- subset(DMR_GF_mLN_pLN_TSS, mean_Diff >= 0.3)
DMR_GF_mLN_pLN_TSS_hyper_all <- subset(DMR_GF_mLN_pLN_TSS, mean_Diff > 0)


#List of Tables
l_DMRs_to_homer <- list(DMR_SPF_mLN_pLN_TSS_hypo_thresh,DMR_SPF_mLN_pLN_TSS_hypo_all,
                        DMR_SPF_mLN_pLN_TSS_hyper_thresh,DMR_SPF_mLN_pLN_TSS_hyper_all,
                        DMR_GF_mLN_pLN_TSS_hypo_thresh,DMR_GF_mLN_pLN_TSS_hypo_all,
                        DMR_GF_mLN_pLN_TSS_hyper_thresh,DMR_GF_mLN_pLN_TSS_hyper_all)
names(l_DMRs_to_homer) <- c("DMR_SPF_mLN_pLN_TSS_hypo_thresh","DMR_SPF_mLN_pLN_TSS_hypo_all",
                            "DMR_SPF_mLN_pLN_TSS_hyper_thresh","DMR_SPF_mLN_pLN_TSS_hyper_all",
                            "DMR_GF_mLN_pLN_TSS_hypo_thresh","DMR_GF_mLN_pLN_TSS_hypo_all",
                            "DMR_GF_mLN_pLN_TSS_hyper_thresh","DMR_GF_mLN_pLN_TSS_hyper_all")

#Empty list for storage
l_DMRs_to_homer_fini <- list()
for(i in seq(length(l_DMRs_to_homer))){
  t_DMR_i <- l_DMRs_to_homer[[i]]
  
  #build table input
  v_chr <- as.character(t_DMR_i$seqnames)
  l_chr <- strsplit(v_chr, "r", fixed = FALSE, perl = FALSE, useBytes = FALSE)
  seqnames_i <- unlist(lapply(l_chr, function(x) x[[2]]))
  start_i = t_DMR_i$start
  end_i = t_DMR_i$end
  id_i = t_DMR_i$peak
  blankVar_i = rep("",length(id_i))
  feature_strand_i = rep("+",length(id_i))
  t_DMR_adapted_i <- data.frame(seqnames = seqnames_i, start = start_i,
                                end = end_i, id = id_i, blankVar = blankVar_i,
                                feature_strand = feature_strand_i)
  
  #Store tables for motif identification
  write.table(t_DMR_adapted_i, file=paste(output_beds_for_motifs, "/", names(l_DMRs_to_homer)[i],".bed", sep = ""),
              quote=F, sep="\t", row.names=F, col.names=F)
  l_DMRs_to_homer_fini[[i]] <- t_DMR_adapted_i
}

#####
#Obtain DMRs evolving around the single CpGs
#####
#Variables
# delta meth required for CpG to be included and locied
CpG_state_thresh <- 0.4
#Function:
# Merges overlapping entries in Granges object
# Extract clusters from Hits object.
extractClustersFromSelfHits <- function(hits)
{
  stopifnot(is(hits, "Hits"))
  stopifnot(queryLength(hits) == subjectLength(hits))
  hits <- union(hits, t(hits))
  qh <- queryHits(hits)
  sh <- subjectHits(hits)
  cid <- seq_len(queryLength(hits))  # cluster ids
  while (TRUE) {
    h <- Hits(qh, cid[sh],
              queryLength(hits), subjectLength(hits))
    cid2 <- pmin(cid, selectHits(h, "first"))
    if (identical(cid2, cid))
      break
    cid <- cid2
  }
  unname(splitAsList(seq_len(queryLength(hits)), cid))
}
# Merge ranges that are "connected" (directly or indirectly)
# via a hit (or several hits) in 'hits'.
mergeConnectedRanges <- function(x, hits)
{
  stopifnot(is(x, "GenomicRanges"))
  stopifnot(is(hits, "Hits"))
  stopifnot(queryLength(hits) == subjectLength(hits))
  stopifnot(queryLength(hits) == length(x))
  clusters <- extractClustersFromSelfHits(hits)
  ans <- range(extractList(x, clusters))
  #if (any(elementLengths(ans) != 1L))
  # stop(wmsg("some connected ranges are not on the same ",
  #          "chromosome and strand, and thus cannot be ",
  #         "merged"))
  ans <- unlist(ans)
  mcols(ans)$revmap <- clusters
  ans
}

#Explanation
# CpGs are widely spaced and under the assumption that methylation dependent TF binding takes place
# the DMRs will be divided along the CpG
# The new loci will be CpG position +/- X bp
l_DMRs_CpG_regions <- list(DMR_SPF_mLN_pLN_TSS_hypo_thresh,DMR_SPF_mLN_pLN_TSS_hypo_all,
                           DMR_SPF_mLN_pLN_TSS_hyper_thresh,DMR_SPF_mLN_pLN_TSS_hyper_all,
                           DMR_GF_mLN_pLN_TSS_hypo_thresh,DMR_GF_mLN_pLN_TSS_hypo_all,
                           DMR_GF_mLN_pLN_TSS_hyper_thresh,DMR_GF_mLN_pLN_TSS_hyper_all)

names(l_DMRs_CpG_regions) <- list("DMR_SPF_mLN_pLN_TSS_hypo_thresh","DMR_SPF_mLN_pLN_TSS_hypo_all",
                                  "DMR_SPF_mLN_pLN_TSS_hyper_thresh","DMR_SPF_mLN_pLN_TSS_hyper_all",
                                  "DMR_GF_mLN_pLN_TSS_hypo_thresh","DMR_GF_mLN_pLN_TSS_hypo_all",
                                  "DMR_GF_mLN_pLN_TSS_hyper_thresh","DMR_GF_mLN_pLN_TSS_hyper_all")

#initialize lists that will contain the list of the CpGs regions per DMR
l_26bp <- list()
l_50bp <- list()
l_100bp <- list()
l_200bp <- list()

#Note: Output of the following loop is a nested list with the content of the input DMR tables as a list of CpG regions per DMR 
#       which is then stored in list annotating the size of the regions
for(i in 1:length(l_DMRs_CpG_regions)){
  # Grab DMRs threshed in DMRs lists
  print(paste("i == ", i, sep = ""))
  DMRs_i <- l_DMRs_CpG_regions[[i]]
  #Obtain common identifier
  Identifier_CpG_table_i <- paste(DMRs_i$seqnames, "_", DMRs_i$start,"_",DMRs_i$end,sep="")
  #Check in CpG list for lists of respective DMR
  #Take only "mLN_SPF_vs_pLN_SPF"
  CpGs_DMR_SPF_mLN_pLN_i <- subset(CpGs_per_DMR, DMR_comp_origin == "mLN_SPF_vs_pLN_SPF")
  split_CpGs_per_DMR_i <- split(CpGs_DMR_SPF_mLN_pLN_i, CpGs_DMR_SPF_mLN_pLN_i$DMR_identifier)
  split_CpGs_per_DMR_overlap_i <- split_CpGs_per_DMR_i[Identifier_CpG_table_i]
  #Initialize lists objects to store results of CpG extension and merge
  ll_CpG_per_DMR_26bp_i <- list()
  ll_CpG_per_DMR_50bp_i <- list()
  ll_CpG_per_DMR_100bp_i <- list()
  ll_CpG_per_DMR_200bp_i <- list()
  #for each DMR go through the CpGs
  for(j in 1:length(Identifier_CpG_table_i)){
    print(paste("j == ", j, sep = ""))
    #grab ith CpG DMR content
    t_CpGs_per_DMR_j <- as.data.frame(split_CpGs_per_DMR_overlap_i[names(split_CpGs_per_DMR_overlap_i[j])])
    t_CpGs_per_DMR_j <- t_CpGs_per_DMR_j[,2:ncol(t_CpGs_per_DMR_j)]
    colnames(t_CpGs_per_DMR_j) <- c("seqname","start","end",
                                    "DMR_identifier","DMR_comp_origin","DMR_num_CpGs",
                                    "X2_pLN_GF","X3a_pLN_GF","X5_mLN_GF",
                                    "X6_mLN_GF","X7_pLN_SPF","X8_pLN_SPF",
                                    "X9_pLN_SPF","X10_mLN_SPF","X11_mLN_SPF",
                                    "X12_mLN_SPF")
    #Print Number of input CpGs
    print(paste("n input CpGs: ", nrow(t_CpGs_per_DMR_j)))
    
    #Calculate the % methylation per CpG
    meth_mLN_SPF <- rowMeans(t_CpGs_per_DMR_j[,c(14:16)])
    meth_pLN_SPF <- rowMeans(t_CpGs_per_DMR_j[,c(11:13)])
    meth_mLN_GF <- rowMeans(t_CpGs_per_DMR_j[,c(9,10)])
    meth_pLN_GF <- rowMeans(t_CpGs_per_DMR_j[,c(7,8)])
    
    #delta methylation for mLNSPF vs pLN SPF
    dmeth_DMR_SPF_mLN_pLN <- meth_mLN_SPF - meth_pLN_SPF
    dmeth_DMR_GF_mLN_pLN <- meth_mLN_GF - meth_pLN_GF
    
    #append delta meth columns
    t_CpGs_per_DMR_j <- cbind(t_CpGs_per_DMR_j, dmeth_DMR_SPF_mLN_pLN, dmeth_DMR_GF_mLN_pLN)
    
    #select CpGs within DMR based on delta meth
    t_CpGs_per_DMR_thresh_j <- subset(t_CpGs_per_DMR_j, abs(dmeth_DMR_SPF_mLN_pLN) > CpG_state_thresh)
    #Some DMRs will have no CpG with a decent dMeth
    if(nrow(t_CpGs_per_DMR_thresh_j) > 0 ){
      #make Granges object
      gr_CpGs_per_DMR_thresh_j <- makeGRangesFromDataFrame(t_CpGs_per_DMR_thresh_j,
                                                           seqnames.field=c("seqnames", "seqname",
                                                                            "chromosome", "chrom",
                                                                            "chr", "chromosome_name",
                                                                            "seqid"),
                                                           start.field="start",
                                                           end.field=c("end"),
                                                           strand.field="strand",
                                                           starts.in.df.are.0based=FALSE,
                                                           keep.extra.columns = TRUE)
      # Granges Objects to be changed
      gr_CpGs_per_DMR_26bp_j <- gr_CpGs_per_DMR_thresh_j
      gr_CpGs_per_DMR_50bp_j <- gr_CpGs_per_DMR_thresh_j
      gr_CpGs_per_DMR_100bp_j <- gr_CpGs_per_DMR_thresh_j
      gr_CpGs_per_DMR_200bp_j <- gr_CpGs_per_DMR_thresh_j
      
      #Extend by:
      # 13 each direction
      start(gr_CpGs_per_DMR_26bp_j) <- start(gr_CpGs_per_DMR_50bp_j) - 13
      end(gr_CpGs_per_DMR_26bp_j) <- end(gr_CpGs_per_DMR_50bp_j) + 13
      hits <- findOverlaps(gr_CpGs_per_DMR_26bp_j)
      ## Subset 'hits' to keep only hits that achieve 60% overlap.
      x <- gr_CpGs_per_DMR_26bp_j[queryHits(hits)]
      y <- gr_CpGs_per_DMR_26bp_j[subjectHits(hits)]
      relative_overlap <- width(pintersect(x, y)) / pmin(width(x), width(y))
      hits <- hits[relative_overlap >= 0.01]
      ## Merge the ranges in 'gr0' that are connected via one or more hits in 'hits'.
      gr_CpGs_per_DMR_26bp_j <- mergeConnectedRanges(gr_CpGs_per_DMR_26bp_j, hits)
      #Print Number of input CpGs
      print(paste("n output 26bp regions: ", length(gr_CpGs_per_DMR_26bp_j)))
      
      #Extend by:
      # 25 each direction
      start(gr_CpGs_per_DMR_50bp_j) <- start(gr_CpGs_per_DMR_50bp_j) - 25
      end(gr_CpGs_per_DMR_50bp_j) <- end(gr_CpGs_per_DMR_50bp_j) + 25
      hits <- findOverlaps(gr_CpGs_per_DMR_50bp_j)
      ## Subset 'hits' to keep only hits that achieve 60% overlap.
      x <- gr_CpGs_per_DMR_50bp_j[queryHits(hits)]
      y <- gr_CpGs_per_DMR_50bp_j[subjectHits(hits)]
      relative_overlap <- width(pintersect(x, y)) / pmin(width(x), width(y))
      hits <- hits[relative_overlap >= 0.01]
      ## Merge the ranges in 'gr0' that are connected via one or more hits in 'hits'.
      gr_CpGs_per_DMR_50bp_j <- mergeConnectedRanges(gr_CpGs_per_DMR_50bp_j, hits)
      #Print Number of input CpGs
      print(paste("n output 50bp regions: ", length(gr_CpGs_per_DMR_50bp_j)))
      
      # 50 each direction
      start(gr_CpGs_per_DMR_100bp_j) <- start(gr_CpGs_per_DMR_100bp_j) - 50
      end(gr_CpGs_per_DMR_100bp_j) <- end(gr_CpGs_per_DMR_100bp_j) + 50
      hits <- findOverlaps(gr_CpGs_per_DMR_100bp_j)
      ## Subset 'hits' to keep only hits that achieve 60% overlap.
      x <- gr_CpGs_per_DMR_100bp_j[queryHits(hits)]
      y <- gr_CpGs_per_DMR_100bp_j[subjectHits(hits)]
      relative_overlap <- width(pintersect(x, y)) / pmin(width(x), width(y))
      hits <- hits[relative_overlap >= 0.01]
      ## Merge the ranges in 'gr0' that are connected via one or more hits in 'hits'.
      gr_CpGs_per_DMR_100bp_j <- mergeConnectedRanges(gr_CpGs_per_DMR_100bp_j, hits)
      #Print Number of input CpGs
      print(paste("n output 100bp regions: ", length(gr_CpGs_per_DMR_100bp_j)))
      
      # 100 each direction
      start(gr_CpGs_per_DMR_200bp_j) <- start(gr_CpGs_per_DMR_200bp_j) - 100
      end(gr_CpGs_per_DMR_200bp_j) <- end(gr_CpGs_per_DMR_200bp_j) + 100
      hits <- findOverlaps(gr_CpGs_per_DMR_200bp_j)
      ## Subset 'hits' to keep only hits that achieve 60% overlap.
      x <- gr_CpGs_per_DMR_200bp_j[queryHits(hits)]
      y <- gr_CpGs_per_DMR_200bp_j[subjectHits(hits)]
      relative_overlap <- width(pintersect(x, y)) / pmin(width(x), width(y))
      hits <- hits[relative_overlap >= 0.01]
      ## Merge the ranges in 'gr0' that are connected via one or more hits in 'hits'.
      gr_CpGs_per_DMR_200bp_j <- mergeConnectedRanges(gr_CpGs_per_DMR_200bp_j, hits)
      #Print Number of input CpGs
      print(paste("n output 200bp regions: ", length(gr_CpGs_per_DMR_200bp_j)))
      
      #store results in list
      ll_CpG_per_DMR_26bp_i[[j]] <- gr_CpGs_per_DMR_26bp_j
      ll_CpG_per_DMR_50bp_i[[j]] <- gr_CpGs_per_DMR_50bp_j
      ll_CpG_per_DMR_100bp_i[[j]] <- gr_CpGs_per_DMR_100bp_j
      ll_CpG_per_DMR_200bp_i[[j]] <- gr_CpGs_per_DMR_200bp_j
    }
  }
  l_26bp[[i]] <- ll_CpG_per_DMR_26bp_i
  l_50bp[[i]] <- ll_CpG_per_DMR_50bp_i
  l_100bp[[i]] <- ll_CpG_per_DMR_100bp_i
  l_200bp[[i]] <- ll_CpG_per_DMR_200bp_i
}

#
l_all_regions <- list(l_26bp,l_50bp,l_100bp,l_200bp)
names(l_all_regions) <- c("26bp","50bp","100bp","200bp")

for(i in 1:length(l_all_regions)){
  print(i)
  #take object containing all data for one size (e.g. 50bp)
  l_CpG_bp_i <- l_all_regions[[i]]
  #take the name for later naming of BED files
  l_CpG_bp_ID_i <- names(l_all_regions)[i]
  print(l_CpG_bp_ID_i)
  for(j in 1:length(l_CpG_bp_i)){
    print(j)
    #take list containing all CpG regions for selected DMRs
    l_CpG_bp_output_j <- l_CpG_bp_i[[j]]
    #take names of the DMR election
    l_CpG_bp_output_ID_j <- names(l_DMRs_CpG_regions)[j]
    print(l_CpG_bp_output_ID_j)
    #get rid of empty rows in list of Granges object
    l_CpG_bp_output_j <- unlist(l_CpG_bp_output_j)
    #instigate empty dataframe
    t_CpG_bp_output_j <- data.frame(matrix(ncol=0,nrow=0))
    for(k in 1:length(l_CpG_bp_output_j)){
      #make a table of each Granges object
      r_CpG_bp_output_j <- as.data.frame(l_CpG_bp_output_j[[k]])
      #print(t_CpG_bp_output_j)
      #append all CpG elements per DMR selection
      t_CpG_bp_output_j <- rbind(t_CpG_bp_output_j, r_CpG_bp_output_j)
      #print(t_CpG_bp_output_j)
      #Generate homer compatible region file
      t_CpG_bp_output_bed_j <- t_CpG_bp_output_j[,c("seqnames","start","end","strand")] 
      #Make chr1 -> 1
      v_chr <- as.character(t_CpG_bp_output_bed_j$seqnames)
      l_chr <- strsplit(v_chr, "r", fixed = FALSE, perl = FALSE, useBytes = FALSE)
      v_chr_integers <- unlist(lapply(l_chr, function(x) x[[2]]))
      t_CpG_bp_output_bed_j$seqnames <- v_chr_integers
      #add identifier
      t_CpG_bp_output_bed_j$id <- paste(t_CpG_bp_output_bed_j$seqname, "_",
                                        t_CpG_bp_output_bed_j$start, "_",
                                        t_CpG_bp_output_bed_j$end, sep = "")
      
      #Generate empty column required by homer
      t_CpG_bp_output_bed_j$blankVar <- NA
      t_CpG_bp_output_bed_j <- t_CpG_bp_output_bed_j[,c("seqnames","start","end","id","blankVar","strand")]
      t_CpG_bp_output_bed_j$blankVar <- c("")
      t_CpG_bp_output_bed_j$strand <- "+"
      #Store Feature BEDs 
      #Write .bed compatible for homer
      print(head(t_CpG_bp_output_bed_j))
      write.table(t_CpG_bp_output_bed_j, file=paste(output_beds_for_motifs,"/CpG_DMRs_40/", l_CpG_bp_ID_i,"_",l_CpG_bp_output_ID_j,".bed", sep = ""),
                  quote=F, sep="\t", row.names=F, col.names=F)
    }
  }
}

