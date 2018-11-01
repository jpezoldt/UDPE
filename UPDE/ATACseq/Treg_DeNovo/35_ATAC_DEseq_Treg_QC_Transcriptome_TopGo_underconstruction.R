# Author: Vincent Gardeux
# Adapted by: Joern Pezoldt
# 12.07.2018
# Function:
# 1) 
# 

#####
#Libraries
#####

require(data.table)
require(ggplot2)
require(ggfortify)
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
library(org.Mm.eg.db)
library(GenomicFeatures)
library(stringr)
library(reshape)
library(calibrate)

#####
#Global variables and paths
#####
#Param
log2FC_RNA = 0.3
log2FC_ATAC = 0.58
padj = 0.05

#TSS region 
# Standard usage of homer 2000 bp upstream
downstream <- 200
upstream <- 2000

# Input required: Set directory
# 1) Counts per peak per replicate
path_input <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_denovo_Treg/homer/Overlap_Group_Merged/Run_1_in_all"
# 2) Common Peak regions
path_common_peak_regions <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_denovo_Treg/peaks/broad/Run_2/ATAC_DeNovoTreg_broad_merged_peaks.bed"
# 3) Output path
path_output <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_denovo_Treg/DESeq2/Run1"
path_output_Motif <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_denovo_Treg/Motif"
# 4) RNAseq DESeq2 analysis
path_RNAseq_DESeq2 <- "/home/pezoldt/NAS2/pezoldt/Analysis/RNAseq/2018_HZI_DeNovoTreg"

# Input required: 
name = "ATAC_DeNovoTreg_all"
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
number_of_samples <- 16

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
#DEseq2 on Peaks
#####
#Define conditons 
condition <- factor(substr(colnames(counts.mat),1,nchar(colnames(counts.mat))-1))
dds <- DESeqDataSetFromMatrix(counts.mat, DataFrame(condition),~condition)
dds.factors <- estimateSizeFactors(dds)

#Perform DESeq2
dds <- DESeq(dds)
res <- results(dds)

#generate dds-resuls objects for relevant comparisons
#mLN_SPF_GF <- results(dds, contrast=c("condition","ATAC_mLNSPF","ATAC_mLNGF"))
#pLN_SPF_GF <- results(dds, contrast=c("condition","ATAC_pLNSPF","ATAC_pLNGF"))
#SPF_mLN_pLN <- results(dds, contrast=c("condition","ATAC_mLNSPF","ATAC_pLNSPF"))
#GF_mLN_pLN <- results(dds, contrast=c("condition","ATAC_mLNGF","ATAC_pLNGF"))

Treg_SPF_GF <- results(dds, contrast=c("condition","ATAC_TregSPF","ATAC_TregGF"))
Tconv_SPF_GF <- results(dds, contrast=c("condition","ATAC_TconSPF","ATAC_TconGF"))
SPF_Treg_Tconv <- results(dds, contrast=c("condition","ATAC_TregSPF","ATAC_TconSPF"))
GF_Treg_Tconv <- results(dds, contrast=c("condition","ATAC_TregGF","ATAC_TconGF"))
Tnaive_SPF_Tconv <- results(dds, contrast=c("condition","ATAC_Tnaive","ATAC_TconSPF"))
Tnaive_SPF_Treg <- results(dds, contrast=c("condition","ATAC_Tnaive","ATAC_TregSPF"))

#QC
#MA
par(mfrow = c(3,2))
plotMA(Treg_SPF_GF, ylim = c(-6, 6), main = "Treg_SPF_GF")
plotMA(Tconv_SPF_GF, ylim = c(-6, 6), main = "Tconv_SPF_GF" )
plotMA(SPF_Treg_Tconv, ylim = c(-6, 6), main = "SPF_Treg_Tconv" )
plotMA(GF_Treg_Tconv, ylim = c(-6, 6), main = "GF_Treg_Tconv" )
plotMA(Tnaive_SPF_Tconv, ylim = c(-8, 6), main = "Tnaive_SPF_Tconv" )
plotMA(Tnaive_SPF_Treg, ylim = c(-8, 6), main = "Tnaive_SPF_Treg" )

#DispEsts
par(mfrow = c(1,1))
plotDispEsts( dds, ylim = c(1e-6, 1e1), main = "All Samples" )


#Hits pValue
par(mfrow = c(3,2))
hist( Treg_SPF_GF$pvalue, breaks=20, col="grey", main = "Treg_SPF_GF"  )
hist( Tconv_SPF_GF$pvalue, breaks=20, col="grey", main = "Tconv_SPF_GF"  )
hist( SPF_Treg_Tconv$pvalue, breaks=20, col="grey", main = "SPF_Treg_Tconv"  )
hist( GF_Treg_Tconv$pvalue, breaks=20, col="grey", main = "GF_Treg_Tconv"  )
hist( Tnaive_SPF_Tconv$pvalue, breaks=20, col="grey", main = "Tnaive_SPF_Tconv"  )
hist( Tnaive_SPF_Treg$pvalue, breaks=20, col="grey", main = "Tnaive_SPF_Treg"  )

par(mfrow = c(1,1))
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

#Get log2FC and pValue for all peaks across all chosen comparisons
DAR_Treg_SPF_GF <- as.data.frame(results(dds, contrast=c("condition","ATAC_TregSPF","ATAC_TregGF"), format = c("DataFrame")))
DAR_Treg_SPF_GF <- DAR_Treg_SPF_GF[,c("log2FoldChange","padj")]
colnames(DAR_Treg_SPF_GF) <- c(paste("log2FC_","Treg_SPF_GF",sep=""),paste("padj_","Treg_SPF_GF",sep=""))

DAR_Tconv_SPF_GF <- as.data.frame(results(dds, contrast=c("condition","ATAC_TconSPF","ATAC_TconGF"), format = c("DataFrame")))
DAR_Tconv_SPF_GF <- DAR_Tconv_SPF_GF[,c("log2FoldChange","padj")]
colnames(DAR_Tconv_SPF_GF) <- c(paste("log2FC_","Tconv_SPF_GF",sep=""),paste("padj_","Tconv_SPF_GF",sep=""))

DAR_SPF_Treg_Tconv <- as.data.frame(results(dds, contrast=c("condition","ATAC_TregSPF","ATAC_TconSPF"), format = c("DataFrame")))
DAR_SPF_Treg_Tconv <- DAR_SPF_Treg_Tconv[,c("log2FoldChange","padj")]
colnames(DAR_SPF_Treg_Tconv) <- c(paste("log2FC_","SPF_Treg_Tconv",sep=""),paste("padj_","SPF_Treg_Tconv",sep=""))

DAR_GF_Treg_Tconv <- as.data.frame(results(dds, contrast=c("condition","ATAC_TregGF","ATAC_TconGF"), format = c("DataFrame")))
DAR_GF_Treg_Tconv <- DAR_GF_Treg_Tconv[,c("log2FoldChange","padj")]
colnames(DAR_GF_Treg_Tconv) <- c(paste("log2FC_","GF_Treg_Tconv",sep=""),paste("padj_","GF_Treg_Tconv",sep=""))

#Table of all conditions
DARs <- cbind(id = rownames(DAR_Treg_SPF_GF), DAR_Treg_SPF_GF, DAR_Tconv_SPF_GF , DAR_SPF_Treg_Tconv, DAR_GF_Treg_Tconv)
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

#####
#Peak/gene location
#####


#All features
mm10_TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene

#TSS sites BioMart
mm10 = useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl") 
TSS.mouse.mm10 = getAnnotation(mart=mm10, featureType="TSS")
Exon.mouse.mm10 = getAnnotation(mart=mm10, featureType="Exon")

#mm10 whole database
#Grange object from peak file
gr_DARs_regions <- toGRanges(DARs_regions, names = DARs_regions$id)

#Annotate Peaks
gr_DARs_regions_anno_TSS <- annotatePeakInBatch(gr_DARs_regions, AnnotationData=TSS.mouse.mm10)

#add gene name
gr_DARs_regions_anno_TSS <- addGeneIDs(annotatedPeak=gr_DARs_regions_anno_TSS, 
                                   feature_id_type="ensembl_gene_id",
                       orgAnn="org.Mm.eg.db", 
                       IDs2Add="symbol")

#Generate dataframe
DARs_features <- as.data.frame(gr_DARs_regions_anno_TSS)

##########
#Process RNA-seq data
##########
#####
#Load expression data
#####
#Set directory to load DESeq2 tables
setwd(path_RNAseq_DESeq2)
names = list.files(pattern="*.csv")
myfiles = lapply(names, read.csv)

#Extract names of files for column labelling
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
#assign list names
names(myfiles) <- names

#get list name and attach to column names
list_all <- lapply(seq_along(myfiles),function(i){
  name_i <- names(myfiles[i])
  #colnames(myfiles[[i]])
  table_i <- cbind(myfiles[[i]][,c("GeneSymbol","EnsemblID","EntrezID","log2FoldChange","padj")],myfiles[[i]][,grepl("RPKMcounts", colnames(myfiles[[i]]))])
  colnames(table_i)[4] <- paste("log2FC_", name_i, sep = "")
  colnames(table_i)[5] <- paste("padj_", name_i, sep = "")
  table_i
}
)

#make table from list
data_all = Reduce(function(...) merge(..., all=T), list_all)
data_all <- data_all[!duplicated(data_all$GeneSymbol),]
data_all <- data_all[!(is.na(data_all$GeneSymbol)),]
data_all[is.na(data_all)] <- 0

#####
#Heatmap selected genes
#####
genes_interest <- c("Foxp3","Il2rb","Il2ra","Bcl2","Ccr9","Cd72","Rora","Sema4a","Tnfrsf1b","Cxcr6","Lactb","Itgb1")
data_all_select <- subset(data_all, GeneSymbol %in% genes_interest)

data_heatmap <- cbind(data_all_select$GeneSymbol, data_all_select[,grepl(c("RPKM"), colnames(data_all_select))])
colnames(data_heatmap)[1] <- "GeneSymbol"
rownames(data_heatmap) <- as.character(data_heatmap$GeneSymbol)
data_heatmap <- data_heatmap[,c(2:ncol(data_heatmap))]

data_heatmap <- log2(data_heatmap)

data_heatmap[data_heatmap == -Inf] <- -1
min(data_heatmap)
data_heatmap_matrix <- data.matrix(data_heatmap)

pheatmap(data_heatmap_matrix, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 10, treeheight_col = 30, show_rownames = TRUE, cluster_cols = TRUE,
         scale = "row", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("blue", "white","red"), space="rgb")(128))

#####
#Get all differentially expressed genes
#####
defining_diff_all <- cbind(data_all$GeneSymbol,
                           data_all[,grepl(c("log2FC"), colnames(data_all))],
                           data_all[,grepl(c("padj"), colnames(data_all))])

n_groups = (ncol(defining_diff_all) - 1) / 2
diff_all <- defining_diff_all[1,]

for(i in 1:n_groups){
  #get table of diff. expressed genes using column indexing via n_groups
  diff_all <- rbind(diff_all, subset(defining_diff_all,
                                     (defining_diff_all[,i + 1] < -log2FC_RNA | defining_diff_all[,i + 1] > log2FC_RNA) &
                                       defining_diff_all[,i + 1 + n_groups] < padj))
}
colnames(diff_all)[1] <- c("GeneSymbol")
#eliminate first row
diff_all <- diff_all[!duplicated(diff_all$GeneSymbol),]
diff_all <- diff_all[!(is.na(diff_all$GeneSymbol)),]
#get expression values for diff expressed genes
diff_all_exp <- subset(data_all, data_all$GeneSymbol %in% as.character(diff_all$GeneSymbol)) 
diff_all_exp <- na.omit(diff_all_exp)

#####
#Heatmap all diff genes
#####
data_all_select <- diff_all_exp

data_heatmap <- cbind(data_all_select$GeneSymbol, data_all_select[,grepl(c("RPKM"), colnames(data_all_select))])
colnames(data_heatmap)[1] <- "GeneSymbol"
rownames(data_heatmap) <- as.character(data_heatmap$GeneSymbol)
data_heatmap <- data_heatmap[,c(2:ncol(data_heatmap))]

data_heatmap <- log2(data_heatmap)

data_heatmap[data_heatmap == -Inf] <- -5
min(data_heatmap)
data_heatmap_matrix <- data.matrix(data_heatmap)

pheatmap(data_heatmap_matrix, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 10, treeheight_col = 30, show_rownames = TRUE, cluster_cols = TRUE,
         scale = "row", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("blue", "white","red"), space="rgb")(128))


expression <- data_all[,c("GeneSymbol",
                          "log2FC_GF_Treg_VS_GF_Tconv","log2FC_SPF_Tconv_VS_GF_Tconv",
                          "log2FC_SPF_Treg_VS_GF_Treg",
                          "log2FC_SPF_Treg_VS_SPF_Tconv",
                          "padj_GF_Treg_VS_GF_Tconv","padj_SPF_Tconv_VS_GF_Tconv",
                          "padj_SPF_Treg_VS_GF_Treg",
                          "padj_SPF_Treg_VS_SPF_Tconv")]

print(paste("Number of expressed genes across comparison:",nrow(expression[!is.na(expression$GeneSymbol),])))

######
#PCA analysis
######
#Get expression data
#library(rgl)
#library(calibrate)

diff_all_exp_RPKM <- diff_all_exp[,4:(ncol(diff_all_exp) - n_groups * 2)]

diff_all_exp_RPKM <- log2(diff_all_exp_RPKM)
diff_all_exp_RPKM[diff_all_exp_RPKM == -Inf] <- -8
min(diff_all_exp_RPKM)

model <- prcomp(t(diff_all_exp_RPKM), scale = TRUE)
pca_per_sample <- model$x
#add categorizer
categorize <- c(1,1,1,1,
                2,2,2,2,
                3,3,3,
                4,4,4,
                5,5,5)
Sample_Name <- rownames(pca_per_sample)
pca_per_sample <- as.data.frame(cbind(pca_per_sample, categorize))
pca_per_sample <- cbind(pca_per_sample, Sample_Name)

#Variance per Sample
Prop_Variance <- model$sdev^2/sum(model$sdev^2) *100
PC_1 = round(Prop_Variance[1])
PC_2 = round(Prop_Variance[2])
PC_3 = round(Prop_Variance[3])
PC_4 = round(Prop_Variance[4])

# Add colored points for groups
#par(mfrow=c(1,2))
with(pca_per_sample, plot(PC1,PC2, type = "p", main="PCA De Novo Treg", cex = 2, pch = 16,
                          xlab = paste("PC1,", "Variance described", sep = " ", PC_1, "%"), 
                          ylab = paste("PC2,", "Variance described", sep = " ", PC_2, "%")))
with(subset(pca_per_sample, categorize == 1), points(PC1, PC2, pch=20, col="red", cex = 2.4))
with(subset(pca_per_sample, categorize == 2), points(PC1, PC2, pch=20, col="orange", cex = 2.4))
with(subset(pca_per_sample, categorize == 3), points(PC1, PC2, pch=20, col="blue", cex = 2.4))
with(subset(pca_per_sample, categorize == 4), points(PC1, PC2, pch=20, col="green", cex = 2.4))
#with(subset(pca_per_sample, categorize == 5), points(PC1, PC2, pch=20, col="yellow", cex = 2.4))
#with(pca_per_sample, textxy(PC1, PC2, labs = Sample_Name))
legend('bottomright', c("Tconv_SPF", "Treg_SPF","Treg_GF","Tconv_GF"), pch = 20, col = c("red","orange","blue","green","yellow"))



#####
#Compare RNA and ATAC
#####
#Comparisons of  all
#Prep RNAseq tables
#for the four relevant comparisons
#generate a list that contains vectors with the UP and DOWN regulated genes
RNAseq_UP <- list()
RNAseq_DOWN <- list()
print("Number of DEGs per comparison")
for(i in 1:4){
  index_log2FC = i + 1
  index_padj = i + 5
  print(colnames(expression[index_log2FC]))
  UP <- as.character(expression[(expression[,index_log2FC] >= log2FC_RNA & !is.na(expression[,index_log2FC])) & 
                                (expression[,index_padj] <= padj & !is.na(expression[,index_padj])),]$GeneSymbol)
  print(length(UP))
  DOWN <- as.character(expression[(expression[,index_log2FC] <= -log2FC_RNA & !is.na(expression[,index_log2FC])) & 
                                  (expression[,index_padj] <= padj & !is.na(expression[,index_padj])),]$GeneSymbol)
  print(length(DOWN))
  RNAseq_UP[[i]] <- UP
  RNAseq_DOWN[[i]] <- DOWN
}
#Name list elements for later utilization
names(RNAseq_UP) <- c("GF_Treg_VS_GF_Tconv","SPF_Tconv_VS_GF_Tconv",
                       "SPF_Treg_VS_GF_Treg","SPF_Treg_VS_SPF_Tconv")
names(RNAseq_DOWN) <- c("GF_Treg_VS_GF_Tconv","SPF_Tconv_VS_GF_Tconv",
                        "SPF_Treg_VS_GF_Treg","SPF_Treg_VS_SPF_Tconv")

#Prep ATACseq tables
#for the four relevant comparisons
#generate a list that contains vectors with the UP and DOWN regulated genes
DARs_features_NA <- DARs_features[!is.na(DARs_features$symbol),]
ATACseq_UP <- list()
ATACseq_DOWN <- list()
k = 5
for(i in 1:4){
  k = k + 2 
  index_log2FC = k
  print(colnames(DARs_features_NA[index_log2FC]))
  index_padj = k + 1
  UP <- as.character(DARs_features_NA[(DARs_features_NA[,index_log2FC] >= log2FC_ATAC & !is.na(DARs_features_NA[,index_log2FC])) & 
                                  (DARs_features_NA[,index_padj] <= padj & !is.na(DARs_features_NA[,index_padj])),]$symbol)
  print(length(UP))
  #print(UP)
  DOWN <- as.character(DARs_features_NA[(DARs_features_NA[,index_log2FC] <= -log2FC_ATAC & !is.na(DARs_features_NA[,index_log2FC])) & 
                                    (DARs_features_NA[,index_padj] <= padj & !is.na(DARs_features_NA[,index_padj])),]$symbol)
  print(length(DOWN))
  #print(DOWN)
  ATACseq_UP[[i]] <- UP
  ATACseq_DOWN[[i]] <- DOWN
}
#Name list elements for later utilization
names(ATACseq_UP) <- c("mLN_SPF_GF","pLN_SPF_GF",
                      "SPF_mLN_pLN","GF_mLN_pLN")
names(ATACseq_DOWN) <- c("mLN_SPF_GF","pLN_SPF_GF",
                        "SPF_mLN_pLN","GF_mLN_pLN")
#reorder lists according to RNAseq
ATACseq_UP <- ATACseq_UP[names(RNAseq_UP)]
ATACseq_DOWN <- ATACseq_DOWN[names(RNAseq_DOWN)]

#Check the overlaps
#empty table for storage
summary_comparison = data.frame(matrix("", ncol = 8, nrow = 0))
#empty lists for storage of GeneNames
Plus_Plus <- list()
Minus_Minus <- list()
Plus_Zero <- list()
Zero_Plus <- list()
Minus_Zero <- list()
Zero_Minus <- list()

for(i in 1:4){
  #Grab corresponding Gene lists
  RNA_UP_i <- RNAseq_UP[[i]]
  RNA_DOWN_i <- RNAseq_DOWN[[i]]
  ATAC_UP_i <- ATACseq_UP[[i]]
  ATAC_DOWN_i <- ATACseq_DOWN[[i]]
  print("index")
  print(i)
  #Track
  print(names(RNAseq_UP)[i])
  print(names(ATACseq_UP)[i])
  #Overlap in line
  Overlap_in_UP <- sort(RNA_UP_i[RNA_UP_i %in% ATAC_UP_i])
  print(length(Overlap_in_UP))
  Overlap_in_DOWN <- sort(RNA_DOWN_i[RNA_DOWN_i %in% ATAC_DOWN_i])
  print(length(Overlap_in_DOWN))
  
  #Overlap out line
  Overlap_out_UP <- RNA_UP_i[RNA_UP_i %in% ATAC_DOWN_i]
  print(length(Overlap_out_UP))
  Overlap_out_DOWN <- RNA_DOWN_i[RNA_DOWN_i %in% ATAC_UP_i]
  print(length(Overlap_out_DOWN))
  
  #RNA UP but no peak
  NonOverlap_UP_RNA <- RNA_UP_i[!RNA_UP_i %in% c(Overlap_in_UP,Overlap_out_UP)]
  print(length(NonOverlap_UP_RNA))
  #ATAC UP but no expression
  NonOverlap_UP_ATAC <- ATAC_UP_i[!ATAC_UP_i %in% c(Overlap_in_UP,Overlap_out_UP)]
  print(length(NonOverlap_UP_ATAC))
  
  #RNA DOWN but no peak
  NonOverlap_DOWN_RNA <- RNA_DOWN_i[!RNA_DOWN_i %in% c(Overlap_in_DOWN,Overlap_out_DOWN)]
  print(length(NonOverlap_DOWN_RNA))
  #ATAC DOWN but no expression
  NonOverlap_DOWN_ATAC <- ATAC_DOWN_i[!ATAC_DOWN_i %in% c(Overlap_in_DOWN,Overlap_out_DOWN)]
  print(length(NonOverlap_DOWN_ATAC))
  
  #Store counts
  comparison <- c(length(Overlap_in_UP),
                   length(Overlap_in_DOWN),
                   length(Overlap_out_DOWN),
                   length(Overlap_out_UP),
                   length(NonOverlap_UP_ATAC),
                   length(NonOverlap_UP_RNA),
                   length(NonOverlap_DOWN_ATAC),
                   length(NonOverlap_DOWN_RNA))
  #compile into table
  summary_comparison <- rbind(summary_comparison, comparison)
  
  #Store GeneSymbols
  Plus_Plus[[i]] <- Overlap_in_UP
  Minus_Minus[[i]] <- Overlap_in_DOWN
  Plus_Zero[[i]] <- NonOverlap_UP_ATAC
  Zero_Plus[[i]] <- NonOverlap_UP_RNA
  Minus_Zero[[i]] <- NonOverlap_DOWN_ATAC
  Zero_Minus[[i]] <- NonOverlap_DOWN_RNA
  
}

#Annotion of +, - and o combinations
# +/+ open/DEG
# -/- closed/suppressed
# o/- noDAR/suppressed
# +/o DAR/noDEG
colnames(summary_comparison) <- c("+/+","-/-",
                                  "+/-","-/+",
                                  "+/o","o/+",
                                  "-/o","o/-")
rownames(summary_comparison) <- names(RNAseq_UP)
print(summary_comparison)

#####
#Get genes with defined peaks for Motif enrichment
#####
# Rational: 
#1) DEGs that don't overlap with a DAR are potentially influenced by TF networks
#   But Promotor regions needs to contain a peak, aka be open
#   get o/+ and o/-
#SPF
No_DAR_UP_SPF <- unique(c(Zero_Plus[[1]]))
No_DAR_DOWN_SPF <- unique(c(Zero_Minus[[1]]))
#GF
No_DAR_UP_GF <- unique(c(Zero_Plus[[2]]))
No_DAR_DOWN_GF <- unique(c(Zero_Minus[[2]]))

# Rational: 
#2) Regions that are open/closed could regulate expression upon activation/inflammation
#   driven by TF networks
#   But Promotor regions needs to contain a peak, aka be open
#   get +/o and -/o
#SPF
Open_NoDEG_SPF <- unique(c(Plus_Zero[[1]]))
Closed_NoDEG_SPF <- unique(c(Minus_Zero[[1]]))
#GF
Open_NoDEG_GF <- unique(c(Plus_Zero[[2]]))
Closed_NoDEG_GF <- unique(c(Minus_Zero[[2]]))

#Which have a peak in the TSS region
downstream <- 200
upstream <- 10000
# or "overlapStart"

#Extract promotor region
#genes <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
prom <- promoters(TSS.mouse.mm10, upstream=upstream, downstream=downstream)

#Get a symbol to TSS and define promoter region
gr_TSS.mouse.mm10_ensembl <- annotatePeakInBatch(TSS.mouse.mm10, AnnotationData=TSS.mouse.mm10)
gr_TSS.mouse.mm10_symbol <- addGeneIDs(annotatedPeak=gr_TSS.mouse.mm10_ensembl, 
                                    feature_id_type="ensembl_gene_id",
                                    orgAnn="org.Mm.eg.db", 
                                    IDs2Add="symbol")
gr_prom <- promoters(gr_TSS.mouse.mm10_symbol, upstream = downstream, downstream = upstream)
#Make table
t_promotors <- as.data.frame(gr_prom)
rownames(t_promotors) <- c()


#Find TSS for identified genes
# 1) background
#all genes expressed but not differentially
expression_NA <- na.omit(expression)
only_expressed <- as.character(subset(expression_NA, (abs(SPF_mLN_pLN_log2FoldChange) < log2FC_RNA & SPF_mLN_pLN_padj > padj) &
                                        (abs(GF_mLN_pLN_log2FoldChange) < log2FC_RNA & GF_mLN_pLN_padj > padj) &
                                        (abs(mLN_SPF_GF_log2FoldChange) < log2FC_RNA & mLN_SPF_GF_padj > padj) &
                                        (abs(pLN_SPF_GF_log2FoldChange) < log2FC_RNA & pLN_SPF_GF_padj > padj))$GeneSymbol)

#all Peaks detected in TSS region but differentially
only_open <- as.character(subset(DARs_features, (abs(log2FC_mLN_SPF_GF) < log2FC_ATAC & padj_mLN_SPF_GF > padj) &
                                        (abs(log2FC_pLN_SPF_GF) < log2FC_ATAC & padj_pLN_SPF_GF > padj) &
                                        (abs(log2FC_SPF_mLN_pLN) < log2FC_ATAC & padj_SPF_mLN_pLN > padj) &
                                        (abs(log2FC_GF_mLN_pLN) < log2FC_ATAC & padj_GF_mLN_pLN > padj))$symbol)

background_genes <- na.omit(unique(only_expressed, only_open))
background_promotors_open <- subset(t_promotors, symbol %in% background_genes)

# 2) Key genes
#   a) list gene modules not regulated by accessibility
TF_DARs <- list(No_DAR_UP_SPF,No_DAR_DOWN_SPF,No_DAR_UP_GF,No_DAR_DOWN_GF)
#   b) list gene modules potentially regulated by activation
TF_Activity <- list(Open_NoDEG_SPF,Closed_NoDEG_SPF,Open_NoDEG_GF,Closed_NoDEG_GF)

#instigate storage list for:
# a) not accessibility
t_Prom_Feat_TF_reg <- list()
for(i in seq(length(TF_DARs))){
  TF_DARs_i <- TF_DARs[[i]]
  t_Prom_Feat_i <- subset(t_promotors, symbol %in% TF_DARs_i)
  t_Prom_Feat_TF_reg[[i]] <- t_Prom_Feat_i
}
names(t_Prom_Feat_TF_reg) <- c("No_DAR_UP_SPF","No_DAR_DOWN_SPF","No_DAR_UP_GF","No_DAR_DOWN_GF")

#instigate storage list for:
# b) Activation dependent
t_Prom_Feat_TF_act <- list()
for(i in seq(length(TF_Activity))){
  TF_Activity_i <- TF_Activity[[i]]
  t_Prom_Feat_i <- subset(t_promotors, symbol %in% TF_Activity_i)
  t_Prom_Feat_TF_act[[i]] <- t_Prom_Feat_i
}
names(t_Prom_Feat_TF_act) <- c("Open_NoDEG_SPF","Closed_NoDEG_SPF","Open_NoDEG_GF","Closed_NoDEG_GF")

#Regions
# a) not accessibility
for(i in seq(length(t_Prom_Feat_TF_reg))){
  t_Prom_Feat_TF_reg_i <- t_Prom_Feat_TF_reg[[i]]
  t_Prom_Feat_TF_reg_i <- t_Prom_Feat_TF_reg_i[,c("feature","symbol","seqnames","start","end")]
  names(t_Prom_Feat_TF_reg_i)[1] <- c("Acc")
  name_i <- names(t_Prom_Feat_TF_reg)[i]
  print(nrow(t_Prom_Feat_TF_reg_i))
  print(name_i)
  rownames(t_Prom_Feat_TF_act_i) <- c()
  write.table(t_Prom_Feat_TF_reg_i, paste(path_output_Motif, "/", name_i,".txt", sep = ""),
              row.names = FALSE,quote=FALSE,sep="\t")
}

# b) Activation dependent
for(i in seq(length(t_Prom_Feat_TF_act))){
  t_Prom_Feat_TF_act_i <- t_Prom_Feat_TF_act[[i]]
  t_Prom_Feat_TF_act_i <- t_Prom_Feat_TF_act_i[,c("feature","symbol","seqnames","start","end")]
  names(t_Prom_Feat_TF_act_i)[1] <- c("Acc")
  name_i <- names(t_Prom_Feat_TF_act)[i]
  print(nrow(t_Prom_Feat_TF_act_i))
  print(name_i)
  rownames(t_Prom_Feat_TF_act_i) <- c()
  write.table(t_Prom_Feat_TF_act_i, paste(path_output_Motif, "/", name_i,".txt", sep = ""),
              row.names = FALSE,quote=FALSE,sep="\t")
}

#Homer Motif only needs Gene identifier
#Background
background_gene_symbol <- background_promotors_open[,c("feature","symbol","seqnames","start","end")]
names(background_gene_symbol)[1] <- c("Acc")
rownames(background_gene_symbol) <- NULL
write.table(background_gene_symbol, paste(path_output_Motif, "/", "background_genes.txt", sep = ""),
            row.names = FALSE,quote=FALSE,sep="\t")



#####
#Get peaks of defined genes for Motif enrichment
#####
# a) not accessibility
#SPF
No_DAR_UP_SPF <- unique(c(Zero_Plus[[1]]))
No_DAR_DOWN_SPF <- unique(c(Zero_Minus[[1]]))
#GF
No_DAR_UP_GF <- unique(c(Zero_Plus[[2]]))
No_DAR_DOWN_GF <- unique(c(Zero_Minus[[2]]))
#list gene modules not regulated by accessibility
TF_DARs <- list(No_DAR_UP_SPF,No_DAR_DOWN_SPF,No_DAR_UP_GF,No_DAR_DOWN_GF)
names(TF_DARs) <- c("No_DAR_UP_SPF","No_DAR_DOWN_SPF","No_DAR_UP_GF","No_DAR_DOWN_GF")
#storing list for regions
for(i in seq(length(TF_DARs))){
  #load genes
  DARs_genes_i <- TF_DARs[[i]]
  #identify peaks for genes
  DARs_features_i <- subset(DARs_features_NA, symbol %in% DARs_genes_i)
  #Use only genes within Promotor regions
  DARs_features_TSS_i <- subset(DARs_features_i, (distancetoFeature <= downstream & distancetoFeature >= -upstream))
  print(names(TF_DARs)[i])
  print(nrow(DARs_features_TSS_i))
  #Generate Granges object
  gr_DARs_features_TSS_i <- toGRanges(DARs_features_TSS_i)
  #make width bigger by 50 each direction
  start(gr_DARs_features_TSS_i) <- start(gr_DARs_features_TSS_i) - 50
  end(gr_DARs_features_TSS_i) <- end(gr_DARs_features_TSS_i) + 50
  #Save generate .bed compatible for homer
  DARs_features_TSS_i <- as.data.frame(gr_DARs_features_TSS_i)
  DARs_features_TSS_bed_i <- DARs_features_TSS_i[,c("seqnames","start","end","id","feature_strand")] 
  DARs_features_TSS_bed_i$blankVar <- NA
  DARs_features_TSS_bed_i <- DARs_features_TSS_bed_i[,c("seqnames","start","end","id","blankVar","feature_strand")] 
  DARs_features_TSS_bed_i$blankVar <- c("")
    
  #Write .bed compatible for homer
  print(head(DARs_features_TSS_bed_i))
  write.table(DARs_features_TSS_bed_i, file=paste(path_output_Motif, "/", names(TF_DARs)[i],".bed", sep = ""),
              quote=F, sep="\t", row.names=F, col.names=F)
}

# b) Activation dependent
#SPF
Open_NoDEG_SPF <- unique(c(Plus_Zero[[1]]))
Closed_NoDEG_SPF <- unique(c(Minus_Zero[[1]]))
#GF
Open_NoDEG_GF <- unique(c(Plus_Zero[[2]]))
Closed_NoDEG_GF <- unique(c(Minus_Zero[[2]]))
#   b) list gene modules potentially regulated by activation
TF_Activity <- list(Open_NoDEG_SPF,Closed_NoDEG_SPF,Open_NoDEG_GF,Closed_NoDEG_GF)
names(TF_Activity) <- c("Open_NoDEG_SPF","Closed_NoDEG_SPF","Open_NoDEG_GF","Closed_NoDEG_GF")
#storing list for regions
for(i in seq(length(TF_Activity))){
  #load genes
  PotAct_genes_i <- TF_Activity[[i]]
  #identify peaks for genes
  PotAct_features_i <- subset(DARs_features_NA, symbol %in% PotAct_genes_i)
  #Use only genes within Promotor regions
  PotAct_features_TSS_i <- subset(PotAct_features_i, (distancetoFeature <= downstream & distancetoFeature >= -upstream))
  print(names(TF_Activity)[i])
  print(nrow(PotAct_features_TSS_i))
  #Generate Granges object
  gr_PotAct_features_TSS_i <- toGRanges(PotAct_features_TSS_i)
  #make width bigger by 50 each direction
  start(gr_PotAct_features_TSS_i) <- start(gr_PotAct_features_TSS_i) - 50
  end(gr_PotAct_features_TSS_i) <- end(gr_PotAct_features_TSS_i) + 50
  #Save generate .bed compatible for homer
  PotAct_features_TSS_i <- as.data.frame(gr_PotAct_features_TSS_i)
  PotAct_features_TSS_bed_i <- PotAct_features_TSS_i[,c("seqnames","start","end","id","feature_strand")] 
  PotAct_features_TSS_bed_i$blankVar <- NA
  PotAct_features_TSS_bed_i <- PotAct_features_TSS_bed_i[,c("seqnames","start","end","id","blankVar","feature_strand")] 
  PotAct_features_TSS_bed_i$blankVar <- c("")
  
  #Write .bed compatible for homer
  #print(head(PotAct_features_TSS_bed_i))
  write.table(PotAct_features_TSS_bed_i, file=paste(path_output_Motif, "/", names(TF_Activity)[i],".bed", sep = ""),
              quote=F, sep="\t", row.names=F, col.names=F)
}

############
#GO analysis
############
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
idfound <- expression$GeneSymbol %in% mappedRkeys(org.Mm.egSYMBOL)
SYMBOL <- toTable(org.Mm.egSYMBOL)
head(SYMBOL)

#Expression DEGs
m <- match(expression$GeneSymbol, SYMBOL$symbol)
GENE_ID <- SYMBOL$gene_id[m]
expression <- cbind(expression, GENE_ID)

#Accessibility DARs
m <- match(DARs_features_NA$symbol, SYMBOL$symbol)
GENE_ID <- SYMBOL$gene_id[m]
DARs_features_NA <- cbind(DARs_features_NA, GENE_ID)

#Define Gene Universe for which the pValue is set to 1.0
geneNames <- unique(c(levels(as.factor(DARs_features_NA$GENE_ID)),levels(as.factor(expression$GENE_ID))))
#genes_pval_background <- rep(1.0, length(geneNames))
#names(genes_pval_background) <- geneNames


#####
#Perform GO analysis single
#####
# a) not accessibility-------------------------------------
#put all gene vectors in one list
# Note: Input required
# Add the gene sets you are interested in
my_gene_groups <- list(Plus_Plus[[1]],Plus_Plus[[2]],
                       Minus_Minus[[1]],Minus_Minus[[2]],
                       Plus_Zero[[1]],Plus_Zero[[2]],
                       Zero_Plus[[1]],Zero_Plus[[2]],
                       Minus_Zero[[1]],Minus_Zero[[2]],
                       Zero_Minus[[1]],Zero_Minus[[2]])

names(my_gene_groups) <- c("SPF_mp_PP","GF_mp_PP",
                          "SPF_mp_MM","GF_mp_MM",
                          "SPF_mp_PZ","GF_mp_PZ",
                          "SPF_mp_ZP","GF_mp_ZP",
                          "SPF_mp_MZ","GF_mp_MZ",
                          "SPF_mp_MZ","GF_mp_MZ")

#Replace GeneSymbols with GENE_IDs and standard p-value of 0.01
#This is only legit if Fisher-Test is used for GO analysis
l_gene_pval <- list()
for(i in seq(length(my_gene_groups))){
  my_gene_group_i <- my_gene_groups[[i]]
  #replace Symbol with GENE_ID
  m <- match(my_gene_group_i, SYMBOL$symbol)
  GENE_ID <- SYMBOL$gene_id[m]
  my_gene_group_i <- GENE_ID

  #add standard p-Value of 0.01
  genes_pval_i <- rep(0.01, length(my_gene_group_i))
  
  #make names vector
  names(genes_pval_i) <- my_gene_group_i
  
  l_gene_pval[[i]] <- genes_pval_i
}


#get one gene list
#l_gene_pval <- l_all_gene_pval[[1]]

#gene lists
l_gene_List <- list()
for(i in 1:length(l_gene_pval)){
  geneList_i <- factor(as.integer(geneNames %in% names(l_gene_pval[[i]])))
  names(geneList_i) <- geneNames
  l_gene_List[[i]] <- geneList_i 
}


List_allRes <- list()
#Do  GO statistics for all gene lists
for(i in 1:length(l_gene_List)){
  #Access the gene lists and p-values for the differentially expressed genes
  geneList <- l_gene_List[[i]]
  pvalue_of_interest <- l_gene_pval[[i]]
  
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
  
  #Fisher's works with weight
  test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
  resultWeight <- getSigGroups(GOdata, test.stat)
  #KS works with elim but not with weight
  test.stat <- new("elimScore", testStatistic = GOKSTest, name = "Fisher test", cutOff = 0.01)
  resultElim <- getSigGroups(GOdata, test.stat)
  #runTest a high level interface for testing Fisher
  resultFis <- runTest(GOdata, statistic = "fisher")
  #Kolmogorov-Smirnov includes p-values
  test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
  resultKS <- getSigGroups(GOdata,  test.stat)
  #runTest a high level interface for testing KS
  elim.ks <- runTest(GOdata, algorithm = "elim", statistic = "ks")
  
  #make table
  allRes <- GenTable(GOdata, classic = resultFis, KS = resultKS, weight = resultWeight,
                     orderBy = "weight", ranksOf = "KS", topNodes = number_nodes)
  #make list of result tables
  List_allRes[[i]] <- allRes
}

####
#Find the Top GOs from the lists
####
#Build table with Top GOs according to "weight"
List_TopGOs <- list()
for(i in 1:length(List_allRes)){
  table_i <- List_allRes[[i]]
  table_tophits <- subset(table_i, weight < 0.0025)
  #print(i)
  #print(nrow(table_tophits))
  List_TopGOs[[i]] <- table_tophits
}

#collect all TopGos
TopGOs_vector <- c()
for(i in 1:length(List_TopGOs)){
  table_i <- List_TopGOs[[i]]
  TopGOs_i <- as.character(table_i$GO.ID)
  TopGOs_vector <- c(TopGOs_vector, TopGOs_i) 
}

#Condense tables and sort by GO.ID
l_topGO <- list()
for(i in 1:length(List_allRes)){
  table_i <- List_allRes[[i]]
  table_i_subset <- subset(table_i, table_i$GO.ID %in% TopGOs_vector)
  table_i_subset_GO_weight <- table_i_subset[,c("GO.ID","weight", "Term")]
  k = i - 1
  colnames(table_i_subset_GO_weight)[2] <- names(my_gene_groups)[i]
  table_i_subset_GO_weight[,2] <- as.numeric(table_i_subset_GO_weight[,2])
  table_i_subset_GO_weight <- table_i_subset_GO_weight[order(table_i_subset_GO_weight$"GO.ID", decreasing=TRUE), ]
  l_topGO[[i]] <- table_i_subset_GO_weight
}
#merger
data_TopGO_weight = Reduce(function(...) merge(..., all=T), l_topGO)
#Delete NAs
data_TopGO_weight <- na.omit(data_TopGO_weight)

####
#Make GO comparison heatmap
####
data_heatmap <- data_TopGO_weight
row.names(data_heatmap) <- paste(data_heatmap$"GO.ID", data_heatmap$Term, sep = "_")
data_heatmap <- data_heatmap[,3:ncol(data_heatmap)]

min(data_heatmap)
data_heatmap <- -log10(data_heatmap)
#Common Minimum at 5
data_heatmap[data_heatmap > 5] <- 5
data_heatmap_matrix <- data.matrix(data_heatmap)
title <- c("GO ATAC/RNA-seq overlap")

#reorder according to Heatmap in main figure
#data_heatmap_matrix <- data_heatmap_matrix[,subset_order]

pheatmap(data_heatmap_matrix, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = FALSE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("white", "grey","deepskyblue2"), space="rgb")(128),
         main = title)






















#####
# PCA ---------------------------------------------------------------------
#####
#Continue here 
condition <- factor(substr(colnames(counts.mat),1,nchar(colnames(counts.mat))-1))

dds <- DESeqDataSetFromMatrix(counts.mat, DataFrame(condition),~condition)
dds.factors <- estimateSizeFactors(dds)
count.norm <- log2(1+counts(dds.factors, normalized=T))
count.norm <- data.frame(count.norm)
#count.norm$id <- rownames(count.norm)


pca <- prcomp(t(count.norm[,-97]), scale. =T)
prop <- summary(pca)$importance[2,c(1,2)]*100
PC <- data.frame(pca$x[,c(1,2)])
#PC$timepoint <- rep(c("t0", "t14"), 48)
PC$libsize <- colSums(counts.mat)
PC$enrichment <- enrich.table$percent
#outliers <- c("LA003_t0", "LA016_t0", "LA016_t14","LA035_t14",
 #             "LA040_t14","LA042_t0", "LA046_t0","LA052_t14")
#out <- PC[which(rownames(PC) %in% outliers),]


ggplot(PC) + geom_point(aes(PC1, PC2, color=libsize), size =2) + 
  #geom_text(data=out, aes(x=out$PC1, y=out$PC2, 
   #                       label=strsplit2(rownames(out), "_")[,1]), nudge_x = 15) +
  #labs(x= paste0("PC1 (", prop[1],"%)"), y= paste0("PC2 (", prop[2],"%)"),
   #    color ="Library size", shape ="Timepoint", title ="PCA on normalized ATACseq counts")+
  theme(plot.title = element_text(hjust=0.5))
  

#ggplot(PC) + geom_point(aes(PC1, PC2, shape=timepoint, color=enrichment),  size=2) + 
 # geom_text(data=out, aes(x=out$PC1, y=out$PC2, 
  #                        label=strsplit2(rownames(out), "_")[,1]), nudge_x = 15) +
  #labs(x= paste0("PC1 (", prop[1],"%)"), y= paste0("PC2 (", prop[2],"%)"),
   #    color ="Enrichment", shape ="Timepoint", title ="PCA on normalized ATACseq counts")+
  #theme(plot.title = element_text(hjust=0.5))
  

#outliers2 <- outliers[-c(2,3,6)]
#count.no.outliers <- count.norm[, -c(which(colnames(count.norm) %in% outliers2)) ]
#pca.no.outliers <- prcomp(t(count.no.outliers))

#PC.no.outliers <- PC[-c(which(rownames(PC) %in% outliers2)), ]
#PC.no.outliers$PC1 <- pca.no.outliers$x[,1]
#PC.no.outliers$PC2 <- pca.no.outliers$x[,2]

#prop2 <- summary(pca.no.outliers)$importance[2,c(1,2)]*100

#out2 <- PC.no.outliers[rownames(PC.no.outliers) %in% c("LA016_t0", "LA016_t14", "LA042_t0"),]

#ggplot(PC.no.outliers) + geom_point(aes(PC1, PC2, shape=timepoint, color=enrichment),  size=2) + 
 # geom_text(data=out2, aes(x=out2$PC1, y=out2$PC2, 
  #                        label=strsplit2(rownames(out2), "_")[,1]), nudge_x = 15) +
  #labs(x= paste0("PC1 (", prop2[1],"%)"), y= paste0("PC2 (", prop2[2],"%)"),
   #    color ="Enrichment", shape ="Timepoint", title ="PCA on normalized ATACseq counts")+
  #theme(plot.title = element_text(hjust=0.5))


#count.no.outliers2 <- count.norm[, -c(which(colnames(count.norm) %in% outliers)) ]
#pca.no.outliers2 <- prcomp(t(count.no.outliers2))

#PC.no.outliers2 <- PC[-c(which(rownames(PC) %in% outliers)), ]
#PC.no.outliers2$PC1 <- pca.no.outliers2$x[,1]
#PC.no.outliers2$PC2 <- pca.no.outliers2$x[,2]

#prop2 <- summary(pca.no.outliers2)$importance[2,c(1,2)]*100

#ggplot(PC.no.outliers2) + geom_point(aes(PC1, PC2, shape=timepoint, color=enrichment),  size=2) + 
 # labs(x= paste0("PC1 (", prop2[1],"%)"), y= paste0("PC2 (", prop2[2],"%)"),
  #     color ="Enrichment", shape ="Timepoint", title ="PCA on normalized ATACseq counts")+
  #theme(plot.title = element_text(hjust=0.5))


#####
# clustering --------------------------------------------------------------
#####

d <- dist(t(count.norm[,1:6]))
clust <- hclust(d, method="ward.D")
plot(clust, cex=0.8)

#d2 <- dist(t(count.no.outliers))
#clust2 <- hclust(d2, method="ward.D")
#plot(clust2, main = "Cluster dendrogram (no outliers)", cex=0.8)

#d3 <- dist(t(count.no.outliers2))
#clust3 <- hclust(d3, method="ward.D")
#plot(clust3, main = "Cluster dendrogram (no outliers)", cex=0.8)

#outliers

# png(filename = "peak_distr1.png", width = 1700, height = 900)
# boxplot(log2(1 + 2*(counts[,-1])), main="Peak distribution before normalisation and filtering", las ="2")
# dev.off()

# png(filename = "peak_distr2.png", width = 1700, height = 900)
# boxplot(counts.norm, main="Peak distribution after library size normalisation", las ="2")
# dev.off()





#Make table




#Overlap in line with hypothesis
Overlap_in_UP <- sort(RNAseq_SPF_mLN_pLN_UP_name[RNAseq_SPF_mLN_pLN_UP_name %in% DARs_features_SPF_mLN_pLN_UP_name])
print(length(Overlap_in_UP))
Overlap_in_DOWN <- sort(RNAseq_SPF_mLN_pLN_DOWN_name[RNAseq_SPF_mLN_pLN_DOWN_name %in% DARs_features_SPF_mLN_pLN_DOWN_name])
print(length(Overlap_in_DOWN))

#Overlap out of line with standard hypothesis
Overlap_out_UP <- RNAseq_SPF_mLN_pLN_UP_name[RNAseq_SPF_mLN_pLN_UP_name %in% DARs_features_SPF_mLN_pLN_DOWN_name]
print(length(Overlap_out_UP))
Overlap_out_DOWN <- RNAseq_SPF_mLN_pLN_DOWN_name[RNAseq_SPF_mLN_pLN_DOWN_name %in% DARs_features_SPF_mLN_pLN_UP_name]
print(length(Overlap_out_DOWN))

#RNA UP but no peak
NonOverlap_UP_RNA <- RNAseq_SPF_mLN_pLN_UP_name[!RNAseq_SPF_mLN_pLN_UP_name %in% c(Overlap_in_UP,Overlap_out_UP)]
print(length(NonOverlap_UP_RNA))
#ATAC UP but no expression
NonOverlap_UP_ATAC <- DARs_features_SPF_mLN_pLN_UP_name[!DARs_features_SPF_mLN_pLN_UP_name %in% c(Overlap_in_UP,Overlap_out_UP)]
print(length(NonOverlap_UP_ATAC))

#RNA DOWN but no peak
NonOverlap_DOWN_RNA <- RNAseq_SPF_mLN_pLN_DOWN_name[!RNAseq_SPF_mLN_pLN_DOWN_name %in% c(Overlap_in_DOWN,Overlap_out_DOWN)]
print(length(NonOverlap_DOWN_RNA))
#ATAC DOWN but no expression
NonOverlap_DOWN_ATAC <- DARs_features_SPF_mLN_pLN_DOWN_name[!DARs_features_SPF_mLN_pLN_DOWN_name %in% c(Overlap_in_DOWN,Overlap_out_DOWN)]
print(length(NonOverlap_DOWN_ATAC))

RNAseq_SPF_mLN_pLN_UP <- subset(expression, SPF_mLN_pLN_log2FoldChange >= log2FC_RNA & SPF_mLN_pLN_padj <= padj)
RNAseq_SPF_mLN_pLN_UP_name <- as.character(RNAseq_SPF_mLN_pLN_UP$GeneSymbol)
print(length(RNAseq_SPF_mLN_pLN_UP_name))

RNAseq_SPF_mLN_pLN_DOWN <- subset(expression, SPF_mLN_pLN_log2FoldChange <= -log2FC_RNA & SPF_mLN_pLN_padj <= padj)
RNAseq_SPF_mLN_pLN_DOWN_name <- as.character(RNAseq_SPF_mLN_pLN_DOWN$GeneSymbol)
print(length(RNAseq_SPF_mLN_pLN_DOWN_name))

DARs_features_SPF_mLN_pLN_UP <- subset(DARs_features, log2FC_SPF_mLN_pLN >= log2FC_ATAC & padj_SPF_mLN_pLN <= padj)
DARs_features_SPF_mLN_pLN_UP <- DARs_features_SPF_mLN_pLN_UP[!is.na(DARs_features_SPF_mLN_pLN_UP$symbol),]
DARs_features_SPF_mLN_pLN_UP_name <- as.character(DARs_features_SPF_mLN_pLN_UP$symbol[!duplicated(DARs_features_SPF_mLN_pLN_UP$symbol)])
print(length(DARs_features_SPF_mLN_pLN_UP_name))

DARs_features_SPF_mLN_pLN_DOWN <- subset(DARs_features, log2FC_SPF_mLN_pLN <= -log2FC_ATAC & padj_SPF_mLN_pLN <= padj)
DARs_features_SPF_mLN_pLN_DOWN <- DARs_features_SPF_mLN_pLN_DOWN[!is.na(DARs_features_SPF_mLN_pLN_DOWN$symbol),]
DARs_features_SPF_mLN_pLN_DOWN_name <- as.character(DARs_features_SPF_mLN_pLN_DOWN$symbol[!duplicated(DARs_features_SPF_mLN_pLN_DOWN$symbol)])
print(length(DARs_features_SPF_mLN_pLN_DOWN_name))


#Overlap in line with hypothesis
Overlap_in_UP <- sort(RNAseq_SPF_mLN_pLN_UP_name[RNAseq_SPF_mLN_pLN_UP_name %in% DARs_features_SPF_mLN_pLN_UP_name])
print(length(Overlap_in_UP))
Overlap_in_DOWN <- sort(RNAseq_SPF_mLN_pLN_DOWN_name[RNAseq_SPF_mLN_pLN_DOWN_name %in% DARs_features_SPF_mLN_pLN_DOWN_name])
print(length(Overlap_in_DOWN))

#Overlap out of line with standard hypothesis
Overlap_out_UP <- RNAseq_SPF_mLN_pLN_UP_name[RNAseq_SPF_mLN_pLN_UP_name %in% DARs_features_SPF_mLN_pLN_DOWN_name]
print(length(Overlap_out_UP))
Overlap_out_DOWN <- RNAseq_SPF_mLN_pLN_DOWN_name[RNAseq_SPF_mLN_pLN_DOWN_name %in% DARs_features_SPF_mLN_pLN_UP_name]
print(length(Overlap_out_DOWN))

#RNA UP but no peak
NonOverlap_UP_RNA <- RNAseq_SPF_mLN_pLN_UP_name[!RNAseq_SPF_mLN_pLN_UP_name %in% c(Overlap_in_UP,Overlap_out_UP)]
print(length(NonOverlap_UP_RNA))
#ATAC UP but no expression
NonOverlap_UP_ATAC <- DARs_features_SPF_mLN_pLN_UP_name[!DARs_features_SPF_mLN_pLN_UP_name %in% c(Overlap_in_UP,Overlap_out_UP)]
print(length(NonOverlap_UP_ATAC))

#RNA DOWN but no peak
NonOverlap_DOWN_RNA <- RNAseq_SPF_mLN_pLN_DOWN_name[!RNAseq_SPF_mLN_pLN_DOWN_name %in% c(Overlap_in_DOWN,Overlap_out_DOWN)]
print(length(NonOverlap_DOWN_RNA))
#ATAC DOWN but no expression
NonOverlap_DOWN_ATAC <- DARs_features_SPF_mLN_pLN_DOWN_name[!DARs_features_SPF_mLN_pLN_DOWN_name %in% c(Overlap_in_DOWN,Overlap_out_DOWN)]
print(length(NonOverlap_DOWN_ATAC))




nrows <- 20; ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
rowRanges <- GRanges(rep(c("chr1", "chr2"), c(5, 15)),
                     IRanges(sample(1000L, 20), width=100),
                     strand=Rle(c("+", "-"), c(12, 8)))
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])
rse0 <- SummarizedExperiment(assays=SimpleList(counts=counts),
                             rowRanges=rowRanges, colData=colData)

rse1 <- shift(rse0, 1)
stopifnot(identical(
  rowRanges(rse1),
  shift(rowRanges(rse0), 1)
))

se2 <- narrow(rse0, start=10, end=-15)
stopifnot(identical(
  rowRanges(se2),
  narrow(rowRanges(rse0), start=10, end=-15)
))

se3 <- resize(rse0, width=75)
stopifnot(identical(
  rowRanges(se3),
  resize(rowRanges(rse0), width=75)
))

se4 <- flank(rse0, width=20)
stopifnot(identical(
  rowRanges(se4),
  flank(rowRanges(rse0), width=20)
))

se5 <- restrict(rse0, start=200, end=700, keep.all.ranges=TRUE)
stopifnot(identical(
  rowRanges(se5),
  restrict(rowRanges(rse0), start=200, end=700, keep.all.ranges=TRUE)
))

