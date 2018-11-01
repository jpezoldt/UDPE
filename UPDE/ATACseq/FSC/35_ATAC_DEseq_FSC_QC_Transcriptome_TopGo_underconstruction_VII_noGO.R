#!/usr/bin/env Rscript

# Author: Joern Pezoldt
# Advice: Maria Litovchenko, Vincent Gardeaux
# 12.08.2018
# Function:
# 1) 
# 

#####
#Libraries & PATHS & Data & Global_Variables
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


#####
#
#####

# Input required: Set directory
# 1) Counts per peak per replicate
path_input <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/homer/Overlap_Group_Merged/Run_3_Outlier_excluded"
# 2) Common Peak regions
path_common_peak_regions <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/peaks/broad/MACS2_broad/Run_3/ATAC_FSC_all_broad_merged_peaks.bed"
# 3) Output path
path_output <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/DESeq2/Run_3_Outlier_excluded"
path_output_Motif <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/Motif/Regions_of_Interest"
# 4) RNAseq DESeq2 analysis
path_RNAseq_DESeq2 <- "/home/pezoldt/NAS2/pezoldt/Analysis/RNAseq/2014_HZI_FSC_endo_DESeq2"

# Input required: 
name = "ATAC_FSC_all"
paste(path_input, "/",name,".txt",sep="")
#read count table
counts <- fread(paste(path_input, "/",name,".txt",sep=""))

#Global variables
#Param
log2FC_RNA = 1.0
padj = 0.05

log2FC_ATAC = 1.0

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
#DEseq2
#####
#Define conditons 
condition <- factor(substr(colnames(counts.mat),1,nchar(colnames(counts.mat))-1))
dds <- DESeqDataSetFromMatrix(counts.mat, DataFrame(condition),~condition)
dds.factors <- estimateSizeFactors(dds)

#Perform DESeq2
dds <- DESeq(dds)
res <- results( dds )

#generate dds-resuls objects for relevant comparisons
mLN_SPF_GF <- results(dds, contrast=c("condition","ATAC_mLNSPF","ATAC_mLNGF"))
pLN_SPF_GF <- results(dds, contrast=c("condition","ATAC_pLNSPF","ATAC_pLNGF"))
SPF_mLN_pLN <- results(dds, contrast=c("condition","ATAC_mLNSPF","ATAC_pLNSPF"))
GF_mLN_pLN <- results(dds, contrast=c("condition","ATAC_mLNGF","ATAC_pLNGF"))

#QC
#MA
#par(mfrow = c(2,2))
#plotMA(mLN_SPF_GF, ylim = c(-6, 6), main = "mLN_SPF_GF")
#plotMA(pLN_SPF_GF, ylim = c(-6, 6), main = "pLN_SPF_GF" )
#plotMA(SPF_mLN_pLN, ylim = c(-6, 6), main = "SPF_mLN_pLN" )
#plotMA(GF_mLN_pLN, ylim = c(-6, 6), main = "GF_mLN_pLN" )

#DispEsts
#plotDispEsts( dds, ylim = c(1e-6, 1e1), main = "All" )


#Hits pValue
#par(mfrow = c(2,2))
#hist( mLN_SPF_GF$pvalue, breaks=20, col="grey", main = "mLN_SPF_GF"  )
#hist( pLN_SPF_GF$pvalue, breaks=20, col="grey", main = "pLN_SPF_GF"  )
#hist( SPF_mLN_pLN$pvalue, breaks=20, col="grey", main = "SPF_mLN_pLN"  )
#hist( GF_mLN_pLN$pvalue, breaks=20, col="grey", main = "GF_mLN_pLN"  )

#par(mfrow = c(1,1))
#plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

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
#Generate dataframe
DARs_features <- as.data.frame(gr_DARs_regions_anno_TSS)

#####
#Process RNA-seq data
#####
#Load expression data
mLN_SPF_GF <- read.csv(paste(path_RNAseq_DESeq2,"/","DESeq2_genes_diffexp_by_fc_mLN-GF_vs_mLN-SPF.csv",sep=""))
#pLNSPF/pLNGF
pLN_SPF_GF <- read.csv(paste(path_RNAseq_DESeq2,"/","DESeq2_genes_diffexp_by_fc_pLN-GF_vs_pLN-SPF.csv",sep=""))
#mLNSPF/pLNSPF
SPF_mLN_pLN <- read.csv(paste(path_RNAseq_DESeq2,"/","DESeq2_genes_diffexp_by_fc_pLN-SPF_vs_mLN-SPF.csv",sep=""))
#mLNGF/pLNGF
GF_mLN_pLN <- read.csv(paste(path_RNAseq_DESeq2,"/","DESeq2_genes_diffexp_by_fc_pLN-GF_vs_mLN-GF.csv",sep=""))

#prep tables for merging
#Perform for every DEseq2 DEG table
#1) Select genes that have an annotated GeneSymbol
#2) Take RPKM values per replica
#3) Rename columns
mLN_SPF_GF_NA <- mLN_SPF_GF[!(is.na(mLN_SPF_GF$GeneSymbol)),]
mLN_SPF_GF_NA_short <- subset(mLN_SPF_GF_NA, select = c("GeneSymbol", "log2FoldChange", "padj",
                                                        "RPKMcounts.4_mLN_GF", "RPKMcounts.5_mLN_GF",
                                                        "RPKMcounts.10_mLN_SPF", "RPKMcounts.12_mLN_SPF"))
colnames(mLN_SPF_GF_NA_short) <- c("GeneSymbol", "mLN_SPF_GF_log2FoldChange", "mLN_SPF_GF_padj",
                                   "mLN_SPF_GF_RPKM_mLN_GF_1", "mLN_SPF_GF_RPKM_mLN_GF_2",
                                   "mLN_SPF_GF_RPKM_mLN_SPF_1", "mLN_SPF_GF_RPKM_mLN_SPF_2")

pLN_SPF_GF_NA <- pLN_SPF_GF[!(is.na(pLN_SPF_GF$GeneSymbol)),]
pLN_SPF_GF_NA_short <- subset(pLN_SPF_GF_NA, select = c("GeneSymbol", "log2FoldChange", "padj",
                                                        "RPKMcounts.1_pLN_GF", "RPKMcounts.3a_pLN_GF",
                                                        "RPKMcounts.7_pLN_SPF", "RPKMcounts.8_pLN_SPF"))
colnames(pLN_SPF_GF_NA_short) <- c("GeneSymbol", "pLN_SPF_GF_log2FoldChange", "pLN_SPF_GF_padj",
                                   "pLN_SPF_GF_RPKM_pLN_GF_1", "pLN_SPF_GF_RPKM_pLN_GF_2",
                                   "pLN_SPF_GF_RPKM_pLN_SPF_1", "pLN_SPF_GF_RPKM_pLN_SPF_2")

SPF_mLN_pLN_NA <- SPF_mLN_pLN[!(is.na(SPF_mLN_pLN$GeneSymbol)),]
SPF_mLN_pLN_NA_short <- subset(SPF_mLN_pLN_NA, select = c("GeneSymbol", "log2FoldChange", "padj",
                                                          "RPKMcounts.7_pLN_SPF", "RPKMcounts.8_pLN_SPF",
                                                          "RPKMcounts.10_mLN_SPF", "RPKMcounts.12_mLN_SPF"))
colnames(SPF_mLN_pLN_NA_short) <- c("GeneSymbol", "SPF_mLN_pLN_log2FoldChange", "SPF_mLN_pLN_padj", 
                                    "SPF_mLN_pLN_RPKM_pLN_SPF_1", "SPF_mLN_pLN_RPKM_pLN_SPF_2",
                                    "SPF_mLN_pLN_RPKM_mLN_SPF_1", "SPF_mLN_pLN_RPKM_mLN_SPF_2")

GF_mLN_pLN_NA <- GF_mLN_pLN[!(is.na(GF_mLN_pLN$GeneSymbol)),]
GF_mLN_pLN_NA_short <- subset(GF_mLN_pLN_NA, select = c("GeneSymbol", "log2FoldChange", "padj", 
                                                        "RPKMcounts.1_pLN_GF", "RPKMcounts.3a_pLN_GF",
                                                        "RPKMcounts.4_mLN_GF", "RPKMcounts.5_mLN_GF"))
colnames(GF_mLN_pLN_NA_short) <- c("GeneSymbol", "GF_mLN_pLN_log2FoldChange", "GF_mLN_pLN_padj", 
                                   "GF_mLN_pLN_RPKM_pLN_GF_1", "GF_mLN_pLN_RPKM_pLN_GF_2",
                                   "GF_mLN_pLN_RPKM_mLN_GF_1", "GF_mLN_pLN_RPKM_mLN_GF_2")

#Merge all four tables, includes NAs
list_all4 = list(mLN_SPF_GF_NA_short, pLN_SPF_GF_NA_short, SPF_mLN_pLN_NA_short, GF_mLN_pLN_NA_short)
overall_expression = Reduce(function(...) merge(..., all=T), list_all4)


#Get columns of interest
expression <- overall_expression[,c("GeneSymbol",
                                    "SPF_mLN_pLN_log2FoldChange","GF_mLN_pLN_log2FoldChange",
                                    "mLN_SPF_GF_log2FoldChange","pLN_SPF_GF_log2FoldChange",
                                    "SPF_mLN_pLN_padj","GF_mLN_pLN_padj",
                                    "mLN_SPF_GF_padj","pLN_SPF_GF_padj")]
nrow(expression[!is.na(expression$GeneSymbol),])
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
names(RNAseq_UP) <- c("SPF_mLN_pLN","GF_mLN_pLN",
                      "mLN_SPF_GF","pLN_SPF_GF")
names(RNAseq_DOWN) <- c("SPF_mLN_pLN","GF_mLN_pLN",
                        "mLN_SPF_GF","pLN_SPF_GF")


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
names(t_Prom_Feat_TF_act) <- c("NoDEG_Open_SPF","NoDEG_Closed_SPF","NoDEG_Open_GF","NoDEG_Closed_GF")

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
DARs_features_TSS <- list()
for(i in seq(length(TF_DARs))){
  #load genes
  DARs_genes_i <- TF_DARs[[1]]
  #identify peaks for genes
  DARs_features_i <- subset(DARs_features_NA, symbol %in% DARs_genes_i)
  #Use only genes within Promotor regions
  DARs_features_TSS_i <- subset(DARs_features_i, (distancetoFeature <= downstream & distancetoFeature >= -upstream))
  print(names(TF_DARs)[i])
  print(nrow(DARs_features_TSS_i))
  #Generate Granges object
  gr_DARs_features_TSS_i <- toGRanges(DARs_features_TSS_i)
  #adjust width
  start(gr_DARs_features_TSS_i) <- start(gr_DARs_features_TSS_i) - 50
  end(gr_DARs_features_TSS_i) <- end(gr_DARs_features_TSS_i) + 50
  #Save generate .bed compatible for homer
  DARs_features_TSS_i <- as.data.frame(gr_DARs_features_TSS_i)
  DARs_features_TSS_bed_i <- DARs_features_TSS_i[,c("seqnames","start","end","id","feature_strand")] 
  #Make chr1 -> 1
  v_chr <- as.character(DARs_features_TSS_bed_i$seqnames)
  l_chr <- strsplit(v_chr, "r", fixed = FALSE, perl = FALSE, useBytes = FALSE)
  v_chr_integers <- unlist(lapply(l_chr, function(x) x[[2]]))
  DARs_features_TSS_bed_i$seqnames <- v_chr_integers
  #Generate empty column required by homer
  DARs_features_TSS_bed_i$blankVar <- NA
  DARs_features_TSS_bed_i <- DARs_features_TSS_bed_i[,c("seqnames","start","end","id","blankVar","feature_strand")] 
  DARs_features_TSS_bed_i$blankVar <- c("")
  #Store Feature BEDs 
  DARs_features_TSS[[i]] <- DARs_features_TSS_bed_i
  #Write .bed compatible for homer
  print(head(DARs_features_TSS_bed_i))
  write.table(DARs_features_TSS_bed_i, file=paste(path_output_Motif, "/", names(TF_DARs)[i],"_Ext50.bed", sep = ""),
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
names(TF_Activity) <- c("NoDEG_Open_SPF","NoDEG_Closed_SPF","NoDEG_Open_GF","NoDEG_Closed_GF")
#storing list for regions
PotAct_features_TSS <- list()
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
  #adjust width
  start(gr_PotAct_features_TSS_i) <- start(gr_PotAct_features_TSS_i) - 50
  end(gr_PotAct_features_TSS_i) <- end(gr_PotAct_features_TSS_i) + 50
  #Save generate .bed compatible for homer
  PotAct_features_TSS_i <- as.data.frame(gr_PotAct_features_TSS_i)
  PotAct_features_TSS_bed_i <- PotAct_features_TSS_i[,c("seqnames","start","end","id","feature_strand")] 
  #Make chr1 -> 1
  v_chr <- as.character(PotAct_features_TSS_bed_i$seqnames)
  l_chr <- strsplit(v_chr, "r", fixed = FALSE, perl = FALSE, useBytes = FALSE)
  v_chr_integers <- unlist(lapply(l_chr, function(x) x[[2]]))
  PotAct_features_TSS_bed_i$seqnames <- v_chr_integers
  #Generate empty column required by homer
  PotAct_features_TSS_bed_i$blankVar <- NA
  PotAct_features_TSS_bed_i <- PotAct_features_TSS_bed_i[,c("seqnames","start","end","id","blankVar","feature_strand")] 
  PotAct_features_TSS_bed_i$blankVar <- c("")
  #Store Feature BEDs 
  PotAct_features_TSS[[i]] <- PotAct_features_TSS_bed_i
  #Write .bed compatible for homer
  #print(head(PotAct_features_TSS_bed_i))
  write.table(PotAct_features_TSS_bed_i, file=paste(path_output_Motif, "/", names(TF_Activity)[i],"_Ext50.bed", sep = ""),
              quote=F, sep="\t", row.names=F, col.names=F)
}

#Average peak region size
hist(PotAct_features_TSS_bed_i$end - PotAct_features_TSS_bed_i$start, freq = NULL)

#####
#Get Background Peaks
#####
#Three types of background are tested
# 1) All genes names with open promotors
# 2) All DARs not used for motif analsis
# 3) All TSS regions of open promotors 
#Homer Motif only needs Gene identifier in the context of TSS/promotor dependent Motif identification
# 1) Background
background_gene_symbol <- background_promotors_open[,c("feature","symbol","seqnames","start","end")]
names(background_gene_symbol)[1] <- c("Acc")
rownames(background_gene_symbol) <- NULL
write.table(background_gene_symbol, paste(path_output_Motif, "/", "background_genes.txt", sep = ""),
            row.names = FALSE,quote=FALSE,sep="\t")

# 2) ####Obtain BED for all peak files for HOMER as background
#Obtain BED for all peak files for HOMER as background
# All Peaks that are not used in any enrichment analysis
PotAct_id <- unlist(lapply(PotAct_features_TSS, function(x){x$id}))
TFs_id <- unlist(lapply(DARs_features_TSS, function(x){x$id}))
background <- DARs_regions[!(DARs_regions$id %in% c(TFs_id,PotAct_id)),]
#Generate Granges object
gr_background <- toGRanges(background)
#adjust width
start(gr_background) <- start(gr_background) - 50
end(gr_background) <- end(gr_background) + 50
#Save generate .bed compatible for homer
background <- as.data.frame(gr_background)
background <- background[,c("seqnames","start","end","id","strand")]
background$blankVar <- NA
colnames(background) <- c("seqnames","start","end","id","feature_strand","blankVar")
background <- background[,c("seqnames","start","end","id","blankVar","feature_strand")] 
background$blankVar <- c("")
#annotate all strand as positive, information is not provided in the context of macs2 and ATAC-seq
background$feature_strand <- rep("+", nrow(background))
#Make chr1 -> 1
v_chr <- as.character(background$seqnames)
l_chr <- strsplit(v_chr, "r", fixed = FALSE, perl = FALSE, useBytes = FALSE)
v_chr_integers <- unlist(lapply(l_chr, function(x) x[[2]]))
background$seqnames <- v_chr_integers
#Export BED file for background
write.table(background, file=paste(path_output_Motif, "/background_all_peaks_minus_input_Ext50",".bed", sep = ""),
            quote=F, sep="\t", row.names=F, col.names=F)

# 3) All TSS regions of open promotors 
#DARs_features contains all peaks
#Only the ones annotated to TSS "carry" gene name in column Symbol 
#Grab TSS that are open
DARs_features_TSS <- subset(DARs_features, (!is.na(DARs_features$symbol)))
background_genes <- unique(DARs_features_TSS[!(DARs_features_TSS$id %in% c(TFs_id,PotAct_id)),]$symbol)
#add gene name
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
listFilters(mouse)
listAttributes(mouse)
ENSEMBLE_to_GeneSymbol <- getBM(attributes=c("ensembl_gene_id", "mgi_symbol"),
                                filters="mgi_symbol",
                                values=background_genes,
                                mart=mouse)
#Grab only background genes from TSS file
#Get Genomic locations of TSS
d_TSS.mouse.mm10 <- as.data.frame(TSS.mouse.mm10)
background_TSS <- subset(d_TSS.mouse.mm10, rownames(d_TSS.mouse.mm10) %in% ENSEMBLE_to_GeneSymbol$ensembl_gene_id)
#Make export bed file
# id is not relevant for homer in the context of background files
background_TSS <- background_TSS[,c("seqnames","start","end","width","strand")]
background_TSS$blankVar <- NA
colnames(background_TSS) <- c("seqnames","start","end","id","feature_strand","blankVar")
background_TSS <- background_TSS[,c("seqnames","start","end","id","blankVar","feature_strand")] 
background_TSS$blankVar <- c("")
#Eliminate all intel with CHR_..._PATCH
background_TSS <- background_TSS[- grep("CHR_", background_TSS$seqnames),]
#chr to integers in seqname
#background_TSS$seqnames <- paste("chr",background_TSS$seqnames,sep="")
write.table(background_TSS, file=paste(path_output_Motif, "/background_all_TSS_minus_input_genes",".bed", sep = ""),
            quote=F, sep="\t", row.names=F, col.names=F)

