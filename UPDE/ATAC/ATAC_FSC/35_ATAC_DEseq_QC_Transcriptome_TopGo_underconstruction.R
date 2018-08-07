# Author: Vincent Gardeux
# Adapted by: Joern Pezoldt
# 12.07.2018
# Function:
# 1) 
# 

#####
#Libraries & PATHS & Data
#####

require(data.table)
require(ggplot2)
require(ggfortify)
require(limma)
require(DESeq2)

# Input required: Set directory
path_input <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/homer/Overlap_Group_Merged/Run_2_Across_All_Replicates"
path_output <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/DESeq2/Run_2_Across_All_Replicates"
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

# !Unclear functionality
#percent_50 <- number_of_samples / 2
#to.filter <- which(rowSums(counts==0)<percent_50) ## 34176
#setkey(counts, id)
#counts.2 <- counts[-to.filter]

# !Unclear functionality
#setkey(counts, id)
#rownames(counts) <- counts[counts$id]$id

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
# batch effect ------------------------------------------------------------
#####
#batch <- fread("BRBseq/brbseq_analysis/complete_batches.txt")
#setkey(batch, name.timepoint)
#batch <- batch[colnames(counts.2)[-1]]
#batch <- na.omit(batch)
#nrow(batch)

#counts.2 <- counts.2[,colnames(counts.2) %in% batch$name.timepoint, with=F]
#var <- list(batch$main.batch,
 #           batch$gender.female,
  #          batch$PC1,
   #         batch$PC2, 
    #        batch$PC3)
#var.names <- c("main batch", "sex", "Genotype PC1","Genotype PC2","Genotype PC3")

#plot_anova_PC(counts.2, var = var, var.names = var.names)
#batch_PC(counts.2, batch$main.batch)

#####
# DEseq2 -------------------------------------------------------------------
#####

#Split table into groups of two
#mLN_SPF vs pLN_SPF
# Input required: Colonization status (SPF, GF) or lymph node type (pLN, mLN)
counts.mat_two_pair <- counts.mat[,grepl(c("SPF"), colnames(counts.mat))]
# Input required: Name of comparison
group_name_1 <- "mLNSPF"
group_name_2 <- "pLNSPF"

#Continue here 
condition <- factor(substr(colnames(counts.mat_two_pair),1,nchar(colnames(counts.mat_two_pair))-1))

#condition <- factor(rep(c("t0", "t14"), ncol))
dds <- DESeqDataSetFromMatrix(counts.mat_two_pair, DataFrame(condition),~condition)
dds.factors <- estimateSizeFactors(dds)
count.norm <- log2(1+counts(dds.factors, normalized=T))
count.norm <- data.frame(count.norm)
count.norm$id <- rownames(count.norm)

#write.table(count.norm, paste(path_output, "/peak_count_norm.txt", sep=""), quote = F, sep="\t", row.names=T)


#mLN_SPF vs mLN_GF
# Input required: Colonization status (SPF, GF) or lymph node type (pLN, mLN)
counts.mat_two_pair <- counts.mat[,grepl(c("mLN"), colnames(counts.mat))]
# Input required: Name of comparison
group_name_1 <- "mLNSPF"
group_name_2 <- "mLNGF"

#Continue here 
condition <- factor(substr(colnames(counts.mat_two_pair),1,nchar(colnames(counts.mat_two_pair))-1))

#condition <- factor(rep(c("t0", "t14"), ncol))
dds <- DESeqDataSetFromMatrix(counts.mat_two_pair, DataFrame(condition),~condition)
dds.factors <- estimateSizeFactors(dds)
count.norm <- log2(1+counts(dds.factors, normalized=T))
count.norm <- data.frame(count.norm)
count.norm$id <- rownames(count.norm)

#write.table(count.norm, paste(path_output, "/peak_count_norm.txt", sep=""), quote = F, sep="\t", row.names=T)

#mLN_GF vs pLN_GF
# Input required: Colonization status (SPF, GF) or lymph node type (pLN, mLN)
counts.mat_two_pair <- counts.mat[,grepl(c("GF"), colnames(counts.mat))]
# Input required: Name of comparison
group_name_1 <- "mLNGF"
group_name_2 <- "pLNGF"

#Continue here 
condition <- factor(substr(colnames(counts.mat_two_pair),1,nchar(colnames(counts.mat_two_pair))-1))

#condition <- factor(rep(c("t0", "t14"), ncol))
dds <- DESeqDataSetFromMatrix(counts.mat_two_pair, DataFrame(condition),~condition)
dds.factors <- estimateSizeFactors(dds)
count.norm <- log2(1+counts(dds.factors, normalized=T))
count.norm <- data.frame(count.norm)
count.norm$id <- rownames(count.norm)

#write.table(count.norm, paste(path_output, "/peak_count_norm.txt", sep=""), quote = F, sep="\t", row.names=T)

#pLN_SPF vs pLN_GF
# Input required: Colonization status (SPF, GF) or lymph node type (pLN, mLN)
counts.mat_two_pair <- counts.mat[,grepl(c("pLN"), colnames(counts.mat))]
# Input required: Name of comparison
group_name_1 <- "pLNSPF"
group_name_2 <- "pLNGF"

#Continue here 
condition <- factor(substr(colnames(counts.mat_two_pair),1,nchar(colnames(counts.mat_two_pair))-1))

#condition <- factor(rep(c("t0", "t14"), ncol))
dds <- DESeqDataSetFromMatrix(counts.mat_two_pair, DataFrame(condition),~condition)
dds.factors <- estimateSizeFactors(dds)
count.norm <- log2(1+counts(dds.factors, normalized=T))
count.norm <- data.frame(count.norm)
count.norm$id <- rownames(count.norm)

#write.table(count.norm, paste(path_output, "/peak_count_norm.txt", sep=""), quote = F, sep="\t", row.names=T)

#Perform DESeq2
dds <- DESeq(dds)
res <- results( dds )

#QC
plotMA( res, ylim = c(-6, 6) )
plotDispEsts( dds, ylim = c(1e-6, 1e1) )
hist( res$pvalue, breaks=20, col="grey" )

#####
#Peak/gene location
#####
library("GenomicFeatures")



#####
#Associate genes with DESeq results
#####
res_table <- as.data.frame(res)
res_table$id <- rownames(res_table)

#Load the genomic regions associated with ids
# Common peak track for Experiment is annotated using homer

# Input required: 
path_regions <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/peaks/broad/ATAC_FSC_all_broad_merged_peaks.bed"
regions <- read.delim(path_regions, header = FALSE)
colnames(regions) <- c("chr","start","end")
res_regions <- cbind(res, regions)



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
