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



###########
#Bsmooth
###########
#load
mLN_SPF_GF <- read.delim2("/home/pezoldt/NAS2/pezoldt/Data/WGBS/WGBS_FSC/BSmooth_DMRs_mLN_SPF_vs_mLN_GF.txt")
pLN_SPF_GF <- read.delim2("/home/pezoldt/NAS2/pezoldt/Data/WGBS/WGBS_FSC/BSmooth_DMRs_pLN_SPF_vs_pLN_GF.txt")
SPF_mLN_pLN <- read.delim2("/home/pezoldt/NAS2/pezoldt/Data/WGBS/WGBS_FSC/BSmooth_DMRs_mLN_SPF_vs_pLN_SPF.txt")
GF_mLN_pLN <- read.delim2("/home/pezoldt/NAS2/pezoldt/Data/WGBS/WGBS_FSC/BSmooth_DMRs_pLN_GF_vs_mLN_GF.txt")
#Crossmapped
SPF_mLN_pLN_Crossmapped <- read.delim2("/home/pezoldt/NAS2/pezoldt/Data/WGBS/WGBS_FSC/BSmooth_DMRs_mLN_SPF_vs_pLN_SPF_CrossMap.txt")
#CpGs per DMR
CpGs_per_DMR <- read.csv("/home/pezoldt/NAS2/pezoldt/Data/WGBS/WGBS_FSC/CpGs_in_DMRs.csv")

#rename group.1 group.2
colnames(mLN_SPF_GF) <- c("X","chr","start","end","n","width","areaStat","meanDiff","meth_mLN_GF",       
                        "meth_mLN_SPF","direction","nearest_tss","distance_to_nearest_tss","nearest_gene_id",
                        "nearest_gene_name","nearest_gene_desc","preceding_tss","following_tss")
colnames(pLN_SPF_GF) <- c("X","chr","start","end","n","width","areaStat","meanDiff","meth_pLN_GF",       
                          "meth_pLN_SPF","direction","nearest_tss","distance_to_nearest_tss","nearest_gene_id",
                          "nearest_gene_name","nearest_gene_desc","preceding_tss","following_tss")
colnames(SPF_mLN_pLN) <- c("X","chr","start","end","n","width","areaStat","meanDiff","meth_pLN_SPF",       
                          "meth_mLN_SPF","direction","nearest_tss","distance_to_nearest_tss","nearest_gene_id",
                          "nearest_gene_name","nearest_gene_desc","preceding_tss","following_tss")
colnames(GF_mLN_pLN) <- c("X","chr","start","end","n","width","areaStat","meanDiff","meth_mLN_GF",       
                          "meth_pLN_GF","direction","nearest_tss","distance_to_nearest_tss","nearest_gene_id",
                          "nearest_gene_name","nearest_gene_desc","preceding_tss","following_tss")
colnames(SPF_mLN_pLN_Crossmapped) <- c("X","chr","start","end","n","width","areaStat","meanDiff","meth_pLN_SPF",       
                           "meth_mLN_SPF","direction","nearest_tss","distance_to_nearest_tss","nearest_gene_id",
                           "nearest_gene_name","nearest_gene_desc","preceding_tss","following_tss")

#add columns indicating the comparison results
mLN_SPF_GF[,"Bsmooth_comp"] <- "mLN_SPF_GF"
pLN_SPF_GF[,"Bsmooth_comp"] <- "pLN_SPF_GF"
SPF_mLN_pLN[,"Bsmooth_comp"] <- "SPF_mLN_pLN"
GF_mLN_pLN[,"Bsmooth_comp"] <- "GF_mLN_pLN"
Bsmooth_list <- list(mLN_SPF_GF,pLN_SPF_GF, SPF_mLN_pLN, GF_mLN_pLN)
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
gr_mLN_SPF_GF <- with(mLN_SPF_GF, GRanges(chr, IRanges(start = start, end = end), 
                                          Gene_Symbol = nearest_gene_name, areastat = areaStat, width_dmr = width, mean_Diff = meanDiff,
                                          mLN_GF = meth_mLN_GF, mLN_SPF = meth_mLN_SPF,
                                          distance_TSS = distance_to_nearest_tss))
gr_pLN_SPF_GF <- with(pLN_SPF_GF, GRanges(chr, IRanges(start = start, end = end), 
                                          Gene_Symbol = nearest_gene_name, areastat = areaStat, width_dmr = width, mean_Diff = meanDiff,
                                          pLN_GF = meth_pLN_GF, pLN_SPF = meth_pLN_SPF,
                                          distance_TSS = distance_to_nearest_tss))
gr_SPF_mLN_pLN <- with(SPF_mLN_pLN, GRanges(chr, IRanges(start = start, end = end), 
                                            Gene_Symbol = nearest_gene_name, areastat = areaStat, width_dmr = width, mean_Diff = meanDiff,
                                            pLN_SPF = meth_pLN_SPF, mLN_SPF = meth_mLN_SPF,
                                            distance_TSS = distance_to_nearest_tss))
gr_GF_mLN_pLN <- with(GF_mLN_pLN, GRanges(chr, IRanges(start = start, end = end), 
                                          Gene_Symbol = nearest_gene_name, areastat = areaStat, width_dmr = width, mean_Diff = meanDiff,
                                          mLN_GF = meth_mLN_GF, pLN_GF = meth_pLN_GF,
                                          distance_TSS = distance_to_nearest_tss))
###Identify overlap location
hits_location <- findOverlaps(gr_SPF_mLN_pLN, gr_GF_mLN_pLN, minoverlap = 50)
#Get index
query_Hits <- queryHits(hits_location)
subject_Hits <- subjectHits(hits_location)
#Get metadata
overlap_location_SPF <- gr_SPF_mLN_pLN[query_Hits]
overlap_location_GF <- gr_GF_mLN_pLN[subject_Hits]
#nonoverlap
#index vector for DMRs not overlapping with location SPF
index_SPF_mLN_pLN <- c(1:nrow(SPF_mLN_pLN))
remove_SPF_mLN_pLN <- query_Hits
non_overlap_SPF <- index_SPF_mLN_pLN[! index_SPF_mLN_pLN %in% remove_SPF_mLN_pLN]
#index vector for DMRs not overlapping with location GF
index_GF_mLN_pLN <- c(1:nrow(GF_mLN_pLN))
remove_GF_mLN_pLN <- query_Hits
non_overlap_GF <- index_GF_mLN_pLN[! index_GF_mLN_pLN %in% remove_GF_mLN_pLN]
#Granges object for non_overlap
nonoverlap_location_SPF <- gr_SPF_mLN_pLN[non_overlap_SPF]
nonoverlap_location_GF <- gr_GF_mLN_pLN[non_overlap_GF]
#make data fram
overlap_location_SPF_mLN_pLN <- as.data.frame(overlap_location_SPF)
nonoverlap_location_SPF_mLN_pLN <- as.data.frame(nonoverlap_location_SPF)
nonoverlap_location_GF_mLN_pLN <- as.data.frame(nonoverlap_location_GF)

###Identify overlap commensal
#hits_commensal <- findOverlaps(gr_mLN_SPF_GF, gr_pLN_SPF_GF, minoverlap = 50)
#Note: no overlap
#make data fram
nonoverlap_commensal_mLN_SPF_GF <- as.data.frame(gr_mLN_SPF_GF)
nonoverlap_commensal_pLN_SPF_GF <- as.data.frame(gr_pLN_SPF_GF)

#merge column chr stop end of each dataframe to generate identifier for DMR CpG average heatmap
o_seq <- with(overlap_location_SPF_mLN_pLN, paste(seqnames, "_", start, "_", end, sep = ""))
non_SPF <- with(nonoverlap_location_SPF_mLN_pLN, paste(seqnames, "_", start, "_", end, sep = ""))
non_GF <- with(nonoverlap_location_GF_mLN_pLN, paste(seqnames, "_", start, "_", end, sep = ""))
non_mLN <- with(nonoverlap_commensal_mLN_SPF_GF, paste(seqnames, "_", start, "_", end, sep = ""))
non_pLN <- with(nonoverlap_commensal_pLN_SPF_GF, paste(seqnames, "_", start, "_", end, sep = ""))
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
dMeth_mLN_SPF_GF <- df_DMR_mean_Bsmooth_CpGmeaned_NA$mLN_SPF - df_DMR_mean_Bsmooth_CpGmeaned_NA$mLN_GF
dMeth_pLN_SPF_GF <- df_DMR_mean_Bsmooth_CpGmeaned_NA$pLN_SPF - df_DMR_mean_Bsmooth_CpGmeaned_NA$pLN_GF
dMeth_SPF_mLN_pLN <- df_DMR_mean_Bsmooth_CpGmeaned_NA$mLN_SPF - df_DMR_mean_Bsmooth_CpGmeaned_NA$pLN_SPF
dMeth_GF_mLN_pLN <- df_DMR_mean_Bsmooth_CpGmeaned_NA$mLN_GF - df_DMR_mean_Bsmooth_CpGmeaned_NA$pLN_GF
df_DMR_mean_Bsmooth_CpGmeaned_dMeth <- cbind(df_DMR_mean_Bsmooth_CpGmeaned_NA, dMeth_mLN_SPF_GF)
df_DMR_mean_Bsmooth_CpGmeaned_dMeth <- cbind(df_DMR_mean_Bsmooth_CpGmeaned_dMeth, dMeth_pLN_SPF_GF)
df_DMR_mean_Bsmooth_CpGmeaned_dMeth <- cbind(df_DMR_mean_Bsmooth_CpGmeaned_dMeth, dMeth_SPF_mLN_pLN)
df_DMR_mean_Bsmooth_CpGmeaned_dMeth <- cbind(df_DMR_mean_Bsmooth_CpGmeaned_dMeth, dMeth_GF_mLN_pLN)
df_DMR_mean_Bsmooth_CpGmeaned_dMeth_TSS <- subset(df_DMR_mean_Bsmooth_CpGmeaned_dMeth, distance_to_nearest_tss < 10000)
subset(df_DMR_mean_Bsmooth_CpGmeaned_dMeth_TSS, nearest_gene_name %in% c("Aldh1a2"))

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


################
#Generate Output for home Motif analysis
################
#SPF_mLN_pLN_Crossmapped
#SPF_mLN_pLN_Cross_TSS <- subset(SPF_mLN_pLN_Crossmapped, distance_to_nearest_tss < 10000)
#Classify DMRs according to location
#####
#Peak/gene location
#####
library("GenomicFeatures")
# @Joern: Continue here
#All features
mm10_TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#TSS sites BioMart
mm10 = useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl") 
TSS.mouse.mm10 = getAnnotation(mart=mm10, featureType="TSS")

#Annotate Peaks SPF-------------------------------
gr_SPF_mLN_pLN_anno_TSS <- annotatePeakInBatch(gr_SPF_mLN_pLN, AnnotationData=TSS.mouse.mm10)

#add gene name
gr_SPF_mLN_pLN_anno_TSS <- addGeneIDs(annotatedPeak=gr_SPF_mLN_pLN_anno_TSS, 
                                       feature_id_type="ensembl_gene_id",
                                       orgAnn="org.Mm.eg.db", 
                                       IDs2Add="symbol")
#Generate dataframe
SPF_mLN_pLN_anno_TSS <- as.data.frame(gr_SPF_mLN_pLN_anno_TSS)

#Annotate Peaks GF-------------------------------
gr_GF_mLN_pLN_anno_TSS <- annotatePeakInBatch(gr_GF_mLN_pLN, AnnotationData=TSS.mouse.mm10)

#add gene name
gr_GF_mLN_pLN_anno_TSS <- addGeneIDs(annotatedPeak=gr_GF_mLN_pLN_anno_TSS, 
                                      feature_id_type="ensembl_gene_id",
                                      orgAnn="org.Mm.eg.db", 
                                      IDs2Add="symbol")
#Generate dataframe
GF_mLN_pLN_anno_TSS <- as.data.frame(gr_GF_mLN_pLN_anno_TSS)

#####
#Export for Motif analysis
#####
### Subset Tables
# SPF
SPF_mLN_pLN_TSS <- SPF_mLN_pLN_anno_TSS[!(is.na(SPF_mLN_pLN_anno_TSS$symbol)),]
SPF_mLN_pLN_TSS_hypo_thresh <- subset(SPF_mLN_pLN_TSS, mean_Diff <= -0.3)
SPF_mLN_pLN_TSS_hypo_all <- subset(SPF_mLN_pLN_TSS, mean_Diff < 0)
SPF_mLN_pLN_TSS_hyper_thresh <- subset(SPF_mLN_pLN_TSS, mean_Diff >= 0.3)
SPF_mLN_pLN_TSS_hyper_all <- subset(SPF_mLN_pLN_TSS, mean_Diff > 0)
# GF
GF_mLN_pLN_TSS <- GF_mLN_pLN_anno_TSS[!(is.na(GF_mLN_pLN_anno_TSS$symbol)),]
GF_mLN_pLN_TSS_hypo_thresh <- subset(GF_mLN_pLN_TSS, mean_Diff <= -0.3)
GF_mLN_pLN_TSS_hypo_all <- subset(GF_mLN_pLN_TSS, mean_Diff < 0)
GF_mLN_pLN_TSS_hyper_thresh <- subset(GF_mLN_pLN_TSS, mean_Diff >= 0.3)
GF_mLN_pLN_TSS_hyper_all <- subset(GF_mLN_pLN_TSS, mean_Diff > 0)


#List of Tables
l_DMRs_to_homer <- list(SPF_mLN_pLN_TSS_hypo_thresh,SPF_mLN_pLN_TSS_hypo_all,
                        SPF_mLN_pLN_TSS_hyper_thresh,SPF_mLN_pLN_TSS_hyper_all,
                        GF_mLN_pLN_TSS_hypo_thresh,GF_mLN_pLN_TSS_hypo_all,
                        GF_mLN_pLN_TSS_hyper_thresh,GF_mLN_pLN_TSS_hyper_all)
names(l_DMRs_to_homer) <- c("SPF_mLN_pLN_TSS_hypo_thresh","SPF_mLN_pLN_TSS_hypo_all",
                            "SPF_mLN_pLN_TSS_hyper_thresh","SPF_mLN_pLN_TSS_hyper_all",
                            "GF_mLN_pLN_TSS_hypo_thresh","GF_mLN_pLN_TSS_hypo_all",
                            "GF_mLN_pLN_TSS_hyper_thresh","GF_mLN_pLN_TSS_hyper_all")

#Empty list for storage
l_DMRs_to_homer_fini <- list()
for(i in seq(length(l_DMRs_to_homer))){
  t_DMR_i <- l_DMRs_to_homer[[i]]
  
  #build table input
  v_chr <- as.character(t_DMR_i$seqnames)
  l_chr <- strsplit(v_chr, "r", fixed = FALSE, perl = FALSE, useBytes = FALSE)
  seqnames_i <- unlist(lapply(l_chr, function(x) x[[2]]))
  #  seqnames_i <- as.character(t_DMR_i$seqnames)
  start_i = t_DMR_i$start
  end_i = t_DMR_i$end
  id_i = t_DMR_i$peak
  blankVar_i = rep("",length(id_i))
  feature_strand_i = rep("+",length(id_i))
  t_DMR_adapted_i <- data.frame(seqnames = seqnames_i, start = start_i,
                              end = end_i, id = id_i, blankVar = blankVar_i,
                              feature_strand = feature_strand_i)

  #Store tables for motif identification
  write.table(t_DMR_adapted_i, file=paste(output_beds_for_motifs, "/", names(l_DMRs_to_homer)[i],"_CHR.bed", sep = ""),
              quote=F, sep="\t", row.names=F, col.names=F)
  l_DMRs_to_homer_fini[[i]] <- t_DMR_adapted_i
}

#####
#Obtain DMRs evolving around the single CpGs
#####
#Variables
# delta meth required for CpG to be included and locied
#CpG_state_thresh <- 0.4
#Function:
# Merges overlapping entries in Granges object
# Extract clusters from Hits object.
#extractClustersFromSelfHits <- function(hits){
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
#mergeConnectedRanges <- function(x, hits){
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
#l_DMRs_CpG_regions <- list(SPF_mLN_pLN_TSS_hypo_thresh,SPF_mLN_pLN_TSS_hypo_all,
 #                          SPF_mLN_pLN_TSS_hyper_thresh,SPF_mLN_pLN_TSS_hyper_all,
  #                         GF_mLN_pLN_TSS_hypo_thresh,GF_mLN_pLN_TSS_hypo_all,
   #                        GF_mLN_pLN_TSS_hyper_thresh,GF_mLN_pLN_TSS_hyper_all)

#names(l_DMRs_CpG_regions) <- list("SPF_mLN_pLN_TSS_hypo_thresh","SPF_mLN_pLN_TSS_hypo_all",
 #                                 "SPF_mLN_pLN_TSS_hyper_thresh","SPF_mLN_pLN_TSS_hyper_all",
  #                                "GF_mLN_pLN_TSS_hypo_thresh","GF_mLN_pLN_TSS_hypo_all",
   #                               "GF_mLN_pLN_TSS_hyper_thresh","GF_mLN_pLN_TSS_hyper_all")

#initialize lists that will contain the list of the CpGs regions per DMR
#l_26bp <- list()
#l_50bp <- list()
#l_100bp <- list()
#l_200bp <- list()

#Note: Output of the following loop is a nested list with the content of the input DMR tables as a list of CpG regions per DMR 
#       which is then stored in list annotating the size of the regions
#for(i in 1:length(l_DMRs_CpG_regions)){
  # Grab DMRs threshed in DMRs lists
  print(paste("i == ", i, sep = ""))
  DMRs_i <- l_DMRs_CpG_regions[[i]]
  #Obtain common identifier
  Identifier_CpG_table_i <- paste(DMRs_i$seqnames, "_", DMRs_i$start,"_",DMRs_i$end,sep="")
  #Check in CpG list for lists of respective DMR
  #Take only "mLN_SPF_vs_pLN_SPF"
  CpGs_SPF_mLN_pLN_i <- subset(CpGs_per_DMR, DMR_comp_origin == "mLN_SPF_vs_pLN_SPF")
  split_CpGs_per_DMR_i <- split(CpGs_SPF_mLN_pLN_i, CpGs_SPF_mLN_pLN_i$DMR_identifier)
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
    dmeth_SPF_mLN_pLN <- meth_mLN_SPF - meth_pLN_SPF
    dmeth_GF_mLN_pLN <- meth_mLN_GF - meth_pLN_GF
    
    #append delta meth columns
    t_CpGs_per_DMR_j <- cbind(t_CpGs_per_DMR_j, dmeth_SPF_mLN_pLN, dmeth_GF_mLN_pLN)
    
    #select CpGs within DMR based on delta meth
    t_CpGs_per_DMR_thresh_j <- subset(t_CpGs_per_DMR_j, abs(dmeth_SPF_mLN_pLN) > CpG_state_thresh)
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
#l_all_regions <- list(l_26bp,l_50bp,l_100bp,l_200bp)
#names(l_all_regions) <- c("26bp","50bp","100bp","200bp")

#for(i in 1:length(l_all_regions)){
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



#################
#Transcriptome Conventional
#################
#load
T_mLN_SPF_GF <- read.csv("/home/pezoldt/NAS2/pezoldt/Data/RNAseq/2014_FSC_mLN_pLN_GF_SPF/DESeq2_genes_diffexp_by_fc_mLN-GF_vs_mLN-SPF.csv")
T_GF_mLN_pLN <- read.csv("/home/pezoldt/NAS2/pezoldt/Data/RNAseq/2014_FSC_mLN_pLN_GF_SPF/DESeq2_genes_diffexp_by_fc_pLN-GF_vs_mLN-GF.csv")
T_SPF_mLN_pLN <- read.csv("/home/pezoldt/NAS2/pezoldt/Data/RNAseq/2014_FSC_mLN_pLN_GF_SPF/DESeq2_genes_diffexp_by_fc_pLN-SPF_vs_mLN-SPF.csv")
T_pLN_SPF_GF <- read.csv("/home/pezoldt/NAS2/pezoldt/Data/RNAseq/2014_FSC_mLN_pLN_GF_SPF/DESeq2_genes_diffexp_by_fc_pLN-GF_vs_pLN-SPF.csv")

mLN_SPF_GF_NA <- T_mLN_SPF_GF[!(is.na(T_mLN_SPF_GF$GeneSymbol)),]
mLN_SPF_GF_NA_short <- subset(mLN_SPF_GF_NA, select = c("GeneSymbol", "log2FoldChange", "padj",
                                                        "RPKMcounts.4_mLN_GF_R1_clip", "RPKMcounts.5_mLN_GF_R1_clip",
                                                        "RPKMcounts.10_mLN_SPF_R1_clip", "RPKMcounts.12_mLN_SPF_R1_clip",
                                                        "meanRPKM.mLN.GF", "meanRPKM.mLN.SPF"))
colnames(mLN_SPF_GF_NA_short) <- c("GeneSymbol", "mLN_SPF_GF_log2FoldChange", "mLN_SPF_GF_padj",
                                   "mLN_SPF_GF_RPKM_mLN_GF_1", "mLN_SPF_GF_RPKM_mLN_GF_2",
                                   "mLN_SPF_GF_RPKM_mLN_SPF_1", "mLN_SPF_GF_RPKM_mLN_SPF_2",
                                   "meanRPKM_mLN_GF", "meanRPKM_mLN_SPF")

pLN_SPF_GF_NA <- T_pLN_SPF_GF[!(is.na(T_pLN_SPF_GF$GeneSymbol)),]
pLN_SPF_GF_NA_short <- subset(pLN_SPF_GF_NA, select = c("GeneSymbol", "log2FoldChange", "padj",
                                                        "RPKMcounts.1_pLN_GF_R1_clip", "RPKMcounts.3a_pLN_GF_R1_clip",
                                                        "RPKMcounts.7_pLN_SPF_R1_clip", "RPKMcounts.8_pLN_SPF_R1_clip",
                                                        "meanRPKM.pLN.GF", "meanRPKM.pLN.SPF"))
colnames(pLN_SPF_GF_NA_short) <- c("GeneSymbol", "pLN_SPF_GF_log2FoldChange", "pLN_SPF_GF_padj",
                                   "pLN_SPF_GF_RPKM_pLN_GF_1", "pLN_SPF_GF_RPKM_pLN_GF_2",
                                   "pLN_SPF_GF_RPKM_pLN_SPF_1", "pLN_SPF_GF_RPKM_pLN_SPF_2",
                                   "meanRPKM_pLN_GF", "meanRPKM_pLN_SPF")

SPF_mLN_pLN_NA <- T_SPF_mLN_pLN[!(is.na(T_SPF_mLN_pLN$GeneSymbol)),]
SPF_mLN_pLN_NA_short <- subset(SPF_mLN_pLN_NA, select = c("GeneSymbol", "log2FoldChange", "padj",
                                                          "RPKMcounts.7_pLN_SPF_R1_clip", "RPKMcounts.8_pLN_SPF_R1_clip",
                                                          "RPKMcounts.10_mLN_SPF_R1_clip", "RPKMcounts.12_mLN_SPF_R1_clip",
                                                          "meanRPKM.pLN.SPF","meanRPKM.mLN.SPF"))
colnames(SPF_mLN_pLN_NA_short) <- c("GeneSymbol", "SPF_mLN_pLN_log2FoldChange", "SPF_mLN_pLN_padj", 
                                    "SPF_mLN_pLN_RPKM_pLN_SPF_1", "SPF_mLN_pLN_RPKM_pLN_SPF_2",
                                    "SPF_mLN_pLN_RPKM_mLN_SPF_1", "SPF_mLN_pLN_RPKM_mLN_SPF_2",
                                    "meanRPKM_pLN_SPF","meanRPKM_mLN_SPF")

GF_mLN_pLN_NA <- T_GF_mLN_pLN[!(is.na(T_GF_mLN_pLN$GeneSymbol)),]
GF_mLN_pLN_NA_short <- subset(GF_mLN_pLN_NA, select = c("GeneSymbol", "log2FoldChange", "padj", 
                                                        "RPKMcounts.1_pLN_GF_R1_clip", "RPKMcounts.3a_pLN_GF_R1_clip",
                                                        "RPKMcounts.4_mLN_GF_R1_clip", "RPKMcounts.5_mLN_GF_R1_clip",
                                                        "meanRPKM.pLN.GF","meanRPKM.mLN.GF"))
colnames(GF_mLN_pLN_NA_short) <- c("GeneSymbol", "GF_mLN_pLN_log2FoldChange", "GF_mLN_pLN_padj", 
                                   "GF_mLN_pLN_RPKM_pLN_GF_1", "GF_mLN_pLN_RPKM_pLN_GF_2",
                                   "GF_mLN_pLN_RPKM_mLN_GF_1", "GF_mLN_pLN_RPKM_mLN_GF_2",
                                   "meanRPKM_pLN_GF","meanRPKM_mLN_GF")

#Merge total data, includes NAs
list_all4 = list(mLN_SPF_GF_NA_short, pLN_SPF_GF_NA_short, SPF_mLN_pLN_NA_short, GF_mLN_pLN_NA_short)
overall_expression = Reduce(function(...) merge(..., all=T), list_all4)

#####
#Thresh expression data
#####
overall_expression_2FC_padj <- subset(overall_expression, 
                                      (mLN_SPF_GF_log2FoldChange > 1.0 & mLN_SPF_GF_padj < 0.05) | 
                                        (mLN_SPF_GF_log2FoldChange < -1.0 & mLN_SPF_GF_padj < 0.05) |                                 
                                        (pLN_SPF_GF_log2FoldChange > 1.0 & pLN_SPF_GF_padj < 0.05) | 
                                        (pLN_SPF_GF_log2FoldChange < -1.0 & pLN_SPF_GF_padj < 0.05) |
                                        (SPF_mLN_pLN_log2FoldChange > 1.0 & SPF_mLN_pLN_padj < 0.05) | 
                                        (SPF_mLN_pLN_log2FoldChange < -1.0 & SPF_mLN_pLN_padj < 0.05) |
                                        (GF_mLN_pLN_log2FoldChange > 1.0 & GF_mLN_pLN_padj < 0.05) | 
                                        (GF_mLN_pLN_log2FoldChange < -1.0 & GF_mLN_pLN_padj < 0.05))


overall_expression_2FC_padj_RPKM <- subset(overall_expression_2FC_padj,
                                           SPF_mLN_pLN_RPKM_pLN_SPF_1 > 1.0 | SPF_mLN_pLN_RPKM_pLN_SPF_2 > 1.0 |
                                             SPF_mLN_pLN_RPKM_mLN_SPF_1 > 1.0 | SPF_mLN_pLN_RPKM_mLN_SPF_2 > 1.0 |
                                             GF_mLN_pLN_RPKM_pLN_GF_1 > 1.0 | GF_mLN_pLN_RPKM_pLN_GF_2 > 1.0 |
                                             GF_mLN_pLN_RPKM_mLN_GF_1 > 1.0 | GF_mLN_pLN_RPKM_mLN_GF_2 > 1.0)

overall_expression_2FC_padj_RPKM <- overall_expression_2FC_padj_RPKM[!(is.na(overall_expression_2FC_padj_RPKM$mLN_SPF_GF_log2FoldChange)),]

#Keep Gene name and expression data of differentially expressed genes
diff_Genes <- overall_expression_2FC_padj_RPKM[,c("GeneSymbol",
                                                  "meanRPKM_mLN_GF","meanRPKM_mLN_SPF",
                                                  "meanRPKM_pLN_GF","meanRPKM_pLN_SPF",
                                                  "mLN_SPF_GF_log2FoldChange","mLN_SPF_GF_padj",
                                                  "pLN_SPF_GF_log2FoldChange","pLN_SPF_GF_padj",
                                                  "SPF_mLN_pLN_log2FoldChange","SPF_mLN_pLN_padj",
                                                  "GF_mLN_pLN_log2FoldChange","GF_mLN_pLN_padj")]

#####
#Transcriptome LowInput
#####
#FRC LI 6000 cells
# Get DEGs
#FRC_LI6000_pLN_mLN <- read.csv("/Users/Pezoldt/PowerFolders/R/2018_ATAC_WGBS_RNA/2017_FSC_LEC_BEC/Data/DESeq2_genes_diffexp_by_fc_mLN_SPF_FRC_vs_pLN_SPF_FRC.csv")
#FRC_LI6000_pLN_mLN <- FRC_LI6000_pLN_mLN[!(is.na(FRC_LI6000_pLN_mLN$GeneSymbol)),]
#log2FC = 1.0
#padjust = 0.05
#UP_mLN_LI6000 <- as.character(subset(FRC_LI6000_pLN_mLN, log2FoldChange < -log2FC & padj < padjust)$GeneSymbol)
#UP_pLN_LI6000  <- as.character(subset(FRC_LI6000_pLN_mLN, log2FoldChange > log2FC & padj < padjust)$GeneSymbol)


##################
#ATACseq
##################
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
#DEseq2 ATAC-seq
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

##################
#Compare DEG vs. DMR
##################
colnames(SPF_mLN_pLN_Crossmapped)
#significant DMRs
DMR_sig_SPF_mLN_pLN <- subset(SPF_mLN_pLN_Crossmapped, n >= 5 & (meanDiff >= 0.40 | meanDiff <= -0.40) &
                                distance_to_nearest_tss <= 10000)
DMR_sig_SPF_mLN_pLN$nearest_gene_name <- as.character(DMR_sig_SPF_mLN_pLN$nearest_gene_name) 
split_DMR_sig_SPF_mLN_pLN <- split(DMR_sig_SPF_mLN_pLN, DMR_sig_SPF_mLN_pLN$nearest_gene_name)
#Quick check if typ of methylation is unidirectional
#Note: Within TSS for mLN vs. pLN that is exclusivly the case

#check <- lapply(split_DMR_sig_SPF_mLN_pLN, function(x){
# length(levels(as.factor(as.character(x$direction))))
#print(n_levels)
#if(n_levels = 1){
#  print(x$nearest_gene_name)
#}
#if(n_levels = 1){
#  print("consistent")
#}else{
#  print("arbitrary")
#}
#})

#Take first DMR per Gene Locus
Meaned_split_DMR_sig_SPF_mLN_pLN <- lapply(split_DMR_sig_SPF_mLN_pLN, function(x){
  mean_per_gene <- mean(x$meanDiff)
  #For adjacent data use first element
  data.frame(GeneSymbol = x$nearest_gene_name[1],
             mean_dMeth_gene = mean_per_gene,
             distance_tss = x$distance_to_nearest_tss[1])
})

Meaned_DMR_sig_SPF_mLN_pLN <- do.call(rbind.data.frame, Meaned_split_DMR_sig_SPF_mLN_pLN)

overlap_DMR_DEG <- subset(diff_Genes, GeneSymbol %in% Meaned_DMR_sig_SPF_mLN_pLN$GeneSymbol)
#Merge DMR and DEG table
merge_DMR_DEG <- merge(diff_Genes, Meaned_DMR_sig_SPF_mLN_pLN,
                       by = "GeneSymbol")
#Note: Not completely consistent with the dogma of Demeth makes expression go high,
#      but OK
plot(merge_DMR_DEG$SPF_mLN_pLN_log2FoldChange, merge_DMR_DEG$mean_dMeth_gene)
mLN_DMR_DEG <- subset(merge_DMR_DEG, mean_dMeth_gene < 0 & SPF_mLN_pLN_log2FoldChange > 0)

####################
#DMR / DAR
####################
#Keep only DARs with annotated GeneSymbol and close to the TSS -10000 +200
DARs_features_TSS <- DARs_features[!(is.na(DARs_features$symbol)),]
DARs_features_TSS <- subset(DARs_features_TSS, shortestDistance <= 10000)
DARs_features_TSS_sig <- subset(DARs_features_TSS, padj_SPF_mLN_pLN < 0.05)
#Merge DMR and DAR Table
DARs_DMRs_TSS <- merge(DARs_features_TSS, Meaned_DMR_sig_SPF_mLN_pLN,
                       by.x = "symbol", by.y = "GeneSymbol")

#####
#FC plot
#####
plot(DARs_DMRs_TSS$log2FC_SPF_mLN_pLN, DARs_DMRs_TSS$mean_dMeth_gene)

#####
#Overlap between DMRs and Peaks/DARs
#####
#Make Granges object of DARs
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
#Make Granges object of DMRs
gr_DMRs_TSS <- makeGRangesFromDataFrame(SPF_mLN_pLN,
                                        seqnames.field=c("seqnames", "seqname",
                                                         "chromosome", "chrom",
                                                         "chr", "chromosome_name",
                                                         "seqid"),
                                        start.field="start",
                                        end.field=c("end"),
                                        strand.field="strand",
                                        starts.in.df.are.0based=FALSE,
                                        keep.extra.columns = TRUE)
#Using DARs as reference
# Note: Check whether CrossMap makes a difference, currently 7 hits in the TSS
hits <- findOverlaps(gr_DARs_TSS, gr_DMRs_TSS, minoverlap = 50)
#check hits
gr_DARs_TSS_overlap <- gr_DARs_TSS[queryHits(hits)]
gr_DMRs_TSS_overlap <- gr_DMRs_TSS[subjectHits(hits)]
DARs_TSS_overlap <- as.data.frame(gr_DARs_TSS_overlap)
DMRs_TSS_overlap <- as.data.frame(gr_DMRs_TSS_overlap)

#Compute percent of overlap and thresh
overlaps <- pintersect(gr_DARs_TSS[queryHits(hits)], gr_DMRs_TSS[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(gr_DMRs_TSS[subjectHits(hits)])
hits <- hits[percentOverlap > 0.05]

################
#
################


