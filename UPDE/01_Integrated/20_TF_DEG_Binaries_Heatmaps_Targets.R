# dissect compiled and curated Output from Homer
# Author: Joern Pezoldt
# Date: 12.10.2018
# Functions:
# 1) Compile known motifs from homer output into one motif list
# 2) Use "homer" to annotate motifs to peaks
# 3) Obtain Peaks that contain motif
# 4) Check to which gene Peaks belong 
# 5) Classify TFs into Core and Unique sets per condition
# 6) Plot expression of genes for Core and unique TF sets
# 7) Plot Binaries of presence of TFs
# 8) Identify DEG TFs



#Libraries
library(UpSetR)
library(pheatmap)
library(biomaRt)
library(pheatmap)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(GenomicFeatures)
library(ChIPpeakAnno)
library(stringr)
library(foreign)
library(GO.db)
library(topGO)
library(org.Mm.eg.db)
library(ggplot2)

#####
#Global & PATHs
#####
table_ID <- "SPF_known_homer_compilation_curated.txt"
sample_ID <- "SPF"
PATH_input <- paste("/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/Motif/Homer_Output/allHomerResults/size_500bp/", sample_ID, sep = "")
table_ID <- "SPF_known_homer_compilation_curated.txt"
# RNAseq DESeq2 analysis
path_RNAseq_DESeq2 <- "/home/pezoldt/NAS2/pezoldt/Data/RNAseq/2017_FSC_LEC_BEC_mLN_pLN_GF_SPF/FSC"
# Note: Input required
PATH_general <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/Motif/Homer_Output/allHomerResults/size_500bp/SPF/"
PATH_motifs <- "/knownResults/"
homer_find_conditions <- "si500_v46_bgTSS_noExt_"
conditions <- c("Closed_DEG_SPF_lowRNAseq",
                "Open_DEG_SPF_lowRNAseq",
                "No_DAR_DOWN_SPF_lowRNAseq",
                "No_DAR_UP_SPF_lowRNAseq",
                "NoDEG_Closed_SPF_lowRNAseq",
                "NoDEG_Open_SPF_lowRNAseq")
PATH_TF_motifs_to_TF_GeneSybmol <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/Motif/Homer_Output/allHomerResults/size_500bp/SPF/SPF_known_homer_compilation_curated.txt" 
PATH_TF_motifs_in_gene <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/Motif/Regions_by_Motif"
PATH_output_RDS <- "/home/pezoldt/NAS2/pezoldt/Analysis/01_Integrated/DARs_DEGs"
#Global variables
log2FC_RNA = 1.0
padj = 0.05

#Databases
#All features
mm10_TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#TSS sites BioMart
mm10 = useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl") 
TSS.mouse.mm10 = getAnnotation(mart=mm10, featureType="TSS")
Exon.mouse.mm10 = getAnnotation(mart=mm10, featureType="Exon")

#Global variables
# score for motifs obtained for mapping motifs to peaks
# Note: distribution of hist(t_TF_binding_sites_i$MotifScore) yields 5 as a robust cutoff
motif_score = 5

#####
### Check for Unique TFs ---------------------------------------
#####
# load
t_homer_TF <- read.table(paste(PATH_input,"/",table_ID, sep=""), sep="\t",dec=".",header=TRUE)


# split according to comparisons
l_homer_TF <- split(t_homer_TF, t_homer_TF$group)
#rename
l_names <- strsplit(names(l_homer_TF), "noExt_")
names(l_homer_TF) <- unlist(sapply(l_names, function(x) {x[2]}))
#Delete global DAR sets
l_homer_TF <- l_homer_TF[c("Closed_DEG_SPF_lowRNAseq","Open_DEG_SPF_lowRNAseq","No_DAR_DOWN_SPF_lowRNAseq","No_DAR_UP_SPF_lowRNAseq",
                           "NoDEG_Closed_SPF_lowRNAseq","NoDEG_Open_SPF_lowRNAseq")]
names(l_homer_TF) <- c("pLN_Open_UP","mLN_Open_UP","pLN_peak_UP","mLN_peak_UP",
                       "pLN_Open_None","mLN_Open_None")

# grab TFs per comparison
l_motifs_per_group <- lapply(l_homer_TF, function(x){
  #x <- l_homer_TF[[1]]                              
  TFs_i <- as.character(unlist(strsplit(as.character(x$Motifs_TF), ";")))
})
#Generate binary Matrix for upset()
t_TFs_binary <- as.data.frame.matrix(table(stack(l_motifs_per_group)))
#Replace 
t_TFs_binary[t_TFs_binary >= 2] <- 1
t_TFs_binary <- t_TFs_binary[,c("pLN_Open_UP","mLN_Open_UP","pLN_peak_UP","mLN_peak_UP","pLN_Open_None","mLN_Open_None")]
# Plot intersection
upset(t_TFs_binary, sets = colnames(t_TFs_binary),
      mainbar.y.label = "TFs ATACseq",
      number.angles = 30, point.size = 5,
      text.scale = 2, keep.order = TRUE)

# Grab intersection of interest
# Note: Only one intersection interesting
Unique_TFs_NoDEG_Open <- subset(t_TFs_binary, 
                                pLN_Open_UP == 0 &
                                  mLN_Open_UP == 0 &
                                  pLN_peak_UP == 1 &
                                  mLN_peak_UP == 0 &
                                  pLN_Open_None == 0 &
                                  mLN_Open_None == 0)

# Write binary table
write.table(t_TFs_binary, paste(PATH_input,"/",sample_ID,"_existence_matrix_known_homer_compile.txt",sep = "") ,dec=".", sep="\t")

#####Make TF binding site heatmap----------------------------------------------
#Get TF and pValue
l_homer_TF_pValue <- list()
for(i in 1:length(l_homer_TF)){
  name_i <- names(l_homer_TF)[i]
  print(name_i)
  t_i <- l_homer_TF[[i]][,c("Motifs_TF","q.value")]
  colnames(t_i) <- c("Motifs_TF", paste("pValue_",name_i,sep=""))
  l_homer_TF_pValue[[i]] <- t_i
  names(l_homer_TF_pValue)[i] <- name_i
}

# Merge tables by Motifs
t_homer_TF_pValue <- Reduce(function(...) merge(..., all=T), l_homer_TF_pValue)

# Condense per TF by meaning accross conditions
l_homer_TF_pValue <- split(t_homer_TF_pValue, t_homer_TF_pValue$Motifs_TF)
# Replace NAs with 1 and mean
l_homer_TF_pValue_mean <- lapply(l_homer_TF_pValue, function(x){
  x[is.na(x)] <- 1
  l_homer_TF_pValue_i <- colMeans(x[,2:ncol(x)])
  l_homer_TF_pValue_i
})

#rowbind list and add TF names
t_homer_TF_pValue_mean <- do.call("rbind", l_homer_TF_pValue_mean)
rownames(t_homer_TF_pValue_mean) <- names(l_homer_TF_pValue)
#eliminate first row as it is not a TF
t_homer_TF_pValue_mean <- t_homer_TF_pValue_mean[2:nrow(t_homer_TF_pValue_mean),]
#eliminate NaN which are derived from the orignial input, that contained the two lists of DARs only
t_homer_TF_pValue_mean <- t_homer_TF_pValue_mean[complete.cases(t_homer_TF_pValue_mean), ]

#-log10 of pValue for plotting
t_homer_TF_pValue_mean <- -log10(t_homer_TF_pValue_mean)
# Replace outlier high pValues with 12
t_homer_TF_pValue_mean[t_homer_TF_pValue_mean == Inf] <- 4
t_homer_TF_pValue_mean <- t_homer_TF_pValue_mean[order(rownames(t_homer_TF_pValue_mean)), ]

output_pheatmap <- pheatmap(t_homer_TF_pValue_mean, cluster_rows = TRUE, legend = TRUE,
                            treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = FALSE,
                            scale = "none", border_color = "black", cellwidth = 10,
                            cellheigth = 10, color = colorRampPalette(c("white","darkolivegreen1","darkolivegreen2","darkolivegreen3","darkolivegreen3","darkolivegreen4","darkgreen","darkgreen"), space="rgb")(128),
                            main = "TF by motif -log10(qValue)",
                            silent = TRUE)


####################
#RNAseq data
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
                  DF[,grepl(c("padj"), colnames(DF))],
                  DF[,grepl(c("log2FC"), colnames(DF))],
                  DF[,c("RPKMcounts.mLN_GF_FRC_1.x","RPKMcounts.mLN_GF_FRC_2.x","RPKMcounts.mLN_GF_FRC_3.x",
                        "RPKMcounts.pLN_GF_FRC_1.x","RPKMcounts.pLN_GF_FRC_2.x","RPKMcounts.pLN_GF_FRC_3.x",
                        "RPKMcounts.mLN_SPF_FRC_1.y","RPKMcounts.mLN_SPF_FRC_2.y","RPKMcounts.mLN_SPF_FRC_3.y",
                        "RPKMcounts.pLN_SFP_FRC_1.x","RPKMcounts.pLN_SFP_FRC_2.x","RPKMcounts.pLN_SFP_FRC_3.x")])
colnames(data_all)[1] <- "GeneSymbol"
#Eliminate columns
data_all <- data_all[,c(1,2,4:7,9:ncol(data_all))]
data_all <- data_all[!(is.na(data_all$GeneSymbol)),]
data_all <- data_all[!duplicated(data_all$GeneSymbol),]
#reorder
data_all <- data_all[,c(1,6:9,2:5,10:ncol(data_all))]
colnames(data_all) <- c("GeneSymbol",
                        "GF_mLN_pLN_log2FoldChange","mLN_SPF_GF_log2FoldChange",
                        "SPF_mLN_pLN_log2FoldChange","pLN_SPF_GF_log2FoldChange",
                        "GF_mLN_pLN_padj","mLN_SPF_GF_padj",
                        "SPF_mLN_pLN_padj","pLN_SPF_GF_padj",
                        "RPKM_mLN_GF_1","RPKM_mLN_GF_2","RPKM_mLN_GF_3",
                        "RPKM_pLN_GF_1","RPKM_pLN_GF_2","RPKM_pLN_GF_3",
                        "RPKM_mLN_SPF_1","RPKM_mLN_SPF_2","RPKM_mLN_SPF_3",
                        "RPKM_pLN_SPF_1","RPKM_pLN_SPF_2","RPKM_pLN_SPF_3")

#Reorder dataframe
data_all <- data_all[,c("GeneSymbol",
                        "SPF_mLN_pLN_log2FoldChange","GF_mLN_pLN_log2FoldChange",
                        "mLN_SPF_GF_log2FoldChange","pLN_SPF_GF_log2FoldChange",
                        "SPF_mLN_pLN_padj","GF_mLN_pLN_padj",
                        "mLN_SPF_GF_padj","pLN_SPF_GF_padj",
                        "RPKM_mLN_GF_1","RPKM_mLN_GF_2","RPKM_mLN_GF_3",
                        "RPKM_pLN_GF_1","RPKM_pLN_GF_2","RPKM_pLN_GF_3",
                        "RPKM_mLN_SPF_1","RPKM_mLN_SPF_2","RPKM_mLN_SPF_3",
                        "RPKM_pLN_SpF_1","RPKM_pLN_SPF_2","RPKM_pLN_SPF_3")]
#store
expression <- data_all

######
#Check expression of TFs
######
#Get rownames in order of pValue heatmap
TFs_by_Motif <- rownames(t_homer_TF_pValue_mean[output_pheatmap$tree_row[["order"]],])
TFs_by_Motif <- unlist(lapply(TFs_by_Motif, function(x){
  unlist(strsplit(x, ";", fixed = TRUE))
}))
TFs_by_Motif <- unique(TFs_by_Motif)
#Get expression values TFs identified
expression_TFs <- subset(expression, GeneSymbol %in% TFs_by_Motif)

#Sort Table
#expressed_TFs_ordered <- TFs_by_Motif[TFs_by_Motif %in% expression_TFs$GeneSymbol]
#expressed_TFs_ordered <- TFs_by_Motif[match(TFs_by_Motif, expressed_TFs_ordered)]
#expression_TFs <- expression_TFs[expressed_TFs_ordered,]
#expressed_TFs_ordered <- expressed_TFs_ordered[!is.na(expressed_TFs_ordered)]
#expression_TFs <- expression_TFs[match(expression_TFs$GeneSymbol, expressed_TFs_ordered),]

#Prep data for heatmap
data_heatmap <- expression_TFs
rownames(data_heatmap) <- as.character(data_heatmap$GeneSymbol)
data_heatmap <- data_heatmap[,10:21]
data_heatmap <- as.matrix(data_heatmap)
data_heatmap <- log2(data_heatmap)
#Find minimum generically
find_min_data <- data_heatmap
find_min_data[find_min_data == -Inf] <- 1000
gen_minimum <- floor(min(find_min_data))
data_heatmap[data_heatmap == -Inf] <- gen_minimum
#Mean replicates
meanRPKM_mLN_GF <- rowMeans(data_heatmap[,1:3])
meanRPKM_pLN_GF <- rowMeans(data_heatmap[,4:6])
meanRPKM_mLN_SPF <- rowMeans(data_heatmap[,7:9])
meanRPKM_pLN_SPF <- rowMeans(data_heatmap[,10:12])
data_heatmap <- cbind(meanRPKM_mLN_GF,meanRPKM_pLN_GF,meanRPKM_mLN_SPF,meanRPKM_pLN_SPF)

pheatmap(data_heatmap, cluster_rows = FALSE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 30, show_rownames = TRUE, cluster_cols = TRUE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("blue", "white","red"), space="rgb")(128))

#####
#Compare TFs
#####
data_heatmap <- expression_TFs
rownames(data_heatmap) <- as.character(data_heatmap$GeneSymbol)
data_heatmap <- data_heatmap[,16:21]

#Calculate pValues
p_Values_expression_TF <- c()
for(i in 1:nrow(data_heatmap)){
  data_heatmap_row_i <- data_heatmap[i,]
  p_Value_i <- t.test(x = data_heatmap_row_i[,1:3],y =  data_heatmap_row_i[,4:6])$p.value
  p_Values_expression_TF[i] <- p_Value_i
}

#FDR adjust pValues
pAdjust_expression_TF <- p.adjust(p_Values_expression_TF, method = "fdr")
names(pAdjust_expression_TF) <- rownames(data_heatmap)

#Select significantly DEGs
sig_expression_TFs <- names(pAdjust_expression_TF[pAdjust_expression_TF < 0.10])

#Select expression data for significant pValues
data_heatmap_sig <- subset(data_heatmap, rownames(data_heatmap) %in% sig_expression_TFs)

#Set groups
group_1 = "mLN"
group_2 = "pLN"
v_classifier <- c(rep(group_1, 3),rep(group_2,3))
#List to store graphs
l_ggplot_images <- list()
for(i in 1:nrow(data_heatmap_sig)){
  TF_expression_i <- data_heatmap_sig[i,]
  print(i)
  name_TF_i <- rownames(TF_expression_i)
  print(name_TF_i)
  #vector of expression values
  Expression_i <- as.numeric(TF_expression_i)
  
  #Jitterplot
  t_ggplot_i <- data.frame(v_classifier = v_classifier, Expression_i = Expression_i)
  p <- ggplot(t_ggplot_i, aes(x=v_classifier, y=Expression_i, shape=v_classifier)) + geom_jitter(position=position_jitter(0.1), cex=3) + 
    scale_color_manual(values = c("blue", "red")) +
    scale_y_continuous(limits = c(floor(min(Expression_i)), ceiling(max(Expression_i)))) + 
    scale_shape_manual(values=c(1,19)) +
    xlab(name_TF_i) +
    ylab("RPKM") +
    theme(legend.position="none")
  l_ggplot_images[[i]] <- p
}
#
multiplot(plotlist = l_ggplot_images, cols = 5)

#TFs but not DEGs
TFs_not_DEG <- TFs_by_Motif[!(TFs_by_Motif %in% rownames(data_heatmap_sig))]

#####
#Functions
#####
#Concatenate motifs -----------------------------------
sampleNames_complete <- paste(homer_find_conditions, conditions, sep = "")
PATH_sampleNames <- paste(PATH_general,sampleNames_complete,PATH_motifs,sep = "")
#' concatenateKnownMotifs
#' loads single motif files from homer (known motifs) and compiles them into one txt file
#' @param sampleNames_complete cahractter vector with path to motif files
#' @param PATH_general string to store all files
#' @param conditions character vector of sample names
#' @return concatenated motif tables within "knownMotifs" and higher folder

concatenateKnownMotifs <- function(PATH_sampleNames,conditions,PATH_general){
  for(i in 1:length(PATH_sampleNames)){
    PATH_sampleNames_i <- PATH_sampleNames[i]
    setwd(PATH_sampleNames_i)
    names_i <- list.files(pattern="*.motif")
    names_i <- paste(PATH_sampleNames_i,names_i,sep = "")
    myfiles_i <- lapply(names_i, read.delim, header = FALSE)
    motif_file_i <- do.call("rbind", myfiles_i)
    motif_file_i[is.na(motif_file_i)] <- ""
    colnames(motif_file_i) <- NULL
    write.table(motif_file_i,paste(PATH_sampleNames_i,conditions[i],"_knownMotif_compiled.txt", sep = ""), sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)
    write.table(motif_file_i,paste(PATH_general,conditions[i],"_knownMotif_compiled.txt", sep = ""), sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)
  }
}

#Run Function
concatenateKnownMotifs(PATH_sampleNames,conditions,PATH_general)


#Identify genes targeted by TFs--------------------------------------------------------------------------
# Explanation Input:
# Input 1: TXT with the peak identifiers, the start and the stop regions (and TF bindings per peak)
# Input 2: TXT with altered peak identifier (e.g. preceding BG), offset from peak center and Motif name
# ToDo's: 
#     1) Align TFs by Motif and not by name
#     2) 
#' pinpointTFMotifstoGenes
#' generates RDS:
#' 1) containing a table that counts the incidence number (motif occurences) for each TF within each gene 
#' @param PATH_sampleNames cahractter vector with path to motif files
#' @param PATH_TF_motifs_to_TF_GeneSybmol character vector with the path to the curated motif file annotating the TFs to the GeneSymbol's
#' @param TSS.mouse.mm10 bimoart annotation for TSS
#' @return l_TF_binding_freq returns table with incidence rates of TF binding per gene and replaces the home dervied annotation of the TF/motif with the curated TF GeneSymbols
PATH_samples_TF_Motif_Gene <- paste(PATH_TF_motifs_in_gene,"/",conditions,sep = "")


######Continue here: tracking of homer ID with TF ID is not functioning
pinpointTFMotifstoGenes <- function(PATH_samples_TF_Motif_Gene, conditions,PATH_TF_motifs_to_TF_GeneSybmol){
  l_TF_binding_freq <- list()
  for(i in 1:length(PATH_samples_TF_Motif_Gene)){
    PATH_samples_TF_Motif_Gene_i <- PATH_samples_TF_Motif_Gene[i]
    conditions_i <- conditions[i]
    print(i)
    #Load TF binding sits
    #t_TF_binding_sites_i <- read.delim(paste(PATH_samples_TF_Motif_Gene_i,"/",conditions_i,"_MotifInstances.txt",sep = ""),stringsAsFactors=FALSE)
    
    #Load Peak location
    t_Peak_location_i <- read.delim(paste(PATH_samples_TF_Motif_Gene_i,"/",conditions_i,"_MotifInstances_Location.txt",sep = ""),stringsAsFactors=FALSE,
                                    header = FALSE)
    #Change Column for TFs
    colnames_i <- as.character(t(as.vector(t_Peak_location_i[1,])))
    colnames_i <- unlist(lapply(strsplit(colnames_i, "\\("), function(x){
      x[[1]]
    }))
    colnames_i <- unlist(lapply(strsplit(colnames_i, "\\/"), function(x){
      x[[1]]
    }))
    colnames_descriptive_i <- colnames_i[1:21]
    #Eliminate SeqBias: columns
    #colnames_i <- unlist(lapply(strsplit(colnames_i, "Bias:"), function(x){
    # x[[1]]
    #}))
    #Load table for annotating homer Motif names to TF GeneSybmols
    TF_motifs_to_TF_GeneSybmol <- read.delim(PATH_TF_motifs_to_TF_GeneSybmol,dec = ".")
    MotifTF_to_GeneSymbolTF <- TF_motifs_to_TF_GeneSybmol[,c("MotifName","Motifs_TF")]
    MotifTF_to_GeneSymbolTF_unique <- MotifTF_to_GeneSymbolTF[!duplicated(MotifTF_to_GeneSymbolTF),]
    #colnames(t_Peak_location_i) <- colnames_i
    #By Name
    v_colnames_new_i_j <- c()
    for(j in 1:length(colnames_i[22:length(colnames_i)])){
      colnames_i_j <-  colnames_i[21+j]
      colnames_new_i_j <- as.character(subset(MotifTF_to_GeneSymbolTF_unique, MotifName  == colnames_i_j)$Motifs_TF)
      colnames_new_i_j <- str_replace_all(string=colnames_new_i_j, pattern=" ", repl="")
      print(colnames_new_i_j)
      if(identical(colnames_new_i_j, character(0))){
        colnames_new_i_j <- paste("NOT_a_TF_",j,sep = "")
      }
      if(is.na(colnames_new_i_j) == FALSE){
        colnames_new_i_j <- colnames_new_i_j
      }
      if(colnames_new_i_j == ""){
        colnames_new_i_j <- paste("NOT_a_TF_",j,sep="")
      }
      v_colnames_new_i_j <- c(v_colnames_new_i_j,colnames_new_i_j)
    }
    
    colnames_total_i <- c(colnames_descriptive_i,v_colnames_new_i_j)
    t_Peak_location_i <- t_Peak_location_i[2:nrow(t_Peak_location_i),]
    colnames(t_Peak_location_i) <- colnames_total_i
    
    colnames(t_Peak_location_i)[1] <- "ID"
    colnames(t_Peak_location_i)[colnames(t_Peak_location_i) == ""] <- "No_TF"
    #t_Peak_location_i <- t_Peak_location_i[,2:ncol(t_Peak_location_i)]
    #Annotate Peak to gene
    gr_Peaks_i <- toGRanges(t_Peak_location_i, names = t_Peak_location_i[,1])
    #Annotate Peaks
    gr_Peaks_TSS_i <- annotatePeakInBatch(gr_Peaks_i, AnnotationData=TSS.mouse.mm10)
    #add gene name
    gr_Peaks_TSS_GeneName_i <- addGeneIDs(annotatedPeak=gr_Peaks_TSS_i, 
                                          feature_id_type="ensembl_gene_id",
                                          orgAnn="org.Mm.eg.db", 
                                          IDs2Add="symbol")
    #Introduces ten new columns
    t_Peaks_TSS_GeneName_i <- as.data.frame(gr_Peaks_TSS_GeneName_i)
    #Rename the column names
    #Rearrange table for counting the TF binding events per TF
    t_Peaks_TSS_GeneName_i <- cbind(t_Peaks_TSS_GeneName_i[,1:22],
                                    t_Peaks_TSS_GeneName_i[,ncol(t_Peaks_TSS_GeneName_i)],
                                    #t_Peaks_TSS_GeneName_i[,(ncol(t_Peaks_TSS_GeneName_i)-12):ncol(t_Peaks_TSS_GeneName_i)],
                                    t_Peaks_TSS_GeneName_i[,23:(ncol(t_Peaks_TSS_GeneName_i)-8)])
    #Determine number of TF binding sites per peak/Gene per TF
    colnames(t_Peaks_TSS_GeneName_i)[23] <- "GeneSymbol"
    colnames(t_Peaks_TSS_GeneName_i)[24:ncol(t_Peaks_TSS_GeneName_i)] <- colnames_total_i[22:length(colnames_total_i)]
    t_Peaks_TSS_GeneName_TFs_i <- t_Peaks_TSS_GeneName_i[,24:ncol(t_Peaks_TSS_GeneName_i)]
    #string split for each TF for each Peak
    numbers_i <- matrix(ncol = ncol(t_Peaks_TSS_GeneName_TFs_i), nrow = 0)
    for(j in 1:nrow(t_Peaks_TSS_GeneName_TFs_i)){
      row_j <- as.character(t_Peaks_TSS_GeneName_TFs_i[j,] )
      n_j <- unlist(lapply(strsplit(row_j, '[(]'),function(x){
        length(x) - 1
      }))
      n_j[n_j == -1] <- 0
      numbers_i <- rbind(numbers_i, n_j)
    }
    numbers_i <- as.data.frame(numbers_i)
    t_Peaks_TSS_GeneName_i[,24:ncol(t_Peaks_TSS_GeneName_i)] <- numbers_i
    print(nrow(t_Peaks_TSS_GeneName_i))
    l_TF_binding_freq[[i]] <- t_Peaks_TSS_GeneName_i
  }
  #returned
  names(l_TF_binding_freq) <- conditions
  l_TF_binding_freq
}

#Run function
l_TF_binding_freq <- pinpointTFMotifstoGenes(PATH_samples_TF_Motif_Gene, conditions,PATH_TF_motifs_to_TF_GeneSybmol)


# Adjust columns -----------------------------------
#' splitMotifGroups
#' Split multiple TFs and assign binding frequencies to all TFs of Motif group
#' @param l_TF_binding_freq list of tables annotating the number of binding motifs per Motif/TF for promotor region
#' @return list of tables for each set of DARs with each column containing the TF as GeneSymbol

splitMotifGroups <- function(l_TF_binding_freq){
  l_TF_binding_freq_cleared <- list()
  for(i in 1:length(l_TF_binding_freq)){
    t_TF_binding_freq_i <- l_TF_binding_freq[[i]]
    rownames_i <- rownames(t_TF_binding_freq_i)
    descriptive_columns_i <- t_TF_binding_freq_i[,c(1:23)]
    TF_columns_i <- t_TF_binding_freq_i[,c(24:ncol(t_TF_binding_freq_i))]
    
    t_reformed_TF_i=data.frame(matrix(ncol=0,nrow=nrow(t_TF_binding_freq_i)))
    for(j in 1:ncol(TF_columns_i)){
      #print(j)
      column_name_j <- colnames(TF_columns_i)[j]
      #print(column_name_j)
      split_column_name_j <- unlist(strsplit(column_name_j, "\\;"))
      #print(split_column_name_j)
      if(length(split_column_name_j) > 1){
        print(unlist(split_column_name_j))
        t_split_TFs <- TF_columns_i[j]
        for(k in 1:(length(split_column_name_j)-1)){
          t_split_TFs <- cbind(t_split_TFs,TF_columns_i[j])
        }
        colnames(t_split_TFs) <- split_column_name_j
      }
      if(length(split_column_name_j) == 1){
        t_split_TFs <- TF_columns_i[j]
      }
      t_reformed_TF_i <- cbind(t_reformed_TF_i,t_split_TFs)
      rownames(t_reformed_TF_i) <- paste(t_TF_binding_freq_i$GeneSymbol, 1:length(t_TF_binding_freq_i$GeneSymbol))
    }
    t_reformed_TF_i <- cbind(t_TF_binding_freq_i$GeneSymbol, t_reformed_TF_i)
    #Eliminate ".1" in TF name
    colnames(t_reformed_TF_i) <- unlist(lapply(strsplit(colnames(t_reformed_TF_i), "\\."), function(x){x[[1]]}))
    l_TF_binding_freq_cleared[[i]] <- t_reformed_TF_i
  }
  names(l_TF_binding_freq_cleared) <- conditions
  l_TF_binding_freq_cleared
}

#Run function
l_TF_binding_freq_cleared <- splitMotifGroups(l_TF_binding_freq)

#####
#Exploratory
#####
#Plot Binding incidence of TFs---------------------------------------
#Generate Heatmap of log2 number of TF binding sites per gene
index = 3
test_TF <- l_TF_binding_freq_cleared[[index]]
data_heatmap <- test_TF
#drops <- c("Seq","GAGA.repeat")
#data_heatmap <- data_heatmap[ , !(names(data_heatmap) %in% drops)]
rownames(data_heatmap) <- paste(data_heatmap$GeneSymbol, 1:length(data_heatmap$GeneSymbol),sep="")
data_heatmap <- data_heatmap[,24:ncol(data_heatmap)]
data_heatmap <- as.matrix(data_heatmap)
#data_heatmap <- log2(data_heatmap)
data_heatmap[data_heatmap == -Inf] <- 0
pheatmap(data_heatmap, scale = "none", cluster_cols = FALSE,
         color = colorRampPalette(c("white", "lightgoldenrod1","chartreuse3","chartreuse4"), space="rgb")(128))

#Plot binding incidence of TFs that are DEGs (script 55_...)-----------------------
DEGs_TFs <- rownames(data_heatmap_sig)
TF_binding_numbers <- l_TF_binding_freq_cleared[[index]]
TF_binding_numbers_DEG_TF <- TF_binding_numbers[,colnames(TF_binding_numbers) %in% DEGs_TFs]
TF_binding_numbers_DEG_TF <- as.matrix(TF_binding_numbers_DEG_TF)
#Keep columns with inflammatory TFs
# ToDo: Check inflammatory response element selection
TF_DEG_inflammatory <- c("Irf1","Nfkb2","Irf3","Irf5","Irf7","Irf8")
TF_DEG_NONinflammatory <- c("Bach1","Cebpa","E2f1","Ebf2","Isl1","Jun","Klf2","Klf3","Klf4","Klf5","Klf7","Xbp1")
TF_binding_numbers_DEG_TF <- TF_binding_numbers_DEG_TF[,colnames(TF_binding_numbers_DEG_TF) %in% TF_DEG_inflammatory]
#delete all rows with only Zeros
TF_binding_numbers_DEG_TF <- TF_binding_numbers_DEG_TF[apply(TF_binding_numbers_DEG_TF[,-1], 1, function(x) !all(x==0)),]

pheatmap(TF_binding_numbers_DEG_TF, scale = "none", cluster_cols = FALSE, show_rownames = TRUE,
         color = colorRampPalette(c("white", "lightgoldenrod1","chartreuse3","chartreuse4"), space="rgb")(128),
         main = names(l_TF_binding_freq[index]))

#####
#Create Binary Map of DEG TF binding in gene sets
#####

# Make binary table comparing vector with colnames -----------------------------------
#' makeBinaryVectorVSColnames
#' Compares character vector against list of tables the colnames
#' @param vector_to_compare character vector
#' @param l_tables names list of tables with colnames to be compared to vector
#' @return binary table indicating presence or absence o
makeBinaryVectorVSColnames <- function(vector_to_compare,l_tables){
  #setup table
  t_binary <- data.frame(matrix(nrow=0, ncol=length(vector_to_compare)))
  for(i in 1:length(l_tables)){
    tables_i <- l_tables[[i]]
    #compare and binarize
    binary_i <- (vector_to_compare %in% colnames(tables_i)) * 1
    #row bind
    t_binary <- rbind(t_binary,binary_i)
  }
  #Name elements
  colnames(t_binary) <- vector_to_compare
  rownames(t_binary) <- names(l_tables)
  #return
  t_binary
}

#Run function
t_binary_noDEG_TFs_in_sets <- makeBinaryVectorVSColnames(TFs_not_DEG,l_TF_binding_freq_cleared)
t_binary_DEG_TFs_in_sets <- makeBinaryVectorVSColnames(DEGs_TFs,l_TF_binding_freq_cleared)

pheatmap(as.matrix(t_binary_noDEG_TFs_in_sets), cluster_rows = FALSE,
         color = colorRampPalette(c("white", "black"), space="rgb")(2),
         treeheight_row = 0, treeheight_col = 0,
         border_color = "black")

#####
#Genes Targeted by unique TFs binding in peak reagion
#####
# TFs DEGs associated Genes
t_binary_DEG_TFs_in_sets <- as.data.frame(t(t_binary_DEG_TFs_in_sets))
TFs_DEG_pLN <- rownames(subset(t_binary_DEG_TFs_in_sets, Closed_DEG_SPF_lowRNAseq == 0 &
                        Open_DEG_SPF_lowRNAseq == 0 &
                        No_DAR_DOWN_SPF_lowRNAseq == 1  & 
                        No_DAR_UP_SPF_lowRNAseq == 0 & 
                        NoDEG_Closed_SPF_lowRNAseq == 0 & 
                        NoDEG_Open_SPF_lowRNAseq == 0))
TFs_DEG_mLN <- rownames(subset(t_binary_DEG_TFs_in_sets, Closed_DEG_SPF_lowRNAseq == 0 &
                        Open_DEG_SPF_lowRNAseq == 0 &
                        No_DAR_DOWN_SPF_lowRNAseq == 0  & 
                        No_DAR_UP_SPF_lowRNAseq == 0 & 
                        NoDEG_Closed_SPF_lowRNAseq == 0 & 
                        NoDEG_Open_SPF_lowRNAseq == 1))

# TFs DEGs associated Genes
t_binary_noDEG_TFs_in_sets <- as.data.frame(t(t_binary_noDEG_TFs_in_sets))
TFs_noDEG_pLN <- rownames(subset(t_binary_noDEG_TFs_in_sets, 
                        Closed_DEG_SPF_lowRNAseq == 0 &
                        Open_DEG_SPF_lowRNAseq == 0 &
                        No_DAR_DOWN_SPF_lowRNAseq == 1  &
                        No_DAR_UP_SPF_lowRNAseq == 0 & 
                        NoDEG_Closed_SPF_lowRNAseq == 0 & 
                        NoDEG_Open_SPF_lowRNAseq == 0))
TFs_noDEG_mLN <- rownames(subset(t_binary_noDEG_TFs_in_sets, 
                        Closed_DEG_SPF_lowRNAseq == 0 &
                        Open_DEG_SPF_lowRNAseq == 0 &
                        No_DAR_DOWN_SPF_lowRNAseq == 0  & 
                        No_DAR_UP_SPF_lowRNAseq == 0 & 
                        NoDEG_Closed_SPF_lowRNAseq == 0 & 
                        NoDEG_Open_SPF_lowRNAseq == 1))

# Check for vector of TFs which Genes are targeted in table -----------------------------------
#' vectorCHECKrownames
#' Compares character vector against list of tables the colnames
#' @param vector_to_compare character vector
#' @param condition_DAR_DEG string defining condition to check
#' @param l_tables names list of tables with colnames to be compared to vector
#' @return named list of vectors containing the genes targeted by the TF

makeBinaryVectorVSColnames <- function(vector_to_compare,condition_DAR_DEG,l_tables){
  #get table of DAR / DEG quadrant 
  t_TFs_for_Genes <- l_tables[[condition_DAR_DEG]]
  l_Genes_for_TFs <- list()
  for(i in 1:length(vector_to_compare)){
    #TF to compare
    TF_i <- vector_to_compare[i]
    print(TF_i)
    #Genes bound by TF_i    
    Genes_i <- as.character(t_TFs_for_Genes[t_TFs_for_Genes[,TF_i] > 0,][,1])
    l_Genes_for_TFs[[i]] <- Genes_i
    print(Genes_i)
    names(l_Genes_for_TFs)[i] <- TF_i
  }
  l_Genes_for_TFs
}

Genes_by_DEG_TF_mLN_active <- makeBinaryVectorVSColnames(TFs_DEG_mLN,"NoDEG_Open_SPF_lowRNAseq",l_TF_binding_freq_cleared)
Genes_by_noDEG_TF_mLN_active <- makeBinaryVectorVSColnames(TFs_noDEG_mLN,"NoDEG_Open_SPF_lowRNAseq",l_TF_binding_freq_cleared)

Genes_by_DEG_TF_pLN_active <- makeBinaryVectorVSColnames(TFs_DEG_pLN,"No_DAR_DOWN_SPF_lowRNAseq",l_TF_binding_freq_cleared)
Genes_by_noDEG_TF_pLN_active <- makeBinaryVectorVSColnames(TFs_noDEG_pLN,"No_DAR_DOWN_SPF_lowRNAseq",l_TF_binding_freq_cleared)
Genes_by_noDEG_TF_pLN_active <- Genes_by_noDEG_TF_pLN_active[lapply(Genes_by_noDEG_TF_pLN_active,length)>0]

#####
#GO per TF
#####
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
#Input data
#####
myfiles <- Genes_by_noDEG_TF_pLN_active
#####
#GeneSymbol to GeneID
#####
for(i in 1:length(myfiles)){
  idfound <- myfiles[[i]] %in% mappedRkeys(org.Mm.egSYMBOL)
  SYMBOL <- toTable(org.Mm.egSYMBOL)
  head(SYMBOL)
  m <- match(myfiles[[i]], SYMBOL$symbol)
  GENE_ID <- SYMBOL$gene_id[m]
  GeneSymbol <- myfiles[[i]]
  myfiles[[i]] <- cbind(GENE_ID, GeneSymbol)
}
#####
#Perform GO analysis single
#####
#determine number of clusters: split by cluster
#link gene IDs with p-values
#perform GO
l_gene_pval <- list()
l_all_gene_pval <- list()
l_gene_List <- list()

#add a pvalue of 0.05 as it is not included in the uitlized Fisher-test
for(i in 1:length(myfiles)){
  names_i <- myfiles[[i]][,1]
  p_value_i <- rep(0.05, length(names_i))
  names(p_value_i) <- names_i
  l_gene_p_val <- p_value_i
  l_gene_pval[[i]] <- l_gene_p_val
}

names(l_gene_pval) <- names(myfiles)

#Gene universe is the entity of expressed genes
geneNames <- expression$GeneSymbol
#change GeneID
idfound <- geneNames %in% mappedRkeys(org.Mm.egSYMBOL)
SYMBOL <- toTable(org.Mm.egSYMBOL)
head(SYMBOL)
m <- match(geneNames, SYMBOL$symbol)
GENE_ID <- SYMBOL$gene_id[m]
geneNames <- GENE_ID

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
  colnames(table_i_subset_GO_weight)[2] <- names(myfiles)[i]
  table_i_subset_GO_weight[,2] <- as.numeric(table_i_subset_GO_weight[,2])
  table_i_subset_GO_weight <- table_i_subset_GO_weight[order(table_i_subset_GO_weight$"GO.ID", decreasing=TRUE), ]
  tail(table_i_subset_GO_weight)
  l_topGO[[i]] <- table_i_subset_GO_weight
}
#merger
data_TopGO_weight = Reduce(function(...) merge(..., all=T), l_topGO)

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
title <- c("Genes_Target_TF_noDEG_pLN_Open")

#reorder according to Heatmap in main figure
#data_heatmap_matrix <- data_heatmap_matrix[,subset_order]
saveRDS(data_heatmap_matrix, file = paste(PATH_output_RDS,"/Genes_Target_TF_noDEG_pLN_Open.Rds"))
data_heatmap_matrix <- readRDS(paste(PATH_output_RDS,"/Genes_Target_TF_noDEG_pLN_Open.Rds",sep=""))
pheatmap(data_heatmap_matrix, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = TRUE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("white", "grey","deepskyblue2"), space="rgb")(128),
         main = title)



