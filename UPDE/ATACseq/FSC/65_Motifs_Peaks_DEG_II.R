# Annotate motifs to cohorts of DARs and DEGs and explore expression

# Author: Joern Pezoldt
# Date: 26.10.2018

# Functions:
# 1) Compile known motifs from homer output into one motif list
# 2) Use "homer" to annotate motifs to peaks
# 3) Obtain Peaks that contain motif
# 4) Check to which gene Peaks belong 
# 5) Classify TFs into Core and Unique sets per condition
# 6) Plot expression of genes for Core and unique TF sets


#Library
library(biomaRt)
library(pheatmap)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(GenomicFeatures)
library(ChIPpeakAnno)
library(stringr)

#####
#PATHs and global variables
#####
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
test_TF <- l_TF_binding_freq[[index]]
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
TF_DEG_inflammatory <- c("Irf1","Nfkb2","Irf1","Irf3","Irf5","Irf7","Irf8")
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

vector_to_compare <- DEGs_TFs
l_tables <- l_TF_binding_freq_cleared

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
t_binary_DEG_TFs_in_sets <- makeBinaryVectorVSColnames(TFs_not_DEG,l_TF_binding_freq_cleared)

pheatmap(as.matrix(t_binary_DEG_TFs_in_sets), cluster_rows = FALSE,
         color = colorRampPalette(c("white", "black"), space="rgb")(2),
         treeheight_row = 0, treeheight_col = 0,
         border_color = black)

