# Process Genomatix Output to obtain lists of TFs
# Author: Joern Pezoldt
# Date: 23.11.2014
# Function:
# 1) Compiles genomatix Output to extract TFs based on Regx
# Note: i) For TFs containing "-" in the official GeneSymbol and for multiple PWM annotation
#        TFs (e.g. Tcfl1, Tcfl2) the character vectors need to be manually curated.
#       ii) Genomatix Output also conains for some TFs a only group name (e.g. Stat)

#Libraries
# Base

#####
#PATH & Global Variables
#####
PATH_input <- "/home/pezoldt/NAS2/pezoldt/Analysis/WGBS/FSC/Genomatix_TF/ANNa/Mouse/To_Use"
PATH_output <- "/home/pezoldt/NAS2/pezoldt/Analysis/WGBS/FSC/Genomatix_TF/ANNa/Mouse/To_Use/Output"

#Global variables
Z_score_overrepresentation <- 2

#####
#Load tables
#####
Genomatix_hyper_all <- read.delim(paste(PATH_input, "/Genomatix_hyper_all.txt",sep=""))
Genomatix_hyper_thresh <- read.delim(paste(PATH_input, "/Genomatix_hyper_thresh.txt",sep=""))
Genomatix_hypo_all <- read.delim(paste(PATH_input, "/Genomatix_hypo_all.txt",sep=""))
Genomatix_hypo_thresh <- read.delim(paste(PATH_input, "/Genomatix_hypo_thresh.txt",sep=""))
Genomatix_dataset_names <- c("Genomatix_hyper_all","Genomatix_hyper_thresh","Genomatix_hypo_all","Genomatix_hypo_thresh")

#####
#First up Character
#####
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

#####
#Extract TFs
#####
#Store tables in list
Genomatix_list <- list(Genomatix_hyper_all,Genomatix_hyper_thresh,Genomatix_hypo_all,Genomatix_hypo_thresh)

#Storage list for TFs
l_TFs_genomatix <- list()
for(i in 1:length(Genomatix_list)){
  Genomatix_i <- Genomatix_list[[i]]
  Genomatix_thresh_i <- subset(Genomatix_i, Z.Score..promoters. >= Z_score_overrepresentation)
  TF_i <- as.character(Genomatix_thresh_i$TF.Matrices)
  TF_i <- unlist(lapply(strsplit(TF_i, "\\$"), function(x){x[[2]]}))
  TF_i <- unlist(lapply(strsplit(TF_i, "\\."), function(x){x[[1]]}))
  TF_i <- unlist(lapply(strsplit(TF_i, "\\_"), function(x){x[[1]]}))
  TF_i <- unlist(lapply(strsplit(TF_i, "\\-"), function(x){x[[1]]}))
  TF_i <- tolower(TF_i)
  TF_i <- firstup(TF_i)
  TF_i <- unique(TF_i)
  TF_i <- TF_i[order(TF_i)]
  l_TFs_genomatix[[i]] <- TF_i
}

names(l_TFs_genomatix) <- Genomatix_dataset_names
#Save vectors as tables
#Certain TF GeneSymbols need to be curated (e.g. "-" in Nkx2-3)
#for(i in 1:length(l_TFs_genomatix)){
#  save_it_i <- l_TFs_genomatix[[i]]
#  name_i <- names(l_TFs_genomatix)[i]
#  write.table(save_it_i,paste(PATH_output,"/Compiled_TFs_",name_i,".txt",sep=""),
#              sep="\t",quote = FALSE,row.names = FALSE, col.names = FALSE)
#}

#Check overlap of identified TFs
# All
Genomatix_hyper_TFs_all <- read.delim(paste(PATH_output, "/Compiled_TFs_Genomatix_hyper_all.txt",sep=""), col.names = "GeneSymbol")
Genomatix_hypo_TFs_all <- read.delim(paste(PATH_output, "/Compiled_TFs_Genomatix_hypo_all.txt",sep=""), col.names = "GeneSymbol")

Genomatix_hyper_TFs_all <- Genomatix_hyper_TFs_all[!duplicated(Genomatix_hyper_TFs_all$GeneSymbol),]
Genomatix_hypo_TFs_all <- Genomatix_hypo_TFs_all[!duplicated(Genomatix_hypo_TFs_all$GeneSymbol),]

write.table(Genomatix_hyper_TFs_all,paste(PATH_output, "/Compiled_TFs_Genomatix_hyper_all_nondup.txt",sep=""), sep = "\t")
write.table(Genomatix_hypo_TFs_all,paste(PATH_output, "/Compiled_TFs_Genomatix_hypo_all_nondup.txt",sep=""), sep = "\t")

#Thresh
Genomatix_hyper_TFs_thresh <- read.delim(paste(PATH_output, "/Compiled_TFs_Genomatix_hyper_thresh.txt",sep=""), col.names = "GeneSymbol")
Genomatix_hypo_TFs_thresh <- read.delim(paste(PATH_output, "/Compiled_TFs_Genomatix_hypo_thresh.txt",sep=""), col.names = "GeneSymbol")

Genomatix_hyper_TFs_thresh <- Genomatix_hyper_TFs_thresh[!duplicated(Genomatix_hyper_TFs_thresh$GeneSymbol),]
Genomatix_hypo_TFs_thresh <- Genomatix_hypo_TFs_thresh[!duplicated(Genomatix_hypo_TFs_thresh$GeneSymbol),]

write.table(Genomatix_hyper_TFs_thresh,paste(PATH_output, "/Compiled_TFs_Genomatix_hyper_thresh_nondup.txt",sep=""), sep = "\t")
write.table(Genomatix_hypo_TFs_thresh,paste(PATH_output, "/Compiled_TFs_Genomatix_hypo_thresh_nondup.txt",sep=""), sep = "\t")


Overlap_TFs_Meth <- Genomatix_hyper_TFs_thresh[Genomatix_hyper_TFs_all %in% Genomatix_hypo_TFs_all]

length(Overlap_TFs_Meth)
length(Genomatix_hyper_TFs_all) - length(Overlap_TFs_Meth)
length(Genomatix_hypo_TFs_all) - length(Overlap_TFs_Meth)

