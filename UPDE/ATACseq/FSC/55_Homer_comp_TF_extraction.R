# dissect compiled and curated Output from Homer
# Author: Joern Pezoldt
# Date: 12.10.2018
# Type of script: Exploratory

#Libraries
library("UpSetR")

#####
#Global & PATHs
#####
table_ID <- "SPF_known_homer_compilation_curated.txt"
sample_ID <- "SPF"
PATH_input <- paste("/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/Motif/Homer_Output/allHomerResults/size_500bp/", sample_ID, sep = "")
table_ID <- "SPF_known_homer_compilation_curated.txt"

### Check for Unique TFs ---------------------------------------
# load
t_homer_TF <- read.table(paste(PATH_input,"/",table_ID, sep=""), sep="\t",dec=".",header=TRUE)


# split according to comparisons
l_homer_TF <- split(t_homer_TF, t_homer_TF$group)
#rename
l_names <- strsplit(names(l_homer_TF), "noExt_")
names(l_homer_TF) <- unlist(sapply(l_names, function(x) {x[2]}))
#Delete global DAR sets
l_homer_TF <- l_homer_TF[c("Closed_DEG_SPF_lowRNAseq","No_DAR_DOWN_SPF_lowRNAseq","No_DAR_UP_SPF_lowRNAseq",
                "NoDEG_Closed_SPF_lowRNAseq","NoDEG_Open_SPF_lowRNAseq","Open_DEG_SPF_lowRNAseq")]
names(l_homer_TF) <- c("pLN_Open_UP","pLN_peak_UP","mLN_peak_UP",
                          "pLN_Open_None","mLN_Open_None","mLN_Open_UP")

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
#ToDo:
# Generate Heatmap of all TFs across conditions plotting the -log10(pvalue)



