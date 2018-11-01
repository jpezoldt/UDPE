## dissect compiled and curated Output from Homer
# Author: Joern Pezoldt
# Date: 28.10.2018
# Type of script: Exploratory

#Libraries
library("UpSetR")


#load table
curated_TF_ATAC <- read.delim("C:/Users/pezoldt/Desktop/Transfer/SPF_known_homer_compilation_curated.txt")



#####
#Global & PATHs
#####
#table_ID <- "SPF_known_homer_compilation_curated.txt"
#sample_ID <- "SPF"
#PATH_input <- paste("/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/Motif/Homer_Output/allHomerResults/size_500bp/", sample_ID, sep = "")

#table_ID <- "SPF_known_homer_compilation_curated.txt"

### Check for Unique TFs ---------------------------------------
# load
#t_homer_TF <- read.table(paste(PATH_input,"/",table_ID, sep=""), sep="\t",dec=".",header=TRUE)
t_homer_TF <- curated_TF_ATAC

# split according to comparisons
l_homer_TF <- split(t_homer_TF, t_homer_TF$group)
#rename
l_names <- strsplit(names(l_homer_TF), "noExt_")
names(l_homer_TF) <- unlist(sapply(l_names, function(x) {x[2]}))
#Delete global DAR sets
l_homer_TF <- l_homer_TF[c("Closed_mLN_DAR_SPF","Open_mLN_DAR_SPF")]

# grab TFs per comparison
l_motifs_per_group <- lapply(l_homer_TF, function(x){
  #x <- l_homer_TF[[1]]                              
  TFs_i <- as.character(unlist(strsplit(as.character(x$Motifs_TF), ";")))
})
#Generate binary Matrix for upset()
t_TFs_binary <- as.data.frame.matrix(table(stack(l_motifs_per_group)))
#Replace 
t_TFs_binary[t_TFs_binary >= 2] <- 1

# Plot intersection
upset(t_TFs_binary, sets = colnames(t_TFs_binary),
      mainbar.y.label = "TFs ATACseq",
      number.angles = 30, point.size = 5,
      text.scale = 2)

# Grab intersection of interest
# Note: Only one intersection interesting
Unique_TFs_Open <- rownames(subset(t_TFs_binary, 
                                  Closed_mLN_DAR_SPF == 0 &
                                  Open_mLN_DAR_SPF == 1))
Unique_TFs_Closed <- rownames(subset(t_TFs_binary, 
                                Closed_mLN_DAR_SPF == 1 &
                                  Open_mLN_DAR_SPF == 0))
Common_TFs_Open_Closed <- rownames(subset(t_TFs_binary, 
                                       Closed_mLN_DAR_SPF == 1 &
                                         Open_mLN_DAR_SPF == 1))
write.table(Unique_TFs_Closed, "C:/Users/pezoldt/Desktop/Transfer/Unique_TFs_Closed.txt",
            row.names = FALSE,quote = FALSE, col.names = FALSE)
write.table(Unique_TFs_Open, "C:/Users/pezoldt/Desktop/Transfer/Unique_TFs_Open.txt",
            row.names = FALSE,quote = FALSE, col.names = FALSE)
write.table(Common_TFs_Open_Closed, "C:/Users/pezoldt/Desktop/Transfer/Common_TFs_Open_Closed.txt",
            row.names = FALSE,quote = FALSE, col.names = FALSE)

# Write binary table
write.table(t_TFs_binary, paste(PATH_input,"/",sample_ID,"_existence_matrix_known_homer_compile.txt",sep = "") ,dec=".", sep="\t")





