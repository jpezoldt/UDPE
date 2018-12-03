# dissect compiled and curated Output from Homer
# Author: Joern Pezoldt
# Date: 12.10.2018
# Type of script: Exploratory

#Libraries
library("UpSetR")
library("pheatmap")

#####
#Global & PATHs
#####
table_ID <- "SPF_known_homer_compilation_curated.txt"
sample_ID <- "SPF"
PATH_input <- paste("/Volumes/Pezoldt_HD_4TB/UPDE_Transfer/20181101/Analysis/ATACseq/ATAC_FSC_all/Motif/Homer_Output/allHomerResults/size_500bp/", sample_ID, sep = "")
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
#Get TF and pValue
l_homer_TF_pValue <- list()
for(i in 1:length(l_homer_TF)){
  name_i <- names(l_homer_TF)[i]
  print(name_i)
  t_i <- l_homer_TF[[i]][,c("Motifs_TF","P.value")]
  colnames(t_i) <- c("Motifs_TF", paste("pValue_",name_i,sep=""))
  l_homer_TF_pValue[[i]] <- t_i
  names(l_homer_TF_pValue)[i] <- name_i
}

# Merge tables by Motifs
t_homer_TF_pValue <- Reduce(function(...) merge(..., all=T), l_homer_TF_pValue)

# Condense per TF by meaning accross conditions
l_homer_TF_pValue <- split(t_homer_TF_pValue, t_homer_TF_pValue$Motifs_TF)
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
pheatmap(t_homer_TF_pValue_mean, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = FALSE,
         scale = "row", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("white", "yellow","darkgreen"), space="rgb")(128),
         main = "TF by motif -log10(pValue)")

