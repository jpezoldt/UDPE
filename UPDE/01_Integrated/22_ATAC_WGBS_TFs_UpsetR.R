# Author: Joern Pezoldt
# Date: 31.11.2018
# Function: 
#    1) Intersect Overlap of TFs identified in WGBS and ATACseq data
#    2) Plot in upsetR

#####
#PATHs & Globals
#####
PATH_WGBS_Genomatix <- "/home/pezoldt/NAS2/pezoldt/Analysis/WGBS/FSC/Genomatix_TF/ANNa/Mouse/To_Use/Output"
PATH_ATAC_Homer <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/Motif/Homer_Output/allHomerResults/size_500bp/SPF/Only_ATAC_Data"
  
  
#####
#Load
#####
# WGBS TFs
Genomatix_hyper_TFs_all <- read.delim(paste(PATH_WGBS_Genomatix, "/Compiled_TFs_Genomatix_hyper_all_nondup.txt",sep=""))
Genomatix_hyper_TFs_all <- as.character(Genomatix_hyper_TFs_all$x)
Genomatix_hypo_TFs_all <- read.delim(paste(PATH_WGBS_Genomatix, "/Compiled_TFs_Genomatix_hypo_all_nondup.txt",sep=""))
Genomatix_hypo_TFs_all <- as.character(Genomatix_hypo_TFs_all$x)

# ATAC TFs
Homer_ATAC_TFs <- read.delim(paste(PATH_ATAC_Homer, "/SPF_existence_matrix_known_homer_compile.txt",sep=""))
Homer_ATAC_open_mLN <- rownames(subset(Homer_ATAC_TFs, 
                                Open_mLN_SPF == 1))
Homer_ATAC_open_pLN <- rownames(subset(Homer_ATAC_TFs, 
                              Open_pLN_SPF == 1))

#Delete global DAR sets
l_TFs_WGBS_ATAC <- list(Genomatix_hyper_TFs_all,
                           Genomatix_hypo_TFs_all,
                           Homer_ATAC_open_mLN,
                           Homer_ATAC_open_pLN)
names(l_TFs_WGBS_ATAC) <- c("demeth_pLN",
                            "demeth_mLN",
                            "open_mLN",
                            "open_pLN")


#Generate binary Matrix for upset()
t_TFs_binary <- as.data.frame.matrix(table(stack(l_TFs_WGBS_ATAC)))
#Replace 
t_TFs_binary[t_TFs_binary >= 2] <- 1

upset(t_TFs_binary, sets = colnames(t_TFs_binary),
      mainbar.y.label = "TFs ATACseq",
      number.angles = 30, point.size = 5,
      text.scale = 2, keep.order = TRUE)
# Check content
Unique_TFs_NoDEG_Open <- subset(t_TFs_binary, 
                                demeth_pLN == 1 &
                                demeth_mLN == 1 &
                                open_mLN == 0 &
                                open_pLN == 0)

# Check content demeth
Unique_TFs_NoDEG_Open <- subset(t_TFs_binary, 
                                demeth_pLN == 1 &
                                  demeth_mLN == 0)

# Check content demeth
Unique_TFs_NoDEG_Open <- subset(t_TFs_binary, 
                                open_mLN == 1 &
                                  open_pLN == 1)
