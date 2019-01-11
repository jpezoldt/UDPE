#!/usr/bin/env Rscript
# FILE: runJoernRun.R ---------------------------------------------------------
#
# DESCRIPTION : Runs dynGENIE3 on expression matrices 
#
# USAGE: 
#
# OPTIONS:  none
# REQUIREMENTS:  
# BUGS: --
# NOTES:  ---
# AUTHOR:  Maria Litovchenko, maria.litovchenko@epfl.ch
# COMPANY:  EPFL, Lausanne, Switzerland
# VERSION:  1
# CREATED:  03.12.2018
# REVISION: 03.12.2018

library(biomaRt)
library(data.table)
library(GENIE3)
source("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/dynGENIE3_R_C_wrapper_v2/dynGENIE3.R")

#args = commandArgs(trailingOnly = T)
#message(args[1])

#####
#PATH
#####
PATH_input_GENIE3 <- "/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/dynGENIE3_R_C_wrapper_v2"
PATH_input_averagedExpression <- "/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/monocle/day0_10_24_56_300/Seurat_selected"
PATH_output <- "/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/dynGENIE3/Seurat_selected_Split_Adv_NonAdv"

# STEP 1: read in the data ----------------------------------------------------
TSexprMatr1 <- read.table(paste(PATH_input_averagedExpression,"/averaged_nonadventi_dynGENIE3.txt",sep=""),
                          header = T, row.names = 1, stringsAsFactor = F)
#TSexprMatr2 <- read.table("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/dynGENIE3_R_C_wrapper_v2/averaged_adventi_dynGENIE3.txt", header = T, row.names = 1, 
 #                         stringsAsFactor = F)


TSexprMatr1 <- TSexprMatr1[apply(TSexprMatr1[,-1], 1, function(x) !all(x==0)),]
#TSexprMatr2 <- TSexprMatr2[apply(TSexprMatr2[,-1], 1, function(x) !all(x==0)),]

#Majority of genes is very lowly expressed log2(average_exp_per_timePoint) < -5
# Thresh for that to reduce noise, average expression > 0.05 
#TSexprMatr1_thresh <- subset(TSexprMatr1, rownames(TSexprMatr1) %in% rownames(TSexprMatr1[rowMeans(TSexprMatr1) >= 0.05,]))
#TSexprMatr2_thresh <- subset(TSexprMatr2, rownames(TSexprMatr2) %in% rownames(TSexprMatr1[rowMeans(TSexprMatr2) >= 0.05,]))
thresh_expression <- 0.05
TSexprMatr1_thresh <- subset(TSexprMatr1, D000 >= thresh_expression |
                               D010 >= thresh_expression |
                               D024 >= thresh_expression |
                               D056 >= thresh_expression |
                               D300 >= thresh_expression)




l_TSexprMatr <- list(as.matrix(TSexprMatr1_thresh))

TSpoints1 <- as.integer(gsub('D', '', colnames(TSexprMatr1)))
#TSpoints2 <- as.integer(gsub('d', '', colnames(TSexprMatr2)))
l_TSpoints <- list(TSpoints1)
#,TSpoints2)

# STEP 2: get all TFs expressed in data set -----------------------------------
# I would suggest you to get all TFs present in dataset to see 
# without bias which are important

# get all TFs from Biomart
while(!exists('ensembl91')) {
  ensembl91 <- useEnsembl("ensembl", version = 91,
                          dataset = "mmusculus_gene_ensembl")
}
egs <- getBM(attributes = c('name_1006', 'external_gene_name'), values = "*",
             mart = ensembl91)
tfGOterms <- paste('transcription, DNA-templated',
                   'transcription factor activity',
                   'transcription factor.*complex',
                   'sequence-specific DNA binding',
                   'transcriptional repressor activity', sep = '|')
egs <- egs[grepl(tfGOterms, egs$name_1006), ]
egs <- egs[!duplicated(egs$external_gene_name), ] 

TFsInTS <- intersect(rownames(TSexprMatr1_thresh), unique(egs$external_gene_name))

# STEP 3: run dynGENIE3 -------------------------------------------------------
setwd(PATH_input_GENIE3)
message(paste('Started simpliest run of dynGENIE3 for', args[1],
              'at', Sys.time()))
res <- dynGENIE3(l_TSexprMatr, l_TSpoints, ncores = 8)


saveRDS(res, paste(PATH_output,"/dynGENIE3_Expr_thresh_",thresh_expression,"_inOne.Rds",sep=""))

#message(paste('\tFinished simpliest run for', args[1], 'at', Sys.time()))

#message(paste('Started run with regulators of dynGENIE3 for', args[1], 
 #             'at', Sys.time()))
res <- dynGENIE3(l_TSexprMatr, l_TSpoints, regulators = TFsInTS, ncores = 8)
saveRDS(res, paste(PATH_output,"/dynGENIE3_Expr_thresh_",thresh_expression,"_inOne_withRegs.Rds",sep=""))
message(paste('\tFinished run with regulators for', args[1], 'at', Sys.time()))

message('DONE!')
