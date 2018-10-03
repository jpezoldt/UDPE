# FILE DESCRIPTION: SCENIC ----------------------------------------------------
#
# DESCRIPTION : runs SCENIC on evolved and not evolved populations
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
# CREATED:  07.08.2018
# REVISION: 07.08.2018

setwd('~/Desktop/BitBucket/BerraErkosar/')
setwd('~/Desktop/BerraErkosar/')
setwd('/home/litovche/BerraErkosar/')

# Installation ----------------------------------------------------------------

#install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/GENIE3_1.2.1.tar.gz", repos=NULL)
#install.packages(c('feather', 'R.utils', 'mixtools'))
#source("https://bioconductor.org/biocLite.R")
#biocLite('GSEABase')
#install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/AUCell_1.2.4.tar.gz", repos=NULL)
#install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/RcisTarget_1.0.2.tar.gz", repos=NULL)
#install.packages("devtools")
#devtools::install_github("aertslab/SCENIC")
# dbFiles <- c("https://resources.aertslab.org/cistarget/databases/drosophila_melanogaster/dm6/flybase_r6.02/mc8nr/gene_based/dm6-5kb-upstream-full-tx-11species.mc8nr.feather")
# download.file(dbFiles, 'data/dm6-5kb-upstream-full-tx-11species.mc8nr.feather')

# Loading of the libraries ----------------------------------------------------
library(biomaRt)
library(data.table)
library(AUCell)
library(GENIE3)
library(RcisTarget)
library(SCENIC)

set.seed(123) # For reproducibility of results

# Functions -------------------------------------------------------------------
runSCENIC_all <- function(countsMat, scenicParams) {
  # STEP 1: Co-expression network 
  # The first step workflow is to infer potential transcription factor targets
  # based on the expression data. To do this use GENIE3. The input are the 
  # expression matrix (filtered), and a list of transcription factors (potential
  # regulators). The output of GENIE3 and a correlation matrix will be used to 
  # create the co-expression modules.
  
  # Gene filter/selection
  countsPerGene <- apply(countsMat, 1, sum)
  summary(countsPerGene)
  # I assume that I have log2 counts, so I will restrict genes to the genes 
  # expressed in more than 6 samples (it removes really little genes)
  countsMat <- countsMat[countsPerGene > 6, ]
  
  # In upcoming steps only the genes that are available on RcisTarget databases 
  # will be used. To save some running time for GENIE3, I ignore the genes that 
  # are not in the databases
  motifRankings <- importRankings(getDatabases(scenicParams)[[1]])
  genesInDatabase <- colnames(getRanking(motifRankings))
  inDB <- which(rownames(countsMat) %in% genesInDatabase)
  message(paste('Genes in DB:', length(inDB)))
  # restrict to genes in DB
  countsMat <- countsMat[inDB, ]
  countsMat <- as.matrix(countsMat)
  
  # GENIE3 can detect both positive and negative associations. In order to 
  # distinguish potential activation from repression, we will split the targets
  # into positive- and negative-correlated targets
  corrMat <- cor(t(countsMat), method = "spearman")
  saveRDS(corrMat, file = getIntName(scenicParams, "corrMat"))
  
  # Step 2 : GENIE3
  # GENIE3 is based on a Random Forest approach, each time it is run the results 
  # will be slightly different. The higher the number of trees used (ntrees), 
  # the lower the variability. 
  runGenie3(countsMat, scenicParams)
  
  # Step 3: SCENIC
  runSCENIC_1_coexNetwork2modules(scenicParams)
  runSCENIC_2_createRegulons(scenicParams)
  runSCENIC_3_scoreCells(scenicParams, countsMat, skipTsne = T)
  
  # Step 4: Binarize the network activity 
  aucellApp <- plotTsne_AUCellApp(scenicParams, countsMat)
  savedSelections <- shiny::runApp(aucellApp)
  
  # Save the modified thresholds:
  newThresholds <- savedSelections$thresholds
  scenicParams@fileNames$int["aucell_thresholds", 1] <- "int/newThresholds.Rds"
  saveRDS(newThresholds, file = getIntName(scenicParams, "aucell_thresholds"))
  saveRDS(scenicParams, file="int/scenicParams.Rds") 
  
  # scenicParams@settings$devType="png"
  runSCENIC_4_aucell_binarize(scenicParams, skipTsne = T)
}

# Inputs ----------------------------------------------------------------------
# counts - normalized 
gutsCounts <- read.table('data/gutsKallistoNcounts.csv', header = T, sep = ',',
                         row.names = 1, stringsAsFactors = F)
larvaeCounts <- read.table('data/gutsKallistoNcounts.csv', header = T, 
                           sep = ',', row.names = 1, stringsAsFactors = F)

# info about samples
gutsINFO <- data.frame(Tissue = 'GUT',
                       Population = gsub('Aceto.|GF.', '', 
                                         colnames(gutsCounts)),
                       Condition = gsub('CTL|SEL', '', 
                                        colnames(gutsCounts)),
                       Rep = gsub('\\D', '', colnames(gutsCounts)))
gutsINFO$Condition <- gsub('\\d', '', gutsINFO$Condition)
gutsINFO$Rep <- as.integer(gutsINFO$Rep)
rownames(gutsINFO) <- colnames(gutsCounts)
larvaeINFO <- data.frame(Tissue = 'Larvae',
                         Population = gsub('Aceto.|GF.', '', 
                                           colnames(gutsCounts)),
                         Condition = gsub('CTL|SEL', '',
                                          colnames(gutsCounts)),
                         Rep = gsub('\\D', '', colnames(gutsCounts)))
larvaeINFO$Condition <- gsub('\\d', '', larvaeINFO$Condition)
larvaeINFO$Rep <- as.integer(larvaeINFO$Rep)
rownames(larvaeINFO) <- colnames(gutsCounts)

# merge together larvae and guts
colnames(gutsCounts) <- paste0('GUT_', colnames(gutsCounts))
colnames(larvaeCounts) <- paste0('Larvae_', colnames(larvaeCounts))
allCounts <-  cbind(gutsCounts, larvaeCounts)
INFO <- rbind(gutsINFO, larvaeINFO)
rownames(INFO) <- colnames(allCounts)

ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "dmelanogaster_gene_ensembl", 
                   host = "jul2015.archive.ensembl.org")
geneNames <- getBM(attributes = c('flybase_gene_id', 'external_gene_name'),
                   filters = 'flybase_gene_id', values = rownames(allCounts),
                   mart = ensembl)
geneNames <- as.data.table(geneNames)
setkey(geneNames, flybase_gene_id)
rownames(allCounts) <- geneNames[rownames(allCounts)]$external_gene_name

# Run SCENIC on larvae, control population, GF --------------------------------
# create directory
dir.create('/home/litovche/BerraErkosar/scenic_larvae_ctrl_GF')
setwd('/home/litovche/BerraErkosar/scenic_larvae_ctrl_GF')
# parameters for SCENIC
scenicOpt <- initializeScenic(org = "dmel", nCores = 8, dbDir = "../data/", 
                              datasetTitle = "SCENIC larvae ctrl GF")
runSCENIC_all(allCounts[, grep('Larvae_CTLGF', colnames(allCounts))], scenicOpt)

# Run SCENIC on larvae, control population, Aceto -----------------------------
# create directory
dir.create('/home/litovche/BerraErkosar/scenic_larvae_ctrl_aceto')
setwd('/home/litovche/BerraErkosar/scenic_larvae_ctrl_aceto')
# parameters for SCENIC
scenicOpt <- initializeScenic(org = "dmel", nCores = 8, dbDir = "../data/", 
                              datasetTitle = "SCENIC larvae ctrl aceto")
runSCENIC_all(allCounts[, grep('Larvae_CTLAceto', colnames(allCounts))],
              scenicOpt)

# Run SCENIC on larvae, evolved population, GF --------------------------------
# create directory
dir.create('/home/litovche/BerraErkosar/scenic_larvae_sel_GF')
setwd('/home/litovche/BerraErkosar/scenic_larvae_sel_GF')
# parameters for SCENIC
scenicOpt <- initializeScenic(org = "dmel", nCores = 8, dbDir = "../data/", 
                              datasetTitle = "SCENIC larvae sel GF")
runSCENIC_all(allCounts[, grep('Larvae_SELGF', colnames(allCounts))],
              scenicOpt)

# Run SCENIC on larvae, evolved population, Aceto -----------------------------
# create directory
dir.create('/home/litovche/BerraErkosar/scenic_larvae_sel_aceto')
setwd('/home/litovche/BerraErkosar/scenic_larvae_sel_aceto')
# parameters for SCENIC
scenicOpt <- initializeScenic(org = "dmel", nCores = 8, dbDir = "../data/", 
                              datasetTitle = "SCENIC larvae sel aceto")
runSCENIC_all(allCounts[, grep('Larvae_SELAceto', colnames(allCounts))], 
              scenicOpt)

# Run SCENIC on gut, control population, GF --------------------------------
# create directory
dir.create('/home/litovche/BerraErkosar/scenic_gut_ctrl_GF')
setwd('/home/litovche/BerraErkosar/scenic_gut_ctrl_GF')
# parameters for SCENIC
scenicOpt <- initializeScenic(org = "dmel", nCores = 8, dbDir = "../data/", 
                              datasetTitle = "SCENIC gut ctrl GF")
runSCENIC_all(allCounts[, grep('GUT_CTLGF', colnames(allCounts))], scenicOpt)

# Run SCENIC on gut, control population, Aceto -----------------------------
# create directory
dir.create('/home/litovche/BerraErkosar/scenic_gut_ctrl_aceto')
setwd('/home/litovche/BerraErkosar/scenic_gut_ctrl_aceto')
# parameters for SCENIC
scenicOpt <- initializeScenic(org = "dmel", nCores = 8, dbDir = "../data/", 
                              datasetTitle = "SCENIC gut ctrl aceto")
runSCENIC_all(allCounts[, grep('GUT_CTLAceto', colnames(allCounts))],
              scenicOpt)

# Run SCENIC on gut, evolved population, GF --------------------------------
# create directory
dir.create('/home/litovche/BerraErkosar/scenic_gut_sel_GF')
setwd('/home/litovche/BerraErkosar/scenic_gut_sel_GF')
# parameters for SCENIC
scenicOpt <- initializeScenic(org = "dmel", nCores = 8, dbDir = "../data/", 
                              datasetTitle = "SCENIC gut sel GF")
runSCENIC_all(allCounts[, grep('GUT_SELGF', colnames(allCounts))],
              scenicOpt)

# Run SCENIC on gut, evolved population, Aceto -----------------------------
# create directory
dir.create('/home/litovche/BerraErkosar/scenic_gut_sel_aceto')
setwd('/home/litovche/BerraErkosar/scenic_gut_sel_aceto')
# parameters for SCENIC
scenicOpt <- initializeScenic(org = "dmel", nCores = 8, dbDir = "../data/", 
                              datasetTitle = "SCENIC gut sel aceto")
runSCENIC_all(allCounts[, grep('GUT_SELAceto', colnames(allCounts))], 
              scenicOpt)