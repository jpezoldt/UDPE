#####
# This script contains all user specific path definitions and required scripts which have to be 
# loaded. 
# User specified part has to be modified to define the path of the home folder of the program.
# It has to be run before all other scripts and adjusted according to the folder structure.
#####
### USER SPECIFIED!!! ####
# Path to home folder containing Data, Filelists, etc subfolders
homepath = "/home/pezoldt/SVFASRAW/mminder/Classifier 2.0/"

#####
### DEPENDENCIES ####
### paths ####
datapath = paste(homepath, "Datasets/", sep="")

### libraries ####
library(glue)

### functions ####
source(paste(datapath, "_general/fulldata.R",sep=""))
source(paste(datapath, "_general/load_util.R",sep=""))
source(paste(datapath, "_general/TPM_util.R",sep=""))
source(paste(homepath, "Datasets/_general/fulldata.R", sep=""))
source(paste(homepath, "Training/functions/testing/test_stats.R",sep=""))
source(paste(homepath, "Training/functions/training/classifier.R",sep=""))
source(paste(homepath, "Training/functions/preliminary/feature_extraction.R",sep=""))
source(paste(homepath, "Training/functions/application/prediction_analysis.R",sep=""))

source(paste(homepath, "generic_methods.R",sep=""))
