#####
### INFO ####
# This is a very basic example of a complete workflow to obtain stem predictions
# from raw data. The example data is a randomly subsampled subset of the 
# Gokce2016 data (GSE82187). This data is present in a transformed format,
# values are the log10 of the transcript per million values plus 1, e.g. 
# log10(TPM+1)
# Before running the script, the userspecs.R script has to be adapted to the 
# specific user and be run.
### HEADER ####
library(glue)

# Define path where the example files are at.
path = glue("{homepath}Examples/")

### CREATING THE FULLDATA OBJECT ####
# At first, we'll create a fulldata object containing the data in a format which
# is compatible with downstream functions. 

# Reading the raw data
rawdata = read.csv(file = glue("{path}01_data.csv"),
                   stringsAsFactors = FALSE, check.names=FALSE)

# The first column contains the gene names
rownames(rawdata) = rawdata[ ,1]
rawdata = rawdata[ ,-1]

# For the fulldata object, genes have to be columns and cells rows, to fit 
# the standard for machine learning techniques
rawdata = data.frame(t(rawdata))

# Reading celltype annotation, adding celltype. Note that this isn't absolutely
# necessary to generate predictions
anno = read.csv(file=glue("{path}01_annotation.csv"),
                stringsAsFactors = FALSE)
celltype = anno$celltype[match(rownames(rawdata), anno$cellname)]

# Creating custom get_TPM function, since the input data doesn't contain counts 
# but log10(TPM + 1) values. Otherwise, if the input data is in the form of a 
# count matrix, the value of get_TPM could just be set to NULL. 
custom_get_TPM = function(x, cellix = T, geneix = T) {
  raw = x$raw
  # Transform back from logarithm
  retransformed = 10^raw[cellix, geneix] - 1
  # Renormalize by a million to correct rounding errors, return
  return(norm_mio_transcr(retransformed))
}

# Create fulldata object
fd = fulldata(raw = rawdata,
              celltype = celltype,
              get_TPM = custom_get_TPM) 


### CREATE PREDICTIONS ####
# Having created the fulldata object, we can use it to generate predictions,
# using an example classifier. 

# Load classifier
classif = readRDS(glue("{path}01_classif.RDS"))

# Prepare data to have format for prediction input
prepped = prepare(classif, fd)

# Predict stemness probabilities
preds = predict(classif, prepped)

# The object preds now contains stemness predictions for each method and each 
# cell, which can be used for further analysis.  


### EXPLORE RESULTS ####
# The function shinyRainPlots generates a Shiny interface, which allows to 
# explore the predictions.
shinyRainPlots(preds, celltype = fd$celltype)









