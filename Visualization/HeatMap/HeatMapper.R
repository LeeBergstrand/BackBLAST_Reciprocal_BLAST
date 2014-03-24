#!/usr/bin/env Rscript 
# Created by: Lee Bergstrand & Erick Cardenas Poire 
# Descript: Converts BLAST TotalBLASTResultss from BackBlast into a heatmap.  
#             
# Requirements: - reshape2 and RColorBrewer modules
#----------------------------------------------------------------------------------------

# Sets the working directory. 
setwd("/Users/lee/Dropbox/R & D/Repositories/BackBLAST-Gene-Cluster-Finder/Visualization/HeatMap/TestData")

# Gets a list files from the working directory.
fileList = list.files(path = getwd(), all.files = FALSE, pattern = "\\.csv$") 

TotalBLASTResults = matrix(data = 0, nrow = 0, ncol = 7)

fileCounter = 1
while (fileCounter <= length(fileList)){
  BackBLASTResults = read.csv(fileList[fileCounter], header = FALSE, stringsAsFactors = TRUE) # Inport data from csv file from BackBLAST.py
  SubjectAccession = sub(".csv$", "", fileList[fileCounter]) # Remove .csv suffix from csv filename. This should be the subject organism accession.
  SubjectAccessionColumn = matrix(data = SubjectAccession, nrow = nrow(BackBLASTResults), ncol = 1) # Makes a column that contains only the organism accession as data.
  BackBLASTResults  = cbind(BackBLASTResults, SubjectAccessionColumn) # Concatenates SubjectAccessionColumn to the BackBLAST results.
  TotalBLASTResults = rbind(TotalBLASTResults, BackBLASTResults) # Concatenates the current csv files BLAST results to the total BLAST results.
  fileCounter = fileCounter + 1
}

colnames(TotalBLASTResults) = c("QuerySeqID", "SubjectSeqID", "PercentIdent", "Evalue", "QueryCoverage", "Bitscore", "TargetOrganism")

library(reshape2)
library(RColorBrewer)

#check this function help to see what to do when multiple values occur
HeatmapMatrix = acast(TotalBLASTResults, QuerySeqID ~ TargetOrganism, value.var = "PercentIdent") # Converts TotalBLASTResults to wide format data matrix.
HeatmapMatrix[is.na(HeatmapMatrix)] = 0 # Replaces NA values with zero for heatmap function.

heatmap(HeatmapMatrix, col = brewer.pal(8,"Blues"), breaks = c(0,30,40,50,60,70,80,90,100), Rowv = NA, Colv = NA)