#!/usr/bin/env Rscript 
# Created by: Lee Bergstrand 
# Descript: Converts BLAST TotalBLASTResultss from BackBlast into a heatmap.  
#             
# Requirements: - reshape2 and RColorBrewer modules
#----------------------------------------------------------------------------------------
library(reshape2)
library(RColorBrewer)
library(ape)
library(gplots)

# Sets the working directory. 
setwd("/Users/lee/Data/SecondBDHAnalysis/InDraft/Results/Filtered")

# Gets a list files from the working directory.
fileList = list.files(path = getwd(), all.files = FALSE, pattern = "\\.csv$") 

TotalBLASTResults = matrix(data = 0, nrow = 0, ncol = 7)

fileCounter = 1
while (fileCounter <= length(fileList)){
  BackBLASTResults = read.csv(fileList[fileCounter], header = FALSE, stringsAsFactors = TRUE) # Inport data from csv file from BackBLAST.py
  SubjectAccession = sub(".csv$", "", fileList[fileCounter]) # Remove .csv suffix from csv filename. This should be the subject organism accession.
  SubjectAccessionColumn = matrix(data = SubjectAccession, nrow = nrow(BackBLASTResults), ncol = 1) # Makes a column that contains only the organism accession as data.
  BackBLASTResults  = cbind(BackBLASTResults, SubjectAccessionColumn) # Concatenates SubjectAccessionColumn to the BackBLAST results.
  BackBLASTResults[[2]] = as.character(BackBLASTResults[[2]]) # Some SubjectSeqIDs in BackBLAST results are purly numeric (example from JGI id numbers). 
  # This throws a warning when you attempt to append this numeric column in BackBLASTResults 
  # to a character column in TotalBLASTResults. The code left typecasts this coloumn to character.
  TotalBLASTResults = rbind(TotalBLASTResults, BackBLASTResults) # Concatenates the current csv files BLAST results to the total BLAST results.
  fileCounter = fileCounter + 1
}
colnames(TotalBLASTResults) = c("QuerySeqID", "SubjectSeqID", "PercentIdent", "Evalue", "QueryCoverage", "Bitscore", "TargetOrganism")

#check this function help to see what to do when multiple values occur
HeatmapMatrix = acast(TotalBLASTResults, QuerySeqID ~ TargetOrganism, value.var = "PercentIdent") # Converts TotalBLASTResults to wide format data matrix.
HeatmapMatrix[is.na(HeatmapMatrix)] = 0 # Replaces NA values with zero for heatmap function.

HeatmapMatrix = HeatmapMatrix[order(rownames(HeatmapMatrix)), ] # Reorders rows by gene name.

ColourPal = brewer.pal(9,"YlGn") # Gets Colour Palete from R colour brewer.
ColourPal[1] = "#F4F5F6" # Swaps lowest colour for off white.
ColourPal = append(ColourPal, "#00311d")
heatmap.2(HeatmapMatrix, Rowv = NA, col = ColourPal, 
          trace="none", xlab = "Genome", ylab = "Steroid Degrading Gene", margins = c(10,9))
