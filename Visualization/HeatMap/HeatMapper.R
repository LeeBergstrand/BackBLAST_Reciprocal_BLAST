#!/usr/bin/env Rscript 
# Created by: Lee Bergstrand & Erick Cardenas Poire 
# Descript: Converts BLAST TotalBLASTResultss from BackBlast into a heatmap.  
#             
# Requirements: - reshape2 and RColorBrewer modules
#----------------------------------------------------------------------------------------
options(warn = 2) 
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
  BackBLASTResults[[2]] = as.character(BackBLASTResults[[2]]) # Some SubjectSeqIDs in BackBLAST results are purly numeric (example from JGI id numbers). 
                                                              # This throws a warning when you attempt to append this numeric column in BackBLASTResults 
                                                              # to a character column in TotalBLASTResults. The code left typecasts this coloumn to character.
  TotalBLASTResults = rbind(TotalBLASTResults, BackBLASTResults) # Concatenates the current csv files BLAST results to the total BLAST results.
  fileCounter = fileCounter + 1
}
colnames(TotalBLASTResults) = c("QuerySeqID", "SubjectSeqID", "PercentIdent", "Evalue", "QueryCoverage", "Bitscore", "TargetOrganism")

library(reshape2)
library(RColorBrewer)
library(ape)

#check this function help to see what to do when multiple values occur
HeatmapMatrix = acast(TotalBLASTResults, QuerySeqID ~ TargetOrganism, value.var = "PercentIdent") # Converts TotalBLASTResults to wide format data matrix.
HeatmapMatrix[is.na(HeatmapMatrix)] = 0 # Replaces NA values with zero for heatmap function.

tree = read.tree("/Users/lee/Desktop/Steriod16S.nwk")

print("Please select a root.")
plot(tree)

tree = root(tree, resolve.root = TRUE, interactive = TRUE)

print("Thanks, continuing...")
plot(tree)

tree$edge.length[which(tree$edge.length == 0)] = 0.00001
tree <- chronopl(tree, lambda = 0.1, tol = 0)
plot(tree)

dendrogram = as.dendrogram(as.hclust.phylo(tree))
plot(dendrogram)

cladeOrder = order.dendrogram(dendrogram)
cladeName =  labels(dendrogram)
cladePosition = data.frame(cladeName, cladeOrder)
print(cladePosition)
cladePosition = cladePosition[order(cladePosition$cladeOrder),]
print(cladePosition)

heatmap(HeatmapMatrix, col = brewer.pal(8,"Blues"), breaks = c(0,30,40,50,60,70,80,90,100), Rowv = NA, Colv = NA)