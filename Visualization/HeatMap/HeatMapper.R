#!/usr/bin/env Rscript 
# Created by: Lee Bergstrand & Erick Cardenas Poire 
# Descript: Converts BLAST TotalBLASTResultss from BackBlast into a heatmap.  
#             
# Requirements: - reshape2 and RColorBrewer modules
#----------------------------------------------------------------------------------------
library(reshape2)
library(RColorBrewer)
library(ape)

# Sets the working directory. 
setwd("/Users/lee/Dropbox/RandD/Repositories/BackBLAST-Gene-Cluster-Finder/Visualization/HeatMap/TestData")
tree = read.tree("/Users/lee/Dropbox/RandD/Repositories/BackBLAST-Gene-Cluster-Finder/Visualization/HeatMap/SteriodBacteria16SAligned2.nwk")

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

# The Tree needs to be rooted so get the user to interactivly pick a root.
print("Please select a root.")
plot(tree)
tree = root(tree, resolve.root = TRUE, interactive = TRUE)
print("Thanks, continuing...")

# The tree must be ultrametric to be able to be passed to hclust.
tree$edge.length[which(tree$edge.length == 0)] = 0.00001
tree <- chronopl(tree, lambda = 0.1, tol = 0)

dendrogram = as.dendrogram(as.hclust.phylo(tree))

cladeOrder = order.dendrogram(dendrogram)
cladeName =  labels(dendrogram)

# Rearanges matrix to match tip ordering of the dendrogram.
cladePosition = data.frame(cladeName, cladeOrder)
cladePosition = cladePosition[order(cladePosition$cladeOrder),]
newOrder = match(cladePosition$cladeName, colnames(HeatmapMatrix))
OrderedHeatmapMatrix = HeatmapMatrix[,newOrder]

OrderedHeatmapMatrix = OrderedHeatmapMatrix[order(rownames(OrderedHeatmapMatrix)), ] # Reorders rows by gene name.

heatmap(OrderedHeatmapMatrix, col = brewer.pal(9,"Greens"), Rowv = NA, Colv = dendrogram)