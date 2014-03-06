#!/usr/bin/env Rscript 
# Created by: Lee Bergstrand & Erick Cardenas Poire 
# Descript: Converts BLAST results BackBlast into a heatmap.  
#             
# Requirements: - reshape2 and RColorBrewer modules
#----------------------------------------------------------------------------------------

# Set working directory. 
setwd("/Users/lee/Dropbox/Repositories/BackBLAST_Reciprocal_BLAST/Visualization/HeatMap/TestData")

# Gets a list files from the working directory.
alist = list.files(path = getwd(), all.files = FALSE, pattern = "\\.csv$") 

result = matrix(data = 0, nrow = 0,ncol = 7)
colnames(result) = c("V1","V2","V3","V4","V5","V6","filename")
dim(result)

counter = 1
counter=2
while (counter <= length(alist)){
  set_name = gsub("( )", "", paste(alist[counter]))
  data = read.csv(set_name,header=FALSE,stringsAsFactors=TRUE)
  filename = matrix(data = alist[counter], nrow = nrow(data), ncol = 1)
  data2  = cbind(data, filename)
  result = rbind(result,data2)
  #print(filename)
  counter = counter + 1
}

library(reshape2)
library(RColorBrewer)

#check this function help to see what to do when multiple values occur
mymatrix = acast(result, V1~filename, value.var="V3")

head(mymatrix)
mymatrix[is.na(mymatrix)] <- 0
head(mymatrix)

heatmap(mymatrix, col = brewer.pal(8,"Blues"),
        breaks=c(0,30,40,50,60,70,80,90,100),

        , Rowv = NA) #, #Colv = NA
