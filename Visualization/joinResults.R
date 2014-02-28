#!/usr/bin/env Rscript
# Created by: Lee Bergstrand and Erick Cardenas Poire
# Descript: 
#             
# Usage: R joinResults.R <resultsCSV1.csv> .... <resultsCSVN.csv>  
# Example: R joinResults.R ./*.csv 
#----------------------------------------------------------------------------------------

# Get file name from the command line arguments.
Args = commandArgs()
ResultFileList = Args[-c(1:5)]
print(ResultFileList)

check <- read.table(ResultFileList[1], header = F, sep = ",") # Read file
cols = ncol(check) # Get number of columns.

print(cols)

#get list of all possible genes 
#counter=3
#old=matrix(NA,ncol=mycols,nrow=0)
#while (counter <= length(alist)){
#  set_name<-gsub("( )", "", paste("c:/R/temp/",alist[counter]))
#  data<-read.delim(set_name,header=F)
#  new=rbind(old,data)
#  old<-new
#  counter=counter+1
#}

#myrows=length(levels(new[,1]))

#rm(counter)
#rm(data)
#rm(old)
#rm(set_name)



#create empty matrix

#result=matrix(data=NA,nrow=myrows,ncol=(length(alist)-2))
#rownames(result)=levels(new[,1])
#colnames(result)=alist[c(3:length(alist))]
#dim(result)
#colnames(result)

#counter=3
#while (counter<=length(alist)){
#  set_name<-gsub("( )", "", paste("c:/R/temp/",alist[counter]))
#  data<-read.delim(set_name,header=F)
#  c=1
#  while (c<=nrow(data)){
#    key=data[c,1]
#    door=data[c,2]
#    rr=alist[counter]
#    result[key,rr]=door
#    c=c+1
#  }
#  counter=counter+1
#}

#result




