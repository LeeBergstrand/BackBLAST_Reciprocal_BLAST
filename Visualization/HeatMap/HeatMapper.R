#add nane of file as a columnt to csv

setwd("/Users/lee/Dropbox/Repositories/BackBLAST_Reciprocal_BLAST/Visualization/TestData/")

alist <-list.files(path=getwd(), all.files=FALSE)
alist

result=matrix(data=0,nrow=0,ncol=7)
colnames(result)=c("V1","V2","V3","V4","V5","V6","filename")
dim(result)

counter=1

while (counter <= length(alist)){
  set_name<-gsub("( )", "", paste(alist[counter]))
  data<-read.csv(set_name,header=FALSE,stringsAsFactors=TRUE)
  filename=matrix(data=alist[counter],nrow=nrow(data),ncol=1)
  data2=cbind(data,filename)
  result=rbind(result,data2)
  print (counter)
  counter=counter+1
}


library(reshape2)
library(RColorBrewer)

#check this function help to see what to do when multiple values occur
mymatrix=acast(result, V1~filename, value.var="V3", fun.aggregate=max)

head(mymatrix)
