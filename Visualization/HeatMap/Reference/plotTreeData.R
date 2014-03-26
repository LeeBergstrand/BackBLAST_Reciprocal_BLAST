plotTreeData<-function(treeFile,matrixFile=NULL,infoFile=NULL,locFile=NULL,outputPDF=NULL,outputPNG=NULL,w,h,matrix.colours=rev(gray(seq(0,1,0.1))),matrix.legend=F,tip.labels=F,offset=0,tip.colour.cex=0.5,legend=T,legend.pos="bottomleft",ancestral.reconstruction=F,boundaries=c(0.5,0.75),cluster=F,locColours=NULL,lwd=1.5,axis=T,axisPos=3,edge.color="black") {

require(ape)

# ladderize tree and extract tip order
t<-read.tree(treeFile)
tl<-ladderize(t)
tips<-tl$edge[,2] # tl[edge,2]
tip.order<-tips[tips<=length(tl$tip.label)]
tip.label.order<-tl$tip.label[tip.order]

# prepare heatmap matrix
if (!is.null(matrixFile)) {
if (is.matrix(matrixFile)) {
x = data.frame(matrixFile)
}
else if (is.data.frame(matrixFile)) {
x = matrixFile
}
else {
x<-read.csv(matrixFile,row.names=1)
}
y.ordered<-x[tip.label.order,]
if (cluster) {
h<-hclust(dist(t(na.omit(y.ordered))))
y.ordered<-y.ordered[,h$order]
}
}

# prepare coloured labels for tree leaves
if (!is.null(locFile)) { 
loc<-read.csv(locFile,row.names=1)
loc1<-as.matrix(loc)[row.names(loc) %in% tl$tip.label,] #vector
tipLabelSet <- character(length(loc1))
names(tipLabelSet) <- names(loc1)
groups<-table(loc1)
n<-length(groups)
groupNames<-names(groups)
if (is.null(locColours)){
colours<-rainbow(n)
}
else{
colours<-locColours
}
for (i in 1:n) {
g<-groupNames[i]
tipLabelSet[loc1==g]<-colours[i]
}
# ancestral reconstruction
if (ancestral.reconstruction) {
ancestral<-ace(loc1,tl,type="discrete")
}
else{
ancestral=NULL
}
}
else{
ancestral=NULL
}

# order additional info
if (!is.null(infoFile)) {
ids<-read.csv(infoFile,row.names=1)
ids.ordered<-ids[rev(tip.label.order),]
}
else {ids.ordered=NULL}

# open PDF for drawing
if (!is.null(outputPDF)) {
pdf(width=w,height=h,file=outputPDF)
}
# open PNG for drawing
if (!is.null(outputPNG)) {
png(width=w,height=h,file=outputPNG)
}

# plot tree
par(fig=c(0,boundaries[1],0,1))
tlp<-plot(tl,no.margin=T,show.tip.label=tip.labels,label.offset=offset,edge.width=lwd,edge.color=edge.color)
if (!is.null(locFile)) { 
tiplabels(col= tipLabelSet[tl$tip.label],pch=16,cex=tip.colour.cex) 
if (ancestral.reconstruction) {
nodelabels(pie=ancestral$lik.anc, cex=0.5, piecol=colours)
}
if (axis) {
axisPhylo(axisPos)
}
}
if (matrix.legend && ncol(y.ordered)<20) { text(labels=colnames(y.ordered),x=rep(tlp$x.lim[2]/2,ncol(y.ordered)),y=c(ncol(y.ordered):1)*10) }
if (legend && !is.null(locFile)) {
legend(legend.pos,legend=groupNames,fill=colours)
}

# plot info
if (!is.null(infoFile) && length(boundaries)>1) {
par(fig=c(boundaries[1],boundaries[2],0,1),new=T)
plot(rep(ncol(ids.ordered),nrow(ids.ordered)),1:nrow(ids.ordered),axes=F,pch="",xlim=c(0.5,ncol(ids.ordered)+0.5))
for (i in 1:ncol(ids.ordered)) {
text(rep(i,nrow(ids.ordered)),nrow(ids.ordered):1,ids.ordered[,i],cex=0.5)
}
}

# plot heatmap
if (!is.null(matrixFile)) {
par(fig=c(boundaries[length(boundaries)],1,0.035,0.965),new=T)
image(as.matrix(t(y.ordered)),col=matrix.colours,axes=F)
}

# close drawing device
if (!is.null(outputPDF) | !is.null(outputPNG)) {
dev.off()
}

# return ordered info and ancestral reconstruction object
return(list(id=ids.ordered,anc=ancestral,mat=as.matrix(t(y.ordered))))
}
