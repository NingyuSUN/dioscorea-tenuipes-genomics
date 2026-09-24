suppressPackageStartupMessages(library(phangorn));a<-commandArgs(TRUE);d<-a[1]
mat<-as.matrix(read.delim(file.path(d,'distance.tsv'),row.names=1,check.names=FALSE));ord<-scan(file.path(d,'cycle.txt'),quiet=TRUE)
net<-neighborNet(as.dist(mat),ord=ord);ss<-net$splits;labs<-attr(ss,'labels');ww<-attr(ss,'weights');if(is.null(labs))labs<-net$tip.label
if(is.null(ww))ww<-sapply(seq_along(ss),function(i)net$edge.length[which(net$splitIndex==i)[1]])
write.table(data.frame(side1=sapply(ss,function(x)paste(labs[x],collapse=',')),weight=ww),file.path(d,'phangorn_same_order_splits.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
