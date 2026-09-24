suppressPackageStartupMessages({library(ape);library(phangorn)})
a<-commandArgs(TRUE);root<-a[1];folders<-list.dirs(root,recursive=FALSE,full.names=TRUE);checks<-list()
for(d in folders){
 f<-file.path(d,'distance.tsv');if(!file.exists(f))next
 mat<-as.matrix(read.delim(f,row.names=1,check.names=FALSE));dm<-as.dist(mat)
 tr<-nj(dm);other<-read.tree(file.path(d,'nj.nwk'));rf<-RF.dist(tr,other,normalize=FALSE,check.labels=TRUE,rooted=FALSE)
 checks[[length(checks)+1]]<-data.frame(panel=basename(d),RF_ape_Biopython=rf,max_patristic_difference=max(abs(cophenetic(tr)[rownames(mat),rownames(mat)]-cophenetic(other)[rownames(mat),rownames(mat)])))
 if(basename(d)%in%c('supplement17','strict17','complete17','one_per_contig17')){
  # The plot is unrooted; arbitrary drawing orientation has no ancestral meaning.
  labs<-tr$tip.label;cols<-ifelse(grepl('^AIC',labs),'#c65a29',ifelse(grepl('^HCO',labs),'#287f9a',ifelse(grepl('^IZU',labs),'#679b39','#794ca1')))
  for(ext in c('pdf','png')){
   dest<-file.path(d,paste0('nj_unrooted.',ext));if(ext=='pdf')pdf(dest,width=12,height=9) else png(dest,width=2400,height=1800,res=200)
   par(mar=c(3,3,4,3));plot(tr,type='unrooted',tip.color=cols,cex=.8,no.margin=FALSE,lab4ut='axial');title(paste(basename(d),'- unrooted Neighbor-Joining'));add.scale.bar(length=.05,cex=.7);mtext('Diploid IBS distance; branch lengths are descriptive, with no ancestral direction.',side=1,line=1,cex=.8);dev.off()
  }
 }
 if(basename(d)=='supplement17'){
  net<-neighborNet(dm);saveRDS(net,file.path(d,'phangorn_network.rds'))
  writeLines(capture.output(str(net)),file.path(d,'phangorn_structure.txt'))
  ss<-net$splits;labs<-attr(ss,'labels');ww<-attr(ss,'weights'); if(is.null(labs))labs<-net$tip.label
  if(is.null(ww))ww<-sapply(seq_along(ss),function(i)net$edge.length[which(net$splitIndex==i)[1]])
  write.table(data.frame(side1=sapply(ss,function(x)paste(labs[x],collapse=',')),weight=ww),file.path(d,'phangorn_splits.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
  pdf(file.path(d,'network_phangorn.pdf'),width=12,height=9);plot(net,cex=.8,main='Independent phangorn Neighbor-Net');dev.off()
 }
}
write.table(do.call(rbind,checks),file.path(root,'cross_implementation.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(root,'R_sessionInfo.txt'))
