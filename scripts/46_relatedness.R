suppressPackageStartupMessages({library(SNPRelate);library(GENESIS);library(GWASTools);library(BiocParallel)})
args<-commandArgs(TRUE); root<-args[1]; old<-args[2]; out<-file.path(root,'relatedness');dir.create(out,showWarnings=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
panels<-c(supplement17=file.path(old,'structure/supplement17/pruned'),strict17=file.path(root,'strict17/pruned'),complete17=file.path(root,'complete17/pruned'))
status<-list()
for(panel in names(panels)){
 p<-panels[[panel]];d<-file.path(out,panel);dir.create(d,showWarnings=FALSE);gf<-file.path(d,'panel.gds')
 if(!file.exists(gf))snpgdsBED2GDS(paste0(p,'.bed'),paste0(p,'.fam'),paste0(p,'.bim'),gf,snpfirstdim=TRUE,cvt.chr='char')
 gds<-snpgdsOpen(gf);sids<-read.gdsn(index.gdsn(gds,'sample.id'));vids<-read.gdsn(index.gdsn(gds,'snp.id'))
 kg<-snpgdsIBDKING(gds,type='KING-robust',autosome.only=FALSE,missing.rate=1,num.thread=2)
 km<-kg$kinship;dimnames(km)<-list(kg$sample.id,kg$sample.id)
 write.table(km,file.path(d,'king_matrix.tsv'),sep='\t',quote=FALSE,col.names=NA)
 # Export decoded genotypes at selected loci for independent BED reader verification.
 z<-snpgdsGetGeno(gds,snp.id=vids[seq_len(min(50,length(vids)))],snpfirstdim=TRUE,verbose=FALSE)
 write.table(z,file.path(d,'first50_genotypes.tsv'),sep='\t',quote=FALSE,row.names=FALSE,col.names=sids,na='NA')
 airs<-list()
 for(th in c(.022,.044,.088)){
  tag<-sprintf('threshold_%03d',round(th*1000));ad<-file.path(d,tag);dir.create(ad,showWarnings=FALSE)
  pa<-tryCatch(pcair(gds,kinobj=km,divobj=km,kin.thresh=th,div.thresh=-.022,snp.include=vids,autosome.only=FALSE,num.cores=2,eigen.cnt=3),error=function(e)e)
  if(inherits(pa,'error')){writeLines(conditionMessage(pa),file.path(ad,'ERROR.txt'));next}
  saveRDS(pa,file.path(ad,'pcair.rds'));write.table(data.frame(sample=rownames(pa$vectors),pa$vectors),file.path(ad,'pcair.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  writeLines(pa$unrels,file.path(ad,'training_set.txt'));airs[[tag]]<-pa
 }
 snpgdsClose(gds)
 reader<-GdsGenotypeReader(gf);gd<-GenotypeData(reader);it<-GenotypeBlockIterator(gd,snpBlock=10000)
 for(tag in names(airs)){
  pa<-airs[[tag]];ad<-file.path(d,tag)
  for(k in 0:2){
   pcs<-pa$vectors[,seq_len(k),drop=FALSE];design<-cbind(1,pcs[pa$unrels,,drop=FALSE]);rank<-qr(design)$rank
   # Algebraic estimability is necessary, not sufficient scientific calibration.
   if(rank<ncol(design)||length(pa$unrels)<=rank){writeLines('Insufficient training rank/residual degrees of freedom',file.path(ad,paste0('PC',k,'_NOT_ESTIMABLE.txt')));next}
   pr<-tryCatch(pcrelate(it,pcs=pcs,training.set=pa$unrels,ibd.probs=FALSE,small.samp.correct=TRUE,BPPARAM=SerialParam()),error=function(e)e)
   if(inherits(pr,'error')){writeLines(conditionMessage(pr),file.path(ad,paste0('PC',k,'_ERROR.txt')));next}
   write.table(pr$kinBtwn,file.path(ad,paste0('PC',k,'_pairs.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
   write.table(pr$kinSelf,file.path(ad,paste0('PC',k,'_self.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
   status[[length(status)+1]]<-data.frame(panel=panel,threshold=tag,PCs=k,training_n=length(pa$unrels),rank=rank,residual_df=length(pa$unrels)-rank)
  }
 }
 close(it)
}
if(length(status))write.table(do.call(rbind,status),file.path(out,'estimability.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
# Missing-data fixture distinguishes the KING denominator used by this installed package.
f<-file.path(out,'fixture.gds');x<-cbind(c(1,1,0,1,0,2,0,1,2,1,NA,0),c(1,1,2,0,1,2,0,1,1,NA,1,0))
if(!file.exists(f))snpgdsCreateGeno(f,x,sample.id=c('fixture_i','fixture_j'),snp.id=1:12,snp.chromosome=rep(1,12),snp.position=1:12,snpfirstdim=TRUE)
g<-snpgdsOpen(f);v<-snpgdsIBDKING(g,autosome.only=FALSE,remove.monosnp=FALSE,missing.rate=1,type='KING-robust');snpgdsClose(g)
writeLines(capture.output(v),file.path(out,'king_fixture.txt'))
writeLines('Completed; estimates exploratory, no pedigree labels assigned.',file.path(out,'COMPLETE'))
