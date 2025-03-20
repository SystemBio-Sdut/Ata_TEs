tq_go_gene <- function(allgog){
  
  
  n1 <- dim(allgog)[1]
  a1 <- allgog[,1]
  a2 <- allgog[,2]
  
  allc <- c()
  
  for(i in 1:n1){
    tmp <- a2[i]
    if(tmp==""){
      next
    }
    tmp1 <- strsplit(tmp,",")[[1]]
    if(length(tmp1)==1){
      allc <- rbind(allc,c(a2[i],a1[i]))
    }
    if(length(tmp1)>1){
      re1 <- rep(a1[i],length(tmp1))
      re2 <- cbind(tmp1,re1)
      allc <- rbind(allc,re2)
    }
  }
  return(allc)
}


Clust.kmeans = function(dat,standard=TRUE,nclust=NULL,
                        max.clust=50,rss.thresh=0.75,iter.max=1000,nstart=10)
  
{
  if(!is.matrix(dat)){
    dat = as.matrix(dat)
  }
  # gene-wise standardization
  if(standard){
    require("preprocessCore")
    dat2 <- normalize.quantiles(dat)
    #dat = t(apply(dat,1,function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)))
  }
  
  
  if(is.null(nclust)){
    for(nclust in 2:max.clust){
      clust.res = kmeans(dat2,centers=nclust,iter.max=iter.max,nstart=nstart)
      rss.ratio = 1-clust.res$tot.withinss/clust.res$totss
      if(rss.ratio >= rss.thresh){
        break
      }
    }
    if(nclust == max.clust & rss.ratio < rss.thresh){
      cat("Warning: the maximum number of clusters reached before the percentage of variance explained exceeds threshold! \n")
    }
  }else{
    clust.res = kmeans(dat2,centers=nclust,iter.max=iter.max,nstart=nstart)
  }
  
  list(nclust=nclust,clust=clust.res$cluster,center=clust.res$centers,totss=clust.res$totss,
       withinss=clust.res$withinss,tot.withinss=clust.res$tot.withinss,betweenss=clust.res$betweenss,
       size=clust.res$size,dat.clust=dat)
}  
useplot_news1 <- function(ret,filename="Figre_Line-type_genic.pdf"){
  
  ii <- ret$index
  ii1 <- c("Intergenic","Distal_intergenic","Genic","Intron","Promoter","Exon","Downstream","ThreeUTR","FiveUTR")
  
  dat <- ret$tmm
  nn <- names(dat)
  names(dat) <- rep("",length(nn))
  datl <- log10(dat)
  
  cairo_pdf(filename, width=7.5, height=4)
  par(mar=c(0, 5.5, 0.5, 0.7),fig=c(0.25,1,0.5,1),mgp=c(3,0.8,0))
  aa <- barplot(datl[1:27], beside = TRUE,
                col = c("#4EEE94"),
                legend.text = names(datl), args.legend=list(xjust=1.5,yjust=6.2),ylim = c(-0.02,6.52),
                xlim=c(0.7,31.6),
                main = " ", font.main = 4,
                sub = " ",
                cex.names = 1.6,axes=F,border=NA)
  box()
  axis(1,aa[,1],labels=F)
  axis(2,seq(0,6,2),seq(0,6,2),las=2,cex.axis=1.0,adj=0.5)
  mtext(bquote("log"["10"]* "(Interaction Size)"),2,line=1.5,cex=1.0)
  
  par(mar=c(0, 0, 0.2, 0),fig=c(0.23,1,0,0.49),mgp=c(3,0.8,0),new=TRUE)
  plot(NA, NA, type="n", axes=FALSE, ylim = c(-0.02,9.02),xlim=c(-15,75),xlab=" ", ylab=" ", main="")
  for(i in 1:9){
    if(i%%2==0){
      rect(3,9-i+0.1,75.7,9-i+1,border=NA,col="#98FB9830")
    }else{
      rect(3,9-i+0.1,75.7,9-i+1,border=NA,col="#FFE4E1")
    }
    text(4,9-i+0.55,ii1[i],pos=2,cex=0.9)
  }
  nf <- c()
  for(j in 1:27){
    nj <- strsplit(nn[j],"-")[[1]]
    nf <- c(nf,paste(ii[as.numeric(nj)],collapse="-"))
    pd <- 9-as.numeric(nj)+0.55
    rect(5.29+(j-1)*2.65-0.25,min(pd),5.29+(j-1)*2.65+0.25,max(pd),border = NA,col="#B3B3B3")
    points(rep(5.28+(j-1)*2.65,length(nj)),pd,pch=20,col="#DAA520",cex=2)
  }
  
  par(mar=c(0, 0, 0, 0),fig=c(0.02,0.23,0.01,0.47),mgp=c(3,0.8,0),new=TRUE)
  aa1 <- barplot(-log10(ret$size), beside = TRUE,horiz =TRUE,
                 col = c("#836FFF"),
                 legend.text = names(datl), args.legend=list(xjust=-100.5,yjust=-600.2),ylim = c(0.3,10.62),
                 xlim=c(-7.1,0.1),angle=45,
                 main = " ", font.main = 4,
                 sub = " ",
                 cex.names = 1.6,axes=F,border=NA)
  box()
  axis(4,aa1[,1],labels=F)
  axis(3,seq(-6,0,2),rev(seq(0,6,2)),las=2,cex.axis=1.0,adj=0.5,las=1,mgp=c(3,0.5,0))
  mtext(bquote("log"["10"]* "(Set Size)"),3,line=1.5,cex=1.0)
  
  dev.off()
  
  
}




ann_combn2 <- function(diff){
  
  index <- c("Intergenic","distal_intergenic","genic","Intron","Promoter","Exon","downstream","threeUTR","fiveUTR")
  
  
  y1 <- diff@detailGenomicAnnotation
  nn1 <- names(y1)
  y1 <- as.matrix(y1)
  res_diff <- tibble::tibble(anno = lapply(1:nrow(y1), function(i) nn1[y1[i,]]))[[1]]
  res_diff1 <- unlist(lapply(1:nrow(y1),function(i) paste(sort(match(res_diff[[i]],index)),collapse="-")))
  
  t1 <- table(res_diff1)
  
  tt <- unique(c(names(t1)))
  ttm <- matrix(0,nrow=1,ncol=length(tt))
  ttm[1,match(names(t1),tt)] <- t1
  colnames(ttm) <- tt
  
  ttm1 <- ttm[,rev(order(colSums(ttm)))]
  
  y11 <- y1[,match(index,nn1)]
  ys <- colSums(y11)
  
  ret <- list(tmm=ttm1,index=index,t1=res_diff1,size=ys)
  return(ret)
}


genic_rep <- function(x1,x2){
  
  sn <- c("Promoter","Intron", "Exon","downstream","fiveUTR","threeUTR")
  alln <- colnames(x2)
  tn <- names(sort(table(x1$tclass),decreasing = T))
  prol <- c()
  for(i in 1:length(sn)){
    ii <- which(x2[,match(sn[i],alln)]==1)
    tmpii <- x1[ii,]
    etm <- c()
    for(j in 1:length(tn)){
      jj <- which(tmpii$tclass==tn[j])
      if(length(jj)<=1){
        etm <- c(etm,0)
      }else{
        etm <- c(etm,sum(tmpii$qend[jj]-tmpii$qstart[jj]))
      }
    }
    prol <- rbind(prol,etm)
  }
  
  rownames(prol) <- sn
  colnames(prol) <- tn
  return(prol)
}




genic_rep_number <- function(x1,x2,x3){
  
  sn <- c("Promoter","Intron", "Exon","downstream","fiveUTR","threeUTR")
  alln <- colnames(x2)
  tn <- names(sort(table(x1$tclass),decreasing = T))
  prol <- c()
  for(i in 1:length(sn)){
    ii <- which(x2[,match(sn[i],alln)]==1)
    tmpii <- x1[ii,]
    tmpiii <- x3[ii,]
    etm <- c()
    for(j in 1:length(tn)){
      jj <- which(tmpii$tclass==tn[j])
      if(length(jj)<=1){
        etm <- c(etm,0)
      }else{
        etm <- c(etm,length(unique(tmpiii$geneId[jj])))
      }
    }
    prol <- rbind(prol,etm)
  }
  
  rownames(prol) <- sn
  colnames(prol) <- tn
  return(prol)
}




Figure_ls_genic_pro <- function(x1,y1,x2,y2){
  
  
  p1 <- sum(y1$qend-y1$qstart)/sum(x1$qend-x1$qstart)
  
  p2 <- sum(y2$qend-y2$qstart)/sum(x2$qend-x2$qstart)
  
  pdf("Figure-ls_genic_pro.pdf",width=2,height=4)
  par(oma=c(0,0,0,0),fig=c(0,1,0,1))
  par(mar=c(2.5,3.5,0.5,0.5))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(0.45,2.51),ylim=c(-0.001,0.22),xlab="",ylab="")
  rect(1-0.4,0,1+0.4,p1,border=NA,col="#4EEE94")
  rect(2-0.4,0,2+0.4,p2,border=NA,col="#FF34B3")
  
  axis(1,c(1,2),c("LINE","SINE"),las=1,cex.axis=0.9,tck=-0.033,mgp=c(2.5,0.3,0),col="black")
  axis(2,seq(0,0.2,0.05),seq(0,20,5),las=1,cex.axis=1,tck=-0.043,mgp=c(2.5,0.6,0),col="black")
  mtext("% of repeats in genic regions (< 3 kb)",2,cex=1.2,line=2.1)
  dev.off()
}




Figure_ls_genic_pro_type <- function(y1,y2){
  
  
  tn1 <- names(sort(table(y1$tclass),decreasing = T))[1:5]
  tn2 <- names(sort(table(y2$tclass),decreasing = T))[1:5]
  
  sum1 <- sum(y1$qend-y1$qstart)
  sum2 <- sum(y2$qend-y2$qstart)
  
  p11 <- c()
  for(i in 1:length(tn1)){
    i1 <- which(y1$tclass==tn1[i])
    p11 <- c(p11,sum(y1$qend[i1]-y1$qstart[i1])/sum1)
  }
  
  p22 <- c()
  for(i in 1:length(tn2)){
    i2 <- which(y2$tclass==tn2[i])
    p22 <- c(p22,sum(y2$qend[i2]-y2$qstart[i2])/sum2)
  }
  
  p111 <- p11[c(1,3,4,5,2)]*100
  tnn1 <- c("L1","L2","I-Jockey","RTE-BovB","Unkown")
  p222 <- p22*100
  tnn2 <- c("tRNA","B2","ID","MIR","Unkown")
  
  
  pdf("Figure-ls_genic_pro_type.pdf",width=4,height=4.2)
  par(oma=c(0,0,0,0),fig=c(0,1,0,1))
  par(mar=c(4.5,3.5,0.5,0.5))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(0.45,10.51),ylim=c(-0.02,4.12),xlab="",ylab="")
  for(i in 1:length(p111)){
    if(i==1){
      rect(i-0.4,0,i+0.4,p111[i]*0.04,border=NA,col="#4EEE94")
      
    }else{
      rect(i-0.4,0,i+0.4,p111[i],border=NA,col="#4EEE94")
    }
  }
  
  for(i in 1:length(p111)){
    if(i==1){
      rect(i+5-0.4,0,i+5+0.4,p222[i]*0.04,border=NA,col="#FF34B3")
      
    }else{
      rect(i+5-0.4,0,i+5+0.4,p222[i],border=NA,col="#FF34B3")
    }
  }
  
  axis(1,seq(1:10),c(tnn1,tnn2),las=1,cex.axis=0.8,tck=-0.033,mgp=c(2.5,0.6,0),col="black",las=2)
  axis(2,c(0,1,2,3,4),c(0,1,2,90,100),las=1,cex.axis=1,tck=-0.023,mgp=c(2.5,0.6,0),col="black")
  segments(-300,2.5,200,2.5,col="white",lwd=48)
  segments(1,2.24,1,2.8,col="#4EEE94",lwd=2,lty=3)
  segments(6,2.24,6,2.8,col="#FF34B3",lwd=2,lty=3)
  mtext("% of Super family",2,cex=1.2,line=2.1)
  rect(2.2,3.4,2.7,3.5,col="#4EEE94",border=NA)
  text(3.6,3.5,"LINE",cex=1)
  rect(2.2+5,3.4,2.7+5,3.5,col="#FF34B3",border=NA)
  text(3.6+5,3.5,"SINE",cex=1)
  dev.off()
}


Figure_rep_gene_pro <- function(lm,sm){
  
  lm1 <- lm[,c(1,3,4,5,2)]/sum(lm)*100
  sm1 <- sm[,c(1:5)]/sum(sm)*100
  
  tnn1 <- c("L1","L2","I-Jockey","RTE-BovB","Unkown")
  tnn2 <- c("tRNA","B2","ID","MIR","Unkown")
  tnc <- c(tnn1,tnn2)
  lsm <- rbind(t(lm1),t(sm1))
  lsm[which(lsm>20)] <- 20
  
  lsm1 <- round(rbind(t(lm1),t(sm1)),2)
  
  ccn <- c("Promoter","Intron","Exon","Downstream","5'UTR","3'UTR")
  
  library(RColorBrewer)
  coll <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(20)
  nd <- as.numeric(lsm) 
  ndc <- cut(nd,20)
  colm <- matrix(coll[match(ndc,levels(ndc))],ncol=6,byrow=F)
  
  pdf("Figure_rep_gene_pro.pdf",height=2.7,width=5.6)
  par(mar=c(0,0,0,0),oma=c(0,0,0,0),bty="n")
  plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(-1,7.1),ylim=c(-1,11.2),
       xaxt="n",yaxt="n",xaxs="i", yaxs="i")
  for(i in 1:10){
    for(j in 1:6){
      rect(j-0.5,11-i,j+0.5,11-i-1,col=paste(colm[i,j],60,sep=""),border="#878787")
      text(j,11-i-0.5,paste(lsm1[i,j],"%",sep=""),cex=0.7,font=2)
    }
    text(0.2,11-i-0.5,tnc[i],cex=0.9,adj=1)
  }
  for(i in 1:10){
    rect(0.5,11-i-0.45,0.3,11-i-0.55,col="grey",border=NA)
  }
  for(i in 1:6){
    text(i,10.5,ccn[i],cex=0.9)
  }
  
  rect(6.6,10,6.7,5,col="#4EEE94",border=NA)
  rect(6.6,5,6.7,0,col="#FF34B3",border=NA)
  text(6.9,7.5,"LINE",srt=-90,cex=1.1)
  text(6.9,2.5,"SINE",srt=-90,cex=1.1)
  dev.off()
  
}


Figure_genic_pro_pie <- function(ret1,ret2){
  
  lsgenic  <- ret1$size[-c(1:3)]
  lsp <- lsgenic/sum(lsgenic)
  ssgenic  <- ret2$size[-c(1:3)]
  ssp <- ssgenic/sum(ssgenic)
  lsp1 <- c(lsp[1:3],sum(lsp[-c(1:3)]))
  names(lsp1) <- c(names(lsp[1:3]),"Other")
  
  ssp1 <- c(ssp[1:3],sum(ssp[-c(1:3)]))
  names(ssp1) <- c(names(ssp[1:3]),"Other")
  
  cols = c("#6495ED","#FF7F50","#A2CD5A","#EE7600")
  pdf("figue-genic-pro-pie.pdf",height=4,width=8.2)
  par(fig=c(0,0.5,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0))
  pie(lsp1, labels=round(lsp*100,2), col=cols)
  # 添加颜色样本标注
  legend("topright", names(lsp1), cex=0.8, fill=cols)
  par(fig=c(0.5,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE)
  pie(ssp1, labels=round(ssp*100,2), col=cols)
  dev.off()
}

Figure_rep_gene_size <- function(ls,ss){
  
  lm1 <- ls[,c(1,3,4,5,2)]/10^3
  sm1 <- ss[,c(1:5)]/10^3
  
  tnn1 <- c("L1","L2","I-Jockey","RTE-BovB","Unkown")
  tnn2 <- c("tRNA","B2","ID","MIR","Unkown")
  tnc <- c(tnn1,tnn2)
  lsm <- rbind(t(lm1),t(sm1))
  #lsm[which(lsm>20)] <- 20
  
  lsm1 <- round(rbind(t(lm1),t(sm1)),2)
  
  ccn <- c("Promoter","Intron","Exon","Downstream","5'UTR","3'UTR")
  
  library(RColorBrewer)
  coll <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(20)
  nd <- as.numeric(lsm) 
  ndc <- cut(nd,20)
  colm <- matrix(coll[match(ndc,levels(ndc))],ncol=6,byrow=F)
  
  pdf("Figure_rep_gene_size.pdf",height=2.7,width=5.6)
  par(mar=c(0,0,0,0),oma=c(0,0,0,0),bty="n")
  plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(-1,7.1),ylim=c(-1,11.2),
       xaxt="n",yaxt="n",xaxs="i", yaxs="i")
  for(i in 1:10){
    for(j in 1:6){
      rect(j-0.5,11-i,j+0.5,11-i-1,col=paste(colm[i,j],60,sep=""),border="#878787")
      text(j,11-i-0.5,lsm1[i,j],cex=0.7,font=2)
    }
    text(0.2,11-i-0.5,tnc[i],cex=0.9,adj=1)
  }
  for(i in 1:10){
    rect(0.5,11-i-0.45,0.3,11-i-0.55,col="grey",border=NA)
  }
  for(i in 1:6){
    text(i,10.5,ccn[i],cex=0.9)
  }
  
  rect(6.6,10,6.7,5,col="#4EEE94",border=NA)
  rect(6.6,5,6.7,0,col="#FF34B3",border=NA)
  text(6.9,7.5,"LINE",srt=-90,cex=1.1)
  text(6.9,2.5,"SINE",srt=-90,cex=1.1)
  dev.off()
  
}


cu_line1 <- function(set1,set2,bin=300,filename="Figure_LINE_enriched_gene.pdf"){
  
  require(drc)
  lcmu <- c()
  lrcmu <- c()
  lll <- c()
  for(i in 1:length(set2)){
    ll <- min(c(length(set1),length(set2[[i]])))
    set21 <- -log10(set1)[1:ll]
    set22 <- -log10(set2[[i]])[1:ll]
    setc <- c(set21,set22)
    setcc <- cut(setc,breaks=bin)
    
    setc1 <- setcc[1:ll]
    setc2 <- setcc[-c(1:ll)]
    
    cul <- cumsum(table(setc1)/length(setc1))
    culc <- cumsum(table(setc2)/length(setc2))
    ll <- names(cul)
    ll1 <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", ll) ),
                 upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", ll) ))
    ll11 <- ll1[,1]+(ll1[,2]-ll1[,1])/2
    lcmu <- rbind(lcmu,cul)
    lrcmu <- rbind(lrcmu,culc)
    lll <- rbind(lll,ll11)
  }
  lcmum <- colMeans(lcmu)
  lrcmum <- colMeans(lrcmu)
  llm <- colMeans(lll)
  zx <- apply(lrcmu,2,sd)
  
  pv <- wilcox.test(lcmum,lrcmum)
  dat <- data.frame(x=llm,y=lrcmum,z=lcmum,s=zx)
  fm <- drm(y ~ x, data = dat, fct = G.3())
  fm1 <- drm(z ~ x, data = dat, fct = G.3())
  fmsd <- drm(zx ~ x, data = dat, fct = G.3())
  
  #Wilcoxon P-value <2.2e-15
  pdf(filename, height=3.5, width=4)
  par(mar=c(3,4, 0.5, 0.5),mgp=c(3,0.8,0),fig=c(0,1,0,1))
  plot(NA, NA, type="n", yaxt="n",xaxt="n", ylim = c(0.38,1.01),xlim=c(0.45,5.18),xlab="", 
       ylab="", main="")
  lines(dat$x,predict(fm),col="grey",lwd=3)
  lines(dat$x,predict(fm1),col="#4EEE94",lwd=3)
  #lines(dat$x,predict(fm)+predict(fmsd),col="black",lwd=2)
  #lines(dat$x,predict(fm)-predict(fmsd),col="black",lwd=2)
  
  axis(2,seq(0.2,1,0.2),c(seq(0.2,0.8,0.2),"1.0"),las=2,cex.axis=1.1,adj=0.5)
  axis(1,seq(0,5,1),seq(0,5,1),las=1,cex.axis=1.1,adj=0.5)
  #mtext( bquote("log"["10"] * "(P-value)"),1,line=2,cex=1.2)
  mtext( expression(paste(plain("log")[10],"(",italic(P),"-value)",sep="")),1,line=2,cex=1.2)
  mtext("Cumulative fraction",2,line=2.7,cex=1.2)
  
  text(3,0.81,expression(paste("Wilcoxon ",italic(P),"-value <2.2e-15",sep="")),cex=1.0)
  rect(1.5,0.59,2.5,0.6,col="grey",border=NA)
  rect(1.5,0.52,2.5,0.53,col="#4EEE94",border=NA)
  text(2.5,0.595,"Random selection",cex=1,pos=4)
  text(2.5,0.525,"LINE-enriched ",cex=1,pos=4)
  dev.off()
  
}

cu_sine1 <- function(set1,set2,bin=300,filename="Figure_LINE_enriched_gene.pdf"){
  
  require(drc)
  lcmu <- c()
  lrcmu <- c()
  lll <- c()
  for(i in 1:length(set2)){
    ll <- min(c(length(set1),length(set2[[i]])))
    set21 <- -log10(set1)[1:ll]
    set22 <- -log10(set2[[i]])[1:ll]
    setc <- c(set21,set22)
    setcc <- cut(setc,breaks=bin)
    
    setc1 <- setcc[1:ll]
    setc2 <- setcc[-c(1:ll)]
    
    cul <- cumsum(table(setc1)/length(setc1))
    culc <- cumsum(table(setc2)/length(setc2))
    ll <- names(cul)
    ll1 <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", ll) ),
                 upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", ll) ))
    ll11 <- ll1[,1]+(ll1[,2]-ll1[,1])/2
    lcmu <- rbind(lcmu,cul)
    lrcmu <- rbind(lrcmu,culc)
    lll <- rbind(lll,ll11)
  }
  lcmum <- colMeans(lcmu)
  lrcmum <- colMeans(lrcmu)
  llm <- colMeans(lll)
  zx <- apply(lrcmu,2,sd)
  
  pv <- wilcox.test(lcmum,lrcmum)
  dat <- data.frame(x=llm,y=lrcmum,z=lcmum,s=zx)
  fm <- drm(y ~ x, data = dat, fct = G.3())
  fm1 <- drm(z ~ x, data = dat, fct = G.3())
  fmsd <- drm(zx ~ x, data = dat, fct = G.3())
  
  #Wilcoxon P-value =5.5e-14
  pdf(filename, height=3.5, width=4)
  par(mar=c(3,4, 0.5, 0.5),mgp=c(3,0.8,0),fig=c(0,1,0,1))
  plot(NA, NA, type="n", yaxt="n",xaxt="n", ylim = c(0.38,1.01),xlim=c(0.38,5.18),xlab="", 
       ylab="", main="")
  lines(dat$x,predict(fm),col="grey",lwd=3)
  lines(dat$x,predict(fm1),col="#FF34B3",lwd=3)
  #lines(dat$x,predict(fm)+predict(fmsd),col="black",lwd=2)
  #lines(dat$x,predict(fm)-predict(fmsd),col="black",lwd=2)
  
  axis(2,seq(0.2,1,0.2),c(seq(0.2,0.8,0.2),"1.0"),las=2,cex.axis=1.1,adj=0.5)
  axis(1,seq(0,5,1),seq(0,5,1),las=1,cex.axis=1.1,adj=0.5)
  #mtext( bquote("log"["10"] * "(P-value)"),1,line=2,cex=1.2)
  mtext( expression(paste(plain("log")[10],"(",italic(P),"-value)",sep="")),1,line=2,cex=1.2)
  mtext("Cumulative fraction",2,line=2.7,cex=1.2)
  
  text(3,0.81,expression(paste("Wilcoxon ",italic(P),"-value =5.5e-14",sep="")),cex=1.0)
  rect(1.5,0.59,2.5,0.6,col="grey",border=NA)
  rect(1.5,0.52,2.5,0.53,col="#FF34B3",border=NA)
  text(2.5,0.595,"Random selection",cex=1,pos=4)
  text(2.5,0.525,"SINE-enriched ",cex=1,pos=4)
  dev.off()
  
}



Figure_LS_GO <- function(gl,gs){

pv <- c(-log10(gl$pvalue[c(1:5,7:9,11:12,14:15)]),-log10(gs$pvalue))
type <- c(gl$Term[c(1:5,7:9,11:12,14:15)],gs$Term)
leve <- c(gl$Level[c(1:5,7:9,11:12,14:15)],gs$Level)
count <- c(gl$Count[c(1:5,7:9,11:12,14:15)],gs$Count)
coll <- c(rep("#4EEE94",12),rep("#FF34B3",12))
id <- c(gl$Description[c(1:5,7:9,11:12,14:15)],gs$Description)
coll1 <- paste(coll,"50",sep="")

pdf("Figure_LS_go_en.pdf",height=4.6,width=7)
par(mar=c(0,0,0,0),oma=c(0,0,0,0),bty="n")
plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(-17.2,17.7),ylim=c(-5,24.1),
     xaxt="n",yaxt="n",xaxs="i", yaxs="i")
for(i in 1:length(pv)){
  rect(0,24-i+0.4,pv[i],24-i-0.4,col=coll[i],border="#878787")
  capitalized_strings <- toupper(substr(id[i], start = 1, stop = 1))
  text(0,24-i,paste0(capitalized_strings, substring(id[i], first = 2,last=100000)),cex=0.8,pos=2)
  rect(-0.3,24-i-0.05,0,24-i+0.05,col="grey",border=NA)
  
  rect(11,24-i+0.45,12,24-i-0.45,col=coll1[i],border="grey")
  text(11.5,24-i+0.02,type[i],cex=0.7)
  rect(12.3,+24-i+0.45,13.3,24-i-0.45,col=coll1[i],border="grey")
  text(12.8,24-i+0.02,leve[i],cex=0.7)
  rect(13.6,+24-i+0.45,14.6,24-i-0.45,col=coll1[i],border="grey")
  text(14.1,24-i+0.02,count[i],cex=0.7)
  
}
text(11.1,-0.8,"GO Type",cex=0.8,srt=-90,pos=4)
text(12.4,-0.8,"GO Level",cex=0.8,srt=-90,pos=4)
text(13.7,-0.8,"GO Count",cex=0.8,srt=-90,pos=4)
text(-7,-2.5,"GO Term",cex=1)
rect(0,-0.7,10,-0.85,col="grey",border=NA)
for(i in 0:4){
  rect(i*2,-0.85,i*2+0.1,-1.3,col="grey",border=NA)
  text(i*2+0.05,-1.9,i*2,cex=0.8)
}
rect(5*2,-0.85,5*2-0.1,-1.3,col="grey",border=NA)
text(5*2-0.05,-1.9,10,cex=0.8)
text(5.5,-3.5,expression(paste(plain("log")[10],"(",italic(P),"-value)",sep="")),cex=1)

rect(-15.3,23.4,-15.1,11.6,col=coll[1],border=NA)
text(-16,17.6,"LINE-enriched genes",cex=1,srt=90,font=1)
rect(-15.3,11.4,-15.1,-0.4,col=coll[13],border=NA)
text(-16,5.5,"SINE-enriched genes",cex=1,srt=90,font=1)
dev.off()

}






Figure_L_GO <- function(gl){
  
  pv <- c(-log10(gl$pvalue[-c(1:5,7:9,11:12,14:15)]))
  type <- c(gl$Term[-c(1:5,7:9,11:12,14:15)])
  leve <- c(gl$Level[-c(1:5,7:9,11:12,14:15)])
  count <- c(gl$Count[-c(1:5,7:9,11:12,14:15)])
  coll <- c(rep("#FF34B3",length(count)))
  id <- c(gl$Description[-c(1:5,7:9,11:12,14:15)])
  coll1 <- paste(coll,"50",sep="")
  
  pdf("Figure_L_go_en.pdf",height=8.6,width=10.5)
  par(mar=c(0,0,0,0),oma=c(0,0,0,0),bty="n")
  plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(-34.2,16.7),ylim=c(-6,57.5),
       xaxt="n",yaxt="n",xaxs="i", yaxs="i")
  for(i in 1:length(pv)){
    rect(0,57-i+0.4,pv[i],57-i-0.4,col=coll[i],border="#878787")
    capitalized_strings <- toupper(substr(id[i], start = 1, stop = 1))
    text(0,57-i,paste0(capitalized_strings, substring(id[i], first = 2,last=100000)),cex=0.8,pos=2)
    rect(-0.3,57-i-0.05,0,57-i+0.05,col="grey",border=NA)
    
    rect(11,57-i+0.45,12,57-i-0.45,col=coll1[i],border="grey")
    text(11.5,57-i+0.02,type[i],cex=0.7)
    rect(12.3,57-i+0.45,13.3,57-i-0.45,col=coll1[i],border="grey")
    text(12.8,57-i+0.02,leve[i],cex=0.7)
    rect(13.6,57-i+0.45,14.6,57-i-0.45,col=coll1[i],border="grey")
    text(14.1,57-i+0.02,count[i],cex=0.7)
    
  }
  text(11.1,-0.8,"GO Type",cex=0.8,srt=-90,pos=4)
  text(12.4,-0.8,"GO Level",cex=0.8,srt=-90,pos=4)
  text(13.7,-0.8,"GO Count",cex=0.8,srt=-90,pos=4)
  text(-7,-2.5,"GO Term",cex=1)
  rect(0,-0.7,10,-0.85,col="grey",border=NA)
  for(i in 0:4){
    rect(i*2,-0.85,i*2+0.1,-1.3,col="grey",border=NA)
    text(i*2+0.05,-1.9,i*2,cex=0.8)
  }
  rect(5*2,-0.85,5*2-0.1,-1.3,col="grey",border=NA)
  text(5*2-0.05,-1.9,10,cex=0.8)
  text(5.5,-3.5,expression(paste(plain("log")[10],"(",italic(P),"-value)",sep="")),cex=1)
  
  #rect(-15.3,23.4,-15.1,11.6,col=coll[1],border=NA)
  #text(-16,17.6,"LINE-enriched genes",cex=1,srt=90,font=1)
  #rect(-15.3,11.4,-15.1,-0.4,col=coll[13],border=NA)
  #text(-16,5.5,"SINE-enriched genes",cex=1,srt=90,font=1)
  dev.off()
  
}


Figure_LS_rep_heatmap <- function(datl,dats){
  
  require(RColorBrewer)
  require(shape)
  
  p3l <- pheatmap(datl,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu"))[3:7])(100),
                  cluster_cols=F,cluster_rows = T,treeheight_row=0,show_colnames = F,show_rownames = F,
                  border_color="grey")
  p3s <- pheatmap(dats,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu"))[4:7])(100),
                  scale="none",cluster_cols=F,cluster_rows = T,treeheight_row=0,show_colnames = F,show_rownames = F,
                  border_color="grey")
  
  datln <- datl[p3l$tree_row$order,]
  datsn <- dats[p3s$tree_row$order,]
  
  llb1 <- c("L1","L2","I-Jockey","RTE-BovB","Unkown")
  ssb1 <- c("tRNA","B2","ID","MIR","Unkown")
  gb1 <- c("Promoter","5'UTR","Exon","Intron","3'UTR","Downstream")
  nl <- dim(datln)[1]
  datls <- rbind(datln,datsn)
  
  coll <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu"))[4:7])(100)
  
  ndc <- cut(datls,100)
  colm <- matrix(coll[match(ndc,levels(ndc))],ncol=30,byrow=F)
  
  datlm <- colm[1:nl,]
  datsm <- colm[-c(1:nl),]
  
  pdf("Figure_LS_heatmap.pdf",height=7.6,width=4.5)
  par(mar=c(0,0,0,0),oma=c(0,0,0,0),bty="n")
  plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(-2,31),ylim=c(-537,6000.1),
       xaxt="n",yaxt="n",xaxs="i", yaxs="i")
  
  for(i in 1:dim(datsm)[1]){
    for(j in 1:30){
      rect(j-0.5,i,j+0.5,i+1,col=datsm[i,j],border=NA)
    }
  }
  for(i in 1:dim(datlm)[1]){
    for(j in 1:30){
      rect(j-0.5,i+4400,j+0.5,i+1+4400,col=datlm[i,j],border=NA)
    }
  }
  for(j in 1:7){
    rect((j-1)*5+1-0.6,1,(j-1)*5+1-0.5,3933,col="grey",border=NA)
    rect((j-1)*5+1-0.6,4401,(j-1)*5+1-0.5,5174,col="grey",border=NA)
  }
  
  for(j in 1:4){
    rect(j+1-0.55,1,j+1-0.5,3933,col="#FF34B3",border=NA)
    rect(j+1-0.55+5,1,j+1-0.5+5,3933,col="#FF34B3",border=NA)
    rect(j+1-0.55+10,1,j+1-0.5+10,3933,col="#FF34B3",border=NA)
    rect(j+1-0.55+15,1,j+1-0.5+15,3933,col="#FF34B3",border=NA)
    rect(j+1-0.55+20,1,j+1-0.5+20,3933,col="#FF34B3",border=NA)
    rect(j+1-0.55+25,1,j+1-0.5+25,3933,col="#FF34B3",border=NA)
    
    rect(j+1-0.55,4401,j+1-0.5,5174,col="#4EEE94",border=NA)
    rect(j+1-0.55+5,4401,j+1-0.5+5,5174,col="#4EEE94",border=NA)
    rect(j+1-0.55+10,4401,j+1-0.5+10,5174,col="#4EEE94",border=NA)
    rect(j+1-0.55+15,4401,j+1-0.5+15,5174,col="#4EEE94",border=NA)
    rect(j+1-0.55+20,4401,j+1-0.5+20,5174,col="#4EEE94",border=NA)
    rect(j+1-0.55+25,4401,j+1-0.5+25,5174,col="#4EEE94",border=NA)
    #rect((j-1)*5+1-0.6,4401,(j-1)*5+1-0.5,5174,col="grey",border=NA)
  }
  for(j in 1:5){
    for(k in 1:6){
      text(j+(k-1)*5,4290,j,cex=0.9)
      text(j+(k-1)*5,-100,j,cex=0.9)
    }
  }
  rect(0.5,3950,30.4,4225,col="#4EEE9450",border="grey")
  rect(0.5,-450,30.4,-175,col="#FF34B350",border="grey")
  poss <- c(3,7,13,20,27)
  for(j in 1:length(poss)){
    text(poss[j],4087.5,paste(c(1:5)[j],": ",llb1[j],sep=""),cex=0.9)
    text(poss[j],-312.5,paste(c(1:5)[j],": ",ssb1[j],sep=""),cex=0.9)
  }
  rect(0.5,5600,5.5,5650,col="black",border=NA)
  rect(5.5,5580,10.5,5670,col="darkgrey",border=NA)
  rect(10.5,5550,16.5,5700,col="darkgrey",border=NA)
  rect(19.5,5550,20.5,5700,col="darkgrey",border=NA)
  rect(20.5,5580,25.5,5670,col="darkgrey",border=NA)
  rect(25.5,5600,30.5,5650,col="black",border=NA)
  rect(5.5,5580,5.7,5850,col="darkgrey",border=NA)
  rect(5.7,5850,7,5820,col="darkgrey",border=NA)
  Arrows(6,5835,7,5835,arr.width=0.25,col="darkgrey",arr.type="triangle",arr.length = 0.3)
  segments(18,5900,16.5,5700,lwd=2,col="darkgrey")
  segments(18,5900,19.5,5700,lwd=2,col="darkgrey")
  posf <- c(3,8,13,18,23,28)
  
  for(i in 1:6){
    text(posf[i],5450,gb1[i],cex=0.9)
    segments(posf[i],5380,posf[i]-2.5,5174,lwd=2,col="grey",lty=2)
    segments(posf[i],5380,posf[i]+2.5,5174,lwd=2,col="grey",lty=2)
  }
  
  text(-0.8,4787.5,"LINE-enriched",srt=90,cex=1)
  text(-0.8,1967,"SINE-enriched",srt=90,cex=1)
  dev.off()
  
}

