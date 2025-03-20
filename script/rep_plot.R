Figure_chr_size <- function(gf1,gf2,kp){
  
  nchr <- length(kp[,1])
  
  g1s <- c()
  g2s <- c()
  for(i in 1:nchr){
    i1 <- which(gf1$qname==kp[i,1])
    i2 <- which(gf2$qname==kp[i,1])
    d1 <- sum(gf1$qend[i1]-gf1$qstart[i1])/kp[i,2]
    d2 <- sum(gf2$qend[i1]-gf2$qstart[i1])/kp[i,2]
    g1s <- c(g1s,d1)
    g2s <- c(g2s,d2)
  }
  
  gs <- g1s + g2s
  
  c1m <- matrix(c(g1s*kp[,2]/10^6,sum((g1s*kp[,2]))/sum(kp[,2])*kp[,2]/10^6),nrow=2,byrow=T)
  c2m <- matrix(c(g2s*kp[,2]/10^6,sum((g2s*kp[,2]))/sum(kp[,2])*kp[,2]/10^6),nrow=2,byrow=T)
  #c1p <- chisq.test(c1m)$p.value
  #c2p <- chisq.test(c2m)$p.value
  
  pdf("Figure-pro_chr_sl.pdf",width=6,height=3)
  par(oma=c(0,0,0,0),fig=c(0,1,0,1))
  par(mar=c(2.5,3.5,0.5,0.5))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(0.3,100.5),ylim=c(-0.001,0.62),xlab="",ylab="")
  for(i in 1:nchr){
    ii <- (i-1)*4+1
    rect(ii,0,ii+1,g1s[i],border=NA,col="#4EEE94")
    rect(ii+1,0,ii+2,g2s[i],border=NA,col="#FF34B3")
    rect(ii+2,0,ii+3,g1s[i]+g2s[i],border=NA,col="#87CEFF")
  }
  axis(1,seq(2,98,4)[seq(1,23,2)],seq(1,25,1)[seq(1,23,2)],las=1,cex.axis=0.9,tck=-0.033,mgp=c(2.5,0.3,0),col="black")
  axis(1,seq(2,98,4)[seq(2,23,2)],seq(1,25,1)[seq(2,23,2)],las=1,cex.axis=0.9,tck=-0.033,mgp=c(2.5,0.3,0),col="black")
  axis(1,seq(2,98,4)[c(24)],c("X"),las=1,cex.axis=0.9,tck=-0.033,mgp=c(2.5,0.3,0),col="black")
  axis(1,seq(2,98,4)[c(25)],c("Y"),las=1,cex.axis=0.9,tck=-0.033,mgp=c(2.5,0.3,0),col="black")
  axis(2,seq(0,0.6,0.2),seq(0,60,20),las=1,cex.axis=1,tck=-0.033,mgp=c(2.5,0.6,0),col="black")
  
  mtext("Chromosome",1,cex=1.2,line=1.4)
  mtext("Proportion (%)",2,cex=1.2,line=2.1)
  coll <- c("#4EEE94","#FF34B3","#87CEFF")
  ll <- c("LINE","SINE","LINE+SINE")
  iii <- seq(15,65,length=3)
  for(i in 1:length(iii)){
    rect(iii[i],0.58,iii[i]+5,0.6,border=NA,col=coll[i])
    text(iii[i]+5,0.59,ll[i],pos=4)
  }
  
  dev.off()
  return(list(linep=g1s,sinep=g2s))
}




Figure_rep_family_LINE <- function(gf,fi){
  
  st <- sort(table(gf$tname),decreasing = T)
  
  
  st1 <- st[which(st>1000)]
  nl <- length(st1)
  snn <- names(st1)
  
  stl <- fi[match(snn,fi$ID),2]
  
  sts <- c()
  
  for(i in 1:nl){
    ii <- which(snn[i]==gf$tname)
    sts <- c(sts,sum(gf$qend[ii]-gf$qstart[ii])/10^6)
  }
  
  pdf("Figure-LINE_family.pdf",width=11,height=9)
  par(oma=c(3,10.2,1.5,0.5),fig=c(0,1/3,0,1))
  par(mar=c(0,0.7,0,0))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(-1,87),ylim=c(-0.001,141),xlab="",ylab="")
  for(i in 1:nl){
    rect(0,141-i-0.4,st1[i]/10^3,141-i+0.4,border=NA,col="#FFBBFF")
    if(i%%2==1){
      axis(2,seq(1,nl,1)[i],rev(snn[seq(1,nl,1)])[i],las=1,cex.axis=0.6,tck=-0.033,mgp=c(2.5,0.6,0),col="black")
    }else{
      axis(2,seq(1,nl,1)[i],rev(snn[seq(1,nl,1)])[i],las=1,cex.axis=0.6,tck=-0.333,mgp=c(2.5,5,0),col="black")
    }
    
  }
  
  #axis(1,seq(0,80,20),seq(0,80,20),las=1,cex.axis=1,tck=-0.033,mgp=c(2.5,0.5,0),col="black")
  axis(1,seq(0,80,20),c(0,expression(paste(20,"×", 10^3,sep="")),
                        expression(paste(40,"×", 10^3,sep="")),expression(paste(60,"×", 10^3,sep="")),
                        expression(paste(80,"×", 10^3,sep=""))),las=1,cex.axis=1,tck=-0.033,mgp=c(2.5,0.55,0),col="black")
  mtext("Copy number",1,cex=1.2,line=1.7)
  mtext("Number of consensus sequences",2,cex=1.2,line=9.3)
  mtext("A",3,cex=1.3,line=0,adj=0)
  
  par(fig=c(1/3,2/3,0,1),new=TRUE)
  par(mar=c(0,0.7,0,0))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(-0.1,18.5),ylim=c(-0.001,141),xlab="",ylab="")
  for(i in 1:nl){
    rect(0,141-i-0.4,stl[i]/10^3,141-i+0.4,border=NA,col="#FF6347")
  }
  
  axis(2,seq(1,nl,1),rep("",nl),las=1,cex.axis=0.6,tck=-0.033,mgp=c(2.5,0.6,0),col="black")
  axis(1,seq(0,16,4),c(0,expression(paste(4,"×", 10^3,sep="")),
                       expression(paste(8,"×", 10^3,sep="")),expression(paste(12,"×", 10^3,sep="")),
                       expression(paste(16,"×", 10^3,sep=""))),las=1,cex.axis=1,tck=-0.033,mgp=c(2.5,0.55,0),col="black")
  #axis(1,seq(0,16,4)+1.7,rep(expression(paste("×", 10^3,sep="")),length(seq(0,16,4))),las=1,cex.axis=1,tck=0,mgp=c(2.5,0.55,0),col="black")
  
  mtext("Length (bp)",1,cex=1.2,line=1.7)
  mtext("B",3,cex=1.3,line=0,adj=0)
  
  par(fig=c(2/3,1,0,1),new=TRUE)
  par(mar=c(0,0.7,0,0))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(-1,105),ylim=c(-0.001,141),xlab="",ylab="")
  for(i in 1:nl){
    rect(0,141-i-0.4,sts[i],141-i+0.4,border=NA,col="#EE9A49")
  }
  axis(2,seq(1,nl,1),rep("",nl),las=1,cex.axis=0.6,tck=-0.033,mgp=c(2.5,0.6,0),col="black")
  axis(1,seq(0,100,20),seq(0,100,20),las=1,cex.axis=1,tck=-0.033,mgp=c(2.5,0.5,0),col="black")
  mtext("Total length in genome (Mb)",1,cex=1.2,line=1.7)
  mtext("C",3,cex=1.3,line=0,adj=0)
  
  dev.off()
  
  
}




Figure_rep_family_SINE <- function(gf,fi){
  
  st <- sort(table(gf$tname),decreasing = T)
  
  
  st1 <- st[which(st>1000)]
  nl <- length(st1)
  snn <- names(st1)
  
  stl <- fi[match(snn,fi$ID),2]
  
  sts <- c()
  
  for(i in 1:nl){
    ii <- which(snn[i]==gf$tname)
    sts <- c(sts,sum(gf$qend[ii]-gf$qstart[ii])/10^6)
  }
  
  pdf("Figure-SINE_family.pdf",width=11,height=7.95)
  par(oma=c(3,10.2,1.5,0.5),fig=c(0,1/3,0,1))
  par(mar=c(0,0.7,0,0))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(-1,267),ylim=c(-0.001,123),xlab="",ylab="")
  for(i in 1:nl){
    rect(0,123-i-0.4,st1[i]/10^3,123-i+0.4,border=NA,col="#FFBBFF")
    if(i%%2==1){
      axis(2,seq(1,nl,1)[i],rev(snn[seq(1,nl,1)])[i],las=1,cex.axis=0.6,tck=-0.033,mgp=c(2.5,0.6,0),col="black")
    }else{
      axis(2,seq(1,nl,1)[i],rev(snn[seq(1,nl,1)])[i],las=1,cex.axis=0.6,tck=-0.333,mgp=c(2.5,5,0),col="black")
    }
    
  }
  
  #axis(1,seq(0,80,20),seq(0,80,20),las=1,cex.axis=1,tck=-0.033,mgp=c(2.5,0.5,0),col="black")
  axis(1,seq(0,250,50),c(0,expression(paste(5,"×", 10^4,sep="")),
                        expression(paste(10,"×", 10^4,sep="")),expression(paste(15,"×", 10^4,sep="")),
                        expression(paste(20,"×", 10^4,sep="")),
                        expression(paste(25,"×", 10^4,sep=""))),las=1,cex.axis=1,tck=-0.033,mgp=c(2.5,0.55,0),col="black")
  mtext("Copy number",1,cex=1.2,line=1.7)
  mtext("Number of consensus sequences",2,cex=1.2,line=9.2)
  mtext("A",3,cex=1.3,line=0,adj=0)
  
  par(fig=c(1/3,2/3,0,1),new=TRUE)
  par(mar=c(0,0.7,0,0))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(-0.1,8.5),ylim=c(-0.001,123),xlab="",ylab="")
  for(i in 1:nl){
    rect(0,123-i-0.4,stl[i]/10^3,123-i+0.4,border=NA,col="#FF6347")
  }
  
  axis(2,seq(1,nl,1),rep("",nl),las=1,cex.axis=0.6,tck=-0.033,mgp=c(2.5,0.6,0),col="black")
  axis(1,seq(0,8,2),c(0,expression(paste(2,"×", 10^3,sep="")),
                       expression(paste(4,"×", 10^3,sep="")),expression(paste(6,"×", 10^3,sep="")),
                       expression(paste(8,"×", 10^3,sep=""))),las=1,cex.axis=1,tck=-0.033,mgp=c(2.5,0.55,0),col="black")
  #axis(1,seq(0,16,4)+1.7,rep(expression(paste("×", 10^3,sep="")),length(seq(0,16,4))),las=1,cex.axis=1,tck=0,mgp=c(2.5,0.55,0),col="black")
  
  mtext("Length (bp)",1,cex=1.2,line=1.7)
  mtext("B",3,cex=1.3,line=0,adj=0)
  
  par(fig=c(2/3,1,0,1),new=TRUE)
  par(mar=c(0,0.7,0,0))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(-1,65),ylim=c(-0.001,123),xlab="",ylab="")
  for(i in 1:nl){
    rect(0,123-i-0.4,sts[i],123-i+0.4,border=NA,col="#EE9A49")
  }
  axis(2,seq(1,nl,1),rep("",nl),las=1,cex.axis=0.6,tck=-0.033,mgp=c(2.5,0.6,0),col="black")
  axis(1,seq(0,100,20),seq(0,100,20),las=1,cex.axis=1,tck=-0.033,mgp=c(2.5,0.5,0),col="black")
  mtext("Total length in genome (Mb)",1,cex=1.2,line=1.7)
  mtext("C",3,cex=1.3,line=0,adj=0)
  
  dev.off()
  
  
}



Figure_chr_number <- function(gf1,gf2,kp){
  
  nchr <- length(kp[,1])
  
  afl <- c();afld <- c();aflh <- c()
  afs <- c();afsd <- c();afsh <- c()
  for(i in 1:nchr){
    ii1 <- unique(gf1$tname[which(gf1$qname==kp[i,1])])
    afl <- c(afl,length(ii1))
    ii11 <- which(unlist(strsplit(ii1,"-"))=="rnd")
    afld <- c(afld,length(ii11))
    aflh <- c(aflh,length(ii1)-length(ii11))
    ii2 <- unique(gf2$tname[which(gf2$qname==kp[i,1])])
    afs <- c(afs,length(ii2))
    ii22 <- which(unlist(strsplit(ii2,"-"))=="rnd")
    afsd <- c(afsd,length(ii22))
    afsh <- c(afsh,length(ii2)-length(ii22))
  }
  
  af <- afl+afs
  
  
  pdf("Figure-number_chr_sl.pdf",width=7,height=8)
  par(oma=c(2,0,0.5,0.5),fig=c(0,1,2/3,1))
  par(mar=c(0.8,4,0,0))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(0.3,100.5),ylim=c(-0.001,760),xlab="",ylab="")
  for(i in 1:nchr){
    ii <- (i-1)*4+1
    rect(ii,0,ii+1,afl[i],border=NA,col="#4EEE94")
    rect(ii+1,0,ii+2,afs[i],border=NA,col="#FF34B3")
    rect(ii+2,0,ii+3,af[i],border=NA,col="#87CEFF")
  }
  axis(1,seq(2,98,4),rep(" ",25),las=1,cex.axis=0.9,tck=-0.033,mgp=c(2.5,0.3,0),col="black")
  axis(2,seq(0,600,200),seq(0,600,200),las=1,cex.axis=1,tck=-0.033,mgp=c(2.5,0.6,0),col="black")
  #mtext("Number of subfamilies",2,cex=1.1,line=2.4)
  
  coll <- c("#4EEE94","#FF34B3","#87CEFF")
  ll <- c("LINE","SINE","LINE+SINE")
  iii <- seq(15,65,length=3)
  for(i in 1:length(iii)){
    rect(iii[i],680,iii[i]+5,700,border=NA,col=coll[i])
    text(iii[i]+5,690,ll[i],pos=4)
  }
  
  par(fig=c(0,1,1/3,2/3),new=TRUE)
  par(mar=c(0.8,4,0,0))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(0.3,100.5),ylim=c(-0.001,520),xlab="",ylab="")
  for(i in 1:nchr){
    ii <- (i-1)*4+1
    rect(ii,0,ii+1,afld[i],border="#4EEE94",col="#4EEE94",density=30,angle=30)
    rect(ii+1,0,ii+2,aflh[i],border="#4EEE94",col="#4EEE94",density=30,angle=-30)
    rect(ii+2,0,ii+3,afld[i]+aflh[i],border=NA,col="#4EEE94")
  }
  axis(1,seq(2,98,4),rep(" ",25),las=1,cex.axis=0.9,tck=-0.033,mgp=c(2.5,0.3,0),col="black")
  axis(2,seq(0,600,200),seq(0,600,200),las=1,cex.axis=1,tck=-0.033,mgp=c(2.5,0.6,0),col="black")
  mtext("Number of consensus sequences",2,cex=1.1,line=2.4)
  
  ll <- c(expression(italic(denovo)),"Homologous","Total")
  dd <- c(30,30,NULL)
  ag <- c(30,-30,30)
  br <- c("#4EEE94","#4EEE94",NA)
  iii <- seq(15,65,length=3)
  for(i in 1:length(iii)){
    rect(iii[i],470,iii[i]+5,485,border=br[i],col="#4EEE94",density=dd[i],angle=ag[i])
    text(iii[i]+5,477.5,ll[i],pos=4)
  }
  
  par(fig=c(0,1,0,1/3),new=TRUE)
  par(mar=c(0.8,4,0,0))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(0.3,100.5),ylim=c(-0.001,520),xlab="",ylab="")
  for(i in 1:nchr){
    ii <- (i-1)*4+1
    rect(ii,0,ii+1,afsd[i],border="#FF34B3",col="#FF34B3",density=30,angle=30)
    rect(ii+1,0,ii+2,afsh[i],border="#FF34B3",col="#FF34B3",density=30,angle=-30)
    rect(ii+2,0,ii+3,afsd[i]+afsh[i],border=NA,col="#FF34B3")
  }
  axis(1,seq(2,98,4),rep(" ",25),las=1,cex.axis=0.9,tck=-0.033,mgp=c(2.5,0.3,0),col="black")
  axis(2,seq(0,600,200),seq(0,600,200),las=1,cex.axis=1,tck=-0.033,mgp=c(2.5,0.6,0),col="black")
  #mtext("Number of categories",2,cex=1.1,line=2.4)
  
  ll <- c(expression(italic(denovo)),"Homologous","Total")
  dd <- c(30,30,NULL)
  ag <- c(30,-30,30)
  br <- c("#FF34B3","#FF34B3",NA)
  iii <- seq(15,65,length=3)
  for(i in 1:length(iii)){
    rect(iii[i],470,iii[i]+5,485,border=br[i],col="#FF34B3",density=dd[i],angle=ag[i])
    text(iii[i]+5,477.5,ll[i],pos=4)
  }
  
  mtext("Chromosome",1,cex=1.2,line=1.4)
  axis(1,seq(2,98,4)[seq(1,23,2)],seq(1,25,1)[seq(1,23,2)],las=1,cex.axis=0.9,tck=-0.033,mgp=c(2.5,0.4,0),col="black")
  axis(1,seq(2,98,4)[seq(2,23,2)],seq(1,25,1)[seq(2,23,2)],las=1,cex.axis=0.9,tck=-0.033,mgp=c(2.5,0.4,0),col="black")
  axis(1,seq(2,98,4)[c(24)],c("X"),las=1,cex.axis=0.9,tck=-0.033,mgp=c(2.5,0.4,0),col="black")
  axis(1,seq(2,98,4)[c(25)],c("Y"),las=1,cex.axis=0.9,tck=-0.033,mgp=c(2.5,0.4,0),col="black")
  dev.off()
}

Figure_density_sd <- function(d1,d2,kp){
  
  nchr <- length(kp[,1])
  
  g1s <- c()
  g2s <- c()
  g1ss <- c()
  g2ss <- c()
  for(i in 1:nchr){
    i1 <- which(d1$Chr==kp[i,1])
    i2 <- which(d2$Chr==kp[i,1])
    dd1 <- sd(d1$Count[i1])
    dd2 <- sd(d2$Count[i2])
    dd11 <- sd(d1$Size[i1]/1000)
    dd22 <- sd(d2$Size[i2]/1000)
    g1s <- c(g1s,dd1)
    g2s <- c(g2s,dd2)
    g1ss <- c(g1ss,dd11)
    g2ss <- c(g2ss,dd22)
  }
  
  
  pdf("Figure-density_chr_sd.pdf",width=5.5,height=5)
  par(oma=c(2,0,0.5,0),fig=c(0,1,0.5,1))
  par(mar=c(0.5,3.5,0,0.5))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(0.3,75.5),ylim=c(-0.001,95),xlab="",ylab="")
  for(i in 1:nchr){
    ii <- (i-1)*3+1
    rect(ii,0,ii+1,g1s[i],border=NA,col="#4EEE94")
    rect(ii+1,0,ii+2,g2s[i],border=NA,col="#FF34B3")
  }
  axis(1,seq(2,75,3),rep(" ",25),las=1,cex.axis=0.9,tck=-0.033,mgp=c(2.5,0.3,0),col="black")
  axis(2,seq(0,90,30),seq(0,90,30),las=1,cex.axis=1,tck=-0.033,mgp=c(2.5,0.6,0),col="black")
  
  mtext("Sd of count",2,cex=1.2,line=2.1)
  
  coll <- c("#4EEE94","#FF34B3")
  ll <- c("LINE","SINE")
  iii <- seq(15,50,length=2)
  for(i in 1:length(iii)){
    rect(iii[i],88,iii[i]+5,91,border=NA,col=coll[i])
    text(iii[i]+5,89.5,ll[i],pos=4)
  }
  
  par(fig=c(0,1,0,0.5),new=TRUE)
  par(mar=c(0.5,3.5,0,0.5))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(0.3,75.5),ylim=c(-0.001,95),xlab="",ylab="")
  for(i in 1:nchr){
    ii <- (i-1)*3+1
    rect(ii,0,ii+1,g1ss[i],border=NA,col="#4EEE94")
    rect(ii+1,0,ii+2,g2ss[i],border=NA,col="#FF34B3")
  }
  
  coll <- c("#4EEE94","#FF34B3")
  ll <- c("LINE","SINE")
  iii <- seq(15,50,length=2)
  for(i in 1:length(iii)){
    rect(iii[i],88,iii[i]+5,91,border=NA,col=coll[i])
    text(iii[i]+5,89.5,ll[i],pos=4)
  }
  
  axis(1,seq(2,75,3)[seq(1,23,2)],seq(1,25,1)[seq(1,23,2)],las=1,cex.axis=0.9,tck=-0.033,mgp=c(2.5,0.3,0),col="black")
  axis(1,seq(2,75,3)[seq(2,23,2)],seq(1,25,1)[seq(2,23,2)],las=1,cex.axis=0.9,tck=-0.033,mgp=c(2.5,0.3,0),col="black")
  axis(1,seq(2,75,3)[c(24)],c("X"),las=1,cex.axis=0.9,tck=-0.033,mgp=c(2.5,0.3,0),col="black")
  axis(1,seq(2,75,3)[c(25)],c("Y"),las=1,cex.axis=0.9,tck=-0.033,mgp=c(2.5,0.3,0),col="black")
  axis(2,seq(0,90,30),seq(0,90,30),las=1,cex.axis=1,tck=-0.033,mgp=c(2.5,0.6,0),col="black")
  mtext("Sd of content",2,cex=1.2,line=2.1)
  mtext("Chromosome",1,cex=1.2,line=1.4)
  dev.off()
}




Figure_age_pro <- function(K1,aout1,aout2,r){
  
  T1 <- as.numeric(K1[,5])/100/(2*r)/10^6
  T1[which(T1<0)] <- 0
  T1[which(is.na(T1))] <- 0
  
  ii <- match(unique(aout1$tname),K1[,2])
  T11 <- T1[ii]
  l11 <- K1[ii,2]
  
  iii <- match(unique(aout2$tname),K1[,2])
  T22 <- T1[iii]
  l22 <- K1[iii,2]
  
  gss1 <- c()
  for(k1 in 1:length(l11)){
    si1 <- which(l11[k1]==aout1$tname)
    gss1 <- c(gss1,sum(aout1$qend[si1]-aout1$qstart[si1])/10^6)
  }
  
  gss2 <- c()
  for(k2 in 1:length(l22)){
    si2 <- which(l22[k2]==aout2$tname)
    gss2 <- c(gss2,sum(aout2$qend[si2]-aout2$qstart[si2])/10^6)
  }
  
  or1 <- order(T11)
  or2 <- order(T22)
  
  ll111 <- l11[or1]
  ll222 <- l22[or2]
  
  T111 <- T11[or1]
  T222 <- T22[or2]
  gss11 <- gss1[or1]
  gss22 <- gss2[or2]
  
  c1 <- cut(T111,breaks=seq(0,95,5))
  c2 <- cut(T222,breaks=seq(0,95,5))
  cl1 <- levels(c1)
  cl2 <- levels(c2)
  gs11h <- c()
  gs22h <- c()
  for(i in 1:length(cl1)){
    gs11h <- c(gs11h,sum(gss11[which(c1==cl1[i])]))
    gs22h <- c(gs22h,sum(gss22[which(c2==cl1[i])]))
  }
  
  pdf("Figure-age.pdf",width=3,height=3)
  par(oma=c(0,0,0,0),fig=c(0,1,0,1))
  par(mar=c(2.7,3.5,0.5,0.5))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(-0.3,20.5),ylim=c(-0.01,9.62),xlab="",ylab="")
  for(i in 1:length(gs11h)){
    rect(i-0.4,0,i+0.4,gs11h[i]/2780*100,border=NA,col="#4EEE94")
    rect(i-0.4,gs11h[i]/2780*100,i+0.4,(gs11h[i]+gs22h[i])/2780*100,border=NA,col="#FF34B3")
  }
  axis(1,seq(0,20,4),seq(0,100,20),las=1,cex.axis=0.8,tck=-0.033,mgp=c(2.5,0.3,0),col="black")
  axis(2,seq(0,9,3),seq(0,9,3),las=1,cex.axis=1,tck=-0.033,mgp=c(2.5,0.6,0),col="black")
  mtext("Age (Mya)",1,cex=1.2,line=1.4)
  mtext("Percent of genome (%)",2,cex=1.2,line=2.1)
  coll <- c("#4EEE94","#FF34B3")
  ll <- c("LINE","SINE")
  iii <- seq(5,13,length=2)
  for(i in 1:length(iii)){
    rect(iii[i],8.6,iii[i]+1,9,border=NA,col=coll[i])
    text(iii[i]+1,8.8,ll[i],pos=4)
  }
  
  dev.off()
}




Figure_dist_pro_heatmap <- function(m1,m2,filename="Figure_dist_pro.pdf",height=12,width=12,
                                    breaks=seq(0,0.6,0.002),fontsize_col=18){
  
  require(pheatmap)
  require(grid)
  require(gridExtra)
  require(ggplot2)
  require(ggplotify)
  require(patchwork)
  require(RColorBrewer)
  
  newnames1 <- rownames(m1)
  newnames <- c("0-200","200-500","500-1000","1000-3000","2000-5000","5000-10000",">10000")
  annotation_col = data.frame(
    Type = factor(rep(c("LINE", "SINE"), each=7))
    
  )
  m <- cbind(m1,m2)
  colnames(m) <- paste("test",1:14,sep="")
  rownames(annotation_col) = colnames(m)
  
  pp1 <- as.ggplot(pheatmap((m),cluster_rows= F,cluster_cols = F,
                            color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(breaks)),
                            labels_col = c(newnames,newnames),labels_row = as.expression(newnames1),fontsize_row=22,
                            annotation_col=annotation_col,fontsize=18,
                            fontsize_col=fontsize_col,angle_col=45,breaks=breaks,display_numbers = TRUE,fontsize_number=17))+
    theme(plot.margin = margin(l = 6,t=10,r = 32))
  ggsave(pp1,file=filename,height=height,width=width)
  
}


