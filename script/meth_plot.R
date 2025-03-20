density_meth_curve <- function(chrn,dd1){
  
  
  pdf("Figure_meth_density.pdf",width=5,height=4)
  par(oma=c(0,0,1,1),fig=c(0,1,0,1))
  par(mar=c(4,3.5,0,0))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(-0.05,1.05),ylim=c(-0.1,7.2),xlab="",ylab="")
  for(j in 1:length(chrn)){
    lines(dd1[[j]],lwd=1,col="#FFBBFF")
  }
  
  axis(1,seq(0,1,0.2),seq(0,1,0.2),cex.axis=1.5,mgp=c(2.5,1,0),col="black")
  axis(2,seq(1,7,1),seq(1,7,1),cex.axis=1.5,mgp=c(2.5,1,0),col="black",las=2)
  mtext("Methylation level",1,cex=1.6,line=2.5)
  mtext("Density",2,cex=1.6,line=2.1)
  dev.off()
}

density_meth_curve_chr <- function(chrn,dd1){
  
  nl <- length(chrn)
  
  pdf("Figure_meth_density_chr.pdf",width=17,height=14)
  par(oma=c(5,6.2,1.5,0.5),mfrow=c(5,5))
  for(i in 1:nl){
    
    par(mar=c(2,2,0,0))
    plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(-0.05,1.05),ylim=c(-0.1,7.2),xlab="",ylab="")
    lines(dd1[[i]],lwd=2,col="#FFBBFF")
    text(0.5,7.3*0.9,chrn[i],cex=2.8)
    if(i >20){
      axis(1,seq(0,1,0.2),seq(0,1,0.2),cex.axis=2.4,tck=-0.053,mgp=c(2.5,2.1,0),col="black")
    }else{
      axis(1,seq(0,1,0.2),labels = F,cex.axis=2.4,tck=-0.053,mgp=c(2.5,0.3,0),col="black")
    }
    if(i==1||i==6||i==11||i==16||i==21){
      axis(2,seq(1,7,1),seq(1,7,1),cex.axis=2.4,las=1,tck=-0.053,mgp=c(2.5,1.8,0),col="black")
    }else{
      axis(2,seq(1,7,1),labels = F,cex.axis=2.4,tck=-0.053,mgp=c(2.5,0.3,0),col="black")
    }
    if(i==23){
      mtext("Methylation level",1,cex=2.1,line=5.1)
    }
    if(i==11){
      mtext("Density",2,cex=2.1,line=5.1)
    }
  }
  
  dev.off()
}

meth_bin_chr <- function(karyotype,meth_bn){
  
  nc <- dim(karyotype)[1]
  pdf("Figure_meth_bin_plot.pdf",width=11,height=15)
  par(oma=c(7,8,1,1),mfrow=c(25,1))
  for(i in 1:nc){
    par(mar=c(0,0,0,0))
    plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(-2,255),ylim=c(-300,5200.2),xlab="",ylab="",bty="n")
    ii1 <- which(meth_bn[,1]==karyotype[i,1])
    methsub <- meth_bn[ii1,]
    cm <- methsub[,4]
    cm[which(cm>5000)] <- 5000
    for(j in 1:length(ii1)){
      if(cm[j]<4000){
        rect(methsub[j,2]/10^6,0,methsub[j,3]/10^6,cm[j],border=NA,col="grey")
      }else{
        rect(methsub[j,2]/10^6,0,methsub[j,3]/10^6,cm[j],border=NA,col="#EE636370")
      }
    }
    axis(2,c(1000,4000),c(1000,4000),cex.axis=2,mgp=c(2.5,1,0),col="black",las=2)
    if(i==13){
      mtext("Distribution of Methylation sites",2,cex=1.6,line=5.7)
    }
    text(karyotype[i,2]/10^6+10,2500,karyotype[i,1],cex=2)
  }
  axis(1,seq(0,250,50),seq(0,250,50),cex.axis=2,mgp=c(2.5,1.5,0),col="black")
  mtext("Chromsome (Mb)",1,cex=1.6,line=4.5)
  dev.off()
}

rep_meth_cor_plot <- function(meth_bn,lbn,sbn){
  
  pdf("Figure_rep_meth_cor.pdf",width=9,height=4)
  par(oma=c(0,0,1,1),fig=c(0,0.5,0,1))
  par(mar=c(3.5,4,0,0))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(2.4,4.3),ylim=c(1.8,6.2),xlab="",ylab="")
  delt <- unique(c(which(meth_bn[,4]==0),which(lbn[,5]==0),which(meth_bn[,4]>20000)))
  tm <- log10(meth_bn[-delt,4])
  tl <- log10(lbn[-delt,5])
  para <- c(6,-0.8)
  CC <- optim(para,EI.sum,x=tm,y=tl,method="Nelder-Mead",control=list(maxit=21000))
  tCC <- EI(CC$par,seq(2.5,4.2,length=300))
  points(tm,tl,col="#ADADAD50",ylim=c(1.8,6.2))
  lines(seq(2.5,4.2,length=300),tCC,col="#4EEE94",lwd=3)
  axis(1,seq(2.5,4,0.5),c(expression(10^"2.5"),expression(10^"3.0"),expression(10^"3.5"),expression(10^"4.0")),cex.axis=1.2,mgp=c(2.5,1,0),col="black")
  axis(2,seq(2,6,1),c(expression(10^2),expression(10^3),expression(10^4),
                      expression(10^5),expression(10^6)),cex.axis=1.2,mgp=c(2.5,0.8,0),col="black",las=2)
  mtext("Distribution of Methylation site",1,cex=1.2,line=2.2)
  mtext("Distribution of LINE content",2,cex=1.2,line=2.5)
  text(2.8,2.5,bquote(italic(r)~"="~.(-0.6628)),cex=1.2,col="#4EEE94",font=2)
  
  par(fig=c(0.5,1,0,1),new=TRUE)
  par(mar=c(3.5,4,0,0))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(2.4,4.3),ylim=c(1.8,6.2),xlab="",ylab="")
  delt <- unique(c(which(meth_bn[,4]==0),which(sbn[,5]==0),which(meth_bn[,4]>20000)))
  tm <- log10(meth_bn[-delt,4])
  tl <- log10(sbn[-delt,5])
  para <- c(6,-0.8)
  CC <- optim(para,EI.sum,x=tm,y=tl,method="Nelder-Mead",control=list(maxit=21000))
  tCC <- EI(CC$par,seq(2.5,4.2,length=300))
  points(tm,tl,col="#ADADAD50",ylim=c(1.8,6.2))
  lines(seq(2.5,4.2,length=300),tCC,col="#FF34B3",lwd=3)
  axis(1,seq(2.5,4,0.5),c(expression(10^"2.5"),expression(10^"3.0"),expression(10^"3.5"),expression(10^"4.0")),cex.axis=1.2,mgp=c(2.5,1,0),col="black")
  axis(2,seq(2,6,1),c(expression(10^2),expression(10^3),expression(10^4),
                      expression(10^5),expression(10^6)),cex.axis=1.2,mgp=c(2.5,0.8,0),col="black",las=2)
  mtext("Distribution of Methylation site",1,cex=1.2,line=2.2)
  mtext("Distribution of SINE content",2,cex=1.2,line=2.5)
  text(2.8,2.5,bquote(italic(r)~"="~.(0.5509)),cex=1.2,col="#FF34B3",font=2)
  dev.off()
}


rep_meth_cor_plot_chr_line <- function(karyotype,meth_bn,lbn){
  
  chrn <- karyotype[,1]
  nl <- length(chrn)
  
  pdf("Figure_rep_meth_cor_chr_line.pdf",width=17,height=14)
  par(oma=c(5,6.2,1.5,0.5),mfrow=c(5,5))
  for(i in 1:nl){
    
    par(mar=c(2,2,0,0))
    plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(2.4,4.3),ylim=c(1.8,6.2),xlab="",ylab="")
    i1s <- meth_bn[which(meth_bn[,1]==chrn[i]),4]
    i2s <- lbn[which(lbn$Chr==chrn[i]),5]
    del <- unique(c(which(i1s==0),which(i2s==0)))
    if(length(del)==0){
      i1sd <- i1s;i2sd <- i2s
    }else{
      i1sd <- i1s[-del];i2sd <- i2s[-del]
    }
    i1sdl <- log10(i1sd);i2sdl <- log10(i2sd);
    points(i1sdl,i2sdl,col="#ADADAD50")
    para <- c(6,-0.8)
    CC <- optim(para,EI.sum,x=i1sdl,y=i2sdl,method="Nelder-Mead",control=list(maxit=21000))
    tCC <- EI(CC$par,seq(2.5,4.2,length=300))
    lines(seq(2.5,4.2,length=300),tCC,col="#4EEE94",lwd=3)
    text(3.35,5.8,chrn[i],cex=2.8)
    if(i >20){
      axis(1,seq(2.5,4,0.5),c(expression(10^"2.5"),expression(10^"3.0"),expression(10^"3.5"),expression(10^"4.0")),cex.axis=2.4,tck=-0.053,mgp=c(2.5,2.4,0),col="black")
    }else{
      axis(1,seq(2.5,4,0.5),labels = F,cex.axis=2.4,tck=-0.053,mgp=c(2.5,0.3,0),col="black")
    }
    if(i==1||i==6||i==11||i==16||i==21){
      axis(2,seq(2,6,1),c(expression(10^2),expression(10^3),expression(10^4),
                          expression(10^5),expression(10^6)),cex.axis=2.4,las=1,tck=-0.053,mgp=c(2.5,1.5,0),col="black")
    }else{
      axis(2,seq(2,6,1),labels = F,cex.axis=2.4,tck=-0.053,mgp=c(2.5,0.3,0),col="black")
    }
    if(i==23){
      mtext("Distribution of Methylation site",1,cex=2.1,line=5.1)
    }
    if(i==11){
      mtext("Distribution of LINE content",2,cex=2.1,line=5.1)
    }
    r1 <- round(cor(i1sdl,i2sdl),4)
    text(3,2.5,bquote(italic(r)~"="~.(r1)),cex=2.2,col="#4EEE94",font=2)
  }
  
  dev.off()
}


rep_meth_cor_plot_chr_sine <- function(karyotype,meth_bn,lbn){
  
  chrn <- karyotype[,1]
  nl <- length(chrn)
  
  pdf("Figure_rep_meth_cor_chr_sine.pdf",width=17,height=14)
  par(oma=c(5,6.2,1.5,0.5),mfrow=c(5,5))
  for(i in 1:nl){
    
    par(mar=c(2,2,0,0))
    plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(2.4,4.3),ylim=c(1.8,6.2),xlab="",ylab="")
    i1s <- meth_bn[which(meth_bn[,1]==chrn[i]),4]
    i2s <- lbn[which(lbn$Chr==chrn[i]),5]
    del <- unique(c(which(i1s==0),which(i2s==0)))
    if(length(del)==0){
      i1sd <- i1s;i2sd <- i2s
    }else{
      i1sd <- i1s[-del];i2sd <- i2s[-del]
    }
    i1sdl <- log10(i1sd);i2sdl <- log10(i2sd);
    points(i1sdl,i2sdl,col="#ADADAD50")
    para <- c(6,-0.8)
    CC <- optim(para,EI.sum,x=i1sdl,y=i2sdl,method="Nelder-Mead",control=list(maxit=21000))
    tCC <- EI(CC$par,seq(2.5,4.2,length=300))
    lines(seq(2.5,4.2,length=300),tCC,col="#FF34B3",lwd=3)
    text(3.35,5.8,chrn[i],cex=2.8)
    if(i >20){
      axis(1,seq(2.5,4,0.5),c(expression(10^"2.5"),expression(10^"3.0"),expression(10^"3.5"),expression(10^"4.0")),cex.axis=2.4,tck=-0.053,mgp=c(2.5,2.4,0),col="black")
    }else{
      axis(1,seq(2.5,4,0.5),labels = F,cex.axis=2.4,tck=-0.053,mgp=c(2.5,0.3,0),col="black")
    }
    if(i==1||i==6||i==11||i==16||i==21){
      axis(2,seq(2,6,1),c(expression(10^2),expression(10^3),expression(10^4),
                          expression(10^5),expression(10^6)),cex.axis=2.4,las=1,tck=-0.053,mgp=c(2.5,1.5,0),col="black")
    }else{
      axis(2,seq(2,6,1),labels = F,cex.axis=2.4,tck=-0.053,mgp=c(2.5,0.3,0),col="black")
    }
    if(i==23){
      mtext("Distribution of Methylation site",1,cex=2.1,line=5.1)
    }
    if(i==11){
      mtext("Distribution of SINE content",2,cex=2.1,line=5.1)
    }
    r1 <- round(cor(i1sdl,i2sdl),4)
    text(3,2.5,bquote(italic(r)~"="~.(r1)),cex=2.2,col="#FF34B3",font=2)
  }
  
  dev.off()
}


Figure_rep_meth_barplot <- function(a2,aout1_mp,aout2_mp){
  
  
  pdf("Figure_rep_meth_barplot.pdf",width=5,height=4)
  par(oma=c(0,0,1,0),fig=c(0,0.5,0,1))
  par(mar=c(2,4,0,1))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(0.5,2.5),ylim=c(-1,32),xlab="",ylab="")
  rect(1-0.3,0,1+0.3,sum(aout1_mp)/dim(a2)[1]*100,col="#4EEE94",border=NA)
  rect(2-0.3,0,2+0.3,sum(aout2_mp)/dim(a2)[1]*100,col="#FF34B3",border=NA)
  axis(1,c(1,2),c("LINE","SINE"),cex.axis=1.3,mgp=c(2.5,0.8,0),col="black")
  axis(2,seq(0,30,10),seq(0,30,10),cex.axis=1.5,mgp=c(2.5,0.8,0),col="black",las=2)
  mtext("Proportion of methylation site (%)",2,cex=1.4,line=2.5)
  
  par(oma=c(0,0,1,0),fig=c(0.5,1,0,1),new=TRUE)
  par(mar=c(2,4.5,0,0.5))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(0.5,2.5),ylim=c(-1,82),xlab="",ylab="")
  rect(1-0.3,0,1+0.3,sum(aout1_mp!=0)/length(aout1_mp)*100,col="#4EEE94",border=NA)
  rect(2-0.3,0,2+0.3,sum(aout2_mp!=0)/length(aout2_mp)*100,col="#FF34B3",border=NA)
  axis(1,c(1,2),c("LINE","SINE"),cex.axis=1.3,mgp=c(2.5,0.8,0),col="black")
  axis(2,seq(0,80,20),seq(0,80,20),cex.axis=1.5,mgp=c(2.5,0.8,0),col="black",las=2)
  mtext("Proportion of methylation",2,cex=1.4,line=4)
  mtext("sequence (%)",2,cex=1.4,line=2.5)
  dev.off()
}



rep_meth_cor_length <- function(a2,aout1,aout2,aout1_mp,aout2_mp){
  
  il <- which(aout1_mp!=0)
  is <- which(aout2_mp!=0)
  lnn <- aout1_mp[il]
  snn <- aout2_mp[is]
  
  ll <- (aout1$qend[il]-aout1$qstart[il])/1000
  sl <- (aout2$qend[is]-aout2$qstart[is])/1000
  
  pdf("Figure_rep_meth_cor_length.pdf",width=9,height=4)
  par(oma=c(0,0,1,1),fig=c(0,0.5,0,1))
  par(mar=c(3.5,4,0,0))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(-0.5,15.4),ylim=c(-5,210),xlab="",ylab="")
  para <- c(2,3)
  CC <- optim(para,EI.sum,x=ll,y=lnn,method="Nelder-Mead",control=list(maxit=21000))
  tCC <- EI(CC$par,seq(min(ll),max(ll),length=300))
  sll <- sample(1:length(ll),10000)
  points(ll[sll],lnn[sll],col="#ADADAD50")
  lines(seq(min(ll),max(ll),length=300),tCC,col="#4EEE94",lwd=3)
  axis(1,seq(0,15,3),seq(0,15,3),cex.axis=1.2,mgp=c(2.5,1,0),col="black")
  axis(2,seq(0,200,50),seq(0,200,50),cex.axis=1.2,mgp=c(2.5,0.8,0),col="black",las=2)
  mtext("Length of methylation LINEs (kb)",1,cex=1.2,line=2.2)
  mtext("Number of methylation sites",2,cex=1.2,line=2.7)
  text(7,210*0.93,bquote(italic(r)~"="~.(0.7723)),cex=1.2,col="#4EEE94",font=2)
  
  par(fig=c(0.5,1,0,1),new=TRUE)
  par(mar=c(3.5,4,0,0))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(-0.5,10.5),ylim=c(-15,610),xlab="",ylab="")
  delt <- unique(c(which(snn>600),which(sl>10)))
  sl1 <- sl[-delt]
  snn1 <- snn[-delt]
  para <- c(2,1)
  CC <- optim(para,EI.sum,x=sl1,y=snn1,method="Nelder-Mead",control=list(maxit=21000))
  tCC <- EI(CC$par,seq(min(sl1),max(sl1),length=300))
  sll <- sample(1:length(sl1),10000)
  points(sl1[sll],snn1[sll],col="#ADADAD50")
  lines(seq(min(sl1),max(sl1),length=300),tCC,col="#FF34B3",lwd=3)
  axis(1,seq(0,10,2),seq(0,10,2),cex.axis=1.2,mgp=c(2.5,1,0),col="black")
  axis(2,seq(0,600,100),seq(0,600,100),cex.axis=1.2,mgp=c(2.5,0.8,0),col="black",las=2)
  mtext("Length of methylation SINEs (kb)",1,cex=1.2,line=2.2)
  mtext("Number of methylation sites",2,cex=1.2,line=2.7)
  text(5,610*0.93,bquote(italic(r)~"="~.(0.5809)),cex=1.2,col="#FF34B3",font=2)
  
  dev.off()
}


Figure_rep_gene_meth_barplot <- function(anno,aout1_mp,aout2_mp){
  
  gli <- which(anno$lan[,1]!="Distal Intergenic")
  gsi <- which(anno$san[,1]!="Distal Intergenic")
  
  glmp <- aout1_mp[gli]
  gsmp <- aout2_mp[gsi]
  
  nglmp <- aout1_mp[-gli]
  ngsmp <- aout2_mp[-gsi]
  
  pdf("Figure_rep_gene_meth_barplot.pdf",width=5,height=4)
  par(oma=c(0,0,1,0),fig=c(0,0.5,0,1))
  par(mar=c(2,4,0,1))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(0.5,2.5),ylim=c(-1,82),xlab="",ylab="")
  rect(1-0.3,0,1+0.3,length(which(glmp!=0))/length(glmp)*100,col="#4EEE94",border=NA)
  rect(2-0.3,0,2+0.3,length(which(gsmp!=0))/length(gsmp)*100,col="#FF34B3",border=NA)
  axis(1,c(1,2),c("LINE","SINE"),cex.axis=1.3,mgp=c(2.5,0.8,0),col="black")
  axis(2,seq(0,80,20),seq(0,80,20),cex.axis=1.5,mgp=c(2.5,0.8,0),col="black",las=2)
  text(1.5,80*0.95,"Neighboring genes",cex=1)
  mtext("Proportion of methylation site (%)",2,cex=1.4,line=2.5)
  
  par(oma=c(0,0,1,0),fig=c(0.5,1,0,1),new=TRUE)
  par(mar=c(2,4,0,1))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(0.5,2.5),ylim=c(-1,82),xlab="",ylab="")
  rect(1-0.3,0,1+0.3,length(which(nglmp!=0))/length(nglmp)*100,col="#4EEE94",border=NA)
  rect(2-0.3,0,2+0.3,length(which(ngsmp!=0))/length(ngsmp)*100,col="#FF34B3",border=NA)
  axis(1,c(1,2),c("LINE","SINE"),cex.axis=1.3,mgp=c(2.5,0.8,0),col="black")
  axis(2,seq(0,80,20),seq(0,80,20),cex.axis=1.5,mgp=c(2.5,0.8,0),col="black",las=2)
  mtext("Proportion of methylation site (%)",2,cex=1.4,line=2.5)
  text(1.5,80*0.95,"Distal intergenic",cex=1)
  dev.off()
}


rep_meth_age_cor <- function(stl1,agel1,ssl1,ages1){
  
  pdf("Figure_rep_meth_age_cor.pdf",width=9,height=4)
  par(oma=c(0,0,1,1),fig=c(0,0.5,0,1))
  par(mar=c(3.5,4,0,0))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(-5,95),ylim=c(-0.02,105),xlab="",ylab="")
  mpf <- stl1[,2]/rowSums(stl1)*100
  points(agel1,mpf,col="#4EEE94")
  para <- c(12,-3)
  CC <- optim(para,EI.sum,x=agel1,y=mpf,method="Nelder-Mead",control=list(maxit=21000))
  tCC <- EI(CC$par,seq(min(agel1),max(agel1),length=300))
  lines(seq(min(agel1),max(agel1),length=300),tCC,col="#4EEE94",lwd=3)
  axis(1,seq(0,80,20),seq(0,80,20),cex.axis=1.2,mgp=c(2.5,1,0),col="black")
  axis(2,seq(0,100,20),seq(0,100,20),cex.axis=1.2,mgp=c(2.5,0.8,0),col="black",las=2)
  mtext("Insert age of LINEs (Mya)",1,cex=1.2,line=2.2)
  mtext("Proportion of methylation sequence (%)",2,cex=1.2,line=2.7)
  text(50,105*0.95,bquote(italic(r)~"="~.(-0.4767)),cex=1.2,col="#4EEE94",font=2)
  
  par(fig=c(0.5,1,0,1),new=TRUE)
  par(mar=c(3.5,4,0,0))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(-5,95),ylim=c(-0.02,105),xlab="",ylab="")
  mpf <- ssl1[,2]/rowSums(ssl1)*100
  points(ages1,mpf,col="#FF34B3")
  para <- c(12,-3)
  CC <- optim(para,EI.sum,x=ages1,y=mpf,method="Nelder-Mead",control=list(maxit=21000))
  tCC <- EI(CC$par,seq(min(ages1),max(ages1),length=300))
  lines(seq(min(ages1),max(ages1),length=300),tCC,col="#FF34B3",lwd=3)
  axis(1,seq(0,80,20),seq(0,80,20),cex.axis=1.2,mgp=c(2.5,1,0),col="black")
  axis(2,seq(0,100,20),seq(0,100,20),cex.axis=1.2,mgp=c(2.5,0.8,0),col="black",las=2)
  mtext("Insert age of SINEs (Mya)",1,cex=1.2,line=2.2)
  mtext("Proportion of methylation sequence (%)",2,cex=1.2,line=2.7)
  text(50,105*0.95,bquote(italic(r)~"="~.(-0.5991)),cex=1.2,col="#FF34B3",font=2)
  dev.off()
}
