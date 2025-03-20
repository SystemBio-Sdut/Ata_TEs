

AB_rep_each <- function(gf1,gf2){
  
  id <- unique(gf2$tclass)
  allt <- list()
  for(i in 1:length(id)){
    tout <- gf2[which(gf2$tclass==id[i]),]
    allt[[i]] <- AB_rep(gf1=gf1,gf2=tout)
  }
  return(list(id=id,allt=allt))
}

lsnorm <- function(lab,sab){
  
  dl <- lab$ar/c(lab$ar+lab$br)
  dl1 <- lab$br/c(lab$ar+lab$br)
  dll <- c(dl,dl1)
  dlz <- (dll-mean(dll))/sd(dll)
  
  dlz[which(dlz>1.5)] <- 1.5
  dlz[which(dlz< -1.5)] <- -1.5
  
  dlza <- dlz[1:25];dlzb <- dlz[-c(1:25)];
  
  ds <- sab$ar/c(sab$ar+sab$br)
  ds1 <- sab$br/c(sab$ar+sab$br)
  dss <- c(ds,ds1)
  dsz <- (dss-mean(dss))/sd(dss)
  
  dsz[which(dsz>1.5)] <- 1.5
  dsz[which(dsz< -1.5)] <- -1.5
  
  dsza <- dsz[1:25];dszb <- dsz[-c(1:25)];
  
  return(list(dlza=dlza,dlzb=dlzb,dsza=dsza,dszb=dszb))
}





AB_rep <- function(gf1,gf2){
  
  library(dplyr)
  cl <- paste("chr",c(1:23,"X","Y"),sep="")
  ar <- br <- as <- bs <-c()
  svc <- list()
  svn <- list()
  for(k in 1:length(cl)){
    tmp <- gf1[which(gf1$V1==cl[k]),]
    tmp1 <- gf2[which(gf2$qname==cl[k]),]
    sv <- c()
    sn1 <- c()
    for( i in 1:dim(tmp)[1]){
      sss <- which(between(tmp1$qend,tmp[i,2],tmp[i,3])==1)
      sv <- c(sv,sum(tmp1$qend[sss]-tmp1$qstart[sss]))
      sn1 <- c(sn1,length(sss))
    }
    svc[[k]] <- sv
    svn[[k]] <- sn1
    ia <- which(tmp$V4=="A")
    ib <- which(tmp$V4=="B")
    as <- c(as,sum(tmp$V3[ia]-tmp$V2[ia])/10^6)
    bs <- c(bs,sum(tmp$V3[ib]-tmp$V2[ib])/10^6)
    ar <- c(ar,sum(sv[ia])/10^6)
    br <- c(br,sum(sv[ib])/10^6)
  }
  return(list(ar=ar,br=br,svn=svn,svc=svc,as=as,bs=bs))
}




Figure_AB_rep_pro <- function(lab,sab){
  
  m1 <- list(a=lab$ar/lab$as,b=lab$br/lab$bs,c=sab$ar/sab$as,d=sab$br/sab$bs)
  
  pdf("Figure_AB_rep_pro.pdf",width=3,height=4)
  par(oma=c(0,0,0,0),fig=c(0,1,0,1))
  par(mar=c(1.5,3.5,0.5,0.5))
  ax <- boxplot(m1,xaxt="n", yaxt="n",xaxs="i", yaxs="i",ylim=c(0.081,0.53),
                col=paste(c("#CD3333","#00B2EE","#CD3333","#00B2EE"),"20",sep=""),pars = list(boxwex = 0.7, staplewex = 0.5, outwex = 0.5),
                boxlwd = 2,border = c("#CD3333","#00B2EE","#CD3333","#00B2EE"))
  rect(0.98,0.44,1.02,0.49,col="grey",border=NA)
  rect(1.98,0.464,2.02,0.49,col="grey",border=NA)
  rect(0.98,0.488,2.02,0.49,col="grey",border=NA)
  text(1.5,0.5,"2.3e-08",cex=1)
  rect(2.98,0.407,3.02,0.47,col="grey",border=NA)
  rect(3.98,0.276,4.02,0.47,col="grey",border=NA)
  rect(2.98,0.468,4.02,0.47,col="grey",border=NA)
  text(3.5,0.48,"6.5e-10",cex=1)
  
  axis(1,1:4,c("LINE","LINE","SINE","SINE"),las=1,cex.axis=0.9,tck=-0.033,mgp=c(2.5,0.45,0),col="black")
  axis(2,seq(0.1,0.5,0.1),seq(0.1,0.5,0.1),las=1,cex.axis=1,tck=-0.033,mgp=c(2.5,0.6,0),col="black")
  mtext("% of compartment",2,cex=1.2,line=2.2)
  coll <- c("#CD3333","#00B2EE")
  ll <- c("A","B")
  iii <- seq(1.3,3.2,length=2)
  for(i in 1:length(iii)){
    rect(iii[i],0.512,iii[i]+0.2,0.522,border=NA,col=coll[i])
    text(iii[i]+0.4,0.518,ll[i],cex=1)
  }
  
  dev.off()
}





Figure_rep_AB_pro <- function(lab,sab){
  
  m1 <- c(sum(lab$ar)/sum((lab$ar+lab$br)),sum(lab$br)/sum((lab$ar+lab$br)),
          sum(sab$ar)/sum((sab$ar+sab$br)),sum(sab$br)/sum((sab$ar+sab$br)))
  m11 <- cumsum(m1[1:2])
  m22 <- cumsum(m1[3:4])
  
  pdf("Figure_rep_AB_pro.pdf",width=2,height=4)
  par(oma=c(0,0,0,0),fig=c(0,1,0,1))
  par(mar=c(2.5,4.5,0.5,0.5))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(0.45,2.51),ylim=c(-0.001,1.1),xlab="",ylab="")
  
  rect(1-0.3,0,1+0.3,m11[1],col="#CD3333",border=NA)
  rect(1-0.3,m11[1],1+0.3,m11[2],col="#00B2EE",border=NA)
  
  rect(2-0.3,0,2+0.3,m22[1],col="#CD3333",border=NA)
  rect(2-0.3,m22[1],2+0.3,m22[2],col="#00B2EE",border=NA)
  
  axis(1,1:2,c("LINE","SINE"),las=1,cex.axis=0.9,tck=-0.043,mgp=c(2.5,0.45,0),col="black")
  axis(2,seq(0,1,0.2),seq(0,1,0.2),las=1,cex.axis=1,tck=-0.043,mgp=c(2.5,0.6,0),col="black")
  mtext("% of repeats",2,cex=1.2,line=3)
  mtext("% of in A or B compartments",2,cex=1.2,line=2)
  coll <- c("#CD3333","#00B2EE")
  ll <- c("A","B")
  iii <- seq(0.8,1.8,length=2)
  for(i in 1:length(iii)){
    rect(iii[i],1.03,iii[i]+0.2,1.06,border=NA,col=coll[i])
    text(iii[i]+0.4,1.0455,ll[i],cex=1)
  }
  
  dev.off()
}










Figure_rep_AB_chr <- function(lschr){
  
  require(shape)
  pdf("Figure_rep_AB_chr.pdf",width=5,height=3)
  par(oma=c(0,0,0,0),fig=c(0,1,0,1))
  par(mar=c(1.5,4,0.5,0.5))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(-2,55),ylim=c(-2.2,2.2),xlab="",ylab="")
  
  rect(-5,-2.5,26.5,3,col="#CD333320",border = NA)
  rect(26.5,-2.5,56,3,col="#00B2EE20",border = NA)
  
  points(1:25,lschr$dlza,col="#4EEE94",pch=19,cex=1.2)
  lines(1:25,lschr$dlza,col="#4EEE9450",lwd=2)
  
  points(28:52,lschr$dlzb,col="#4EEE94",pch=19,cex=1.2)
  lines(28:52,lschr$dlzb,col="#4EEE9450",lwd=2)
  
  points(1:25,lschr$dsza,col="#FF34B3",pch=19,cex=1.2)
  lines(1:25,lschr$dsza,col="#FF34B350",lwd=2)
  
  points(28:52,lschr$dszb,col="#FF34B3",pch=19,cex=1.2)
  lines(28:52,lschr$dszb,col="#FF34B350",lwd=2)
  
  text(1,-1.9,"Chr1",cex=0.8);text(24,-1.9,"ChrY",cex=0.8)
  Arrows(4, -1.9, 20.5, -1.9,col="grey",lwd=3,arr.type="triangle",arr.width = 0.2,arr.length=0.2)
  text(2+27,-1.9,"Chr1",cex=0.8);text(24+27,-1.9,"ChrY",cex=0.8)
  Arrows(5+27, -1.9, 20.5+27, -1.9,col="grey",lwd=3,arr.type="triangle",arr.width = 0.2,arr.length=0.2)
  segments(-100,0,100,0,col="grey",lwd=2)
  axis(2,seq(-2,2,0.5),seq(-2,2,0.5),las=1,cex.axis=1,tck=-0.033,mgp=c(2.5,0.7,0),col="black")
  axis(1,c(12.5,28+12.5),c("A","B"),las=1,cex.axis=1.2,tck=0,mgp=c(2.5,0.35,0),col="black")
  mtext("Relative content",2,cex=1.2,line=2.5)
  
  coll <- c("#4EEE94","#FF34B3")
  ll <- c("LINE","SINE")
  iii <- c(10.5,38.5)
  for(i in 1:length(iii)){
    points(iii[i],1.95,pch=19,cex=1.6,col=coll[i])
    text(iii[i]+3.8,1.98,ll[i],cex=1)
  }
  Arrows(19.5, 1.64, 21.5, 1.54,col="grey",lwd=1.5,arr.type="triangle",arr.width = 0.2,arr.length=0.2)
  text(16.7,1.64,"Chr22",cex=0.7)
  Arrows(28.5, 1.54, 26.5, 1.34,col="grey",lwd=1.5,arr.type="triangle",arr.width = 0.2,arr.length=0.2)
  text(31,1.54,"ChrY",cex=0.7)
  dev.off()
}





Figure_rep_ABeach_pro_line <- function(labe,sabe){
  
  ii1 <- unique(labe$id)[c(1,4,2,6,5,7,8,3)]
  ii1n <- c("L1","L2","I-Jockey","CR1","RTE-X","RTE-BovB","Dong-R4","Unknown")
  lta <- c();ltb <- c()
  for(i in 1:length(ii1)){
    ii1i <- which(ii1[i]==labe$id)
    tmpa <- labe$allt[[ii1i]]$ar/labe$allt[[ii1i]]$as
    tmpb <- labe$allt[[ii1i]]$br/labe$allt[[ii1i]]$bs
    lta <- cbind(lta,tmpa)
    ltb <- cbind(ltb,tmpb)
  }
  
  ii2 <- unique(sabe$id)[c(1:3,5,7,8,9,10,6)]
  ii2n <- c("tRNA","B2","ID","MIR","Alu","5S-Deu-L2","tRNA-RTE","tRNA-Deu","Unknown")
  sta <- c();stb <- c()
  for(i in 1:length(ii2)){
    ii2i <- which(ii2[i]==sabe$id)
    tmpa <- sabe$allt[[ii2i]]$ar/sabe$allt[[ii2i]]$as
    tmpb <- sabe$allt[[ii2i]]$br/sabe$allt[[ii2i]]$bs
    sta <- cbind(sta,tmpa)
    stb <- cbind(stb,tmpb)
  }
  
  labee <- list()
  for(i in 1:8){
    labee[[i]] <- cbind(lta[,i],ltb[,i])*100
  }
  
  sabee <- list()
  for(i in 1:9){
    sabee[[i]] <- cbind(sta[,i],stb[,i])*100
  }
  
  pdf("Figure_ABeach_rep_pro_line.pdf",width=10,height=3)
  liney <- rbind(c(0,50),c(0,0.25),c(0,0.24),c(0,0.03),c(0,0.012),c(0,0.50),c(0,0.004),c(0,0.8))
  ll1 <- c(10,0.05,0.06,0.01,0.003,0.1,0.001,0.2)
  par(oma=c(0,2,1,0))
  par(fig=c(0,1,0,1),mfrow=c(1,8))
  par(mar=c(2,2.5,1,0.5))
  for(i in 1:8){
    ax <- boxplot(labee[[i]],xaxt="n", yaxt="n",xaxs="i", yaxs="i",ylim=c(0,max(labee[[i]],na.rm=T)*1.15),
                  col=paste(c("#CD3333","#00B2EE","#CD3333","#00B2EE"),"20",sep=""),pars = list(boxwex = 0.7, staplewex = 0.5, outwex = 0.5),
                  boxlwd = 2,border = c("#CD3333","#00B2EE","#CD3333","#00B2EE"))
    axis(1,1.5,ii1n[i],las=1,cex.axis=1.2,tck=0,mgp=c(2.5,0.5,0),col="black",font=2)
    axis(2,seq(liney[i,1],liney[i,2],ll1[i]),seq(liney[i,1],liney[i,2],ll1[i]),las=1,cex.axis=1.05,tck=-0.043,mgp=c(2.5,0.6,0),col="black")
    
    if(i==1){
      mtext("% of compartment",2,cex=1,line=2.2)
    }
  }
  dev.off()
  
  pdf("Figure_ABeach_rep_pro_sine.pdf",width=11.2,height=3)
  
  liney <- rbind(c(0,40),c(0,0.8),c(0,0.3),c(0,0.4),c(0,0.004),c(0,0.004),c(0,0.003),c(0,0.002),c(0,0.016))
  ll1 <- c(10,0.2,0.1,0.1,0.001,0.001,0.001,0.001,0.004)
  par(oma=c(0,2,1,0))
  par(fig=c(0,1,0,1),mfrow=c(1,9))
  par(mar=c(2,2.5,1,0.5))
  for(i in 1:9){
    ax <- boxplot(sabee[[i]],xaxt="n", yaxt="n",xaxs="i", yaxs="i",ylim=c(0,max(sabee[[i]],na.rm=T)*1.15),
                  col=paste(c("#CD3333","#00B2EE","#CD3333","#00B2EE"),"20",sep=""),pars = list(boxwex = 0.7, staplewex = 0.5, outwex = 0.5),
                  boxlwd = 2,border = c("#CD3333","#00B2EE","#CD3333","#00B2EE"))
    axis(1,1.5,ii2n[i],las=1,cex.axis=1.2,tck=0,mgp=c(2.5,0.5,0),col="black",font=2)
    axis(2,seq(liney[i,1],liney[i,2],ll1[i]),seq(liney[i,1],liney[i,2],ll1[i]),las=1,cex.axis=1.05,tck=-0.043,mgp=c(2.5,0.6,0),col="black")
    
    if(i==1){
      mtext("% of compartment",2,cex=1,line=2.2)
    }
  }
  dev.off()
}