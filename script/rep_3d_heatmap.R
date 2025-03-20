snp.sq1 <- function(i, j, color, polygon.par) {
  #if(i == j) return(); # pas de plot d'un SNP avec lui même
  if(j < i) { tmp <- j; j <- i; i <- tmp } # i < j
  d <- (j-i)/2
  cx <- i+d
  cy <- d
  do.call(polygon, c(list(x = cx+c(-1,0,1,0)/2, y = cy + c(0,1,0,-1)/2, col = color), polygon.par))
  #if(write.ld.val) text(cx, cy, ld, cex = cex.ld)
}


sq_heatmap <- function(smat,pc,lab,sab,ii,kp,llb,filepre="Figure_ls_cor_chr"){
  
  library(RColorBrewer)
  n <- dim(smat)[1]
  coll <- colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(200)
  nd <- as.numeric(unlist(smat))
  nd[which(nd>0.75)] <- 0.75
  nd[which(nd< -0.75)] <- -0.75
  ndc <- cut(nd,200)
  colm <- matrix(coll[match(ndc,levels(ndc))],ncol=n,byrow=F)
  
  pdf(paste(filepre,"_",llb[ii],".pdf",sep=""),height=6.4,width=5)
  par(mar=c(0,0,0,0),oma=c(2,0,0,0),fig=c(0,1,0.94,1))
  plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(1,5),ylim=c(1,4),
       xaxt="n",yaxt="n",xaxs="i", yaxs="i",bty="n")
  text(3,2.3,"Correlation matrix",cex=1.2)
  par(mar=c(0.5,3.3,0,1),fig=c(0,1,0.3,0.94),new=TRUE)
  plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(0,n+0.6),ylim=c(0,n+0.6),
       xaxt="n",yaxt="n",xaxs="i", yaxs="i")
  for(i in 1:n){
    for(j in 1:n){
      #rect(j-0.5,i-0.5,j+0.5,i+0.5,col=colm[i,j],border=NA)
      do.call(rect, c(list(xleft=j-0.5,ybottom=i-0.5,xright=j+0.5,ytop=i+0.5,col = colm[i,j],border=NA)))
    }
  }
  par(mar=c(0.5,0,0,0),fig=c(0,0.15,0.3,0.94),new=TRUE)
  plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(1,5),ylim=c(1,4),
       xaxt="n",yaxt="n",xaxs="i", yaxs="i",bty="n")
  sc <- round(seq(-0.7,0.7,0.1),2)
  scl <- seq(1,3.99,length=length(sc))[-c(1,15)]
  for(k in 1:length(scl)){
    if(k%%2!=0){
      rect(2.95,scl[k],3.35,scl[k]+0.01,border=NA,col="black")
      text(2.1,scl[k]+0.005,sc[k+1],cex=0.8)
    }
  }
  
  cs <- seq(1,4,length=201)
  for(k in 1:200){
    rect(3.3,cs[k],4.1,cs[k+1],border = NA,col=coll[k])
  }
  
  labt <- as.numeric(scale(lab$svc[[ii]]))
  sabt <- as.numeric(scale(sab$svc[[ii]]))
  
  par(mar=c(0.5,3.2,0,1),fig=c(0,1,0.2,0.3),new=TRUE)
  plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(0.5,n+0.5),ylim=c(-2.2,2.2),
       xaxt="n",yaxt="n",xaxs="i", yaxs="i")
  for(i in 1:n){
    rect(i-0.5,0,i+0.5,labt[i],col="#4EEE94",border=NA)
  }
  axis(2,c(-2,0,2),c(-2,0,2),las=1,cex.axis=0.8,tck=-0.073,mgp=c(2.5,0.4,0),col="black")
  mtext("LINE size",2,cex=0.8,line=1.5)
  mtext("NORMD",2,cex=0.8,line=2.3)
  par(mar=c(0.5,3.2,0,1),fig=c(0,1,0.1,0.2),new=TRUE)
  plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(0.5,n+0.5),ylim=c(-2.2,2.2),
       xaxt="n",yaxt="n",xaxs="i", yaxs="i")
  for(i in 1:n){
    rect(i-0.5,0,i+0.5,sabt[i],col="#FF34B3",border=NA)
  }
  axis(2,c(-2,0,2),c(-2,0,2),las=1,cex.axis=0.8,tck=-0.073,mgp=c(2.5,0.4,0),col="black")
  mtext("SINE size",2,cex=0.8,line=1.5)
  mtext("NORMD",2,cex=0.8,line=2.3)
  par(mar=c(0.5,3.2,0,1),fig=c(0,1,0,0.1),new=TRUE)
  plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(0.5,n+0.5),ylim=c(-0.042,0.042),
       xaxt="n",yaxt="n",xaxs="i", yaxs="i")
  pct <- pc[which(pc[,1]==paste("chr",llb[ii],sep="")),5]
  for(i in 1:n){
    if(pct[i]==0){
      next
    }
    if(pct[i] >0){
      rect(i-0.5,0,i+0.5,pct[i],col="#CD3333",border=NA) 
    }else{
      rect(i-0.5,0,i+0.5,pct[i],col="#00B2EE",border=NA) 
    }
  }
  if(ii==23||ii==25){
    ll <- seq(1,kp[ii],floor(kp[[ii]])/4-0.5)
  }else{
    ll <- seq(1,kp[ii],floor(kp[[ii]])/4-1)
  }
  mtext("score",2,cex=0.8,line=1.5)
  mtext("PC",2,cex=0.8,line=2.3)
  axis(2,c(-0.04,0,0.04),c(-0.04,0,0.04),las=1,cex.axis=0.8,tck=-0.073,mgp=c(2.5,0.4,0),col="black")
  axis(1,seq(1,n,length=5),ll,las=1,cex.axis=1,tck=-0.093,mgp=c(2.5,0.5,0),col="black")
  mtext(paste("Chromsome",llb[ii]," (Mb)",sep=""),1,cex=1,line=1.4)
  dev.off()
  
}

rep_3d_pro <- function(lab,sab,mat,ii){
  
  labt <- as.numeric(scale(lab$svc[[ii]]))
  sabt <- as.numeric(scale(sab$svc[[ii]]))
  
  i1 <- which(labt>2*sabt)
  i2 <- which(sabt>2*labt)
  
  tmpl <- as.matrix(mat[i1,i1])
  lmm <- tmpl[lower.tri(tmpl)]
  lpi <- which(lmm > 0)
  lni <- which(lmm < 0)
  lp <- lmm[lpi];ln <- lmm[lni]
  lp_pro <- length(lpi)/length(lmm);ln_pro <- length(lni)/length(lmm)
  
  tmps <- as.matrix(mat[i2,i2])
  smm <- tmps[lower.tri(tmps)]
  spi <- which(smm > 0)
  sni <- which(smm < 0)
  sp <- smm[spi];sn <- smm[sni]
  sp_pro <- length(spi)/length(smm);sn_pro <- length(sni)/length(smm)
  
  tmpls <- as.numeric(as.matrix(mat[i1,i2]))
  lspi <- which(tmpls > 0)
  lsni <- which(tmpls < 0)
  lsp <- tmpls[lspi];lsn <- tmpls[lsni]
  lsp_pro <- length(lspi)/length(tmpls);lsn_pro <- length(lsni)/length(tmpls)
  
  return(list(lp=lp,ln=ln,sp=sp,sn=sn,lsp=lsp,lsn=lsn,lp_pro=lp_pro,ln_pro=ln_pro,
              sp_pro=sp_pro,sp_pro=sp_pro,sn_pro=sn_pro,lsp_pro=lsp_pro,lsn_pro=lsn_pro))
}

tad_rep <- function(gf1,gf2){
  
  library(dplyr)
  cl <- paste("chr",c(1:23,"X","Y"),sep="")
  svc <- list()
  svn <- list()
  fullp <- list();partp <- list()
  for(k in 1:length(cl)){
    tmp <- gf1[which(gf1$V1==cl[k]),]
    tmp1 <- gf2[which(gf2$qname==cl[k]),]
    sv <- c()
    sn1 <- c()
    full <- c();part <- c()
    for( i in 1:dim(tmp)[1]){
      sss1 <- which(between(tmp1$qend,tmp[i,2],tmp[i,3])==1)
      sss2 <- which(between(tmp1$qstart,tmp[i,2],tmp[i,3])==1)
      sss <- unique(c(sss1,sss2))
      ssl <- length(which(table(c(sss1,sss2))==2))
      full <- c(full,ssl);part <- c(part,length(sss)-ssl)
      sv <- c(sv,sum(tmp1$qend[sss]-tmp1$qstart[sss]))
      sn1 <- c(sn1,length(sss))
    }
    svc[[k]] <- sv
    svn[[k]] <- sn1
    fullp[[k]] <- full
    partp[[k]] <- part
  }
  return(list(svn=svn,svc=svc,fullp=fullp,partp=partp))
}

loop_rep <- function(gf1,gf2){
  
  library(dplyr)
  cl <- paste("chr",c(1:23,"X","Y"),sep="")
  svc <- list()
  svn <- list()
  fullp <- list();partp <- list()
  for(k in 1:length(cl)){
    tmp <- gf1[which(gf1[,1]==cl[k]),]
    tmp1 <- gf2[which(gf2$qname==cl[k]),]
    sv <- c()
    sn1 <- c()
    full <- c();part <- c()
    for( i in 1:dim(tmp)[1]){
      sss1 <- which(between(tmp1$qend,tmp[i,2],tmp[i,3])==1)
      sss2 <- which(between(tmp1$qstart,tmp[i,2],tmp[i,3])==1)
      sss <- unique(c(sss1,sss2))
      ssl <- length(which(table(c(sss1,sss2))==2))
      full <- c(full,ssl);part <- c(part,length(sss)-ssl)
      sv <- c(sv,sum(tmp1$qend[sss]-tmp1$qstart[sss]))
      sn1 <- c(sn1,length(sss))
    }
    svc[[k]] <- sv
    svn[[k]] <- sn1
    fullp[[k]] <- full
    partp[[k]] <- part
  }
  return(list(svn=svn,svc=svc,fullp=fullp,partp=partp))
}




rep_interaction_boxplot <- function(pall){
  
  pdf("Figure_rep_interboxplot.pdf",width=5,height=4)
  par(oma=c(0,0,0,0),fig=c(0,1,0,1))
  par(mar=c(1.8,3.5,0.5,0.5))
  ax <- boxplot(pall,xaxt="n", yaxt="n",xaxs="i", yaxs="i",ylim=c(0.081,0.93),
                col=paste(c("#4EEE94","#4EEE94","#FF34B3","#FF34B3","#808080","#808080"),"20",sep=""),pars = list(boxwex = 0.7, staplewex = 0.5, outwex = 0.5),
                boxlwd = 2,border = c("#4EEE94","#4EEE94","#FF34B3","#FF34B3","#808080","#808080"))
  
  axis(1,1:6,rep(c("Pos","Neg"),3),las=1,cex.axis=0.9,tck=-0.023,mgp=c(2.5,0.45,0),col="black")
  axis(2,seq(0.1,0.9,0.2),seq(0.1,0.9,0.2),las=1,cex.axis=1,tck=-0.023,mgp=c(2.5,0.5,0),col="black")
  mtext("% of interaction",2,cex=1.2,line=2.2)
  
  
  rect(0.8,0.85,2.2,0.86,border=NA,col="#4EEE94")
  text(1.5,0.885,"LINE-LINE",cex=1)
  rect(2.8,0.85,4.2,0.86,border=NA,col="#FF34B3")
  text(3.5,0.885,"SINE-SINE",cex=1)
  rect(4.8,0.85,6.2,0.86,border=NA,col="#808080")
  text(5.5,0.885,"LINE-SINE",cex=1)
  
  dev.off()
}


sq_commat1_1 <- function(lab,sab,tad,loop,gene1,llb,kp,filepre="Figure_density_tad_loop1_1.pdf"){
  
  
  i1 <- rev(seq(0,1,length=8))
  i2 <- list()
  k <- 1
  for(i in 1:(length(i1)-1)){
    td <- seq(i1[i],i1[i+1],length=7)
    pii <- c()
    pii1 <- c()
    if(i%%2!=0){
      for(j in 1:(length(td)-1)){
        pii <- rbind(pii,c(0,0.5,td[j],td[j+1]))
      }
      i2[[k]] <- pii
      k <- k + 1
      for(j in 1:(length(td)-1)){
        pii1 <- rbind(pii1,c(0.5,1,td[j],td[j+1]))
      }
      i2[[k]] <- pii1
      k <- k + 1
    }
    if(i%%2==0){
      for(j in 1:(length(td)-1)){
        pii <- rbind(pii,c(0,0.5,td[j],td[j+1]))
      }
      i2[[k]] <- pii
      k <- k + 1
      for(j in 1:(length(td)-1)){
        pii1 <- rbind(pii1,c(0.5,1,td[j],td[j+1]))
      }
      i2[[k]] <- pii1
      k <- k + 1
    }
  }
  
  pdf(filepre,height=29.7,width=24)
  par(oma=c(0,0.5,4.5,0.5))
  for(ii in 1:13){
    chri <- paste("chr",llb[ii],sep="")
    labt <- as.numeric(scale(lab$svc[[ii]]))
    sabt <- as.numeric(scale(sab$svc[[ii]]))
    tadt <- tad[which(tad$V1==chri),]
    loopt <- loop[which(loop$V1==chri),]
    genet <- gene1[which(gene1$V1==chri),]
    if(ii==1){
      par(fig=c(i2[[ii]][1,1],i2[[ii]][1,2],i2[[ii]][1,4],i2[[ii]][1,3]),mar=c(0.5,1,0,1))#1
      n <- length(labt)
      plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(0.5,n+0.5),ylim=c(-2.2,2.2),
           xaxt="n",yaxt="n",xaxs="i", yaxs="i")
      for(i in 1:n){
        rect(i-0.5,0,i+0.5,labt[i],col="#4EEE94",border=NA)
      }
      mtext(paste("Chromsome",llb[ii],sep=""),3,cex=3.2,line=1.2)
      par(fig=c(i2[[ii]][2,1],i2[[ii]][2,2],i2[[ii]][2,4],i2[[ii]][2,3]),mar=c(0.5,1,0,1),new=TRUE)#2
      n <- length(sabt)
      plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(0.5,n+0.5),ylim=c(-2.2,2.2),
           xaxt="n",yaxt="n",xaxs="i", yaxs="i")
      for(i in 1:n){
        rect(i-0.5,0,i+0.5,sabt[i],col="#FF34B3",border=NA)
      }
      par(fig=c(i2[[ii]][3,1],i2[[ii]][3,2],i2[[ii]][3,4],i2[[ii]][3,3]),mar=c(0.5,1,0,1),new=TRUE)#3
      n <- dim(tadt)[1]
      plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(-0.1,kp[ii]+0.5),ylim=c(-2.2,2.2),
           xaxt="n",yaxt="n",xaxs="i", yaxs="i")
      for(i in 1:n){
        rect(tadt[i,2]/10^6,-2.19,tadt[i,3]/10^6,2.19,col="#FFA500",border=NA)
      }
      par(fig=c(i2[[ii]][4,1],i2[[ii]][4,2],i2[[ii]][4,4],i2[[ii]][4,3]),mar=c(0.5,1,0,1),new=TRUE)#4
      n <- dim(loopt)[1]
      plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(-0.1,kp[ii]+0.5),ylim=c(-2.2,2.2),
           xaxt="n",yaxt="n",xaxs="i", yaxs="i")
      for(i in 1:n){
        rect(loopt[i,3]/10^6,-2.19,loopt[i,5]/10^6,2.19,col="#00F5FF",border=NA)
      }
      par(fig=c(i2[[ii]][5,1],i2[[ii]][5,2],i2[[ii]][5,4],i2[[ii]][5,3]),mar=c(0.5,1,0,1),new=TRUE)#5
      n <- dim(genet)[1]
      plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(-0.1,kp[ii]+0.5),ylim=c(-2.2,2.2),
           xaxt="n",yaxt="n",xaxs="i", yaxs="i")
      for(i in 1:n){
        rect(genet[i,4]/10^6,-2.19,genet[i,5]/10^6,2.19,col="grey",border=NA)
      }
    }else{
      par(fig=c(i2[[ii]][1,1],i2[[ii]][1,2],i2[[ii]][1,4],i2[[ii]][1,3]),mar=c(0.5,1,0,1),new=TRUE)#1
      n <- length(labt)
      plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(0.5,n+0.5),ylim=c(-2.2,2.2),
           xaxt="n",yaxt="n",xaxs="i", yaxs="i")
      for(i in 1:n){
        rect(i-0.5,0,i+0.5,labt[i],col="#4EEE94",border=NA)
      }
      mtext(paste("Chromsome",llb[ii],sep=""),3,cex=3.2,line=1.2)
      par(fig=c(i2[[ii]][2,1],i2[[ii]][2,2],i2[[ii]][2,4],i2[[ii]][2,3]),mar=c(0.5,1,0,1),new=TRUE)#2
      n <- length(sabt)
      plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(0.5,n+0.5),ylim=c(-2.2,2.2),
           xaxt="n",yaxt="n",xaxs="i", yaxs="i")
      for(i in 1:n){
        rect(i-0.5,0,i+0.5,sabt[i],col="#FF34B3",border=NA)
      }
      par(fig=c(i2[[ii]][3,1],i2[[ii]][3,2],i2[[ii]][3,4],i2[[ii]][3,3]),mar=c(0.5,1,0,1),new=TRUE)#3
      n <- dim(tadt)[1]
      plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(-0.1,kp[ii]+0.5),ylim=c(-2.2,2.2),
           xaxt="n",yaxt="n",xaxs="i", yaxs="i")
      for(i in 1:n){
        rect(tadt[i,2]/10^6,-3,tadt[i,3]/10^6,3,col="#FFA500",border=NA)
      }
      par(fig=c(i2[[ii]][4,1],i2[[ii]][4,2],i2[[ii]][4,4],i2[[ii]][4,3]),mar=c(0.5,1,0,1),new=TRUE)#4
      n <- dim(loopt)[1]
      plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(-0.1,kp[ii]+0.5),ylim=c(-2.2,2.2),
           xaxt="n",yaxt="n",xaxs="i", yaxs="i")
      for(i in 1:n){
        rect(loopt[i,3]/10^6,-2.19,loopt[i,5]/10^6,2.19,col="#00F5FF",border=NA)
      }
      par(fig=c(i2[[ii]][5,1],i2[[ii]][5,2],i2[[ii]][5,4],i2[[ii]][5,3]),mar=c(0.5,1,0,1),new=TRUE)#5
      n <- dim(genet)[1]
      plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(-0.1,kp[ii]+0.5),ylim=c(-2.2,2.2),
           xaxt="n",yaxt="n",xaxs="i", yaxs="i")
      for(i in 1:n){
        rect(genet[i,4]/10^6,-2.19,genet[i,5]/10^6,2.19,col="grey",border=NA)
      }
    }
  } 
  par(fig=c(i2[[14]][1,1],i2[[14]][1,2],i2[[14]][1,4],i2[[14]][1,3]),mar=c(0.5,1,0,1),new=TRUE,bty="n")#1
  plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(0,10+0.5),ylim=c(-2.2,2.2),
       xaxt="n",yaxt="n",xaxs="i", yaxs="i",frame=F)
  rect(0,-0.5,3,0.5,col="#4EEE94",border=NA)
  text(3.5,0,"NORMD LINE size",adj=0,cex=3)
  par(fig=c(i2[[14]][2,1],i2[[14]][2,2],i2[[14]][2,4],i2[[14]][2,3]),mar=c(0.5,1,0,1),new=TRUE,bty="n")#2
  plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(0,10+0.5),ylim=c(-2.2,2.2),
       xaxt="n",yaxt="n",xaxs="i", yaxs="i",frame=F)
  rect(0,-0.5,3,0.5,col="#FF34B3",border=NA)
  text(3.5,0,"NORMD SINE size",adj=0,cex=3)
  par(fig=c(i2[[14]][3,1],i2[[14]][3,2],i2[[14]][3,4],i2[[14]][3,3]),mar=c(0.5,1,0,1),new=TRUE,bty="n")#3
  plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(0,10+0.5),ylim=c(-2.2,2.2),
       xaxt="n",yaxt="n",xaxs="i", yaxs="i",frame=F)
  rect(0,-0.5,3,0.5,col="#FFA500",border=NA)
  text(3.5,0,"TAD",adj=0,cex=3)
  par(fig=c(i2[[14]][4,1],i2[[14]][4,2],i2[[14]][4,4],i2[[14]][4,3]),mar=c(0.5,1,0,1),new=TRUE,bty="n")#4
  plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(0,10+0.5),ylim=c(-2.2,2.2),
       xaxt="n",yaxt="n",xaxs="i", yaxs="i",frame=F)
  rect(0,-0.5,3,0.5,col="#00F5FF",border=NA)
  text(3.5,0,"Loop",adj=0,cex=3)
  par(fig=c(i2[[14]][5,1],i2[[14]][5,2],i2[[14]][5,4],i2[[14]][5,3]),mar=c(0.5,1,0,1),new=TRUE,bty="n")#5
  plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(0,10+0.5),ylim=c(-2.2,2.2),
       xaxt="n",yaxt="n",xaxs="i", yaxs="i",frame=F)
  rect(0,-0.5,3,0.5,col="grey",border=NA)
  text(3.5,0,"Gene",adj=0,cex=3)
  dev.off()
}

sq_commat1_2 <- function(lab,sab,tad,loop,gene1,llb,kp,filepre="Figure_density_tad_loop1_2.pdf"){
  
  
  i1 <- rev(seq(0,1,length=7))
  i2 <- list()
  k <- 1
  for(i in 1:(length(i1)-1)){
    td <- seq(i1[i],i1[i+1],length=7)
    pii <- c()
    pii1 <- c()
    if(i%%2!=0){
      for(j in 1:(length(td)-1)){
        pii <- rbind(pii,c(0,0.5,td[j],td[j+1]))
      }
      i2[[k]] <- pii
      k <- k + 1
      for(j in 1:(length(td)-1)){
        pii1 <- rbind(pii1,c(0.5,1,td[j],td[j+1]))
      }
      i2[[k]] <- pii1
      k <- k + 1
    }
    if(i%%2==0){
      for(j in 1:(length(td)-1)){
        pii <- rbind(pii,c(0,0.5,td[j],td[j+1]))
      }
      i2[[k]] <- pii
      k <- k + 1
      for(j in 1:(length(td)-1)){
        pii1 <- rbind(pii1,c(0.5,1,td[j],td[j+1]))
      }
      i2[[k]] <- pii1
      k <- k + 1
    }
  }
  
  pdf(filepre,height=27.41538,width=24)
  par(oma=c(0,0.5,4.5,0.5))
  for(ii in 1:12){
    chri <- paste("chr",llb[ii+13],sep="")
    labt <- as.numeric(scale(lab$svc[[ii+13]]))
    sabt <- as.numeric(scale(sab$svc[[ii+13]]))
    tadt <- tad[which(tad$V1==chri),]
    loopt <- loop[which(loop$V1==chri),]
    genet <- gene1[which(gene1$V1==chri),]
    if(ii==1){
      par(fig=c(i2[[ii]][1,1],i2[[ii]][1,2],i2[[ii]][1,4],i2[[ii]][1,3]),mar=c(0.5,1,0,1))#1
      n <- length(labt)
      plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(0.5,n+0.5),ylim=c(-2.2,2.2),
           xaxt="n",yaxt="n",xaxs="i", yaxs="i")
      for(i in 1:n){
        rect(i-0.5,0,i+0.5,labt[i],col="#4EEE94",border=NA)
      }
      mtext(paste("Chromsome",llb[ii],sep=""),3,cex=3.2,line=1.2)
      par(fig=c(i2[[ii]][2,1],i2[[ii]][2,2],i2[[ii]][2,4],i2[[ii]][2,3]),mar=c(0.5,1,0,1),new=TRUE)#2
      n <- length(sabt)
      plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(0.5,n+0.5),ylim=c(-2.2,2.2),
           xaxt="n",yaxt="n",xaxs="i", yaxs="i")
      for(i in 1:n){
        rect(i-0.5,0,i+0.5,sabt[i],col="#FF34B3",border=NA)
      }
      par(fig=c(i2[[ii]][3,1],i2[[ii]][3,2],i2[[ii]][3,4],i2[[ii]][3,3]),mar=c(0.5,1,0,1),new=TRUE)#3
      n <- dim(tadt)[1]
      plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(-0.1,kp[ii+13]+0.5),ylim=c(-2.2,2.2),
           xaxt="n",yaxt="n",xaxs="i", yaxs="i")
      for(i in 1:n){
        rect(tadt[i,2]/10^6,-2.19,tadt[i,3]/10^6,2.19,col="#FFA500",border=NA)
      }
      par(fig=c(i2[[ii]][4,1],i2[[ii]][4,2],i2[[ii]][4,4],i2[[ii]][4,3]),mar=c(0.5,1,0,1),new=TRUE)#4
      n <- dim(loopt)[1]
      plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(-0.1,kp[ii+13]+0.5),ylim=c(-2.2,2.2),
           xaxt="n",yaxt="n",xaxs="i", yaxs="i")
      for(i in 1:n){
        rect(loopt[i,3]/10^6,-2.19,loopt[i,5]/10^6,2.19,col="#00F5FF",border=NA)
      }
      par(fig=c(i2[[ii]][5,1],i2[[ii]][5,2],i2[[ii]][5,4],i2[[ii]][5,3]),mar=c(0.5,1,0,1),new=TRUE)#5
      n <- dim(genet)[1]
      plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(-0.1,kp[ii+13]+0.5),ylim=c(-2.2,2.2),
           xaxt="n",yaxt="n",xaxs="i", yaxs="i")
      for(i in 1:n){
        rect(genet[i,4]/10^6,-2.19,genet[i,5]/10^6,2.19,col="grey",border=NA)
      }
    }else{
      par(fig=c(i2[[ii]][1,1],i2[[ii]][1,2],i2[[ii]][1,4],i2[[ii]][1,3]),mar=c(0.5,1,0,1),new=TRUE)#1
      n <- length(labt)
      plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(0.5,n+0.5),ylim=c(-2.2,2.2),
           xaxt="n",yaxt="n",xaxs="i", yaxs="i")
      for(i in 1:n){
        rect(i-0.5,0,i+0.5,labt[i],col="#4EEE94",border=NA)
      }
      mtext(paste("Chromsome",llb[ii],sep=""),3,cex=3.2,line=1.2)
      par(fig=c(i2[[ii]][2,1],i2[[ii]][2,2],i2[[ii]][2,4],i2[[ii]][2,3]),mar=c(0.5,1,0,1),new=TRUE)#2
      n <- length(sabt)
      plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(0.5,n+0.5),ylim=c(-2.2,2.2),
           xaxt="n",yaxt="n",xaxs="i", yaxs="i")
      for(i in 1:n){
        rect(i-0.5,0,i+0.5,sabt[i],col="#FF34B3",border=NA)
      }
      par(fig=c(i2[[ii]][3,1],i2[[ii]][3,2],i2[[ii]][3,4],i2[[ii]][3,3]),mar=c(0.5,1,0,1),new=TRUE)#3
      n <- dim(tadt)[1]
      plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(-0.1,kp[ii+13]+0.5),ylim=c(-2.2,2.2),
           xaxt="n",yaxt="n",xaxs="i", yaxs="i")
      for(i in 1:n){
        rect(tadt[i,2]/10^6,-3,tadt[i,3]/10^6,3,col="#FFA500",border=NA)
      }
      par(fig=c(i2[[ii]][4,1],i2[[ii]][4,2],i2[[ii]][4,4],i2[[ii]][4,3]),mar=c(0.5,1,0,1),new=TRUE)#4
      n <- dim(loopt)[1]
      plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(-0.1,kp[ii+13]+0.5),ylim=c(-2.2,2.2),
           xaxt="n",yaxt="n",xaxs="i", yaxs="i")
      for(i in 1:n){
        rect(loopt[i,3]/10^6,-2.19,loopt[i,5]/10^6,2.19,col="#00F5FF",border=NA)
      }
      par(fig=c(i2[[ii]][5,1],i2[[ii]][5,2],i2[[ii]][5,4],i2[[ii]][5,3]),mar=c(0.5,1,0,1),new=TRUE)#5
      n <- dim(genet)[1]
      plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(-0.1,kp[ii+13]+0.5),ylim=c(-2.2,2.2),
           xaxt="n",yaxt="n",xaxs="i", yaxs="i")
      for(i in 1:n){
        rect(genet[i,4]/10^6,-2.19,genet[i,5]/10^6,2.19,col="grey",border=NA)
      }
    }
  }
  dev.off()
}



Figure_tad_loop_pro <- function(ltadst,stadst,llopst1d,soopst1d,lsum,ssum){
  
  tadlsum <- sum(unlist(ltadst$svc))/10^6/lsum
  tadssum <- sum(unlist(stadst$svc))/10^6/ssum
  
  
  plsum <- sum(unlist(lloopst1d$svc))/10^6/lsum
  pssum <- sum(unlist(sloopst1d$svc))/10^6/ssum
  
  require(shape)
  pdf("Figure_tad_loop_pro.pdf",width=3.5,height=3.5)
  par(oma=c(0,0,0,0),fig=c(0,0.5,0,1))
  par(mar=c(2,3,0.5,0.5))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(0.5,2.5),ylim=c(-0.1,21),xlab="",ylab="")
  
  rect(0.6,0,1.4,tadlsum*100,col="#4EEE94",border = NA)
  rect(1.6,0,2.4,tadssum*100,col="#FF34B3",border = NA)
  
  axis(2,seq(0,20,5),seq(0,20,5),las=1,cex.axis=1,tck=-0.053,mgp=c(2.5,0.6,0),col="black")
  axis(1,c(1,2),c("LINE","SINE"),las=1,cex.axis=1,tck=-0.043,mgp=c(2.5,0.35,0),col="black")
  mtext("% of TAD",2,cex=1.2,line=1.7)
  par(oma=c(0,0,0,0),fig=c(0.5,1,0,1),new=TRUE)
  par(mar=c(2,3,0.5,0.5))
  plot(NA,NA,xaxt="n", yaxt="n",xaxs="i", yaxs="i",xlim=c(0.5,2.5),ylim=c(-0.1,21),xlab="",ylab="")
  
  rect(0.6,0,1.4,plsum*100,col="#4EEE94",border = NA)
  rect(1.6,0,2.4,pssum*100,col="#FF34B3",border = NA)
  
  axis(2,seq(0,20,5),seq(0,20,5),las=1,cex.axis=1,tck=-0.053,mgp=c(2.5,0.6,0),col="black")
  axis(1,c(1,2),c("LINE","SINE"),las=1,cex.axis=1,tck=-0.043,mgp=c(2.5,0.35,0),col="black")
  mtext("% of loop",2,cex=1.2,line=1.7)
  
  dev.off()
}



sq_heatmap_interval_tran <- function(smat,lab,sab,pc,ii,tad,loop,gene1,interval=c(10,40),xti=c(10,20,30,40),llb,
                                     filepre="Figure_hetamap_tran",polygon.par = list(border = NA)){
  
  library(RColorBrewer)
  n <- dim(smat)[1]
  coll <- colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(200)
  nd <- as.numeric(unlist(smat))
  nd[which(nd>0.75)] <- 0.75
  nd[which(nd< -0.75)] <- -0.75
  ndc <- cut(nd,200)
  colm <- matrix(coll[match(ndc,levels(ndc))],ncol=n,byrow=F)
  
  si <- interval*10^6/100000
  comi <- colm[(si[1]:si[2]),(si[1]:si[2])]
  
  n <- dim(comi)[1]
  max.dist <- Inf
  depth <-  n
  positions <- 1:n
  
  graph.depth <- 0
  for (i in seq(1, n - 1)) for (j in seq(i + 1, min(n, i + 
                                                    depth))) {
    if (positions[j] - positions[i] > max.dist) 
      next
    graph.depth <- max(graph.depth, (j - i)/2)
  }
  
  
  pdf(paste(filepre,"_",llb[ii],"_",interval[1],"-",interval[2],"mb",".pdf",sep=""),height=6,width=5.5)
  par(oma=c(2.5,1.3,0,0.8))
  par(mar=c(0,0,0,0),fig=c(0,1,0.5,1))
  plot(0,0,pch=16,col="#757575",xlab=" ",ylab=" ",xlim=c(12,n-11),ylim=c(0,graph.depth+1),
       type = "n",xaxt = "n", yaxt = "n", bty = "n")
  
  for (i in seq(1, n)) for (j in seq(i, min(n, i + depth))) {
    snp.sq1(i, j,comi[i,j],polygon.par)
  }
  #text(2.6*n/3,graph.depth/1.7,paste("chr",llb[ii],": ",interval[1],"-",interval[2],"mb",sep=""),cex=1)
  labt <- as.numeric(scale(lab$svc[[ii]]))[(si[1]:si[2])]
  sabt <- as.numeric(scale(sab$svc[[ii]]))[(si[1]:si[2])]
  par(mar=c(0,0,0,0),fig=c(0,1,0.41666,0.51),new=TRUE)
  plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(0,n+1),ylim=c(-2.2,2.2),
       xaxt="n",yaxt="n",xaxs="i", yaxs="i")
  for(i in 1:n){
    rect(i-0.5,0,i+0.5,labt[i],col="#4EEE94",border=NA)
  }
  #axis(2,c(-2,0,2),c(-2,0,2),las=1,cex.axis=1,tck=-0.073,mgp=c(2.5,0.5,0),col="black")
  mtext("LINE",2,cex=0.8,line=0.3)
  par(mar=c(0,0,0,0),fig=c(0,1,0.3333333,0.41666),new=TRUE)
  plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(0.5,n+0.5),ylim=c(-2.2,2.2),
       xaxt="n",yaxt="n",xaxs="i", yaxs="i")
  for(i in 1:n){
    rect(i-0.5,0,i+0.5,sabt[i],col="#FF34B3",border=NA)
  }
  #axis(2,c(-2,0,2),c(-2,0,2),las=1,cex.axis=1,tck=-0.073,mgp=c(2.5,0.5,0),col="black")
  mtext("SINE",2,cex=0.8,line=0.3)
  par(mar=c(0,0,0,0),fig=c(0,1,0.25,0.3333333),new=TRUE)
  plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(0.5,n+0.5),ylim=c(-0.042,0.042),
       xaxt="n",yaxt="n",xaxs="i", yaxs="i")
  pct <- pc[which(pc[,1]==paste("chr",llb[ii],sep="")),5][(si[1]:si[2])]
  for(i in 1:n){
    if(pct[i]==0){
      next
    }
    if(pct[i] >0){
      rect(i-0.5,0,i+0.5,pct[i],col="#CD3333",border=NA) 
    }else{
      rect(i-0.5,0,i+0.5,pct[i],col="#00B2EE",border=NA) 
    }
  }
  
  mtext("PC score",2,cex=0.8,line=0.3)
  
  dd <- interval*10^6
  par(mar=c(0,0,0,0),fig=c(0,1,0.1666667,0.25),new=TRUE)
  plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=dd,ylim=c(-2,2),
       xaxt="n",yaxt="n",xaxs="i", yaxs="i")
  tadt <- tad[which(tad$V1==paste("chr",llb[ii],sep="")),]
  tadts <- tadt[which(between(tadt$V3,dd[1],dd[2])==1),]
  ntad <- dim(tadts)[1]
  for(i in 1:ntad){
    rect(tadts[i,2],-2.19,tadts[i,3],2.19,col="#FFA500",border=NA)
  }
  
  mtext("TAD",2,cex=0.8,line=0.3)
  
  par(mar=c(0,0,0,0),fig=c(0,1,0.08333333,0.1666667),new=TRUE)
  plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=dd,ylim=c(-2,2),
       xaxt="n",yaxt="n",xaxs="i", yaxs="i")
  loopt <- loop[which(loop$V1==paste("chr",llb[ii],sep="")),]
  loopts <- loopt[which(between(loop$V3,dd[1],dd[2])==1),]
  nloop <- dim(loopts)[1]
  for(i in 1:n){
    rect(loopts[i,3],-2.19,loopts[i,5],2.19,col="#00F5FF",border=NA)
  }
  
  mtext("loop",2,cex=0.8,line=0.3)
  
  par(mar=c(0,0,0,0),fig=c(0,1,0,0.08333333),new=TRUE)
  plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=dd,ylim=c(-2,2),
       xaxt="n",yaxt="n",xaxs="i", yaxs="i")
  genet <- gene1[which(gene1$V1==paste("chr",llb[ii],sep="")),]
  genets <- genet[which(between(genet$V5,dd[1],dd[2])==1),]
  ngene <- dim(genets)[1]
  for(i in 1:ngene){
    rect(genets[i,4],-2.19,genets[i,5],2.19,col="grey",border=NA)
  }
  
  mtext("gene",2,cex=0.8,line=0.3)
  axis(1,seq(dd[1],dd[2],length=length(xti)),xti,las=1,cex.axis=1,tck=-0.13,mgp=c(2.5,0.5,0),col="black")
  mtext(paste("Chromsome",llb[ii]," (Mb)",sep=""),1,cex=1,line=1.4)
  dev.off()
  
  # 
  # mtext(paste("Chromsome",llb[ii],":",interval[1],"-",interval[2],"Mb",sep=""),1,cex=1,line=1.4)
  # dev.off()
  
}


