ciro.plot1 <- function(ciro,xc=400,yc=400,r=360,xlab=NULL,ylab=NULL,
                       col.out=1,print.chr.lab=T,scale=T,filename="Figure1-1.pdf"){
  
  
  cairo_pdf(filename, 7, 7,family ="Diverda Sans Com")
  par(mar=c(0, 0, 0, 0))
  plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab=xlab, ylab=ylab, main="")
  chr.po <- ciro$chr.po
  chr.num <- nrow(chr.po)
  nnchr <- c(1:23,"X","Y")
  for (chr.i in c(1:chr.num)){
    w1 <- as.numeric(chr.po[chr.i,2])
    w2 <- as.numeric(chr.po[chr.i,3])
    draw.arc.s(xc, yc, r, w1, w2, col="#63B8FF", lwd=7)
    draw.arc.s(xc, yc, r=r-105, w1, w2, col="#63B8FF", lwd=5)
    draw.arc.s(xc, yc, r=r-210, w1, w2, col="#63B8FF", lwd=4)
    draw.arc.s(xc, yc, r=r-315, w1, w2, col="grey", lwd=1)
    
    if (print.chr.lab){
      w.m <- w1+(w2-w1)/2
      r.j <- r/20+30
      chr.t <- chr.po[chr.i,1]
      draw.text.rt(xc, yc, r+r.j, w.m, nnchr[chr.i], cex=1.0,LG=chr.i)
    } 
    if(chr.i!=25){
      if (scale){
        v1 <- as.numeric(chr.po[chr.i,2])
        v2 <- as.numeric(chr.po[chr.i,3])
        v3 <- as.numeric(chr.po[chr.i,6])
        v4 <- as.numeric(chr.po[chr.i,7])
        
        total.num <- as.numeric(chr.po[nrow(chr.po),5])
        do.scale.cir(xc=xc, yc=yc, the.r=r+8, total.num=total.num, 
                     col="#63B8FF", lwd=1,
                     V1=v1, V2=v2, V3=v3, V4=v4)
      }
    }
    draw.arc.s(xc, yc, r=r-7, w1, w2, col="grey", lwd=1)
    draw.arc.s(xc, yc, r=r-112, w1, w2, col="grey", lwd=1)
    draw.arc.s(xc, yc, r=r-217, w1, w2, col="grey", lwd=1)
    
  }
  
  
  r1  <- c(10,115,220)
  r2  <- c(95,200,305)
  collall <- c("#EE7942","#FF34B3","#7A67EE")
  for(col.v in 1:3){
    res <- ciro$L[[col.v]]
    mapping <- as.numeric(res[,4])
    mapping[(which(mapping>400))] <- 400
    poss1 <- as.numeric(res[,2])/10^6
    poss2 <- as.numeric(res[,3])/10^6
    dat.min <- min(mapping, na.rm=T)
    dat.max <- max(mapping, na.rm=T)
    my.R1 <- r-r1[col.v]
    my.R2 <- r-r2[col.v]
    
    for (chr.i in 1:chr.num){
      chr.s <- ciro$chr.n[chr.i] 
      ti   <- which(chr.s==res[,1])
      posti1 <- poss1[ti]
      posti2 <- poss2[ti]
      v1 <- as.numeric(chr.po[chr.i,2])
      v2 <- as.numeric(chr.po[chr.i,3])
      v3 <- as.numeric(chr.po[chr.i,6])
      v4 <- as.numeric(chr.po[chr.i,7])
      nchrp <- length(posti1)
      for (i in 1:nchrp){
        
        my.v1 <- as.numeric(mapping[ti][i])  
        vv1 <- scale.v(my.v1, my.R1, my.R2, dat.min, dat.max)
        po1 <- posti1[i]
        po2 <- posti2[i]
        w1 <- scale.v(po1, v1, v2, v3, v4)
        w2 <- scale.v(po2, v1, v2, v3, v4)
        draw.arc.pg(xc,yc,w1=w1, w2=w2, r1=my.R1, r2=vv1, col=collall[col.v], border=NA, lwd=0.01)
        #draw.line(xc, yc, w=w, l1=vv, l2=vv0, col=collall[col.v], lwd=1, lend=1)
        #draw.point(xc, yc, w1=w, r1=vv, col=collc[[col.v]][ti][i],pch=3,cex=0.4)
      }
    }
  }
  rect(385,415,415,420,border=NA,col="#EE7942")
  text(401,427,"LINE",cex=0.7,col="#EE7942")
  rect(385,392,415,397,border=NA,col="#FF34B3")
  text(401,405,"SINE",cex=0.7,col="#FF34B3")
  rect(385,370,415,375,border=NA,col="#7A67EE")
  text(401,383,"LINE+SINE",cex=0.7,col="#7A67EE")
  dev.off()
}

ciro.plot <- function(ciro,xc=400,yc=400,r=360,xlab=NULL,ylab=NULL,
                      col.out=1,print.chr.lab=T,scale=T,filename="Figure1-1.pdf"){
  
  
  cairo_pdf(filename, 7, 7,family ="Diverda Sans Com")
  par(mar=c(0, 0, 0, 0))
  plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab=xlab, ylab=ylab, main="")
  chr.po <- ciro$chr.po
  chr.num <- nrow(chr.po)
  nnchr <- c(1:23,"X","Y")
  for (chr.i in c(1:chr.num)){
    w1 <- as.numeric(chr.po[chr.i,2])
    w2 <- as.numeric(chr.po[chr.i,3])
    draw.arc.s(xc, yc, r, w1, w2, col="#63B8FF", lwd=7)
    draw.arc.s(xc, yc, r=r-75, w1, w2, col="#63B8FF", lwd=5)
    draw.arc.s(xc, yc, r=r-155, w1, w2, col="#63B8FF", lwd=4)
    draw.arc.s(xc, yc, r=r-230, w1, w2, col="#63B8FF", lwd=4)
    draw.arc.s(xc, yc, r=r-310, w1, w2, col="grey", lwd=2)
    
    if (print.chr.lab){
      w.m <- w1+(w2-w1)/2
      r.j <- r/20+30
      chr.t <- chr.po[chr.i,1]
      draw.text.rt(xc, yc, r+r.j, w.m, nnchr[chr.i], cex=1.0,LG=chr.i)
    } 
    if(chr.i!=25){
      if (scale){
        v1 <- as.numeric(chr.po[chr.i,2])
        v2 <- as.numeric(chr.po[chr.i,3])
        v3 <- as.numeric(chr.po[chr.i,6])
        v4 <- as.numeric(chr.po[chr.i,7])
        
        total.num <- as.numeric(chr.po[nrow(chr.po),5])
        do.scale.cir(xc=xc, yc=yc, the.r=r+8, total.num=total.num, 
                     col="#63B8FF", lwd=1,
                     V1=v1, V2=v2, V3=v3, V4=v4)
      }
    }
    draw.arc.s(xc, yc, r=r-7, w1, w2, col="grey", lwd=1)
    draw.arc.s(xc, yc, r=r-81, w1, w2, col="grey", lwd=1)
    draw.arc.s(xc, yc, r=r-161, w1, w2, col="grey", lwd=1)
    draw.arc.s(xc, yc, r=r-236, w1, w2, col="grey", lwd=1)
    #draw.line3(xc, yc, w1, w2, r1=r-10, r2=r-10, col="blue", lwd=lwd)
  }
  
  pos1 <- ciro$L$Start+(ciro$L$End - ciro$L$Start)/2
  pos2 <- ciro$S$Start+(ciro$S$End - ciro$S$Start)/2
  posss <- cbind(pos1,pos2,pos1,pos2)
  pvall <- cbind(ciro$L$Count,ciro$S$Count,ciro$L$Size,ciro$S$Size)
  
  pvall[which(pvall[,1]>400),1] <- 400
  r1  <- c(10,84,164,239)
  r2  <- c(65,145,220,295)
  collall <- c("#4EEE94","#FF34B3","#4EEE94","#FF34B3")
  for(col.v in 1:4){
    
    mapping <- pvall[,col.v]
    poss <- posss[,col.v]
    dat.min <- min(mapping, na.rm=T)
    dat.max <- max(mapping, na.rm=T)
    my.R1 <- r-r1[col.v]
    my.R2 <- r-r2[col.v]
    
    for (chr.i in 1:chr.num){
      chr.s <- ciro$chr.n[chr.i] 
      ti   <- which(chr.s==ciro$L$Chr)
      posti <- poss[ti]
      v1 <- as.numeric(chr.po[chr.i,2])
      v2 <- as.numeric(chr.po[chr.i,3])
      v3 <- as.numeric(chr.po[chr.i,6])
      v4 <- as.numeric(chr.po[chr.i,7])
      nchrp <- length(posti)
      for (i in 1:(nchrp-1)){
        
        my.v    <- as.numeric(mapping[ti][i])  
        my.v0    <- as.numeric(mapping[ti][i+1])
        vv   <- scale.v(my.v, my.R1, my.R2, dat.min, dat.max)
        vv0   <- scale.v(my.v0, my.R1, my.R2, dat.min, dat.max)
        po      <- posti[i]
        w  <- scale.v(po, v1, v2, v3, v4)
        draw.line(xc, yc, w=w, l1=vv, l2=vv0, col=collall[col.v], lwd=1, lend=1)
        #draw.point(xc, yc, w1=w, r1=vv, col=collc[[col.v]][ti][i],pch=3,cex=0.4)
      }
    }
  }
  
  dev.off()
}






ciro.dat.g <- function(ll,ss,kp){
  
  nchrll <- as.numeric(kp[,2])/10^6
  seg.num    <- length(nchrll)
  
  gap.angle.size <- 1
  seg.angle.from <- 276
  
  seg.full.l  <- sum(nchrll)
  cir.angle.r <- (360-seg.num*gap.angle.size)/seg.full.l
  
  out.s     <- c()
  l.old     <- 0
  gap.angle <- 0
  
  for (i in 1:seg.num){
    seg.n <- paste("Chr",i,sep="")
    len   <- cumsum(nchrll)[i]
    w1    <- cir.angle.r*l.old + gap.angle
    w2    <- cir.angle.r*len   + gap.angle
    out.s     <- rbind(out.s, c(seg.n, w1+seg.angle.from, w2+seg.angle.from, l.old, len, 0, nchrll[i]))
    gap.angle <- gap.angle + gap.angle.size
    l.old     <- len;
  }
  
  ciro <- list(chr.po=out.s,L=ll,S=ss,chr.n=kp[,1])
  return(ciro)
}


ciro.dat.g1 <- function(lb,sb,lsb,kp){
  
  nchrll <- as.numeric(kp[,2])/10^6
  seg.num    <- length(nchrll)
  
  gap.angle.size <- 1
  seg.angle.from <- 276
  
  seg.full.l  <- sum(nchrll)
  cir.angle.r <- (360-seg.num*gap.angle.size)/seg.full.l
  
  out.s     <- c()
  l.old     <- 0
  gap.angle <- 0
  
  for (i in 1:seg.num){
    seg.n <- paste("Chr",i,sep="")
    len   <- cumsum(nchrll)[i]
    w1    <- cir.angle.r*l.old + gap.angle
    w2    <- cir.angle.r*len   + gap.angle
    out.s     <- rbind(out.s, c(seg.n, w1+seg.angle.from, w2+seg.angle.from, l.old, len, 0, nchrll[i]))
    gap.angle <- gap.angle + gap.angle.size
    l.old     <- len;
  }
  ll <- list()
  ll[[1]] <- lb; ll[[2]] <- sb; ll[[3]] <- lsb
  ciro <- list(chr.po=out.s,L=ll,chr.n=kp[,1])
  return(ciro)
}




draw.text.rt <- function(xc, yc, r, w, n, col="black", cex=1, side="out",LG){
  
  w     <- w%%360;
  the.o <- w;
  
  the.w <- 360-w;
  w     <- w/360*2*pi;
  x     <- xc+r*cos(w);
  y     <- yc-r*sin(w);
  
  
  num2  <- 32;
  
  if (side=="out"){
    if (the.w <= 90 ){
      the.pos <- 4;
    } else if (the.w > 90 & the.w <= 180) {
      the.w <- the.w ;
      the.pos <- 2;
    } else if (the.w > 180 & the.w <= 270){
      the.w <- the.w%%180;
      the.pos <- 2;
    } else if (the.w > 270 & the.w <= 360){
      the.w <- the.w +180;
      the.pos <- 4;
    }
    
    if (the.pos==2){
      x <- x+num2;
    }
    if (the.pos==4){
      x <- x-num2;
    }
  } 
  #if(LG==22){
   # text(x=x+8, y+3, adj=0, offset=1, labels=n, srt=the.w+270, 
   #      pos=the.pos, col=col, cex=cex);
  #}else if(LG==23){
  #  text(x=x+7, y+2, adj=0, offset=1, labels=n, srt=the.w+270, 
   #      pos=the.pos, col=col, cex=cex);
  #}else{
    text(x, y, adj=0, offset=1, labels=n, srt=the.w+270, 
         pos=the.pos, col=col, cex=cex);
  #}

  
}


do.scale.cir <-function (xc=xc, yc=yc, the.r=the.r, total.num=total.num, 
                         col="blue", lwd=0.001,
                         V1=V1, V2=V2, V3=V3, V4=V4){
  
  sum.po  <- as.numeric(total.num);         # total number
  scale.w <- 30#sum.po/150;                         # one scale size
  scale.l <- nchar(as.integer(scale.w))-1;       # length of the number
  scale.i <- as.integer(scale.w/(10^scale.l));   # the first digital
  scale.d <- scale.i * 10^scale.l;               # the first digital, then 0
  scale.m <- 1 * 10^scale.l;                     # the unit of the scale   
  #draw.arc.s(xc, yc, w1=V1, w2=V2, r=the.r, col=col, lwd=lwd);
  
  start.p  <- 0;
  w        <- scale.v(as.numeric(start.p), V1, V2, V3, V4);
  scale.n  <- as.integer(as.numeric(V4)/scale.d+1);
  
  for (j in 1:scale.n){
    po.s <- as.integer(start.p/scale.m);
    
    draw.line(xc,    yc,   w, the.r-10, the.r+3, col=col, lwd=lwd);
    #xx draw.line1(xc,    yc,   w, the.r-10, the.r+3, col=col, lwd=lwd);
    draw.text.rt.axes(xc, yc,  the.r+6, w, po.s*10, col=1, cex=0.5);
    
    start.p <- start.p + scale.d/2;
    w        <- scale.v(as.numeric(start.p), V1, V2, V3, V4);
    
    if (w <= V2){  
      draw.line(xc, yc, w, the.r-10, the.r+1, col=col, lwd=lwd);
    }
    
    start.p <- start.p + scale.d/2;
    w        <- scale.v(as.numeric(start.p), V1, V2, V3, V4);
  }
  
}

draw.line1 <- function (xc, yc, w, l1, l2, col=col, lwd=lwd, lend=1) {
  w  <- (w/360)*2*pi;
  x1 <- xc+l1*cos(w);
  y1 <- yc-l1*sin(w);
  x2 <- xc+l2*cos(w);
  y2 <- yc-l2*sin(w);
  #segments(x1, y1, x2, y2, col=col, lwd=lwd, lend=lend);
  return(c(x2,y2))
}


draw.text.rt.axes <- function(xc, yc, r, w, n, col="black", cex=1, side="out"){
  
  w     <- w%%360;
  the.o <- w;
  
  the.w <- 360-w;
  w     <- w/360*2*pi;
  x     <- xc+r*cos(w);
  y     <- yc-r*sin(w);
  
  
  num2  <- 23;
  
  if (side=="out"){
    if (the.w <= 90 ){
      the.pos <- 4;
    } else if (the.w > 90 & the.w <= 180) {
      the.w <- the.w +180;
      the.pos <- 2;
    } else if (the.w > 180 & the.w <= 270){
      the.w <- the.w%%180;
      the.pos <- 2;
    } else if (the.w > 270 & the.w <= 360){
      the.w <- the.w +180+180;
      the.pos <- 4;
    }
    
    if (the.pos==2){
      x <- x+num2;
    }
    if (the.pos==4){
      x <- x-num2;
    }
  } 
  
  text(x, y, adj=0, offset=1, labels=n, srt=the.w, 
       pos=the.pos, col=col, cex=cex);
}



draw.arc.s <- function (xc, yc, r, w1, w2, col="lightblue", lwd=1, lend=1){
  
  ang.d <- abs(w1-w2);
  pix.n <- ang.d * 5;
  if (pix.n < 2){
    pix.n <- 2;
  }
  
  ang.seq <- rev(seq(w1,w2,length.out=pix.n));
  ang.seq <- ang.seq/360*2*pi;
  
  fan.i.x <- xc + cos(ang.seq) * r;
  fan.i.y <- yc - sin(ang.seq) * r;
  lines(fan.i.x, fan.i.y, col=col, lwd=lwd, type="l", lend=lend);
}



scale.v <- function(v, a, b, min.v, max.v) {
  v <- v-min.v; 
  v <- v/(max.v-min.v); 
  v <- v*(b-a);  
  v+a
}

draw.line <- function (xc, yc, w, l1, l2, col=col, lwd=lwd, lend=1) {
  w  <- (w/360)*2*pi;
  x1 <- xc+l1*cos(w);
  y1 <- yc-l1*sin(w);
  x2 <- xc+l2*cos(w);
  y2 <- yc-l2*sin(w);
  segments(x1, y1, x2, y2, col=col, lwd=lwd, lend=lend);
}

draw.arc.pg <- function (xc, yc, 
                         w1, w2, r1, r2, col="lightblue", border="lightblue", lwd=0.01
){
  
  ang.d <- abs(w1-w2);
  pix.n <- ang.d * 10;
  if (pix.n < 10){
    pix.n <- 10;
  }
  
  ang.seq <- rev(seq(w1,w2,length.out=pix.n));
  ang.seq <- ang.seq/360*2*pi;
  
  fan.i.x <- xc + cos(ang.seq) * r1;
  fan.i.y <- yc - sin(ang.seq) * r1;
  
  
  fan.o.x <- xc + cos(ang.seq) * r2;
  fan.o.y <- yc - sin(ang.seq) * r2;
  
  polygon(c(rev(fan.i.x), fan.o.x ), c(rev(fan.i.y), fan.o.y), 
          fillOddEven=F, border=border, col=col, lwd=lwd, lend=1)
  
}




draw.point <- function (xc, yc, w1, r1, col=col,pch=1,cex=1){
  
  ang.seq <-w1/360*2*pi
  fan.i.x <- xc + cos(ang.seq) * r1
  fan.i.y <- yc - sin(ang.seq) * r1
  
  points(fan.i.x, fan.i.y, col=col, pch=pch,cex=cex)
  
}



draw.arc.st <- function (xc, yc, r, w1, w2, col="lightblue", lwd=1, lend=1){
  
  ang.d <- abs(w1-w2);
  pix.n <- ang.d * 5;
  if (pix.n < 2){
    pix.n <- 2;
  }
  
  ang.seq <- rev(seq(w1,w2,length.out=pix.n));
  ang.seq <- ang.seq/360*2*pi;
  
  fan.i.x <- xc + cos(ang.seq) * r;
  fan.i.y <- yc - sin(ang.seq) * r;
  lines(fan.i.x, fan.i.y, col=col, lwd=lwd, lty=1);
}







draw.line3 <- function (xc, yc, w1, w2, r1, r2, col=col, lwd=lwd,lty=1){
  theangle1 <- w1;
  theangle2 <- w2;
  l1        <- r1;
  l2        <- r2;
  
  theangle1 <- (theangle1/360)*2*pi;
  x1        <- xc+l1*cos(theangle1);
  y1        <- yc-l1*sin(theangle1);
  
  theangle2 <- (theangle2/360)*2*pi;
  x2        <- xc+l2*cos(theangle2);
  y2        <- yc-l2*sin(theangle2);
  
  segments(x1, y1, x2, y2, col=col, lwd=lwd,lty=lty,lend="butt");
}






draw.line4 <- function (xc, yc, w1, w2, r1, r2, col=col, lwd=lwd,lty=1){
  theangle1 <- w1;
  theangle2 <- w2;
  l1        <- r1;
  l2        <- r2;
  
  theangle1 <- (theangle1/360)*2*pi;
  x1        <- xc+l1*cos(theangle1);
  y1        <- yc-l1*sin(theangle1);
  
  theangle2 <- (theangle2/360)*2*pi;
  x2        <- xc+l2*cos(theangle2);
  y2        <- yc-l2*sin(theangle2);
  
  segments(x1, y1, x2, y2, col=col, lwd=lwd);
}


draw.link <- function(xc, yc, r, w1, w2, col=col, lwd=lwd) {
  # for translocation
  w3  <- (w1+w2)/2;
  w1  <- w1/360*2*pi;
  w2  <- w2/360*2*pi;
  w3  <- w3/360*2*pi;
  x0  <- xc+r*cos(w1);
  y0  <- yc-r*sin(w1);
  x1  <- xc+r*cos(w2);
  y1  <- yc-r*sin(w2);
  x <- c(x0,xc,xc,x1);
  y <- c(y0,yc,yc,y1);
  points(bezierCurve(x,y,60), type="l", col=col, lwd=lwd, lend="butt")
}

draw.link.pg <- function(xc, yc, r, w1.1, w1.2, w2.1, w2.2, col=col, lwd=lwd) {
  #####################################################
  w1 <- w1.1;
  w2 <- w2.2;
  w3  <- (w1+w2)/2;
  w1  <- w1/360*2*pi;
  w2  <- w2/360*2*pi;
  w3  <- w3/360*2*pi;
  x0  <- xc+r*cos(w1);
  y0  <- yc-r*sin(w1);
  x1  <- xc+r*cos(w2);
  y1  <- yc-r*sin(w2);
  x <- c(x0,xc,xc,x1);
  y <- c(y0,yc,yc,y1);
  bc1 <- bezierCurve(x,y,60);
  
  ang.d <- abs(w1.1-w1.2);
  pix.n <- ang.d * 10;
  if (pix.n < 10){
    pix.n <- 10;
  }
  
  ang.seq <- rev(seq(w1.1,w1.2,length.out=pix.n));
  ang.seq <- ang.seq/360*2*pi;
  
  fan.1.x <- xc + cos(ang.seq) * r;
  fan.1.y <- yc - sin(ang.seq) * r;
  
  ######################################################
  w1 <- w1.2;
  w2 <- w2.1;
  w3  <- (w1+w2)/2;
  w1  <- w1/360*2*pi;
  w2  <- w2/360*2*pi;
  w3  <- w3/360*2*pi;
  x0  <- xc+r*cos(w1);
  y0  <- yc-r*sin(w1);
  x1  <- xc+r*cos(w2);
  y1  <- yc-r*sin(w2);
  x <- c(x0,xc,xc,x1);
  y <- c(y0,yc,yc,y1);
  bc2 <- bezierCurve(x,y,60);
  
  ang.d <- abs(w2.1-w2.2);
  pix.n <- ang.d * 10;
  if (pix.n < 10){
    pix.n <- 10;
  }
  
  ang.seq <- rev(seq(w2.1,w2.2,length.out=pix.n));
  ang.seq <- ang.seq/360*2*pi;
  
  fan.2.x <- xc + cos(ang.seq) * r;
  fan.2.y <- yc - sin(ang.seq) * r;
  
  polygon(c(bc1$x, fan.2.x, rev(bc2$x), rev(fan.1.x)), 
          c(bc1$y, fan.2.y, rev(bc2$y), rev(fan.1.y)), 
          fillOddEven=T, border=col, col=col); 
}



bezierCurve <-function(x, y, n=10)  {	
  outx <- NULL	
  outy <- NULL 	
  i <- 1	
  for (t in seq(0, 1, length.out=n))		{		
    b <- bez(x, y, t)		
    outx[i] <- b$x		
    outy[i] <- b$y 		
    i <- i+1		
  } 	
  return (list(x=outx, y=outy))	
}


bez <-function(x, y, t)	{	
  outx <- 0	
  outy <- 0	
  n <- length(x)-1	
  for (i in 0:n)		{		
    outx <- outx + choose(n, i)*((1-t)^(n-i))*t^i*x[i+1]		
    outy <- outy + choose(n, i)*((1-t)^(n-i))*t^i*y[i+1]		
  } 	
  return (list(x=outx, y=outy))	
}



