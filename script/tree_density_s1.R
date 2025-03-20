library(ggtree)
library(ggplot2)
library(TreeAndLeaf)
library(RedeR)
library(igraph)
library(dendextend)

iss <- list()

iss[[1]] <- c(113,52,13,47,61,11,115,14,70,9,126,79,56,62,71,124,41,32,77,122,37,34,30,112,26,89,46,117,1,
              68,2,88,22,7,99,15,8,59,5,109,125,76,17,48,12,90,96,80,94)
iss[[2]] <- c(19,39,110,6,84,42,25,23,31,29,55,35,72,95,69,36,20,24,51,83,57,53,16,78,91,3,33,45,40,28,18,4,
              21,43,65,85,64)
iss[[3]] <- c(50,105,67,101,100,73,102,119,87,104,60,118,107,121,111)
iss[[4]] <- c(66,108,97,114,63,92,75,86,74,58,123,93) 
iss[[5]] <- c(103,49,82,54,106,120,116,81,98,44,10,38,27)

group <- rep(0,126)
group[iss[[1]]] <- "tRNAA_aAlb";group[iss[[2]]] <- "tRNAB_aAlb";
group[iss[[3]]] <- "tRNAC_aAlb";group[iss[[4]]] <- "tRNAD_aAlb";
group[iss[[5]]] <- "tRNAE_aAlb";






ndatl1 <-data.frame(ID=sname[which(ress2$Class==st[1])],Group=group,CP=round(log10(ress2[which(ress2$Class==st[1]),3]),2))


group.df <- data.frame(label=sname[which(ress2$Class==st[1])], Group=group, Size=round(log10(ress2[which(ress2$Class==st[1]),3]),2))
tree_s11 <- dplyr::full_join(tree_s1, group.df, by='label')

#--- Generate a ggtree with 'daylight' layout
pal <- c("#FF7F50","#BCEE68","#00BFFF","#FF6EB4","#836FFF")
ggt <- ggtree(tree_s11, layout = 'daylight', branch.length='none')

#--- Plot the ggtree
ggt + geom_tippoint(aes(color=Group,size=Size)) +  
  scale_color_manual(values=pal) + scale_y_reverse()

tall1 <- treeAndLeaf(ggt)

#--- Map attributes to the tree-and-leaf
#Note: 'refcol = 1' indicates that 'dat' col 1 will be used as mapping IDs
tall1 <- att.mapv(g = tall1, dat = group.df, refcol = 1)

#--- Set graph attributes using the 'att.setv' wrapper function
tall1 <- att.setv(g = tall1, from = "Group", to = "nodeColor",
                  cols = pal)
tall1 <- att.setv(g = tall1, from = "Size", to = "nodeSize",
                  xlim = c(10, 50, 10),nquant =3)

#--- Set graph attributes using 'att.addv' and 'att.adde' functions
tall1 <- att.addv(tall1, "nodeFontSize", value = 1)
tall1 <- att.addv(tall1, "nodeLineWidth", value = 0)
tall1 <- att.addv(tall1, "nodeColor", value = "black", index=!V(tall1)$isLeaf)
tall1 <- att.adde(tall1, "edgeWidth", value = 3)
tall1 <- att.adde(tall1, "edgeColor", value = "black")
#--- Call RedeR application
rdp <- RedPort()
calld(rdp)
resetd(rdp)
#--- Send the tree-and-leaf to the interactive R/Java interface
addGraph(obj = rdp, g = tall1, gzoom=50)

#--- Select inner nodes, preventing them from relaxing
selectNodes(rdp, V(tall1)$name[!V(tall1)$isLeaf], anchor=TRUE)

#--- Call 'relax' to fine-tune the leaf nodes
relax(rdp, p1=25, p2=100, p3=5, p5=1, p8=5, ps=TRUE)

#--- Add legends
addLegend.color(obj = rdp, tall1, title = "Family", 
                position = "topright",vertical=T,dxborder=60,dyborder=40,ftsize=13,bend=0.5,size=20)
addLegend.size(obj = rdp, tall1, title = "Copy number", 
               position = "topright", 
               vertical=T,dxborder=60,dyborder=200,ftsize=13,intersp=10)





s1s <- list()
s1s[[1]] <- ress2[which(ress2$Class==st[1]),][iss[[1]],]
s1s[[2]] <- ress2[which(ress2$Class==st[1]),][iss[[2]],]
s1s[[3]] <- ress2[which(ress2$Class==st[1]),][iss[[3]],]
s1s[[4]] <- ress2[which(ress2$Class==st[1]),][iss[[4]],]
s1s[[5]] <- ress2[which(ress2$Class==st[1]),][iss[[5]],]



K1 <- as.matrix(read.csv("./data/sesame1.divsum",sep="\t"))
mr1 <- read.csv("./data/mutation rates.csv",header=T)
r <- mean(mr1$Average.yearly.mutation.rate..m_yearly.[which(mr1$Taxonomic.group=="Mammal")])


Figure_age_pro_s1 <- function(K1,aout1,r,L1s){
  
  T1 <- as.numeric(K1[,5])/100/(2*r)/10^6
  T1[which(T1<0)] <- 0
  T1[which(is.na(T1))] <- 0
  
  ii1 <- match(L1s[[1]][,1],K1[,2]);ii2 <- match(L1s[[2]][,1],K1[,2]);ii3 <- match(L1s[[3]][,1],K1[,2]);
  ii4 <- match(L1s[[4]][,1],K1[,2]);ii5 <- match(L1s[[5]][,1],K1[,2])
  
  T11 <- T1[ii1];T22 <- T1[ii2];T33 <- T1[ii3];T44 <- T1[ii4];T55 <- T1[ii5];
  
  gss2 <- list()
  for(k1 in 1:length(L1s)){
    tml <- L1s[[k1]][,1]
    gss1 <- c()
    for(k2 in 1:length(tml)){
      si1 <- which(tml[k2]==aout1$tname)
      gss1 <- c(gss1,sum(aout1$qend[si1]-aout1$qstart[si1])/10^6)
    }
    gss2[[k1]] <- gss1
  }
  
  
  or1 <- order(T11);or2 <- order(T22);or3 <- order(T33);or4 <- order(T44);or5 <- order(T55);
  T111 <- T11[or1];T222 <- T22[or2];T333 <- T33[or3];T444 <- T44[or4];T555 <- T55[or5];
  gss11 <- gss2[[1]][or1];gss22 <- gss2[[2]][or2];gss33 <- gss2[[3]][or3];
  gss44 <- gss2[[4]][or4];gss55 <- gss2[[5]][or5];
  
  c1 <- cut(T111,breaks=seq(0,95,5));c2 <- cut(T222,breaks=seq(0,95,5));c3 <- cut(T333,breaks=seq(0,95,5))
  c4 <- cut(T444,breaks=seq(0,95,5));c5 <- cut(T555,breaks=seq(0,95,5));
  
  cl1 <- levels(c1);cl2 <- levels(c2);cl3 <- levels(c3);cl4 <- levels(c4);cl5 <- levels(c5);
  gs11h <- c();gs22h <- c();gs33h <- c();gs44h <- c();gs55h <- c();
  for(i in 1:length(cl1)){
    gs11h <- c(gs11h,sum(gss11[which(c1==cl1[i])]));gs22h <- c(gs22h,sum(gss22[which(c2==cl1[i])]))
    gs33h <- c(gs33h,sum(gss33[which(c3==cl1[i])]));gs44h <- c(gs44h,sum(gss44[which(c4==cl1[i])]))
    gs55h <- c(gs55h,sum(gss55[which(c5==cl1[i])]));
  }
  gsh <- list()
  gsh[[1]] <- gs11h;gsh[[2]] <- gs22h;gsh[[3]] <- gs33h;gsh[[4]] <- gs44h;gsh[[5]] <- gs55h;
  
  ln <- c("tRNAA_aAlb","tRNAB_aAlb","tRNAC_aAlb","tRNAD_aAlb","tRNAE_aAlb")
  
  
  
  pdf("Figure-age_tRNA.pdf",width=3,height=2.5)
  par(oma=c(0,0,0,0),fig=c(0,1,0,1))
  par(mar=c(2.7,3.7,0.5,0.5))
  pal <- c("#FF7F50","#BCEE68","#00BFFF","#FF6EB4","#836FFF","#00E5EE")
  for(i in 1:length(L1s)){
    x1 <- barplot(gsh[[i]]/sum(unlist(gsh))*100,las=2,col=pal[i],mgp=c(2,0.8,0),border=NA)
    box()
    if(i==1||i==2||i==5){
      text(17,max(gsh[[i]]/sum(unlist(gsh))*100)*0.92,ln[i],cex=1,col=pal[i])
    }else{
      text(8,max(gsh[[i]]/sum(unlist(gsh))*100)*0.85,ln[i],cex=1,col=pal[i])
    }
    axis(1,seq(0,22,length=6),seq(0,100,20),las=1,cex.axis=1,tck=-0.033,mgp=c(2.5,0.4,0),col="black")
    #axis(2,seq(0,9,3),seq(0,9,3),las=1,cex.axis=1,tck=-0.033,mgp=c(2.5,0.6,0),col="black")
    mtext("Age (Mya)",1,cex=1,line=1.3)
    mtext("Percent of tRNA (%)",2,cex=1,line=2.6)
  }
  
  dev.off()
}

Figure_age_pro_s1(K1=K1,aout1=aout2,r=r,L1s=s1s)



pdf("Figure_s1_st.pdf",width=7,height=3.5)
par(oma=c(0,0,0,0),fig=c(0,1/3,0,1))
par(mar=c(4.2,2.7,0.5,1))
sn <- c("tRNAA_aAlb","tRNAB_aAlb","tRNAC_aAlb","tRNAD_aAlb","tRNAE_aAlb")
pal <- c("#FF7F50","#BCEE68","#00BFFF","#FF6EB4","#836FFF","#00E5EE")
s1 <- ress2[which(ress2$Class==st[1]),]
plot(NA,NA,xlim=c(0.4,5.6),ylim=c(-1,85),xlab="",ylab="",xaxt="n", yaxt="n",xaxs="i", yaxs="i")
ssn <- list()
ssc <- list()
for(i in 1:length(iss)){
  sdo <- sum((unlist(strsplit(s1[iss[[i]],1],"-"))=="rnd"))
  rect(i-0.4,0,i+0.4,sdo,col=pal[i],density=30,angle=45)
  rect(i-0.4,sdo,i+0.4,length(iss[[i]]),col=pal[i],density=30,angle=-45)
  ssn[[i]] <- s1[iss[[i]],2]
  ssc[[i]] <- s1[iss[[i]],3]
  
}
rect(1.5,75,2.5,80,density=30,angle=45)
text(3.6,77.5,expression(italic("De novo")),cex=0.8)
rect(1.5,65,2.5,70,density=30,angle=-45)
text(3.8,67.5,"Homology",cex=0.8)
axis(1,seq(1,5,1),label=F,las=1,cex.axis=0.9,tck=-0.053,mgp=c(2.5,0.3,0),col="black")
text(x = 1:5, y = par("usr")[3] - 14.5,labels = sn,xpd = NA,srt = 45,cex = 0.8)
axis(2,seq(0,80,20),seq(0,80,20),las=1,cex.axis=0.9,tck=-0.053,mgp=c(2.5,0.5,0),col="black")
mtext("Number of sequence",2,cex=1,line=1.5)

par(oma=c(0,0,0,0),fig=c(1/3,2/3,0,1),new=TRUE)
par(mar=c(4.2,2.5,0.5,0.5))
a <- boxplot(ssn,col=pal,xaxt="n", yaxt="n",xaxs="i", yaxs="i",ylim=c(0,6200))
axis(1,seq(1,5,1),label=F,las=1,cex.axis=0.9,tck=-0.053,mgp=c(2.5,0.3,0),col="black")
text(x = 1:5, y = par("usr")[3] - 1054.5,labels = sn,xpd = NA,srt = 45,cex = 0.8)
axis(2,seq(0,6,2)*1000,seq(0,6,2),las=1,cex.axis=0.9,tck=-0.053,mgp=c(2.5,0.6,0),col="black")
mtext("Length of sequence (kb)",2,cex=1,line=1.5)

par(oma=c(0,0,0,0),fig=c(2/3,1,0,1),new=TRUE)
par(mar=c(4.2,2.5,0.5,0.5))
a <- boxplot(ssc,col=pal,xaxt="n", yaxt="n",xaxs="i", yaxs="i",ylim=c(0,56000))
axis(1,seq(1,5,1),label=F,las=1,cex.axis=0.9,tck=-0.053,mgp=c(2.5,0.3,0),col="black")
text(x = 1:5, y = par("usr")[3] - 9624.5,labels = sn,xpd = NA,srt = 45,cex = 0.8)
axis(2,seq(0,5,1)*10000,seq(0,5,1),las=1,cex.axis=0.9,tck=-0.053,mgp=c(2.5,0.6,0),col="black")
mtext("Copy number of sequence",2,cex=1,line=1.5)
mtext(expression(paste("×", 10^4,sep="")),3,cex=0.7,line=-1,adj=-0.18)


dev.off()


