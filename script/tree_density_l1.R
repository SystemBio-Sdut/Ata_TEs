
library(ggtree)
library(ggplot2)
library(TreeAndLeaf)
library(RedeR)
library(igraph)
library(dendextend)

ill <- list()

ill[[1]] <- c(31,40,18,56,19,5,29,2,28,1,96,79,113,94,85,73,65,50,124,63,91,190,89,175,161,117,168,172,169,
              68,81,51,42,45,33,41,35,7,59,123,9,147,22,90,20,62,34,131,163,121,159,151,120,158,76,155,
              135,103,101,162, 139,160,148,184,112,165,132,110,138,134,157,100,146,189,171,191,187,154,178,185,
              181,80,87,36)
ill[[2]] <- c(108,115,118,130,77,188,144,39,102,152,106,177,129,136,54,104,128,179,48,98)
ill[[3]] <- c(14,64,58,27,10,8,44,43,25,21,6,141,16,69,12,170,4,11,24,3,32,97,30,37,26,38,23,13,15,17,111)
ill[[4]] <- c(84,86,47,109,61,119,107,93,105,127,149,166,53,137,145,46,99,143,95,140,174,114,122,75,55,92,125,52,
              72,67,164,180,176,156,173,183,116,186,88,57,167,49,70,82,83,60,153,182,66)
ill[[5]] <- c(142,78,133,74) 
ill[[6]] <- c(150,126,71)

group <- rep(0,191)
group[ill[[1]]] <- "LINE1A_aAlb";group[ill[[2]]] <- "LINE1B_aAlb";
group[ill[[3]]] <- "LINE1C_aAlb";group[ill[[4]]] <- "LINE1D_aAlb";
group[ill[[5]]] <- "LINE1E_aAlb";group[ill[[6]]] <- "LINE1F_aAlb";






ndatl1 <-data.frame(ID=fname[which(resl2$Class==lt[1])],Group=group,CP=round(log10(resl2[which(resl2$Class==lt[1]),3]),2))


group.df <- data.frame(label=fname[which(resl2$Class==lt[1])], Group=group, Size=round(log10(resl2[which(resl2$Class==lt[1]),3]),2))
tree_l11 <- dplyr::full_join(tree_l1, group.df, by='label')

#--- Generate a ggtree with 'daylight' layout
pal <- c("#FF7F50","#BCEE68","#00BFFF","#FF6EB4","#836FFF","#00E5EE")
ggt <- ggtree(tree_l11, layout = 'daylight', branch.length='none')

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





L1s <- list()
L1s[[1]] <- resl2[which(resl2$Class==lt[1]),][ill[[1]],]
L1s[[2]] <- resl2[which(resl2$Class==lt[1]),][ill[[2]],]
L1s[[3]] <- resl2[which(resl2$Class==lt[1]),][ill[[3]],]
L1s[[4]] <- resl2[which(resl2$Class==lt[1]),][ill[[4]],]
L1s[[5]] <- resl2[which(resl2$Class==lt[1]),][ill[[5]],]
L1s[[6]] <- resl2[which(resl2$Class==lt[1]),][ill[[6]],]




K1 <- as.matrix(read.csv("./data/sesame1.divsum",sep="\t"))
mr1 <- read.csv("./data/mutation rates.csv",header=T)
r <- mean(mr1$Average.yearly.mutation.rate..m_yearly.[which(mr1$Taxonomic.group=="Mammal")])


Figure_age_pro_l1 <- function(K1,aout1,r,L1s){
  
  T1 <- as.numeric(K1[,5])/100/(2*r)/10^6
  T1[which(T1<0)] <- 0
  T1[which(is.na(T1))] <- 0
  
  ii1 <- match(L1s[[1]][,1],K1[,2]);ii2 <- match(L1s[[2]][,1],K1[,2]);ii3 <- match(L1s[[3]][,1],K1[,2]);
  ii4 <- match(L1s[[4]][,1],K1[,2]);ii5 <- match(L1s[[5]][,1],K1[,2]);ii6 <- match(L1s[[6]][,1],K1[,2]);
  
  T11 <- T1[ii1];T22 <- T1[ii2];T33 <- T1[ii3];T44 <- T1[ii4];T55 <- T1[ii5];T66 <- T1[ii6];
  
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
  
  
  or1 <- order(T11);or2 <- order(T22);or3 <- order(T33);or4 <- order(T44);or5 <- order(T55);or6 <- order(T66)
  T111 <- T11[or1];T222 <- T22[or2];T333 <- T33[or3];T444 <- T44[or4];T555 <- T55[or5];T666 <- T66[or6]
  gss11 <- gss2[[1]][or1];gss22 <- gss2[[2]][or2];gss33 <- gss2[[3]][or3];
  gss44 <- gss2[[4]][or4];gss55 <- gss2[[5]][or5];gss66 <- gss2[[6]][or6];
  
  c1 <- cut(T111,breaks=seq(0,95,5));c2 <- cut(T222,breaks=seq(0,95,5));c3 <- cut(T333,breaks=seq(0,95,5))
  c4 <- cut(T444,breaks=seq(0,95,5));c5 <- cut(T555,breaks=seq(0,95,5));c6 <- cut(T666,breaks=seq(0,95,5))
  
  cl1 <- levels(c1);cl2 <- levels(c2);cl3 <- levels(c3);cl4 <- levels(c4);cl5 <- levels(c5);cl6 <- levels(c6);
  gs11h <- c();gs22h <- c();gs33h <- c();gs44h <- c();gs55h <- c();gs66h <- c()
  for(i in 1:length(cl1)){
    gs11h <- c(gs11h,sum(gss11[which(c1==cl1[i])]));gs22h <- c(gs22h,sum(gss22[which(c2==cl1[i])]))
    gs33h <- c(gs33h,sum(gss33[which(c3==cl1[i])]));gs44h <- c(gs44h,sum(gss44[which(c4==cl1[i])]))
    gs55h <- c(gs55h,sum(gss55[which(c5==cl1[i])]));gs66h <- c(gs66h,sum(gss66[which(c6==cl1[i])]))
  }
  gsh <- list()
  gsh[[1]] <- gs11h;gsh[[2]] <- gs22h;gsh[[3]] <- gs33h;gsh[[4]] <- gs44h;gsh[[5]] <- gs55h;gsh[[6]] <- gs66h;
  
  ln <- c("LINE1A_aAlb","LINE1B_aAlb","LINE1C_aAlb","LINE1D_aAlb","LINE1E_aAlb","LINE1F_aAlb")
  
  
  
  pdf("Figure-age_l1.pdf",width=3,height=2.5)
  par(oma=c(0,0,0,0),fig=c(0,1,0,1))
  par(mar=c(2.7,3.7,0.5,0.5))
  pal <- c("#FF7F50","#BCEE68","#00BFFF","#FF6EB4","#836FFF","#00E5EE")
  for(i in 1:length(L1s)){
    x1 <- barplot(gsh[[i]]/sum(unlist(gsh))*100,las=2,col=pal[i],mgp=c(2,0.8,0),border=NA)
    box()
    if(i==1||i==3){
      text(16,max(gsh[[i]]/sum(unlist(gsh))*100)*0.85,ln[i],cex=1,col=pal[i])
    }else{
      text(8,max(gsh[[i]]/sum(unlist(gsh))*100)*0.85,ln[i],cex=1,col=pal[i])
    }
    axis(1,seq(0,22,length=6),seq(0,100,20),las=1,cex.axis=1,tck=-0.033,mgp=c(2.5,0.4,0),col="black")
    #axis(2,seq(0,9,3),seq(0,9,3),las=1,cex.axis=1,tck=-0.033,mgp=c(2.5,0.6,0),col="black")
    mtext("Age (Mya)",1,cex=1,line=1.3)
    mtext("Percent of LINE1 (%)",2,cex=1,line=2.6)
  }
  
  dev.off()
}

Figure_age_pro_l1(K1=K1,aout1=aout1,r=r,L1s=L1s)





pdf("Figure_l1_st.pdf",width=7,height=3.5)
par(oma=c(0,0,0,0),fig=c(0,1/3,0,1))
par(mar=c(4.2,2.7,0.5,1))
ln <- c("LINE1A_aAlb","LINE1B_aAlb","LINE1C_aAlb","LINE1D_aAlb","LINE1E_aAlb","LINE1F_aAlb")
pal <- c("#FF7F50","#BCEE68","#00BFFF","#FF6EB4","#836FFF","#00E5EE")
l1 <- resl2[which(resl2$Class==lt[1]),]
plot(NA,NA,xlim=c(0.4,6.6),ylim=c(-1,85),xlab="",ylab="",xaxt="n", yaxt="n",xaxs="i", yaxs="i")
lln <- list()
llc <- list()
for(i in 1:length(ill)){
  sdo <- sum((unlist(strsplit(l1[ill[[i]],1],"-"))=="rnd"))
  rect(i-0.4,0,i+0.4,sdo,col=pal[i],density=30,angle=45)
  rect(i-0.4,sdo,i+0.4,length(ill[[i]]),col=pal[i],density=30,angle=-45)
  lln[[i]] <- l1[ill[[i]],2]
  llc[[i]] <- l1[ill[[i]],3]
  
}
rect(2.5,75,3.5,80,density=30,angle=45)
text(4.7,77.5,expression(italic("De novo")),cex=0.8)
rect(2.5,65,3.5,70,density=30,angle=-45)
text(4.9,67.5,"Homology",cex=0.8)
axis(1,seq(1,6,1),label=F,las=1,cex.axis=0.9,tck=-0.053,mgp=c(2.5,0.3,0),col="black")
text(x = 1:6, y = par("usr")[3] - 14.5,labels = ln,xpd = NA,srt = 45,cex = 0.8)
axis(2,seq(0,80,20),seq(0,80,20),las=1,cex.axis=0.9,tck=-0.053,mgp=c(2.5,0.5,0),col="black")
mtext("Number of sequence",2,cex=1,line=1.5)

par(oma=c(0,0,0,0),fig=c(1/3,2/3,0,1),new=TRUE)
par(mar=c(4.2,2.5,0.5,0.5))
a <- boxplot(lln,col=pal,xaxt="n", yaxt="n",xaxs="i", yaxs="i",ylim=c(0,13000))
axis(1,seq(1,6,1),label=F,las=1,cex.axis=0.9,tck=-0.053,mgp=c(2.5,0.3,0),col="black")
text(x = 1:6, y = par("usr")[3] - 2224.5,labels = ln,xpd = NA,srt = 45,cex = 0.8)
axis(2,seq(0,12,3)*1000,seq(0,12,3),las=1,cex.axis=0.9,tck=-0.053,mgp=c(2.5,0.6,0),col="black")
mtext("Length of sequence (kb)",2,cex=1,line=1.5)

par(oma=c(0,0,0,0),fig=c(2/3,1,0,1),new=TRUE)
par(mar=c(4.2,2.5,0.5,0.5))
a <- boxplot(llc,col=pal,xaxt="n", yaxt="n",xaxs="i", yaxs="i",ylim=c(0,56000))
axis(1,seq(1,6,1),label=F,las=1,cex.axis=0.9,tck=-0.053,mgp=c(2.5,0.3,0),col="black")
text(x = 1:6, y = par("usr")[3] - 9624.5,labels = ln,xpd = NA,srt = 45,cex = 0.8)
axis(2,seq(0,5,1)*10000,seq(0,5,1),las=1,cex.axis=0.9,tck=-0.053,mgp=c(2.5,0.6,0),col="black")
mtext("Copy number of sequence",2,cex=1,line=1.5)
mtext(expression(paste("×", 10^4,sep="")),3,cex=0.7,line=-1,adj=-0.18)
dev.off()






