
load("../data/LINE.RData")
load("../data/SINE.RData")

source("rep_3d.R")

a <- read.table("afa_100kb.ev_gc.txt")

kp <- read.table("genome.length")[,2]/10^6

lab <- AB_rep(gf1=a,gf2=aout1)
sab <- AB_rep(gf1=a,gf2=aout2)

Figure_AB_rep_pro(lab=lab,sab=sab)
Figure_rep_AB_pro(lab=lab,sab=sab)

lschr <- lsnorm(lab=lab,sab=sab)
Figure_rep_AB_chr(lschr=lschr)

labe <- AB_rep_each(gf1=a,gf2=aout1)
sabe <- AB_rep_each(gf1=a,gf2=aout2)


Figure_rep_ABeach_pro_line(labe=labe,sabe=sabe)



#hic_100k = fanc.load("afa_contact_map.hic@100000")
#ab = fanc.ABCompartmentMatrix.from_hic(hic_100k)
#ab_chr18 = ab.matrix(('chr18', 'chr18'))
#import pandas as pd
#import numpy as np
#dataframe = pd.DataFrame(ab_chr18.data)
#dataframe_chr18.to_csv("chr18_mat.csv")
load("../data/LINE.RData")
load("../data/SINE.RData")


source("rep_3d.R")

a <- read.table("./tmp_d/afa_100kb.ev_gc.txt")
kpf <- read.table("../../Old/Ahed/genome.length")
kp <- kpf[,2]/10^6
lab <- AB_rep(gf1=a,gf2=aout1)
sab <- AB_rep(gf1=a,gf2=aout2)

source("rep_3d_heatmap.R")

llb <- c(1:23,"X","Y")
for(i in 1:25){
  ab <- read.csv(paste("./mat/chr",llb[i],"_mat.csv",sep=""),header=F)[-1,-1]
  sq_heatmap(smat=ab,pc=a,lab=lab,sab=sab,ii=i,kp=kp,llb=llb,filepre="Figulbre_ls_cor_chr")
}

pall <- c()
for(i in 1:25){
  mat <- read.csv(paste("./mat/chr",llb[i],"_mat.csv",sep=""),header=F)[-1,-1]
  rtmp <- rep_3d_pro(lab=lab,sab=sab,mat=mat,ii=i)
  tmppro <- c(rtmp$lp_pro,rtmp$ln_pro,rtmp$sp_pro,rtmp$sn_pro,
              rtmp$lsp_pro,rtmp$lsn_pro)
  pall <- rbind(pall,tmppro)
}

rep_interaction_boxplot(pall=pall)

gene <- read.table("Ata.gff")
gene1 <- gene[which(gene$V3=="gene"),]
rm(gene)

tad <- read.table("hic_100k_corrected_boundaries.bed")

ltadst <- tad_rep(gf1=tad,gf2=aout1)
stadst <- tad_rep(gf1=tad,gf2=aout2)

loop <- read.table("afa_loops.bedgraph")

lloopst1x <- loop_rep(gf1=loop[,1:3],gf2=aout1)
lloopst1y <- loop_rep(gf1=loop[,4:6],gf2=aout1)
lloopst1d <- loop_rep(gf1=loop[,c(1,3,5)],gf2=aout1)

sloopst1x <- loop_rep(gf1=loop[,1:3],gf2=aout2)
sloopst1y <- loop_rep(gf1=loop[,4:6],gf2=aout2)
sloopst1d <- loop_rep(gf1=loop[,c(1,3,5)],gf2=aout2)

sq_commat1_1(lab=lab,sab=sab,tad=tad,loop=loop,gene1=gene1,llb=llb,kp=kp,filepre="Figure_density_tad_loop1_1.pdf")
sq_commat1_2(lab=lab,sab=sab,tad=tad,loop=loop,gene1=gene1,llb=llb,kp=kp,filepre="Figure_density_tad_loop1_2.pdf")


lsum <- sum(aout1$qend-aout1$qstart)/10^6
ssum <- sum(aout2$qend-aout2$qstart)/10^6

Figure_tad_loop_pro(ltadst=ltadst,stadst=stadst,llopst1d=llopst1d,soopst1d=soopst1d,
                    lsum=lsum,ssum=ssum)



chr2 <- read.csv("./mat/chr2_mat.csv",header=F)[-1,-1]
sq_heatmap_interval_tran(smat=chr2,pc=a,lab=lab,sab=sab,ii=2,tad=tad,loop=loop,gene1=gene1,
                         interval=c(10,40),xti=c(10,20,30,40),llb=llb,
                         filepre="Figure_hetamap_tran")

chr3 <- read.csv("./mat/chr3_mat.csv",header=F)[-1,-1]
sq_heatmap_interval_tran(smat=chr1,pc=a,lab=lab,sab=sab,ii=3,tad=tad,loop=loop,gene1=gene1,
                         interval=c(120,150),xti=c(120,130,140,150),llb=llb,
                         filepre="Figure_hetamap_tran")

