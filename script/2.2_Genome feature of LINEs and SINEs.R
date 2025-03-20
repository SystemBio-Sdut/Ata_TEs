library(repeatR)
library(seqinr)

source("rep_sup.R")
#afa_chr_ngap1.fa.out annotates the result for the genome repeat sequence
ahap <- read_rm("data/afa_chr_ngap1.fa.out")

a <- read.table("data/family_info2.txt",header=T)

una <- read.fasta("data/unkown.fasta.lib")
kp <- read.table("data/genome.length")

unan <- c()
for(i in 1:length(una)){
  tmp1 <- attributes(una[[i]])$Annot
  tmp11 <- strsplit(tmp1,": ")[[1]]
  tmp12 <- strsplit(tmp11[[2]],"|",fixed=T)[[1]][1]
  unan <- c(unan,tmp12)
}

unan[which(unan=="reverse")] <- 0
unan[which(unan=="forward")] <- 0
unan[which(unan=="unclear")] <- 0
unan[which(unan=="reverse")] <- 0
unan1 <- unan[-which(unan=="0")]
unaa <- unique(unlist(strsplit(names(una),"|",fixed=T)))[-2]
unaa1 <- unaa[-which(unan=="0")]

for(i in 1:length(unaa1)){
  ii <- which(ahap$tname==unaa1[i])
  ahap$tclass[ii] <- unan1[i]
}

a[match(unaa1,a[,1]),4] <- unan1


at <- table(a$Class)

atn <- names(at)

lengs <- c()
for(i in 1:length(atn)){
  
  index <- which(ahap$tclass==atn[i])
  ahapi <- ahap[index,]
  tl <- sum(ahapi$qend-ahapi$qstart)/1000000
  lengs <- c(lengs,tl)
}

atn1 <- c("LINE","LINE/CR1","LINE/Dong-R4","LINE/I-Jockey","LINE/L1","LINE/L1-Tx1","LINE/L2","LINE/RTE-BovB","LINE/RTE-X")
atn2 <- c("SINE/5S","SINE/5S-Deu-L2","SINE/Alu","SINE/B2","SINE/B4","SINE/ID","SINE/MIR","SINE/tRNA",
          "SINE/tRNA-Deu","SINE/tRNA-RTE","SINE")


lengs1 <- c()
for(i in 1:length(atn1)){
  
  index <- which(ahap$tclass==atn1[i])
  ahapi <- ahap[index,]
  tl <- sum(ahapi$qend-ahapi$qstart)/1000000
  lengs1 <- c(lengs1,tl)
}


lengs2 <- c()
for(i in 1:length(atn2)){
  
  index <- which(ahap$tclass==atn2[i])
  ahapi <- ahap[index,]
  tl <- sum(ahapi$qend-ahapi$qstart)/1000000
  lengs2 <- c(lengs2,tl)
}


names(lengs2) <- atn2
names(lengs1) <- atn1

slengs1 <- sort(lengs1,decreasing = T)
slengs2 <- sort(lengs2,decreasing = T)

s1p <- slengs1/2650*100;s2p <- slengs2/2650*100


sln1 <- names(s1p)
sln2 <- names(s2p)

r1s <- c()
all1 <- c()
for( i in 1:length(sln1)){
  aii1 <- which(ahap$tclass==sln1[i])
  t1mp <- ahap$tname[aii1]
  tcp1 <- length(t1mp)
  tt1 <- unique(t1mp)
  tt1t <- unlist(strsplit(tt1,"-"))
  dovo1 <- length(which(tt1t=="rnd"))
  homo1 <- length(tt1)-dovo1
  ml1 <- mean(a[match(tt1,a$ID),2])
  all1 <- c(all1,aii1)
  r1s <- rbind(r1s,c(sln1[i],tcp1,length(tt1),dovo1,homo1,ml1))
}


r2s <- c()
all2 <- c()
for( i in 1:length(sln2)){
  aii2 <- which(ahap$tclass==sln2[i])
  t2mp <- ahap$tname[aii2]
  tcp2 <- length(t2mp)
  tt2 <- unique(t2mp)
  tt2t <- unlist(strsplit(tt2,"-"))
  dovo2 <- length(which(tt2t=="rnd"))
  homo2 <- length(tt2)-dovo2
  ml2 <- mean(a[match(tt2,a$ID),2])
  all2 <- c(all2,aii2)
  r2s <- rbind(r2s,c(sln2[i],tcp2,length(tt2),dovo2,homo2,ml2))
}


tab1 <- cbind(rbind(r1s,r2s),c(s1p,s2p),c(slengs1,slengs2))
#write.csv(tab1,file="sine_line_st.csv",row.names = F,quote = T)

all1su <- sort(unique(all1),decreasing = F)
all2su <- sort(unique(all2),decreasing = F)

aout1 <- ahap[all1su,]
aout2 <- ahap[all2su,]





lds <- rep_density(karyotype=kp,gff=aout1,bin=250000)
sds <- rep_density(karyotype=kp,gff=aout2,bin=250000)


lnl <- rep_nl(karyotype=kp,gff=aout1)
snl <- rep_nl(karyotype=kp,gff=aout2)

lrs <- rep_size(karyotype=kp,gff=aout1,bin=250000)
srs <- rep_size(karyotype=kp,gff=aout2,bin=250000)

ldist <- rep_dist(karyotype=kp,gff=aout1,bin=c(0,200,500,1000,2000,5000,10000),dist_tre=100,repn=4)
sdist <- rep_dist(karyotype=kp,gff=aout2,bin=c(0,200,500,1000,2000,5000,10000),dist_tre=100,repn=4)
lsdist <- rep_dist(karyotype=kp,gff=rbind(aout1,aout2),bin=c(0,200,500,1000,2000,5000,10000),dist_tre=100,repn=4)

lb <- c()
sb <- c()
lsb <- c()

for(i in 1:25){
  ld1 <- rep_density_single(ldist$nbgff[[i]],bin=500000)$rett
  sd1 <- rep_density_single(sdist$nbgff[[i]],bin=500000)$rett
  lsd1 <- rep_density_single(lsdist$nbgff[[i]],bin=500000)$rett
  lb <- rbind(lb,as.data.frame(cbind(rep(kp[i,1],dim(ld1)[1]),ld1)))
  sb <- rbind(sb,as.data.frame(cbind(rep(kp[i,1],dim(sd1)[1]),sd1)))
  lsb <- rbind(lsb,as.data.frame(cbind(rep(kp[i,1],dim(lsd1)[1]),lsd1)))
}

aa <- rep_density_single(ldist$nbgff[[i]],bin=500000)$rett
#del <- which(aa$rett[,3]==0)
#plot(aa$rett[,3])

plot(0,0,ylim=c(0,600),xlim=c(0,2e+08))

rett1 <- aa$rett
for(i in 1:dim(rett1)[1]){
  rect(rett1[i,1],0,rett1[i,2],rett1[i,3],border=NA,col="red")
}


########提取BED########
options(scipen=999)
lpos <- aout1[,c(5,6,7)]
write.table(lpos,file="lpos.bed",quote = F,row.names = F,col.names = F,sep="\t")
spos <- aout2[,c(5,6,7)]
write.table(spos,file="spos.bed",quote = F,row.names = F,col.names = F,sep="\t")
#seqkit subseq --bed lpos.bed -u 0 -d 0 -o LINE.fa afa_chr_ngap1.fa
#seqkit subseq --bed spos.bed -u 0 -d 0 -o SINE.fa afa_chr_ngap1.fa


library(ChIPseeker)

gff <- "./data/Ata.gff"

txdb <- makeTxDbFromGFF(file = gff,format="gff3")

lanno <- annotatePeak("./data/lpos.bed", tssRegion=c(-3000, 3000), TxDb=txdb,flankDistance = 5000,
                      assignGenomicAnnotation = TRUE,addFlankGeneInfo = TRUE)




sanno <- annotatePeak("./data/spos.bed", tssRegion=c(-3000, 3000), TxDb=txdb,flankDistance = 5000,
                      assignGenomicAnnotation = TRUE,addFlankGeneInfo = TRUE)
sanno


source("rep_plot1.R")
lanno1 <- as.data.frame(lanno@anno@elementMetadata)
lanno2 <- as.data.frame(lanno@detailGenomicAnnotation)

del <- which(lanno1$annotation=="Distal Intergenic")

lann1 <- lanno1[-del,]
lann2 <- lanno2[-del,]
aout11 <- aout1[-del,]

lst1 <- genic_rep(x1=aout11,x2=lann2)
lst1g <- genic_rep_number(x1=aout11,x2=lann2,x3=lann1)



sanno1 <- as.data.frame(sanno@anno@elementMetadata)
sanno2 <- as.data.frame(sanno@detailGenomicAnnotation)
del1 <- which(sanno1$annotation=="Distal Intergenic")

sann1 <- sanno1[-del1,]
sann2 <- sanno2[-del1,]
aout22 <- aout2[-del1,]

sst1 <- genic_rep(x1=aout22,x2=sann2)
sst1g <- genic_rep_number(x1=aout22,x2=sann2,x3=sann1)


lgt <- unlist(lapply(strsplit(lann1[,1],split=" "),function(x){unlist(x)[1]}))

sgt <- unlist(lapply(strsplit(sann1[,1],split=" "),function(x){unlist(x)[1]}))

allg <- sort(unique(c(lann1$geneId,sann1$geneId)))

allg1 <- sort(unique(c(lann1$geneId)))

allg2 <- sort(unique(c(sann1$geneId)))

com1 <- intersect(allg1,allg2)
lg_sp <- allg1[-match(com1,allg1)]
sg_sp <- allg2[-match(com1,allg2)]

library(clusterProfiler)
library(GO.db)

ann <- read.csv("./data/merge.anno",sep="\t",header=F)
an_go <- ann[,c(1,3)]
an_go1 <- an_go[!is.na(an_go[,2]),]
goterms <- as.matrix(Term(GOTERM))
goterms1 <- cbind(rownames(goterms),goterms)
rownames(goterms1) <- NULL

gogene <- tq_go_gene(allgog=an_go1)

#go
gene_gol <- paste(lg_sp,"1",sep=".")
go_enl <- enricher(gene_gol,TERM2GENE=gogene,TERM2NAME=goterms1,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
go_xxl <- go_enl@result

lrp <- list()
for(i in 1:30){
  
  ts1 <- sample(1:32231,length(lg_sp),replace = F)
  tsr1 <- enricher(an_go[ts1,1],TERM2GENE=gogene,TERM2NAME=goterms1,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
  lrp[[i]] <- tsr1@result$pvalue
}


gene_gos <- paste(sg_sp,"1",sep=".")
go_ens <- enricher(gene_gos,TERM2GENE=gogene,TERM2NAME=goterms1,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
go_xxs <- go_ens@result

srp <- list()
for(i in 1:30){
  
  ts1 <- sample(1:32231,length(sg_sp),replace = F)
  tsr1 <- enricher(an_go[ts1,1],TERM2GENE=gogene,TERM2NAME=goterms1,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
  srp[[i]] <- tsr1@result$pvalue
}

library("GOxploreR")

gi1 <- which(go_xxl$p.adjust<0.01)
gi2 <- which(go_xxs$p.adjust<0.05)

gol <- go_xxl[gi1,]
gos <- go_xxs[gi2,]

gol1  <- cbind(Term=Ontology(gol[,1]),gol)
gos1  <- cbind(Term=Ontology(gos[,1]),gos)

gol11 <- gol1[-which(is.na(gol1[,1])),]
gos11 <- gos1

gol12 <- cbind(Level=rep(NA,dim(gol11)[1]),gol11)
gos12 <- cbind(Level=rep(NA,dim(gos11)[1]),gos11)

gol12$Level[which(gol12$Term=="BP")] <- GOTermBPOnLevel(gol12$ID[which(gol12$Term=="BP")])[,2]
gol12$Level[which(gol12$Term=="MF")] <- GOTermMFOnLevel(gol12$ID[which(gol12$Term=="MF")])[,2]
gol12$Level[which(gol12$Term=="CC")] <- GOTermCCOnLevel(gol12$ID[which(gol12$Term=="CC")])[,2]


gos12$Level[which(gos12$Term=="BP")] <- GOTermBPOnLevel(gos12$ID[which(gos12$Term=="BP")])[,2]
gos12$Level[which(gos12$Term=="MF")] <- GOTermMFOnLevel(gos12$ID[which(gos12$Term=="MF")])[,2]




tn1 <- c("LINE/L1","LINE/L2","LINE/I-Jockey","LINE/RTE-BovB","LINE")
tn2 <- c("SINE/tRNA","SINE/B2","SINE/ID","SINE/MIR","SINE")
gt <- c("Promoter","5'","Exon","Intron","3'","Downstream")
ng1 <- length(lg_sp)
ng2 <- length(sg_sp)

tmpgc1 <- c()
for(i in 1:ng1){
  
  gene <- lg_sp[i]
  ii1 <- which(gene==lann1$geneId)
  if(is.na(ii1[1])){
    tmp11 <- rep(0,30)
  }else{
    tmp11 <- c()
    aouts <- aout11[ii1,]
    for(k in 1:length(gt)){
      ii11 <- which(lgt[ii1]==gt[k])
      if(is.na(ii11[1])){
        tmp111 <- rep(0,5)
      }else{
        tmp111 <- c()
        aoutss <- aouts[ii11,]
        for(kk in 1:length(tn1)){
          ii111 <- which(aoutss$tclass==tn1[kk])
          if(is.na(ii111[1])){
            tmp111 <- c(tmp111,0)
          }else{
            tmp111 <- c(tmp111,sum(aoutss$qend[ii111]-aoutss$qstart[ii111]))
          }
        }
      }
      tmp11 <- c(tmp11,tmp111)
    }
  }
  
  tmpgc1 <- rbind(tmpgc1,tmp11)
}


tmpgc2 <- c()
for(i in 1:ng2){
  gene <- sg_sp[i]
  ii2 <- which(gene==sann1$geneId)
  if(is.na(ii2[1])){
    tmp22 <- rep(0,30)
  }else{
    tmp22 <- c()
    aoutl <- aout22[ii2,]
    for(k in 1:length(gt)){
      ii22 <- which(sgt[ii2]==gt[k])
      if(is.na(ii22[1])){
        tmp222 <- rep(0,5)
      }else{
        tmp222 <- c()
        aoutll <- aoutl[ii22,]
        for(kk in 1:length(tn2)){
          ii222 <- which(aoutll$tclass==tn2[kk])
          if(is.na(ii222[1])){
            tmp222 <- c(tmp222,0)
          }else{
            tmp222 <- c(tmp222,sum(aoutll$qend[ii222]-aoutll$qstart[ii222]))
          }
        }
      }
      tmp22 <- c(tmp22,tmp222)
    }
  }
  tmpgc2 <- rbind(tmpgc2,tmp22)
}


library("preprocessCore")
datl <- log10(tmpgc1+1)
dats <- log10(tmpgc2+1)
rownames(datl) <- lg_sp
rownames(dats) <- sg_sp
