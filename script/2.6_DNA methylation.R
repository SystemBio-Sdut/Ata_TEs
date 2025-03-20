#load("hfi_meth.RData")
source("meth_sup.R")
source("meth_plot.R")
library(dplyr)



chrn <- paste("chr",c(1:23,"X","Y"),sep="")
dd1 <- list()
for(i in 1:length(chrn)){
 dd1[[i]] <-  density(a2[which(a2[,1]==chrn[i]),4])
}

density_meth_curve(chrn=chrn,dd1=dd1)

density_meth_curve_chr(chrn=chrn,dd1=dd1)

karyotype <- read.table("../../Old/Ahed/genome.length")

meth_bn <- meth_size(karyotype=karyotype,gff=a2,bin=250000)

meth_bin_chr(karyotype=karyotype,meth_bn=meth_bn)


load("../data/LINE.RData")
load("../data/SINE.RData")

lbn <- rep_size(karyotype=karyotype,gff=aout1,bin=250000)
sbn <- rep_size(karyotype=karyotype,gff=aout2,bin=250000)

rep_meth_cor_plot(meth_bn=meth_bn,lbn=lbn,sbn=sbn)

rep_meth_cor_plot_chr_line(karyotype=karyotype,meth_bn=meth_bn,lbn=lbn)
rep_meth_cor_plot_chr_sine(karyotype=karyotype,meth_bn=meth_bn,lbn=sbn)

#a2i <- list()
#for(chrn in 1:dim(karyotype)[1]){
#  
#  a2i[[chrn]] <- which(a2[,1]==karyotype[chrn,1])
#}


#aout1_m <- c()
#for(i in 1:dim(aout1)[1]){
#  li <- a2i[[match(aout1$qname[i],karyotype[,1])]]
#  aout1_m <- c(aout1_m,sum(between(a2[li,3],aout1$qstart[i],aout1$qend[i])))
#}

#aout1_mp <- rep_site_cal(gff1=aout1,gff2=a2,karyotype=karyotype,t=100)
#save(aout1_mp,file="line_meth.RData")
#aout2_m <- c()
#for(i in 1:dim(aout2)[1]){
#   si <- which(a2[,1]==aout2$qname[i])
#  aout2_m <- c(aout2_m,sum(between(a2[si,3],aout2$qstart[i],aout2$qend[i])))
#}
#aout2_mp <- rep_site_cal(gff1=aout2,gff2=a2,karyotype=karyotype,t=100)
#save(aout2_mp,file="sine_meth.RData")

load("line_meth.RData")
load("sine_meth.RData")

Figure_rep_meth_barplot(a2=a2,aout1_mp=aout1_mp,aout2_mp=aout2_mp)

rep_meth_cor_length(a2=a2,aout1=aout1,aout2=aout2,aout1_mp=aout1_mp,aout2_mp=aout2_mp)

load("anno.RData")

Figure_rep_gene_meth_barplot(anno=anno,aout1_mp=aout1_mp,aout2_mp=aout2_mp)
  
  

mr1 <- read.csv("../data/mutation rates.csv",header=T)
r <- mean(mr1$Average.yearly.mutation.rate..m_yearly.[which(mr1$Taxonomic.group=="Mammal")])


itt <- as.matrix(read.csv("../data/sesame1.divsum",sep="\t"))
ilu <- unique(match(aout1$tname,itt[,2]))
isu <- unique(match(aout2$tname,itt[,2]))
agel <- as.numeric(itt[ilu,5])/100/(2*r)/10^6
ages <- as.numeric(itt[isu,5])/100/(2*r)/10^6

ilun <- itt[ilu,2]
isun <- itt[isu,2]

tl1 <- c()
for(i in 1:length(ilun)){
  tmpl1 <- aout1_mp[ilun[i]==aout1$tname]
  tl1 <- rbind(tl1,c(sum(tmpl1==0),sum(tmpl1!=0)))
}

stl1 <- tl1[which(rowSums(tl1)>100),]
agel1 <- agel[which(rowSums(tl1)>100)]

sl1 <- c()
for(i in 1:length(isun)){
  tmps1 <- aout2_mp[isun[i]==aout2$tname]
  sl1 <- rbind(sl1,c(sum(tmps1==0),sum(tmps1!=0)))
}


ssl1 <- sl1[which(rowSums(sl1)>100),]
ages1 <- ages[which(rowSums(sl1)>100)]


rep_meth_age_cor(stl1=stl1,agel1=agel1,ssl1=ssl1,ages1=ages1)




  
