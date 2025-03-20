

load("data/SINE.RData")
library(seqinr)

ress <- read.table("data/family_sine_info1.txt",header=T)

ress1 <- ress[which(as.numeric(ress$Length)>100),]
ress2 <- ress1[which(as.numeric(ress1$Number) > 10),]


sfa <- read.fasta("data/family_sine.fasta")

st <- unique(ress2$Class)

for(i in 1:length(st)){
  
  tmp_fasta <- sfa[match(ress2$ID[which(ress2$Class==st[i])],names(sfa))]
  write.fasta(sequences = tmp_fasta,names=ress2$ID[which(ress2$Class==st[i])],
              file.out = paste("./data/tmp/s_",i,".fasta",sep=""))
}


sname <- ress2$ID 
sname[which(ress2$Class=="SINE/tRNA")] <- paste("tRNA-",1:length(which(ress2$Class=="SINE/tRNA")),"_aAlb",sep="")
sname[which(ress2$Class=="SINE/B2")] <- paste("B2-",1:length(which(ress2$Class=="SINE/B2")),"_aAlb",sep="")
sname[which(ress2$Class=="SINE/MIR")] <- paste("MIR-",1:length(which(ress2$Class=="SINE/MIR")),"_aAlb",sep="")
sname[which(ress2$Class=="SINE/ID")] <- paste("ID-",1:length(which(ress2$Class=="SINE/ID")),"_aAlb",sep="")
sname[which(ress2$Class=="SINE/5S-Deu-L2")] <- paste("5S-Deu-L2-",1:length(which(ress2$Class=="SINE/5S-Deu-L2")),"_aAlb",sep="")
sname[which(ress2$Class=="SINE/tRNA-RTE")] <- paste("tRNA-RTE-",1:length(which(ress2$Class=="SINE/tRNA-RTE")),"_aAlb",sep="")
sname[which(ress2$Class=="SINE/Alu")] <- paste("Alu-",1:length(which(ress2$Class=="SINE/Alu")),"_aAlb",sep="")
sname[which(ress2$Class=="SINE")] <- paste("UnS-",1:length(which(ress2$Class=="SINE")),"_aAlb",sep="")
sname[which(ress2$Class=="SINE/tRNA-Deu")] <- paste("SINE/tRNA-Deu-",1:length(which(ress2$Class=="SINE/tRNA-Deu")),"_aAlb",sep="")
#FISH
tq <- c(2,8,56,157)
write.fasta(sequences = sfa[match(ress2$ID[tq],names(sfa))],names=sname[tq],
            file.out = "FISH_sine.fasta")



library("msa")
library("ape")
library("ggtree")

seqs <- readAAStringSet("data/tmp/s_con.fasta")
names(seqs) <- c(st[1:7],paste("SINE/",sname[which(ress2$Class==st[8])],sep=""),st[9])

align_s <- msa(seqs,method=c("ClustalW"),verbose=T)
align_st <- msaConvert(align_s, type="seqinr::alignment")
d_s <- dist.alignment(align_st, "identity") 
tree_s <- njs(d_s)

css <- list(Known=names(seqs)[c(1:7,12)],Unknown=names(seqs)[-c(1:7,12)])
tree_snew <- groupOTU(tree_s,.node=css)

p2 <- ggtree(tree_snew, layout='circular', ladderize=FALSE, size=0.8, branch.length="none",col=c("#171717"))+
  geom_tiplab2(size=4.5,aes(color=group), show.legend=FALSE,hjust=-0.1)+
  scale_color_manual(values=c("#00B2EE","#EE6363"))+
  geom_tippoint(size=3,aes(color=group))+geom_nodepoint(color="#1E90FF", alpha=1/4, size=2)+
  geom_treescale(linesize=1)+theme(legend.title=element_text(face="bold"))+xlim(0, 10)+
  guides(color = guide_legend(title = "Type")) 
p2
ggsave(p2, file="Figure_sine-tree.pdf", height=7, width=8)

#s1
seqs1 <- readAAStringSet("data/tmp/s_1.fasta")
names(seqs1) <- sname[which(ress2$Class==st[1])]
align_s1 <- msa(seqs1,method=c("ClustalW"),verbose=T)
align_st1 <- msaConvert(align_s1, type="seqinr::alignment")
d_s1 <- dist.alignment(align_st1, "identity") 
tree_s1 <- njs(d_s1)

ps1 <- ggtree(tree_s1, layout='daylight', ladderize=T, size=0.8, branch.length="none",col=c("#171717"))+
  geom_tiplab(size=1,hjust=-0.15)+geom_tippoint(size=1,col="#EE1289")+geom_nodepoint(color="#1E90FF", alpha=1/4, size=1)#+
geom_treescale(linesize=1)+theme(legend.title=element_text(face="bold"))
ps1
ggsave(ps1, file="Figure_sine_s1-tree.pdf", height=7, width=7)


#s2
seqs2 <- readAAStringSet("data/tmp/s_2.fasta")
names(seqs2) <- sname[which(ress2$Class==st[2])]
align_s2 <- msa(seqs2,method=c("ClustalW"),verbose=T)
align_st2 <- msaConvert(align_s2, type="seqinr::alignment")
d_s2 <- dist.alignment(align_st2, "identity") 
tree_s2 <- njs(d_s2)

ps2 <- ggtree(tree_s2, layout='daylight', ladderize=FALSE, size=0.8, branch.length="none",col=c("#171717"))+
  geom_tiplab(size=2,hjust=-0.15)+geom_tippoint(size=1,col="#EE1289")+
  geom_nodepoint(color="#1E90FF", alpha=1/4, size=1)+
  xlim(-1.2,1)+ylim(-1.2,0.9)

ps2
ggsave(ps2, file="Figure_sine_s2-tree.pdf", height=6, width=6)

#s3
seqs3 <- readAAStringSet("data/tmp/s_3.fasta")
names(seqs3) <- sname[which(ress2$Class==st[3])]
align_s3 <- msa(seqs3,method=c("ClustalW"),verbose=T)
align_st3 <- msaConvert(align_s3, type="seqinr::alignment")
d_s3 <- dist.alignment(align_st3, "identity") 
tree_s3 <- njs(d_s3)

ps3 <- ggtree(tree_s3, layout='daylight', ladderize=FALSE, size=0.8, branch.length="none",col=c("#171717"))+
  geom_tiplab(size=2,hjust=-0.15)+geom_tippoint(size=1,col="#EE1289")+
  geom_nodepoint(color="#1E90FF", alpha=1/4, size=1)+
  xlim(-2,13)+ylim(-7.5,5)

ps3
ggsave(ps3, file="Figure_sine_s3-tree.pdf", height=6, width=6)

#s4
seqs4 <- readAAStringSet("data/tmp/s_4.fasta")
names(seqs4) <- sname[which(ress2$Class==st[4])]
align_s4 <- msa(seqs4,method=c("ClustalW"),verbose=T)
align_st4 <- msaConvert(align_s4, type="seqinr::alignment")
d_s4 <- dist.alignment(align_st4, "identity") 
tree_s4 <- njs(d_s4)

ps4 <- ggtree(tree_s4, layout='daylight', ladderize=FALSE, size=0.8, branch.length="none",col=c("#171717"))+
  geom_tiplab(size=2,hjust=-0.15)+geom_tippoint(size=1,col="#EE1289")+
  geom_nodepoint(color="#1E90FF", alpha=1/4, size=1)+
  xlim(-3.5,3.5)+ylim(-4,3.8)

ps4
ggsave(ps4, file="Figure_sine_s4-tree.pdf", height=6, width=6)


#s8
seqs8 <- readAAStringSet("data/tmp/s_8.fasta")
names(seqs8) <- sname[which(ress2$Class==st[8])]
align_s8 <- msa(seqs8,method=c("ClustalW"),verbose=T)
align_st8 <- msaConvert(align_s8, type="seqinr::alignment")
d_s8 <- dist.alignment(align_st8, "identity") 
tree_s8 <- njs(d_s8)

ps8 <- ggtree(tree_s8, layout='daylight', ladderize=FALSE, size=0.8, branch.length="none",col=c("#171717"))+
  geom_tiplab(size=2,hjust=-0.15)+geom_tippoint(size=1,col="#EE1289")+
  geom_nodepoint(color="#1E90FF", alpha=1/4, size=1)+
  xlim(-1.8,1.8)+ylim(-1.8,1.5)

ps8
ggsave(ps8, file="Figure_sine_s8-tree.pdf", height=6, width=6)


for(i in 1:length(ress2$Length)){
  
  tmp_fasta <- sfa[match(ress2$ID[i],names(sfa))]
  write.fasta(sequences = tmp_fasta,names=ress2$ID[i],
              file.out = paste("./data/tmp/sine_sub/",ress2$ID[i],".fasta",sep=""))
}

options(scipen=200)
for(i in 1:length(ress2$Length)){
  
  subf <- aout2[which(ress2$ID[i]==aout2$tname),c(5:7)]
  write.table(subf,file=paste("./data/tmp/sine_sub/",ress2$ID[i],".bed",sep=""),row.names = F,col.names = F,
              quote = F,sep="\t")
}

#for i in $(ls *.bed)
#do
#    seqkit subseq --bed $i -u 0 -d 0 -o $i.fa.seq ../../afa_chr_ngap1.fa 
#done

#for i in $(ls *.fasta)
#do
#    cat $i ${i%.*}.bed.fa.seq > ${i%.*}.cb 
#done


#for i in $(ls *.cb)
#do
#   mafft --reorder --thread -1 $i > $i.maf &
#done

#for i in $(ls *.maf)
#do
#   snp-sites $i -v > ${i%.*}.vcf
#   vcftools --vcf ${i%.*}.vcf --freq --out ${i%.*}
#done


f <- read.table("./data/tmp/sin_cb.frq")
#pro <- c()
#pro1 <- c()
#for(j in 1:dim(f)[1]){
#  a <- read.table(as.character(f[j,1]),fill=T,sep=" ",header=T)
#  tmpp <- 0
#  for(i in 1:dim(a)[1]){
#    atmp <- strsplit(as.character(a[i,]),"\t")[[1]][-c(1:4)]
#    aat <- c()
#    for(k in 1:length(atmp)){
#      atmp1 <- strsplit(atmp[k],":")[[1]]
#      aat <- rbind(aat,atmp1)
#    }
#    if(any(aat[,1]=="*")){
#      ii <- aat[which(aat[,1]=="*"),2]
#      if(as.numeric(ii)<0.85){
tf1 <- as.numeric(aat[-which(aat[,1]=="*"),2])
tfp <- tf1/sum(tf1)
if(max(tfp)<0.9){
  tmpp <- tmpp + 1
}
}
#    }else{
#      tf1 <- as.numeric(aat[,2])
#      tfp <- tf1/sum(tf1)
#      if(max(tfp)<0.9){
#        tmpp <- tmpp + 1
#      }
#    }

#  }
#  pro <- c(pro,tmpp)
#  pro1 <- c(pro1,dim(a)[1])
#}

#pr <- cbind(pro,pro1)
#write.table(pr,file="sine_pro.txt") 

pros <- read.table("data/tmp/sine_pro.txt",header=T)
all_ms <- unlist(strsplit(f[,1],".cb.frq"))

nress2 <- ress2[match(all_ms,ress2$ID),]
poyrs <- pros[,1]/pros[,2]

ms <- c();ss <- c()
for(i in 1:length(st)){
  r1s <- poyrs[which(nress2$Class==st[i])]
  ms <- c(ms,mean(r1s))
  ss <- c(ss,sd(r1s))
}


