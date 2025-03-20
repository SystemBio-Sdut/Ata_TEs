
meth_size <- function(karyotype,gff,bin=250000){
  
  library(dplyr)
  # select only interested feature
  all_chrs <- karyotype[,1]
  
  
  # make output data frame
  density <- c()
  
  # loop against each chromosome
  for (chr in all_chrs) {
    chr_end <- karyotype[karyotype[,1] == chr, 2]
    # count, but only consider the start of each gene
    gff1 <- gff[gff[,1] == chr,]
    tmp <- data.frame(table(cut(gff1[, 3],breaks = c(seq(0, chr_end, bin),chr_end))))
    
    # make outcome tidy
    tmp <-tidyr::separate(tmp, 1,into = c("Start", "End"),sep = ",")
    tmp$Start <- as.numeric(gsub("\\(", "", tmp$Start)) + 1
    tmp$End <-as.numeric(gsub("\\]", "", tmp$End))
    
    # add chr info
    tmp$Chr <- chr
    
    # reorder columns
    tmp <- tmp[c(4, 1:3)]
    colnames(tmp) <- c("Chr", "Start", "End", "Count")
    
    # modify last end_pos to chr_end
    tmp[nrow(tmp), "End"] <- chr_end
    density <- rbind(density, tmp)
  }
  
  return(density)
}

rep_size <- function(karyotype,gff,bin=250000){
  
  library(dplyr)
  # select only interested feature
  all_chrs <- karyotype[,1]
  
  
  # make output data frame
  density <- c()
  
  # loop against each chromosome
  for (chr in all_chrs) {
    chr_end <- karyotype[karyotype[,1] == chr, 2]
    # count, but only consider the start of each gene
    gff1 <- gff[gff$qname == chr,]
    tmp <- data.frame(table(cut(gff1[, 6],
                                breaks = c(
                                  seq(0, chr_end, bin),
                                  chr_end
                                ))))
    
    # make outcome tidy
    tmp <-tidyr::separate(tmp,
                          1,
                          into = c("Start", "End"),
                          sep = ",")
    tmp$Start <- as.numeric(gsub("\\(", "", tmp$Start)) + 1
    tmp$End <-as.numeric(gsub("\\]", "", tmp$End))
    
    # add chr info
    tmp$Chr <- chr
    
    # reorder columns
    tmp <- tmp[c(4, 1:3)]
    colnames(tmp) <- c("Chr", "Start", "End", "Count")
    
    # modify last end_pos to chr_end
    tmp[nrow(tmp), "End"] <- chr_end
    sv <- c()
    for( i in 1:dim(tmp)[1]){
      sss <- which(between(gff1$qend,tmp[i,2],tmp[i,3])==1)
      sv <- c(sv,sum(gff1$tend[sss]-gff1$tstart[sss]))
    }
    tmp <- cbind(tmp,Size=sv)
    # combine
    density <- rbind(density, tmp)
  }
  
  return(density)
}

EI <- function(para,x){
  
  A <- para[1]+para[2]*x
  return(A)
}


EI.sum <- function(para,x,y){
  B <- sum((y-EI(para,x))^2)
  return(B)
}


rep_site_cal <- function(gff1,gff2,karyotype,t=100){
  
  require(parallel)
  
  mt <- gff2
  
  a2i <- list()
  for(chrn in 1:dim(karyotype)[1]){
    
    a2i[[chrn]] <- which(mt[,1]==karyotype[chrn,1])
  }
  
  calc_m <- function(at,mt,a2i,karyotype) {
    
    nat <- dim(at)[1]
    ssum <- c()
    for(i in 1:nat){
      li <- a2i[[match(at$qname[i],karyotype[,1])]]
      ssum <- c(ssum,sum(dplyr::between(mt[li,3],at$qstart[i],at$qend[i])))
    }
    
    return(ssum)
  }
  
  chunks <- split(gff1, cut(seq_along(gff1[,1]), t, labels = FALSE))
  
  cl <- makeCluster(t)  # 这里使用 4 个线程，根据您的计算机配置可以适当调整
  
  # 在并行计算环境中计算平均值
  results <- parLapply(cl, chunks, calc_m,mt=mt,a2i=a2i,karyotype=karyotype)
  
  # 关闭并行计算环境
  stopCluster(cl)
  
  res <- as.numeric(unlist(results))
  return(res)
}