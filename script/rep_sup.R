read_rm <- function (file, tibble = FALSE, keep_order = TRUE, include_secondary = FALSE) 
{
  lines <- readLines(file)
  lines <- str_trim(lines[4:length(lines)], "left")
  raw_data <- data.frame(str_split(lines, "\\s+", simplify = TRUE))
  by_strand <- split(raw_data, raw_data[, 9])
  comp_ali <- by_strand[["C"]][c(1:11, 14, 13, 12, 15, 16)]
  col_names <- c("score", "p_sub", "p_del", "p_ins", "qname", 
                 "qstart", "qend", "qextend", "complement", "tname", "tclass", 
                 "tstart", "tend", "textend", "ID", "ali_type")
  names(comp_ali) <- col_names
  names(by_strand[["+"]]) <- col_names
  res <- rbind.data.frame(comp_ali, by_strand[["+"]])
  if (keep_order) {
    res <- res[order(as.numeric(rownames(res))), ]
  }
  res$ali_type <- ifelse(res$ali_type == "*", "secondary", 
                         "primary")
  if (!include_secondary) {
    res <- subset(res, ali_type == "primary")
  }
  res <- .numerify(res, c(1:4, 6, 7, 12:13))
  res <- .de_paren(res, c("textend", "qextend"))
  if (keep_order) {
  }
  class(res) <- c("repeat_table", "data.frame")
  rownames(res) <- NULL
  res
}

.de_paren <- function(df, column_indices) {
  for(i in column_indices){
    df[[i]] <-  as.numeric(gsub("\\(|\\)", "", df[[i]]))
  }
  df
}


.numerify <- function(df, column_indices){
  for(i in column_indices){
    df[[i]] <- as.numeric(as.character(df[[i]]))
  }
  df
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



rep_density_single <- function(gff,bin=250000){
  
  
  library(tidyverse)
  library(purrr)
  
  chr_end <- max(gff[,6])
  tmp <- data.frame(table(cut(gff[, 6],
                              breaks = c(
                                seq(0, chr_end, bin),
                                chr_end
                              ))))
  
  tmp <-tidyr::separate(tmp,
                        1,
                        into = c("Start", "End"),
                        sep = ",")
  tmp$Start <- as.numeric(gsub("\\(", "", tmp$Start)) + 1
  tmp$End <-as.numeric(gsub("\\]", "", tmp$End))
  
  # modify last end_pos to chr_end
  tmp[nrow(tmp), "End"] <- chr_end
  
  
  #combine
  tmp1 <- cbind(tmp,1:dim(tmp)[1])
  tmp2 <- tmp1[-which(tmp1[,3]==0),]
  
  df <- tmp2[,4] %>% tibble(dat = .)
  a <- df %>%group_by(seq_id = cumsum(dat != lag(dat) + 1 | is.na(dat != lag(dat) + 1))) %>%nest()
  b <- sapply(a$data, "[[","dat")
  
  rett <- c()
  for( i in 1:length(b)){
    tc <- c(min(tmp1[b[[i]],1]),max(tmp1[b[[i]],2]),sum(tmp1[b[[i]],3]))
    rett <- rbind(rett,tc)
  }
  res <- list(tmp=tmp,rett=rett)
  return(res)
}



rep_gc <- function(karyotype,gff,egc,bin=250000){
  
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
      sv <- c(sv,mean(egc[sss,2]))
    }
    tmp <- cbind(tmp,Size=sv)
    # combine
    density <- rbind(density, tmp)
  }
  
  return(density)
}





rep_dist <- function(karyotype,gff,bin=c(0,200,500,1000,2000,5000,10000),dist_tre=10,repn=4){
  
  library(dplyr)
  library(tidyverse)
  library(purrr)
  # select only interested feature
  all_chrs <- karyotype[,1]
  
  
  # make output data frame
  distall <- list()
  distdis <- c()
  distdisp <- c()
  nblock <- list()
  nbgff <- list()
  # loop against each chromosome
  i <- 1
  for (chr in all_chrs) {
    chr_end <- karyotype[karyotype[,1] == chr, 2]
    # count, but only consider the start of each gene
    gff1 <- gff[gff$qname == chr,]
    nrep <- dim(gff1)[1]
    dvs <- gff1$qstart[2:nrep]-gff1$qend[1:(nrep-1)]
    dvs[which(dvs<0)] <- 0
    
    distall[[i]] <- dvs
    dvfre <- table(cut(dvs,breaks=c(bin,max(dvs))))
    dvfrep <- dvfre/sum(dvfre)
    distdis <- rbind(distdis,dvfre)
    distdisp <- rbind(distdisp,dvfrep)
    
    st <- which(dvs<dist_tre)
    
    df <- st %>% tibble(dat = .)
    a <- df %>%group_by(seq_id = cumsum(dat != lag(dat) + 1 | is.na(dat != lag(dat) + 1))) %>%nest()
    
    b <- sapply(a$data, "[[","dat")
    
    c <- keep(b,function(x) length(x) >repn)
    
    if(length(c)>0){
      nbgft <- c()
      cb <- c()
      for(k in 1:length(c)){
        cb <- c(cb,gff1$qend[max(c[[k]])]-gff1$qstart[min(c[[k]])-1])
        nbgft <- rbind(nbgft,gff1[(min(c[[k]])-1):max(c[[k]]),])
      }
      nbgff[[i]] <- nbgft
      nblock[[i]] <- cb
    }else{
      nbgff[[i]] <- NULL
      nblock[[i]] <- NULL
    }
    
    i <- i + 1
  }
  rownames(distdis) <- rownames(distdisp) <- karyotype[,1]
  final <- list(distall=distall,distdis=distdis,distdisp=distdisp,nblock=nblock,nbgff=nbgff)
  return(final)
}


rep_density <- function(karyotype,gff,bin=250000){
  
  
  # select only interested feature
  all_chrs <- karyotype[,1]
  
  
  # make output data frame
  density <- c()
  
  # loop against each chromosome
  for (chr in all_chrs) {
    chr_end <- karyotype[karyotype[,1] == chr, 2]
    # count, but only consider the start of each gene
    tmp <- data.frame(table(cut(gff[gff$qname == chr, 6],
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
    
    # combine
    density <- rbind(density, tmp)
  }
  
  return(density)
}





rep_nl <- function(karyotype,gff){
  
  
  
  # select only interested feature
  all_chrs <- karyotype[,1]
  
  
  # loop against each chromosome
  lc <- c()
  lm <- c()
  for (chr in all_chrs) {
    ti <- which(gff$qname==chr)
    
    lc <- c(lc,length(ti))
    lm <- c(lm,sum((gff$qend-gff$qstart)[ti])/10^6)
  }
  pr <- lm/(karyotype[,2]/10^6)
  list(lc=lc,lm=lm,pr=pr)
}

GOlevel<-function (onto, level){
  {
    GOChildren <- function (ontoTab, terms) {
      children <- ontoTab[ontoTab[,2] %in% terms,1]
    }
    
    if (!((onto == "MF") || (onto == "BP") || (onto == "CC"))) {
      on.exit(cat("You must enter MF, BP or CC as ontology"))
      GOLevel <- NULL
    }
    else {
      switch(onto, MF = {
        topNode <- "GO:0003674"
        chTab<-toTable(GOMFCHILDREN)
      }, BP = {
        topNode <- "GO:0008150"
        chTab<-toTable(GOBPCHILDREN)
      }, CC = {
        topNode <- "GO:0005575"
        chTab<-toTable(GOCCCHILDREN)
      })
      
      GOLevel <- topNode
      if (level < 1) {
        on.exit(cat("You must enter a level greater than 1"))
        GOLevel <- NULL
      }
      else {
        if (level==1) {
          GOLevel <- topNode
        }
        else{
          GOLevel<-topNode
          for (i in (1:(level - 1))){
            # childLevel <- chTab[chTab[,2] %in% GOLevel,1]
            childLevel <- GOChildren(chTab, GOLevel)
            GOLevel <- childLevel
            # GOLevel <- GOLevel[!is.na(GOLevel)] # No se si cal?
          }
        }
      }
      
      return(unique(GOLevel))
    }
  }
}