## This code is part of the study:
## Stochastic dispersal rather than deterministic selection explains the spatio-temporal distribution of soil bacteria in a temperate grassland
## doi: 10.3389/fmicb.2020.01391
## ©Franz-Sebastian Krah
## 02 - 27 - 2019

#' @title raup_crick_abu_par
#' @param com community matrix (spXsite)
#' @param reps number of bootstraps
#' @param ncore number of cores (serial: ncore = 1; parallel > 1)
#' @param classic_metric standardizes the metric to range from -1 to 1
#' @param split_ties adds half of the number of null observations that are equal to the observed number of shared species to the calculation- this is highly recommended
#' @details Parallelized version of the Raup-Crick algorithm for "abundance" data (Stegen et al. 2013).
#' Previous code loops over each pairwise community combination;
#' here we randomize the full community matrix and compute Bray-Curtis for the 
#' full matrix and then conduct subsequent Raup-Crick calculations as in Stegen.
#' This increaes computational speed. Further here implemented as multi-core version.
#' The code is acurate and fast (see paper, supplement section D)
#' @author Franz-Sebastian Krah

options(warn = -1) # Turn off warning
## 
#
suppressMessages(library(optparse))
# 
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="dt.txt",
                help="Otu table [default %default]"),
	make_option(c("-r", "--repeats"), type="numeric", default="999",
                help="Repeat number [default %default]"),
	make_option(c("-n", "--ncore"), type="numeric", default="10",
                help="Core number [default %default]"),			
    make_option(c("-o", "--output"), type="character", default="result",
                help="Output RC [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
}

comun1=read.csv(opts$input,head=T,stringsAsFactors=F,row.names=1,sep = "\t")
comun=t(comun1)



raup_crick_abu_par <- function(com, reps, ncore, classic_metric=FALSE, split_ties=TRUE){
  

  require("parallel")
  require("doSNOW")
  require("vegan")
  pb <- txtProgressBar(max =reps, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  cl <- makeCluster(ncore)
  registerDoSNOW(cl)
  
  bray.rand <- foreach(randomize = 1:reps, 
    .options.snow = opts,
    .packages = c("vegan", "picante")) %dopar% {
      
      
      null.dist <- com*0
      
      for(i in 1:nrow(com)){

        com.pa <- (com>0)*1
        gamma<-ncol(com)
        occur<-apply(com>0, MARGIN=2, FUN=sum)
        abundance<-apply(com, MARGIN=2, FUN=sum)
        com1 <- rep(0,gamma)
        
        com1[sample(1:gamma, sum(com.pa[i,]), replace=FALSE, prob=occur)]<-1
        com1.samp.sp = sample(which(com1>0), (sum(com[i,])-sum(com1)),
          replace=TRUE,prob=abundance[which(com1>0)]);
        com1.samp.sp = cbind(com1.samp.sp,1)
        com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2],com1.samp.sp[,1],FUN=sum))
        colnames(com1.sp.counts) = 'counts'
        com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts))
        com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts
        x <- com1
        null.dist[i,] <- x
        rm('com1.samp.sp','com1.sp.counts')
      }
      as.matrix(vegdist(null.dist, "bray"))
    }
  stopCluster(cl)
  
  ## Calculate beta-diversity for obs metacommunity
  bray.obs <- as.matrix(vegdist(com, "bray"))
  
  ##how many null observations is the observed value tied with?
  null_bray_curtis <- bray.rand
  num_exact_matching_in_null <- lapply(null_bray_curtis, function(x) x==bray.obs)
  num_exact_matching_in_null <- apply(simplify2array(num_exact_matching_in_null), 1:2, sum)
  
  ##how many null values are smaller than the observed *dissimilarity*?
  num_less_than_in_null <- lapply(null_bray_curtis, function(x) (x<bray.obs)*1)
  num_less_than_in_null <- apply(simplify2array(num_less_than_in_null), 1:2, sum)
  
  
  rc = (num_less_than_in_null)/reps; # rc;
  
  if(split_ties){
    
    rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
  };
  
  
  if(!classic_metric){
    
    ##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
    
    rc = (rc-.5)*2
  };
  
  return(rc)
  
}

RC=raup_crick_abu_par(comun,opts$repeats, opts$ncore)
write.csv(as.matrix(RC),paste(opts$output,"_RC.csv",sep = ""))







