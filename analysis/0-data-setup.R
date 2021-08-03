## divide kemmeren data into CV splits
## determine definition of strong effect
rm(list=ls())
set.seed(1234)
load("../data/Kemmeren.RData")

## load arguments
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=1){
  stop("Must have exactly one command line argument.")
} else{
  source(args[1])
}



## change names
obs <- data$obs
ko <- data$int

## strong interventional effects (sie)
## defined as intervention resulting in expression
## outside the range of what was observed in observational data
Fn <- apply(obs,2,ecdf)
dat_norm <- ko
for(ii in 1:ncol(dat_norm)){
  dat_norm[,ii] <- Fn[[ii]](ko[,ii])
}
ko_sie <- dat_norm == 0 | dat_norm == 1
##hist(rowSums(ko_sie))
##hist(colSums(ko_sie))


#### alternative definition of sie from pnas article
ko_sie_pnas <- ko

## for intervention on X_j, replace X_j expression
## with mean expression accross all samples
## this prevents intervention the sie of ko ii being gene ii
for(ii in 1:nrow(ko_sie_pnas)){
  ko_sie_pnas[ii,rownames(ko_sie_pnas)[ii]] <- mean(ko_sie_pnas[,rownames(ko_sie_pnas)[ii]])
}
## for each gene, find minimum and maximum expression in ko
## replace min with -1 and max with 1, all other 0
for(ii in 1:ncol(ko_sie_pnas)){
  ix_min <- which.min(ko_sie_pnas[,ii])
  ix_max <- which.max(ko_sie_pnas[,ii])
  ko_sie_pnas[,ii] <- 0
  ko_sie_pnas[ix_min,ii] <- -1
  ko_sie_pnas[ix_max,ii] <- 1
}
## SIE condition i in PNAS
## if a ko was not successful (if X_j ko and X_j does not
## have lowest expression among all obs/int data on X_j ko)
## then remove all sie assigned to this X_j
success <- rep(NA,nrow(ko_sie_pnas))
for(ii in 1:nrow(ko_sie_pnas)){
  gname <- rownames(ko_sie_pnas)[ii]
  gexp <- ko[ii,gname]
  ## if knockout was not successful
  ## do not assign sie
  if(!all(ko[,gname] >= gexp)){
    ko_sie_pnas[ii,] <- 0
  }
}



## run with 300 genes, about 20x less data. useful for testing
if(FAST){
  ix_col <- sample(colnames(obs),300)
  obs <- obs[,colnames(obs) %in% ix_col]
  ko <- ko[rownames(ko) %in% ix_col,colnames(ko) %in% ix_col]
  ko_sie <- ko_sie[rownames(ko_sie) %in% ix_col,colnames(ko_sie) %in% ix_col]
  ko_sie_pnas <- ko_sie_pnas[rownames(ko_sie_pnas) %in% ix_col,
                             colnames(ko_sie_pnas) %in% ix_col]
}

## create folds for cross validation
ix <- sample(1:3,nrow(ko),replace=TRUE)

save(obs,ko,ko_sie,ko_sie_pnas,ix,file="0-data-setup.RData")
