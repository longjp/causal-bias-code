load("0-data-setup.RData")
source("predictEST.R")
source("funcs.R")

## load arguments
args = commandArgs(trailingOnly=TRUE)
# use laptop (fast) or hpc (slow) arguments
if(length(args)==0){
  stop("Must have exactly at least one command line argument.")
} else{
  source(args[1])
}
# set the seed
if(length(args)==2){
  rseed <- args[2]
} else {
  rseed <- 1234
}

pred <- CVModel(predictEST,obs,ko,ix,rseed=rseed,mc.cores=mc.cores)
save(pred,file=paste0("1-fit_",rseed,".RData"))