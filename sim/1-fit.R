## do not run directly, run 1-lsf.R instead. see README
load("0-params.RData")
source("funcs.R")

## load arguments
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=3){
  stop("Must have exactly three command line argument.")
} else{
  source(args[1])
  ii <- as.numeric(args[2])
  rseed <- as.numeric(args[3])
}

set.seed(rseed)
rseeds <- sample(1:10^8,N)
res <- mclapply(1:N,FUN=RunSim,
                params=params_hpc[[ii]],rseeds=rseeds,
                mc.cores=mc.cores)
res <- lapply(res,function(x){x[[1]]})
cnames <- names(res[[1]])
res <- matrix(unlist(res),ncol=length(res[[1]]),byrow=TRUE)
colnames(res) <- cnames
save(res,file=paste0("1-",ii,".RData"))