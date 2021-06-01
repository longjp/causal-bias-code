### generate lsf submission files
### each has a different seed for the boostrap sample
rm(list=ls())
set.seed(12345)

## load arguments for full run or quick run
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=1){
  stop("Must have exactly one command line argument.")
} else{
  source(args[1])
}

load("0-params.RData")

## Generate text for .lsf files
GenLSF <- function(ii,rseed){
  text1 <- '#BSUB -J 1-'
  text2 <- '
  #BSUB -W 00:10
  #BSUB -q short 
  #BSUB -n 28
  #BSUB -M 130
  #BSUB -R rusage[mem=130]
  #BSUB -N
  #BSUB -B
  #BSUB -u jplong@mdanderson.org
  Rscript 1-fit.R args_hpc.R '
  text3 <- ' > 1-'
  text4 <- '.out'
  out <- paste0(text1,
                ii,
                text2,
                ii,
                " ",
                rseed,
                text3,
                ii,
                text4)
  fname <- paste0("1-",ii,"_",rseed,".lsf")
  return(list(fname=fname,out=out))
}

## creates .lsf files
rseeds <- sample(1:10^8,length(params_hpc))
for(ii in 1:length(params_hpc)){
    lsf <- GenLSF(ii,rseeds[ii])
    cat(lsf[[2]],file=lsf[[1]])
    ## submit files for running on server
    ## only don't submit if testing this script
    if(!FAST){
      system(paste0("bsub < ",lsf[[1]]))
    }
}