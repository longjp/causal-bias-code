### generate lsf submission files
### each has a different seed for the boostrap sample
rm(list=ls())
set.seed(1234)

## load arguments for full run or quick run
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=1){
  stop("Must have exactly one command line argument.")
} else{
  source(args[1])
}

## Generate text for .lsf files
GenLSF <- function(rseed,method){
  text1 <- '#BSUB -J 1-fit_'
  text2 <- '
#BSUB -W 01:00
#BSUB -q short 
#BSUB -n 28
#BSUB -M 20
#BSUB -R rusage[mem=20]
#BSUB -N
#BSUB -B
#BSUB -u jplong@mdanderson.org
Rscript 1-fit.R args_hpc.R '
  text3 <- ' > 1-fit_'
  text4 <- '.out'
  out <- paste0(text1,
                rseed,
                text2,
                rseed,
                text3,
                rseed,
                text4)
  fname <- paste0("1-fit_",rseed,".lsf")
  return(list(fname=fname,out=out))
}

## creates .lsf files
rseeds <- sample(1:10^8,Nboot)
for(ii in 1:Nboot){
    lsf <- GenLSF(rseeds[ii])
    cat(lsf[[2]],file=lsf[[1]])
    ## submit files for running on server
    ## only don't submit if testing this script
    if(!FAST){
      system(paste0("bsub < ",lsf[[1]]))
    }
}


