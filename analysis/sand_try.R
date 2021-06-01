##
rm(list=ls())
library(parallel)

5/0


?try

ii <- 22
ii %% 8

Compute <- function(ii){
  if((ii %% 3)==0){
    a <- colMeans(rnorm(10))
  } else {
    a <- "np"
  }
  return(a)
}

out <- mclapply(1:10,Compute,mc.cores=3)
Compute(3)
  
