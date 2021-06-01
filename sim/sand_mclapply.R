# - why does run sim have RunSim <- function(...,params){
#   - why not RunSim <- function(ii,params){
#     - test this out
#     - also get simulations to set seed to entirely reproducible
#     - improve output of .Rmd simulation code
#     
#     - currently kemmeren data will not produce same output because seed is set
#     and then parallel processes begin. figure out how to make this work

rm(list=ls())
library(parallel)
set.seed(1234)


RunSim <- function(ii,A,B,rseeds){
  set.seed(rseeds[ii])
  print(paste0("run: ",ii," B=",B))
  #Sys.sleep(runif(1))
  return(runif(1))
}



runif(1)
A <- 30
B <- 50
mc.cores <- 2
N <- 100
rseeds <- sample(1:10^8,N)
out <- mclapply(1:N,RunSim,A=A,B=B,rseeds=rseeds)

unlist(out)
