library(parallel)
## causalDantzig was removed from current version of InvariantCausalPrediction package. why?
## here we install (if uncomment two lines) / load old version 
old_lib <- "/Users/jplong/Library/R/3.5/library_old/"
##library(remotes)
##install_version("InvariantCausalPrediction", version = "0.7-1", lib = old_lib)
library("InvariantCausalPrediction",lib.loc=old_lib)


n <- 100
X <- matrix(rnorm(n*1000),nrow=n)
Y <- rnorm(n)
EnvInd <- c(rep(1,n/2),rep(2,n/2))
fit <- try(causalDantzig(X,Y,EnvInd),silent=TRUE)
class(fit)



n <- 100
X <- matrix(rnorm(n*10),nrow=n)
Y <- rnorm(n)
EnvInd <- c(rep(1,n/2),rep(2,n/2))
causalDantzig(X,Y,EnvInd)
fit <- try(causalDantzig(X,Y,EnvInd),silent=TRUE)
class(fit)


?try

