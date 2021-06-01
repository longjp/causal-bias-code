## test iv
rm(list=ls())
setwd("../analysis")
load("0-data-setup.RData")
ls()
dim(ko)
dim()

## set up data
exp <- c(rep("obs",nrow(obs)),rep("int",nrow(ko)))
X <- rbind(obs,ko)
m1 <- matrix(0,nrow=nrow(obs),ncol=nrow(ko))
m2 <- diag(nrow(ko))
E <- rbind(m1,m2)




Xsd <- apply(X,2,sd)
hist(Xsd,breaks=100)



## make predicted X
Xobs <- X[exp=="obs",]
alpha0 <- colMeans(Xobs)
sd0 <- apply(Xobs,2,sd)
X <- t(t(X) - alpha0)
lambda <- 4*sd0
## OLS fit    
Dhat <- t(E)%*%X
nk <- pmax(colSums(E),1)
Dhat <- Dhat/nk
## Lasso shrinkage
Dsign <- 2*(Dhat>0)-1
sub <- (matrix(1/nk,ncol=1)%*%matrix(lambda,nrow=1))/2
Dhat <- Dsign*pmax(abs(Dhat) - sub,0)
Xhat <- E%*%Dhat

Xsd_hat <- apply(Xhat,2,sd)


plot(Xsd,Xsd_hat)
abline(a=0,b=1)
hist(Xsd_hat/Xsd)



trim_sd <- function(x){
  quants <- quantile(x,c(0.01,0.99))
  ix <- (x>quants[1]) & (x<quants[2])
  return(sd(x[ix]))
}

Xsd_trim <- apply(X,2,trim_sd)
Xsd_hat_trim <- apply(Xhat,2,trim_sd)

par(mfcol=c(1,2))
hist(Xsd_hat/Xsd)
hist(Xsd_hat_trim/Xsd_trim)

dev.off()
hist((Xsd_hat/Xsd)^2)


ix <- which.max(Xsd_hat/Xsd)
Xsd_hat[ix]
Xsd[ix]
cols <- c("obs"=1,"int"=2)
plot(X[,ix],Xhat[,ix],col=cols[exp])
abline(a=0,b=1)
