rm(list=ls())
library(tidyr)
library(kableExtra)
library(ggplot2)
library(GGally)
source('funcs.R')
load("0-params.RData")
set.seed(1)

ix <- which(sim_params_des$ntop==320 & sim_params_des$p0==1)
sim_params_des[ix,]

params <- list()
for(ii in 1:length(ix)){
  params[[ii]] <- params_hpc[[ix[ii]]]
}

dat <- list()
for(ii in 1:length(ix)){
  dat[[ii]] <- GenerateData(params[[ii]]$pp,
                            params[[ii]]$p0,
                            params[[ii]]$nobs,
                            params[[ii]]$nko,
                            params[[ii]]$ko_strength,
                            params[[ii]]$np,
                            params[[ii]]$S2,
                            params[[ii]]$st1,
                            params[[ii]]$st2,
                            params[[ii]]$ntop)
}


## generate data for plotting
Xplot <- list()
for(ii in 1:length(dat)){
  X <- dat[[ii]]$X
  #A <- dat$Abeta$A
  #beta <- dat$Abeta$beta
  pp <- 6400
  p0 <- sim_params_des$p0[ix[ii]]
  E <- dat[[ii]]$ER$E
  ExpInd <- 1*(rowSums(abs(E))!=0)
  colnames(X) <- c(paste0("X",1:(pp/2)),
                   "Y",
                   paste0("X",
                          (pp/2+1):(pp)))
  ix2 <- pp/2+1
  Y <- X[,ix2]
  X <- X[,-ix2]
  X <- cbind(X1=X[,1],X[,(pp/2-p0+1):(pp/2),drop=FALSE],Y,X[,(pp/2+1):(pp/2+p0),drop=FALSE])
  ExpIndRecode <- c("Obs","KO")
  X <- data.frame(X,Type=ExpIndRecode[ExpInd+1],
                  Strength=sim_params_des$st[ix[ii]],
                  stringsAsFactors=FALSE)
  X <- X[nrow(X):1,]
  Xplot[[ii]] <- X
}


Xplot <- do.call(rbind,Xplot)


Xplot_new <- pivot_longer(Xplot,which(grepl("X",colnames(Xplot))))
p <- ggplot(Xplot_new,aes(x=value,y=Y,color=Type)) +
  geom_point(alpha=0.3) + labs(x="x",y="y") +
  facet_grid(vars(Strength),vars(name))
pdf("../ms-bias/figs/sim_scatter.pdf",width=10,height=5)
plot(p)
dev.off()

