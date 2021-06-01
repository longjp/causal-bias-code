rm(list=ls())
load("0-data-setup.RData")
library(ggplot2)
set.seed(1234)

all_data <- rbind(ko,obs)

cols <- c(rep(rgb(0, 0, 0.9,alpha=0.3),nrow(ko)),
          rep(rgb(0.9, 0, 0,alpha=0.3),nrow(obs)))

cl <- c(rep("ko",nrow(ko)),rep("obs",nrow(obs)))


dat <- data.frame(all_data,"Environment"=cl)
xname <- c("YLL019C")
yname <- c("YLL020C")
p <- ggplot(dat,aes_string(x=xname,y=yname,color="Environment")) +
  geom_point(alpha=0.3)
pdf("../ms-bias/figs/gene_pair_ex1.pdf",width=5,height=4)
plot(p)
dev.off()

xname <- "YKL043W"
yname <- "YKL044W"
p <- ggplot(dat,aes_string(x=xname,y=yname,color="Environment")) +
  geom_point(alpha=0.3)
pdf("../ms-bias/figs/gene_pair_ex2.pdf",width=5,height=4)
plot(p)
dev.off()



# 
# 
# ## check top buhlmann candidates
# xgenes <- c("YMR104C","YPL273W","YCL040W","YLL019C","YMR186W")
# ygenes <- c("YMR103C","YMR321C","YCL042W","YLL020C","YPL240C")
# ii <- 0
# 
# 
# 
# ii <- ii + 1
# xname <- xgenes[ii]
# yname <- ygenes[ii]
# pdf("buhl_candidate.pdf",width=6,height=4)
# plot(all_data[,xname],all_data[,yname],
#      col=cols,xlab=xname,ylab=yname)
# abline(h=range(obs[,yname]),col='grey')
# abline(v=range(obs[,xname]),col='grey')
# points(ko[xname,xname],ko[xname,yname],col='orange',cex=1.3,lwd=2,pch=3)
# if(yname %in% rownames(ko)){
#   points(ko[yname,xname],ko[yname,yname],col='green',cex=1.3,lwd=2,pch=3)
# }
# dev.off()
# 
# 
# 
# 
# 
# ## fit ICP model to these data
# 
# library(InvariantCausalPrediction)
# 
# 
# 
# ii <- 0
# 
# 
# xgenes
# ygenes
# 
# ii <- ii + 1
# ii
# xgenes[ii]
# ygenes[ii]
# ## original preds
# Y <- all_data[,ygenes[ii]]
# X <- all_data[,xgenes[ii],drop=FALSE]
# ExpInd <- c(rep(1,nrow(ko)),rep(2,nrow(obs)))
# ICP(X,Y,ExpInd=ExpInd)
# ## now reverse
# X <- all_data[,ygenes[ii],drop=FALSE]
# Y <- all_data[,xgenes[ii]]
# ExpInd <- c(rep(1,nrow(ko)),rep(2,nrow(obs)))
# ICP(X,Y,ExpInd=ExpInd)
# 
# 
# 
# ###
# xgenes %in% rownames(ko)
# ygenes %in% rownames(ko)
# 
# 
