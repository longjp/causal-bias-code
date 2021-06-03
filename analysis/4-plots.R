rm(list=ls())
load("0-data-setup.RData")
library(ggplot2)
set.seed(1234)

all_data <- rbind(ko,obs)

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