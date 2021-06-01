## find set of 4 or so tightly correlated genes
## use knockouts to see if can find causal structure
## among them


## summary of findings
## 1) found tightly correlated group of 6 mrna
## 2) none of these mrna is causally influcing others (when ko 1, does not effect any others)
##     I think ICP would return null if put all 6 in model, but if put just 2 then might assign causality
## 3) can search for kos which significantly effected all 6 and are decently correlated with them
## 4) there are several and these are presumably causes of them
## QUESTION: is YGL246C a good causal candidate for other 6 because some kos only effect
##           other 6 and no YGL246C suggesting that 6 are downstream from YGL246C
## TODO: clean up this note

rm(list=ls())
load("0-data-setup.RData")
library(ggplot2)
library(GGally)
ls()

rownames(ko) <- make.names(rownames(ko))
colnames(ko) <- make.names(colnames(ko))
rownames(ko_sie) <- make.names(rownames(ko_sie))
colnames(ko_sie) <- make.names(colnames(ko_sie))
rownames(obs) <- make.names(rownames(obs))
colnames(obs) <- make.names(colnames(obs))


dim(obs)

obs_sub <- obs[,rownames(ko)]
obs_corr <- cor(obs_sub)
dim(obs_corr)
sum(abs(obs_corr)>0.8)


## find set of genes which are highly correlated with each other
gene <- which(colSums(abs(obs_corr)>0.9)>5)[1]
gr <- rownames(obs_corr)[obs_corr[,gene]>0.9]

## make scatterplot matrix. it does not appear that any of
## these genes are regulating each other. the knockout for
## gene X does not produce large effects for the other genes
oko <- rbind(ko,obs)
oko <- as.data.frame(oko)
oko$cl <- as.factor(c(rep("ko",nrow(ko)),rep("obs",nrow(obs))))
ix <- which(colnames(oko) %in% gr)
ggpairs(oko,columns=ix,mapping=aes(colour=cl,alpha=0.3))


## find genes which have a large causal impact on all of these genes
## and add a few to scatterplot
causes <- rownames(ko_sie)[rowSums(ko_sie[,gr])==length(gr)]
gr_causes <- c(gr,causes[c(1,3)])
ix <- which(colnames(oko) %in% gr_causes)
ggpairs(oko,columns=ix,mapping=aes(colour=cl,alpha=0.3))
### TODO: for YAL021C on other genes, there are a lot of genes
###     effecting Y, but not a lot effecting X. is this evidence that X->Y


## now find genes with a large causal impact and which are 
## highly correlated with targets
temp <- cor(obs[,c(causes,gr)])
temp <- temp[causes,gr]
causes_corr <- rownames(temp)[rowSums(abs(temp) > 0.5)>=5]
causes_corr

gr_causes <- c(gr,causes_corr[3:4])
gr_causes

## genes with some causes which are strongly correlated with result
ix <- which(colnames(oko) %in% gr_causes)
ggpairs(oko,columns=ix,upper=list("continuous"="points"),
        mapping=aes(colour=cl,alpha=0.3))
gr_causes
causes_corr
## the two genes "YGL078C" "YGL246C" exert significant effects on the other genes
## but do not exert significant effects on each other
## they are decently correlated with other genes, but probably wouldn't pass any prescreening





#### LATER CODE
### view particular candidate
all_data <- rbind(ko,obs)
cols <- c(rep(rgb(0, 0, 0.9,alpha=0.3),nrow(ko)),
          rep(rgb(0.9, 0, 0,alpha=0.3),nrow(obs)))
## check specific candidate
xname <- "YGL037C"
yname <- "YHR087W"
plot(all_data[,xname],all_data[,yname],
     col=cols,xlab=xname,ylab=yname)
abline(h=range(obs[,yname]),col='grey')
abline(v=range(obs[,xname]),col='grey')
points(ko[xname,xname],ko[xname,yname],col='orange',cex=1.3,lwd=2,pch=3)
if(yname %in% rownames(ko)){
  points(ko[yname,xname],ko[yname,yname],col='green',cex=1.3,lwd=2,pch=3)
}


