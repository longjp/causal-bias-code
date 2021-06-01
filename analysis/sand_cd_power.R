## cdantzig has low power when choose many predictors
## maybe can be improved with better environment allocation
## this would be good because can do less aggressive lasso prescreen
## Question: is regularized causal dantzig effective here
rm(list=ls())
library(InvariantCausalPrediction)
load("0-data-setup.RData")

ls()
ICP


xname <- "YPL273W"
yname <- "YMR321C"

fold <- ix[which(rownames(ko)==xname)]
ko_train <- ko[ix!=fold,]
ko_test_names <- rownames(ko)[ix==fold]





dat <- rbind(obs,ko_train)
sum(colnames(dat)==yname)
ExpInd <- c(rep("obs",nrow(obs)),rep("ko",nrow(ko_train)))


##ix <- rownames(dat) != yname
##exp_sub <- exp[ix]
##X <- dat[ix,]


Y <- dat[,yname]
X <- dat[,colnames(dat)!=yname]

fit <- ICP(X,Y,ExpInd)
fit




cvfit <- cv.glmnet(X, Y)
a <- coef(cvfit, s = "lambda.1se")
a <- as.matrix(a)
a <- rownames(a)[a[,1]!=0]
to_use <- colnames(X)[colnames(X) %in% a]
if(length(to_use)>0){
  X_sub <- X[,to_use,drop=FALSE]
  ### run CD on selected
  ## if system is singular, say nothing is significant
  fit <- try(causalDantzig(X_sub,Y,ExpInd=ExpInd),silent=TRUE)
  if(class(fit)!="try-error"){
    pvals <- fit$p_value
  } else {
    pvals <- NULL
  }
  X_sub <- NULL
} else{
  pvals <- NULL
}
pvals
min(pvals)


## now try regularized causal dantzig
## memory explosion, does not work
##fit <- causalDantzig(X,Y,ExpInd,regularization=TRUE,mc.cores=1)

