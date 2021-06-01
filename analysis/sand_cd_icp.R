## problem: causalDantzig did identify causal predictors well
##    when we used glmnet lasso variable
##    selection (CV with tuning parameter min1se)
## solution: glmnet is selecting too many variable, resulting
## in nothing being significant now have it choose top X (~5)
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




##cvfit <- cv.glmnet(X, Y)


coefs <- coef(glmnet(X, Y))
coefs <- coefs[-1,] ## remove first row
coefs_nonzero <- colSums(coefs!=0)
coefs <- coefs[,which.min(abs(coefs_nonzero-5))[1]]
to_use <- names(coefs)[coefs!=0]

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


out <- vector("list",ncol(dat))



## process output
names(out) <- colnames(dat)

## make results matrix
ko_pred <- matrix(1,
                  nrow=length(ko_test_names),
                  ncol=ncol(dat))
rownames(ko_pred) <- ko_test_names
colnames(ko_pred) <- colnames(dat)
for(ii in 1:length(out)){
  if(!is.null(out[[ii]])){
    p_sub <- out[[ii]][names(out[[ii]]) %in% rownames(ko_pred)]
    if(length(p_sub)!=0){
      ko_pred[names(p_sub),ii] <- p_sub
    }
  }
}


testmat <- matrix(1:4,nrow=2)
colnames(testmat) <- c("c","d")
rownames(testmat) <- c("a","b")
testmat
testmat[c("b","a"),1] <- c(10,5)
testmat
