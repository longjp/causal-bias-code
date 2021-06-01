CVModel <- function(predictX,obs,ko,ix,rseed=rseed,mc.cores=4){
  folds <- unique(ix)
  out_predictX <- list()
  for(ii in 1:length(folds)){
    ko_train <- ko[ix!=ii,]
    ko_test_names <- rownames(ko)[ix==ii]
    out_predictX[[ii]] <- predictX(obs,ko_train,ko_test_names,rseed,mc.cores)
  }
  # ## extract predictions, order to match ko
  # ko_pred <- lapply(out_predictX,function(x){x[[1]]})
  # ko_pred <- do.call(rbind,ko_pred)
  # ko_pred <- ko_pred[rownames(ko),]
  # ## extract diagnostic / tuning info returned by method, do not post process
  # out <- lapply(out_predictX,function(x){x[[2]]})
  return(out_predictX)
}

TSLSreg <- function(X,Y,E){
  ExpInd <- rowSums(E)
  Xobs <- X[ExpInd==0,,drop=FALSE]
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
  ## unregularized second stage
  Xhati <- cbind(1,Xhat)
  coefs <- lm.fit(Xhati,Y)$coefficient[-1]
  return(coefs)
}