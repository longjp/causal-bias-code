library(parallel)
library(glmnet)
library(AER)
## causalDantzig was removed from current version of InvariantCausalPrediction package. 
## here we install (if uncomment two lines) / load old version 
#library(remotes)
#install_version("InvariantCausalPrediction", version = "0.7-1")
library(InvariantCausalPrediction)

predictEST <- function(obs,ko_train,ko_test_names,rseed,mc.cores=4){
  set.seed(rseed)
  dat <- rbind(obs,ko_train)
  exp <- c(rep("obs",nrow(obs)),rep("int",nrow(ko_train)))
  ## create E matrix
  m1 <- matrix(0,nrow=nrow(obs),ncol=nrow(ko_train))
  m2 <- diag(nrow(ko_train))
  E <- rbind(m1,m2)
  rownames(E) <- rownames(dat)
  ComputeFit <- function(ii,rseeds){
    set.seed(rseeds[ii])
    if(ii %% 100==0){
      print(paste0("run: ",ii,"/",N))
    }
    ## make gene ii the response, remove Y ko if present
    yname <- colnames(ko_train)[ii]
    ix <- rownames(dat) != yname
    exp_sub <- exp[ix]
    X <- dat[ix,]
    E <- E[ix,]
    Y <- X[,ii]  
    X <- X[,-ii]
    ## bootstrap sample training data
    ix <- sample(1:length(Y),replace=TRUE)
    exp_sub <- exp_sub[ix]
    Y <- Y[ix]
    X <- X[ix,]
    E <- E[ix,]
    ## run lasso to select top predictors
    coefs <- coef(glmnet(X, Y))
    coefs_nonzero <- colSums(coefs!=0)
    coefs <- coefs[,which.min(abs(coefs_nonzero-5))[1]]
    ## create l1 estimate
    l1_est <- coefs[-1]
    l1_est <- l1_est[l1_est!=0]
    ## randomly permute coefficients to get l1 random
    l1r_est <- l1_est
    l1r_est[] <- sample(unname(l1r_est))
    ## determine X subset
    ix <- coefs[-1]!=0
    if(sum(ix)>0){
      Xsub <- X[,ix,drop=FALSE]
    } else {
      ## select random columns if lasso does not return anything
      Xsub <- X[,sample(ncol(X),4)] 
    }
    ## CD coeff
    fit <- try(causalDantzig(Xsub,Y,ExpInd=exp_sub))
    if(class(fit)!="try-error"){
      fn <- names(coef(fit))
      cd_est <- coef(fit)[,1]
      names(cd_est) <- fn
    } else {
      cd_est <- NULL
    }
    ## ICP
    fit <- ICP(Xsub,Y,exp_sub,
               maxNoVariables=100,
               maxNoVariablesSimult=100,
               stopIfEmpty=TRUE,
               showAcceptedSets=FALSE,
               showCompletion=FALSE)
    icp_est <- fit$maximinCoefficients
    ## instrumental variables
    Xsubi <- cbind(1,Xsub)
    iv_est <- ivreg.fit(Xsubi,Y,E)$coefficients[-1]
    ## regularized IV2
    iv2_est <- TSLSreg(Xsub,Y,E)
    Xsub <- NULL
    X <- NULL
    Y <- NULL
    E <- NULL
    gc()
    tout <- list(l1_est=l1_est,
                 l1r_est=l1r_est,
                 cd_est=cd_est,
                 icp_est=icp_est,
                 iv_est=iv_est,
                 iv2_est=iv2_est)
    print(paste0("================ yname:",yname))
    print("tout is:")
    print(tout)
    return(tout)
  }
  N <- ncol(ko_train)
  ## each iteration of mclapply has its own random starting seed
  ## this makes computations exactly reproducible
  rseeds <- sample(1:10^8,N)
  out <- mclapply(1:N,ComputeFit,rseeds=rseeds,mc.cores=mc.cores)
  names(out) <- colnames(ko_train)[1:N]
  return(out)
}