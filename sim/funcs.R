## functions for generating simulation data
library(InvariantCausalPrediction)
library(magic)

GenGraph <- function(p,probs){
  if(length(probs)!=3){
    stop("probs should be length 3")
  }
  vals <- sample(c(-1,1,0),
                 size=p^2,
                 prob=probs,
                 replace=TRUE)
  A <- matrix(vals,nrow=p)
  A[lower.tri(A,diag=TRUE)] <- 0
  return(A)
}

CreateA <- function(pp,p0,st1=1,st2=1,ntop=20){
  ## top node connects with about ntop
  ## nodes, scaling prob1 with p1,p0 ensures
  ## invertibility of A
  prob1 <- min(ntop/(2*pp),1/2)
  prob0 <- 1-2*prob1
  A <- GenGraph(pp+1,c(prob1,prob1,prob0))
  ## eliminate all causes of Y
  A[,pp/2+1] <- 0
  ## eliminate all effects of Y
  A[pp/2+1,] <- 0
  ## make p0/2 causes of Y
  ix <- (pp/2-p0+1):(pp/2)
  A[ix,pp/2+1] <- st1*(2*rbinom(length(ix),prob=0.5,size=1)-1)
  ## make p0/2 effects of Y
  ix <- (pp/2+2):(pp/2+p0+1)
  A[pp/2+1,ix] <- st2*(2*rbinom(length(ix),prob=0.5,size=1)-1)
  ## scale A so variance is stable
  A <- t(t(A) / (sqrt(colSums(A*A) + 1)))
  ## extract causal coefficients on Y from A and remove Y on Y cause
  beta <- A[,pp/2+1]
  beta <- beta[-(pp/2+1)]
  return(list(A=A,beta=beta))
}

GenerateER <- function(pp,p0,nobs,nko,ko_strength,np){
  ko <- diag(np)
  ko <- ko[rep(1:nrow(ko),each=nko),]
  E <- rbind(matrix(0,nrow=nobs,ncol=np),ko)
  R <- ko_strength*diag(pp-2*p0)
  Rm <- matrix(0,nrow=nrow(R),ncol=2*p0+1)
  R <- cbind(R[,1:(ncol(R)/2)],Rm,R[,(ncol(R)/2+1):ncol(R)])
  R <- R[sample(1:nrow(R),np),]
  return(list(E=E,R=R))
}

CreateCov <- function(p,rho){
  Sigma <- matrix(1:p,nrow=p,ncol=p)
  Sigma <- abs(t(t(Sigma) - (1:p)))
  Sigma <- rho^Sigma
  return(Sigma)
}

## computes square root of covariance matrix
## useful to precompute rather than repeating
## each time simulation is run
ComputeS2 <- function(pp,p0,rho1,rho0){
  Sigma1 <- CreateCov(pp/2-p0,rho1)
  Sigma0 <- CreateCov(2*p0+1,rho0)
  s <- svd(Sigma1)
  R1 <- t(s$v %*% (t(s$u) * sqrt(pmax(s$d, 0))))
  s <- svd(Sigma0)
  R0 <- t(s$v %*% (t(s$u) * sqrt(pmax(s$d, 0))))
  return(adiag(R1,R0,R1))
}


GenerateData <- function(pp,p0,nobs,nko,ko_strength,np,S2,st1,st2,ntop){
  Abeta <- CreateA(pp,p0,st1,st2,ntop)
  A <- Abeta$A
  beta <- Abeta$beta
  ER <- GenerateER(pp,p0,nobs,nko,ko_strength,np)
  E <- ER$E
  R <-ER$R
  deltaX <- matrix(rnorm(nrow(E)*(pp+1)),nrow=nrow(E))%*%S2
  IA <- diag(nrow(A))-A
  X <- t(solve(t(IA),t(E%*%R + deltaX)))
  return(list(X=X,Abeta=Abeta,ER=ER))
}

RunSim <- function(jj,params,rseeds){
  set.seed(rseeds[jj])
  Ntrue <- params$p0
  ## generate data
  dat <- GenerateData(params$pp,
                      params$p0,
                      params$nobs,
                      params$nko,
                      params$ko_strength,
                      params$np,
                      params$S2,
                      params$st1,
                      params$st2,
                      params$ntop)
  X <- dat$X
  beta <- dat$Abeta$beta
  E <- dat$ER$E
  ix <- params$pp/2 + 1
  Y <- X[,ix]
  X <- X[,-ix]
  ## fit lasso to regularize
  coefs <- coef(glmnet(X, Y))
  coefs_nonzero <- colSums(coefs!=0)
  coefs <- coefs[,which.min(abs(coefs_nonzero-5))[1]]
  #fitgm <- cv.glmnet(X,Y)
  #coefs <- coef(fitgm)
  ## evaluate lasso
  l1_est <- coefs[-1]
  L1 <- sum(beta[order(abs(l1_est),
                       decreasing=TRUE)[1:Ntrue]]!=0)
  ix <- coefs[-1]!=0
  nused <- sum(ix)
  ## determine X subset
  ix <- coefs[-1]!=0
  if(sum(ix)>0){
    Xsub <- X[,ix,drop=FALSE]
  } else {
    ## select random columns if lasso does not return anything
    Xsub <- X[,sample(ncol(X),4)] 
  }
  ## create random order of L1 chosen coeffs
  ## used to break ties for methods / random guess when
  ## methods fail
  rand_coef <- rep(0,ncol(Xsub))
  rand_coef[ix] <- sample(which(ix))
  ## L1 random - randomly select among L1 selected variables
  L1R <- sum(beta[order(abs(rand_coef),
                        decreasing=TRUE)[1:Ntrue]]!=0)
  ## fit CD
  ExpInd <- 1*(rowSums(abs(E))!=0)
  cd_est <- rep(0,ncol(X))
  fit <- try(causalDantzig(Xsub,Y,ExpInd=ExpInd),silent=TRUE)
  if(class(fit)!="try-error"){
    cd_est[ix] <- coef(fit)[,1]
  } 
  CD <- sum(beta[order(abs(cd_est),rand_coef,
                       decreasing=TRUE)[1:Ntrue]]!=0)
  ## ICP
  fit <- ICP(Xsub,Y,ExpInd,
             maxNoVariables=10,
             maxNoVariablesSimult=10,
             stopIfEmpty=TRUE,
             showAcceptedSets=FALSE,
             showCompletion=FALSE)
  ICPest <- fit$maximinCoefficients
  IC <- rep(0,ncol(X))
  if(!is.null(ICPest)){
    IC[ix] <- ICPest
  }
  IC <- sum(beta[order(abs(IC),rand_coef,
                       decreasing=TRUE)[1:Ntrue]]!=0)
  dat <- NULL
  gc()
  return(list(c(L1=L1,L1R=L1R,CD=CD,ICP=IC),nused))
}