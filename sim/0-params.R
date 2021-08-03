## stores / generates parameters for simulations
source("funcs.R")

## check structure of parameters
CheckParams <- function(params){
  params_names <- c("pp","p0","nobs","nko","np",
                    "ko_strength","rho1","rho0",
                    "st1","st2","S2","ntop")
  if(mean(params_names %in% names(params))!=1){
    print(names(params))
    print(params_names)
    stop("params does not contain all necessary parameters")
  }
  if(mean(names(params) %in% params_names)!=1){
    print(names(params))
    print(params_names)
    stop("params contains unused parameters")
  }
  if(params$np > params$pp - 2*params$p0){
    stop("np must be less than or equal to pp - 2p0")
  }
  if(params$pp %% 2 == 1){
    stop("pp must be even")
  }
  if(params$pp <= 2*params$p0){
    stop("pp must be greater than 2*p0")
  }
  return(params)  
}

## updates parameters and checks validity
## regenerates S2 if necessary 
## ie. whenever (pp,p0,rho1,rho0) changes
UpdateParams <- function(params,update=NULL){
  containsS2 <- "S2" %in% names(params)
  ## determine if need to regenerate S2, only if
  ## updates to select parameters or S2 no currently computed
  if(class(update)=="list"){
    ix <- names(update)[names(update) %in% c("pp","p0","rho1","rho0")]
    regen <- !identical(update[ix],params[ix]) | !containsS2
  } else {
    regen <- !containsS2
  }
  if("S2" %in% names(update)){
    stop("S2 may not be updated. must be derived from other parameters")
  }
  if(!is.null(update)){
    ## update each parameter value
    for(ii in 1:length(update)){
      params[[names(update)[ii]]] <- update[[ii]]
    }
  }
  ## recompute S2 if necessary or originally missing
  if(regen){
    params$S2 <- ComputeS2(params$pp,params$p0,params$rho1,params$rho0)
  }
  print(regen)
  return(CheckParams(params))
}



####
#### parameters for quick demo
####
params1 <- list(pp=10, ## size of ancestor / descendent sets
                p0=1, ## size of parents / children
                nobs=2000, ## observations per condition
                nko=1, ## number of ko per gene knocked out (1 for kemmeren)
                np=8, ## number of genes knocked out
                ko_strength=-40, ## knockout shift strength
                rho1=0, # covariance for ancestors / descendants
                rho0=0, # covariance direct relatives
                st1=1, ## st1/st2 control relative correlation of y with parents/children
                st2=1,
                ntop=8) ## expected number of connections for node 1 (higher=more dense network)
## make parameter values
params1 <- UpdateParams(params1)
params2 <- UpdateParams(params1,list(st1=1,st2=10))
params3 <- UpdateParams(params1,list(pp=100,p0=1,
                                     nobs=500,nko=1,np=98)) ## higher dimensional
params_quick <- list(params1=params1,
                     params2=params2,
                     params3=params3)


####
#### full simulation parameters
####
pp <- 6400
ntop <- c(20,40,80,160,320)
st <- c("strong","weak")
p0 <- c(1,2)
sim_params_des <- data.frame(ntop=rep(ntop,each=length(st)),
                             st=rep(st,length(pp)))
sim_params_des <- cbind(sim_params_des[rep(1:nrow(sim_params_des),length(p0)),],
                        p0=rep(p0,each=nrow(sim_params_des)))

st1_code <- c(1,1)
names(st1_code) <- c("strong","weak")
st2_code <- c(1,1000)
names(st2_code) <- c("strong","weak")
sim_params_des$st1 <- st1_code[sim_params_des$st]
sim_params_des$st2 <- st2_code[sim_params_des$st]


params_hpc <- vector("list",nrow(sim_params_des))
for(ii in 1:nrow(sim_params_des)){
  print(paste0("=== generating params ",ii," / ",nrow(sim_params_des)))
  print(sim_params_des[ii,])
  params_hpc[[ii]] <- UpdateParams(params1,
                                   list(pp=pp,
                                        p0=sim_params_des$p0[ii],
                                        np=pp/4,
                                        nobs=300,
                                        nko=1,
                                        rho1=0,
                                        rho0=0,
                                        st1=sim_params_des$st1[ii],
                                        st2=sim_params_des$st2[ii],
                                        ntop=sim_params_des$ntop[ii]))
  params1 <- params_hpc[[ii]] ## so do not have to regenerate S2 for every update
}

save(params_quick,params_hpc,sim_params_des,file="0-params.RData")
