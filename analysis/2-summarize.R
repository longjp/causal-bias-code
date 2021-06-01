## organizes output of pipeline. outputs
##          res  :  list of matrices with rankings of each cause effect pair
##     top_pred  :  list of data frames with top rankings and reversed rankings
##                  useful for assessing bias in methods
## top_pred_sym  :  list of data frames with top predictions for pairs (i,j)
##                  where both i and j are knocked out

rm(list=ls())
library(Matrix)

load("0-data-setup.RData")

## get results files
f <- list.files()
f <- f[grepl("1-fit",f) & grepl(".RData",f)]

recode <- 1:ncol(ko)
names(recode) <- colnames(ko)

## create structure for res
load(f[1]) ## loads pred, results from model fits
res <- vector("list",length(pred[[1]][[1]]))
names(res) <- names(pred[[1]][[1]])
for(ii in 1:length(res)){
  res[[ii]] <- vector("list",length(f))
}

## load predictions into sparse matrices
for(jj in 1:length(f)){
  print(paste0("run: ",jj))
  load(f[jj]) ## loads pred, results from model fit
  for(method_ix in 1:length(res)){
    print(paste0("  method: ",method_ix))
    aa <- vector("list",length(pred))
    for(cv_ix in 1:length(aa)){
      print(paste0("    cv_ix: ",cv_ix))
      vals <- unlist(lapply(pred[[cv_ix]],function(x){x[[method_ix]]}))
      vals_names <- strsplit(names(vals),".",fixed=TRUE)
      koeffect <- vapply(vals_names,function(x){x[1]},c(""))
      kocause <- vapply(vals_names,function(x){x[2]},c(""))
      kocause_ix <- recode[kocause]
      koeffect_ix <- recode[koeffect]
      aa[[cv_ix]] <- sparseMatrix(kocause_ix,koeffect_ix,x=vals,
                                  dims=c(ncol(ko),ncol(ko)),
                                  dimnames=list(colnames(ko),colnames(ko)))
      ## do not make predictions on training data
      aa[[cv_ix]][rownames(ko)[ix!=cv_ix],] <- 0
    }
    for(ii in 2:length(aa)){
      ## for genes which are not knocked out only make predictions
      ## using first fold 
      ix_rows <- !(rownames(aa[[ii]]) %in% rownames(ko))
      aa[[ii]][ix_rows,] <- 0
      aa[[1]] <- aa[[1]] + aa[[ii]]
    }
    aa <- aa[[1]]
    res[[method_ix]][[jj]] <- aa
  }
}

## visual check
head(koeffect_ix,20)
head(kocause_ix,20)
head(vals,20)
res[["cd_est"]][[1]][1:8,1:8]
res[["cd_est"]][[2]][1:8,1:8]


## distributions of estimates
names(res)
for(ii in 1:length(res)){
  print(names(res)[ii])
  print(summary(as.vector(res[[ii]][[1]])))
}

## threshold estimates
# thres <- 2.5
# for(ii in 1:length(res)){
#   for(jj in 1:length(res[[ii]])){
#     res[[ii]][[jj]][res[[ii]][[jj]] > thres] <- thres
#     res[[ii]][[jj]][res[[ii]][[jj]] < -thres] <- -thres
#   }
# }


## rank then average ranks across simulation runs
for(ii in 1:length(res)){
  print(paste0("averaging performance for: ",ii," / ",length(res)))
  res[[ii]][[1]] <- apply(abs(res[[ii]][[1]]),2,rank)
  for(jj in 2:length(res[[ii]])){
    res[[ii]][[1]] <- res[[ii]][[1]] + apply(abs(res[[ii]][[jj]]),2,rank)
  }
  res[[ii]] <- -res[[ii]][[1]]/length(res[[ii]]) ## negative of average rank, so -large is good
}

## find top effects for each method and rank of flipped effect
## e.g. if X->Y is top prediction followed-up on, what is rank of Y->X
ko_sie_full <- matrix(NA,ncol=ncol(res[[1]]),nrow=nrow(res[[1]]),
                      dimnames=list(rownames(res[[1]]),colnames(res[[1]])))
ko_sie_full[rownames(ko_sie),] <- ko_sie
top_pred <- list()
top_pred_sym <- list()

for(ii in 1:length(res)){
  print(paste0("finding top effects for: ",ii," / ",length(res)))
  res_df <- data.frame(cause=rep(rownames(res[[ii]]),ncol(res[[ii]])),
                       effect=rep(colnames(res[[ii]]),each=nrow(res[[ii]])),
                       pred=as.vector(res[[ii]]),
                       res=as.vector(ko_sie_full))
  res_df <- res_df[order(res_df$pred),]
  res_df$follow_up <- res_df$cause %in% rownames(ko)
  ## compute value of effect->cause prediction and result (TRUE,FALSE,NA)
  ce_names <- paste0(res_df$cause,res_df$effect)
  vals <- res_df$pred
  names(vals) <- ce_names
  res_df$pred_flip <- vals[paste0(res_df$effect,res_df$cause)]
  vals <- res_df$res
  names(vals) <- ce_names
  res_df$res_flip <- vals[paste0(res_df$effect,res_df$cause)]
  ## compute rank of original prediction and flipped prediction
  res_df$rank <- rank(res_df$pred)
  res_df$rank_flip <- rank(res_df$pred_flip)
  head(res_df)
  
  ## output followed up results
  res_df_sub <- res_df[res_df$follow_up,c("cause","effect","res","rank",
                                          "res_flip","rank_flip")]
  ## save top predictions
  top_pred[[ii]] <- res_df_sub[1:5000,]
  ix <- res_df$follow_up & !is.na(res_df$res_flip)
  res_df_sub <- res_df[ix,
                       c("cause","effect","res","rank",
                                          "res_flip","rank_flip")]
  top_pred_sym[[ii]] <- res_df_sub[1:5000,]
}
names(top_pred) <- names(res)
names(top_pred_sym) <- names(res)
save(res,top_pred,top_pred_sym,file="2-summarize.RData")
