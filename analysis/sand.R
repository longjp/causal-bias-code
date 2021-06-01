rm(list=ls())
load("1b-predictINS.RData")
load("0-data-setup.RData")

out <- predINS$out

N <- sum(vapply(out,function(x){dim(x)[1]},c(0)))
pq <- dim(out[[1]])[2:3]
?array
out_merge <- array(NA_integer_,dim=c(N,pq))
inx <- 1
for(ii in 1:length(out)){
  out_merge[inx:(inx+dim(out[[ii]])[1]-1),,] <- out[[ii]]
  inx <- inx + dim(out[[ii]])[1]
}

## match order
dimnames(out_merge)[[1]] <- unlist(lapply(out,function(x){dimnames(x)[[1]]}))
out_merge <- out_merge[rownames(ko_sie),,]
out <- out_merge

## slightly tuned thresholds
pa_frac <- out[,,2] / out[,,1]
ko_pred <- pa_frac<0.25 & out[,,1] > 20

sum(ko_pred)
pred_score <- 1*ko_sie[ko_pred]
binom.test(sum(pred_score),length(pred_score),
                 p=mean(ko_sie))

