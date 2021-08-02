## runs cross validation on models
CVModel <- function(predictX,obs,ko,ix,rseed=rseed,mc.cores=4){
  folds <- unique(ix)
  out_predictX <- list()
  for(ii in 1:length(folds)){
    ko_train <- ko[ix!=ii,]
    ko_test_names <- rownames(ko)[ix==ii]
    out_predictX[[ii]] <- predictX(obs,ko_train,ko_test_names,rseed,mc.cores)
  }
  return(out_predictX)
}