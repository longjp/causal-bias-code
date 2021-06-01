rm(list=ls())
library(tidyr)
library(kableExtra)
load("0-params.RData")

## find .RData files output by hpc
fs <- list.files()
fs <- fs[grepl(".RData",fs) & grepl("1-",fs,fixed=TRUE)]
ix <- as.numeric(gsub(".RData","",gsub("1-","",fs),fixed=TRUE))
fs <- fs[order(ix)]

## create output matrix
load(fs[1])
cnames <- colnames(res)
out <- matrix(NA_real_,
              nrow=length(fs),
              ncol=length(cnames))
colnames(out) <- cnames

for(ii in 1:length(fs)){
  load(fs[ii])
  out[ii,] <- colMeans(res)
}

## merge with simulation parameters
res <- cbind(sim_params_des,out)
res$st1 <- NULL
res$st2 <- NULL

res

## prepare to output
res <- pivot_wider(res,
                   names_from="st",
                   values_from=colnames(res)[4:ncol(res)])
res <- cbind(res[,1:2],
             res[,grepl("strong",colnames(res))],
             res[,grepl("weak",colnames(res))])
colnames(res) <- vapply(strsplit(colnames(res),"_"),function(x){x[1]},c("howdy"))


to_use <- c("ntop","p0","L1","L1R","CD","ICP")
res <- res[,colnames(res) %in% to_use]
colnames(res) <- gsub(".1","",colnames(res),fixed=TRUE)


res$p0 <- NULL
names(res)[1] <- "$n_t$"
head_len <- (ncol(res)-1)/2
## good table advice https://haozhu233.github.io/kableExtra/awesome_table_in_pdf.pdf
caption <- "Simulations results. For the $p_0=1$ rows, the best performance is 1 and the worst performance is 0. For the $p_0=2$ rows, the best performance is 2 and the worst performance is 0. \\label{tab:high-dim}"
tab <- kbl(res,booktabs=TRUE,format="latex",escape=FALSE,caption=caption) %>%
  add_header_above(c(" " = 1, "Strong" = head_len, "Weak" = head_len)) %>%
  pack_rows("$p_0=1$",1,nrow(res)/2,escape=FALSE) %>%
  pack_rows("$p_0=2$",nrow(res)/2+1,nrow(res),escape=FALSE)
save_kable(tab,file="../ms-bias/figs/sim_results.tex")

