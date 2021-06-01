rm(list=ls())
library(ggplot2)
library(Matrix)
library(kableExtra)

load("0-data-setup.RData")
load("2-summarize.RData")

methods_to_use <- c("cd_est","icp_est","l1_est","l1r_est")
mnames <- c("CD","ICP","L1","L1R")
names(mnames) <- methods_to_use


## print table with top predictions
## make table of top top predictions
for(ii in 1:length(methods_to_use)){
    tab_out <- top_pred[[methods_to_use[ii]]][1:10,]
    colnames(tab_out) <- gsub("_","-",colnames(tab_out))
     caption <- paste0(mnames[ii]," rankings  \\label{tab:",mnames[ii],"-ranks}")
    fout <- paste0(mnames[ii],"_results.tex")
    tab_out <- kbl(tab_out,row.names=FALSE,booktabs=TRUE,format="latex",
                   escape=FALSE,caption=caption)
    save_kable(tab_out,file=paste0("../ms-bias/figs/",fout))
}
  

## analysis of binary classifier
for(ii in 1:length(top_pred)){
  print("==============")
  print(names(top_pred)[ii])
  print("==============")
  ## number of strong effects predicted
  ntop <- 20
  ## 95% ci 
  print(binom.test(sum(top_pred[[ii]]$res[1:ntop]),ntop,
                   p=mean(ko_sie)))
}


## ANALYSIS1: all top_pred
## ROC based performance assessment
for(ii in 1:length(top_pred)){
  top_pred[[ii]]$true_pos <- cumsum(top_pred[[ii]]$res)
  top_pred[[ii]]$false_pos <- cumsum(1-top_pred[[ii]]$res)
  top_pred[[ii]]$method <- names(top_pred)[ii]
  ## need to add (0,0) to make ROC curve start at (0,0)
  top_pred[[ii]] <- rbind(as.list(rep(0,ncol(top_pred[[ii]]))),top_pred[[ii]])
  top_pred[[ii]]$method[1] <- names(top_pred)[ii]
}

res_df_sub <- do.call(rbind,top_pred)
res_df_sub <- res_df_sub[res_df_sub$method %in% methods_to_use,]
res_df_sub$method <- mnames[res_df_sub$method]



p <- ggplot(res_df_sub,aes(x=false_pos,y=true_pos,color=method)) +
  geom_path() + labs(x="False Positive",y="True Positive") +
  theme(legend.position = c(0.1,0.8)) +
  geom_abline(intercept=0,slope=mean(ko_sie)) +
  labs(colour="Estimator")
p <- p + coord_cartesian(xlim=c(0,200),ylim=c(0,100))
p
p_zoom <- p + coord_cartesian(xlim=c(0,10),ylim=c(0,15))
p_zoom


pdf("../ms-bias/figs/roc.pdf",width=5,height=4)
plot(p)
dev.off()

pdf("../ms-bias/figs/roc_zoom.pdf",width=5,height=4)
plot(p_zoom)
dev.off()



## ANALYSIS2: all top_pred_sym
## ROC based performance assessment
for(ii in 1:length(top_pred_sym)){
  top_pred_sym[[ii]]$true_pos <- cumsum(top_pred_sym[[ii]]$res)
  top_pred_sym[[ii]]$false_pos <- cumsum(1-top_pred_sym[[ii]]$res)
  top_pred_sym[[ii]]$method <- names(top_pred_sym)[ii]
  ## need to add (0,0) to make ROC curve start at (0,0)
  top_pred_sym[[ii]] <- rbind(as.list(rep(0,ncol(top_pred_sym[[ii]]))),top_pred_sym[[ii]])
  top_pred_sym[[ii]]$method[1] <- names(top_pred_sym)[ii]
}

res_df_sub <- do.call(rbind,top_pred_sym)
res_df_sub <- res_df_sub[res_df_sub$method %in% methods_to_use,]
res_df_sub$method <- mnames[res_df_sub$method]


## get random guessing
ko_sie_sub <- ko_sie[,colnames(ko_sie) %in% rownames(ko_sie)]

p <- ggplot(res_df_sub,aes(x=false_pos,y=true_pos,color=method)) +
  geom_path() + labs(x="False Positive",y="True Positive") +
  theme(legend.position = c(0.1,0.8)) +
  geom_abline(intercept=0,slope=mean(ko_sie_sub)) +
  labs(colour="Estimator")
p <- p + coord_cartesian(xlim=c(0,200),ylim=c(0,40))
p
p_zoom <- p + coord_cartesian(xlim=c(0,10),ylim=c(0,5))
p_zoom


pdf("../ms-bias/figs/roc_sym.pdf",width=5,height=4)
plot(p)
dev.off()

pdf("../ms-bias/figs/roc_zoom_sym.pdf",width=5,height=4)
plot(p_zoom)
dev.off()

