### assess whether models with sie's as
### predictors produce more invariant residual distributions


rm(list=ls())
load("0-data-setup.RData")
library(ggplot2)
library(glmnet)

ls()

rownames(ko) <- make.names(rownames(ko))
colnames(ko) <- make.names(colnames(ko))
rownames(ko_sie) <- make.names(rownames(ko_sie))
colnames(ko_sie) <- make.names(colnames(ko_sie))
rownames(obs) <- make.names(rownames(obs))
colnames(obs) <- make.names(colnames(obs))


class(obs)


oko <- rbind(obs,ko)
oko <- as.data.frame(oko)
oko$cl <- c(rep("obs",nrow(obs)),rep("ko",nrow(ko)))

head(colnames(oko),100)

ii <- 30
sum(ko_sie[,ii])
colnames(ko_sie)[ii]

## ko dist has larger variance
## not invariant ko vs obs
ggplot(oko,aes(x=YAL001C,fill=cl)) +
  geom_density(alpha=0.3)


## find significant effects on gene of interest,
## run regression, compute residuals, 
## compare to dist w/o regression
## result: variance is tighter, but possibly less invariance
sig_eff <- rownames(ko_sie)[ko_sie[,ii]]
f <- as.formula(paste0(paste0(colnames(ko_sie)[ii],"~"),
                    paste0(sig_eff,collapse="+"),collapse=""))
fit <- lm(f,data=oko)
oko$residuals <- fit$residuals
p1 <- ggplot(oko,aes(x=YAL001C,fill=cl)) +
  geom_density(alpha=0.3) +
  theme(legend.position = c(0.8, 0.8))
p2 <- ggplot(oko,aes(x=residuals,fill=cl)) +
  geom_density(alpha=0.3) +
  theme(legend.position = c(0.8, 0.8))
gridExtra::grid.arrange(p1,p2,nrow=1)


## option 2: lasso to select top genes
## result: much, much more invariant

## conclusion: result is disappointing
## because ignoring causal issues results
## in most invariant residual distribution
## but: only expect a small number of interventions
##     to affect target. so by ignoring causal issues
##     should obtain model with is good (invariant) for
##     most knockouts. therefore perhaps not suprising
##
x <- rbind(obs,ko)
y <- x[,ii]
x <- x[,-ii]
glm_fit <- glmnet(x,y)
cvfit <- cv.glmnet(x, y)
plot(cvfit)
preds <- predict(cvfit, newx = x, s = "lambda.min")
oko$residuals_glm <- y - preds

p1 <- ggplot(oko,aes(x=YAL001C,fill=cl)) +
  geom_density(alpha=0.3) +
  theme(legend.position = c(0.8, 0.8))
p2 <- ggplot(oko,aes(x=residuals,fill=cl)) +
  geom_density(alpha=0.3) +
  theme(legend.position = c(0.8, 0.8))
p3 <- ggplot(oko,aes(x=residuals_glm,fill=cl)) +
  geom_density(alpha=0.3) +
  theme(legend.position = c(0.8, 0.8))
gridExtra::grid.arrange(p1,p2,p3,nrow=1)

