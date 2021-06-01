rm(list=ls())
library(glmnet)
?glmnet


# Gaussian
x <- matrix(rnorm(100 * 20), 100, 20)
y <- rnorm(100)
fit1 <- glmnet(x, y)
print(fit1)
a <- coef(fit1)
dim(a)
coef(fit1, s = 0.01)  # extract coefficients at a single value of lambda
predict(fit1, newx = x[1:10, ], s = c(0.01, 0.005))  # make predictions


cvob1 <- cv.glmnet(x, y)
set.seed(1010)
n = 1000
p = 100
nzc = trunc(p/10)
x = matrix(rnorm(n * p), n, p)
beta = rnorm(nzc)
fx = x[, seq(nzc)] %*% beta
eps = rnorm(n) * 5
y = drop(fx + eps)
px = exp(fx)
px = px/(1 + px)
ly = rbinom(n = length(px), prob = px, size = 1)

set.seed(1011)
cvob1 = cv.glmnet(x, y)
plot(cvob1)
coef(cvob1)

a <- c("a"=0,"b"=0)
a
a[a!=0]
class(cvob1)
a <- coef(cvob1)
class(a)
a
a <- as.matrix(a)
class(a)
dim(a)
co <- as.matrix(coef(cvob1))[,1][-1]
