## ----setup, include=TRUE-------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(hmlasso)
library(MASS)

## ----data----------------------------------------------------------------
# setting
n <- 2000 # number of samples
p <- 20 # number of dimensions
r <- 0.5 # correlation levels
eps <- 1 # noise level

# construct covariance matrix
mu <- rep(0, length=p)
Sigma <- matrix(r, nrow=p, ncol=p) # correlation among variables is 0
diag(Sigma) <- 1 # variance of each variable is 1

# generate training data
set.seed(0)
X <- mvrnorm(n, mu, Sigma)
b <- matrix(0, nrow = p, ncol = 1)
b[1:10,] <- (-1)^(1:10) * (1:10)
y <- X %*% b + eps * rnorm(n)
colnames(X) <- paste0("V", 1:p)
rownames(b) <- paste0("V", 1:p)

# introduce missing data
X_incompl <- X
mrc <- runif(p, min=0.1, max=0.9) # missing rate of each column
for (j in 1:p) {
  mr <- sample(1:nrow(X), round(nrow(X) * mrc[j]))
  X_incompl[mr, j] <- NA
}
image(t(apply(apply(X_incompl, c(1,2), function(v){ifelse(is.na(v), 1, 0)}), 2, rev)),
      col=grey(11:0/12), main="missing pattern") # visualize (black:na, white:fill)

# generate test data
n_ts <- 10000
X_ts <- mvrnorm(n_ts, mu, Sigma)
y_ts <- X_ts %*% b + eps * rnorm(n_ts)
colnames(X_ts) <- paste0("V", 1:p)

# setup cv fold
nfolds <- 5
foldid1 <- sample(rep(1:nfolds, (n %/% nfolds)), replace=FALSE)
foldid2 <- sample(1:nfolds, (n %% nfolds), replace=FALSE)
foldid <- c(foldid1, foldid2)

## ----meanimp-------------------------------------------------------------
# MeanImp
cv_fit <- cv.hmlasso(X_incompl, y, nlambda=50, lambda.min.ratio=1e-1,
                     foldid=foldid, impute_method="mean",
                     direct_prediction=FALSE, positify="mean")
plot(cv_fit) # plot lambda-MSE
plot(cv_fit$fit) # plot solution path
cv_fit$cve[cv_fit$lambda.min.index] # minimum cross-valiation error
cv_fit$fit$beta[, cv_fit$lambda.min.index] # estimated coefficient
pred <- predict(cv_fit$fit, X_ts) # predict
sqrt(sum((cv_fit$fit$beta[, cv_fit$lambda.min.index] - as.vector(b))^2)) # estimation error
sqrt(sum((pred[, cv_fit$lambda.min.index] - y_ts)^2) / n_ts) # prediction error

## ----cocolasso-----------------------------------------------------------
# CoCoLasso
cv_fit <- cv.hmlasso(X_incompl, y, nlambda=50, lambda.min.ratio=1e-1, 
                        foldid=foldid, direct_prediction=TRUE, 
                        positify="admm_max", weight_power = 0)
plot(cv_fit) # plot lambda-MSE
plot(cv_fit$fit) # plot solution path
cv_fit$cve[cv_fit$lambda.min.index] # minimum cross-valiation error
cv_fit$fit$beta[, cv_fit$lambda.min.index] # estimated coefficient
pred <- predict(cv_fit$fit, X_ts) # predict
sqrt(sum((cv_fit$fit$beta[, cv_fit$lambda.min.index] - as.vector(b))^2)) # estimation error
sqrt(sum((pred[, cv_fit$lambda.min.index] - y_ts)^2) / n_ts) # prediction error

## ----hmlasso-------------------------------------------------------------
# HMLasso
cv_fit <- cv.hmlasso(X_incompl, y, nlambda=50, lambda.min.ratio=1e-1, 
                        foldid=foldid, direct_prediction=TRUE, 
                        positify="admm_frob", weight_power = 1)
plot(cv_fit) # plot lambda-MSE
plot(cv_fit$fit) # plot solution path
cv_fit$cve[cv_fit$lambda.min.index] # minimum cross-valiation error
cv_fit$fit$beta[, cv_fit$lambda.min.index] # estimated coefficient
pred <- predict(cv_fit$fit, X_ts) # predict
sqrt(sum((cv_fit$fit$beta[, cv_fit$lambda.min.index] - as.vector(b))^2)) # estimation error
sqrt(sum((pred[, cv_fit$lambda.min.index] - y_ts)^2) / n_ts) # prediction error

