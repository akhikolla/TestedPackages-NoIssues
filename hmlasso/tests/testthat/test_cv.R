library(MASS)

context("cross validation")

test_that("complete case", {
  X <- matrix(c(14, 25, 4, 7, 1, 2, 23, 11, 30, 18, 19, 26, 10, 22, 21, 27, 9, 5, 13, 17, 24, 12, 16, 29, 15, 20, 28, 3, 8, 6),
              nrow=10, ncol=3)
  b <- matrix(c(-1, 1, 0), nrow=3, ncol=1)
  e <- matrix(0.01 * (-1)^(1:10), nrow=10, ncol=1)
  y <- X %*% b + e
  cv_fit <- cv.hmlasso(X, as.numeric(y))
  fit <- hmlasso(X, as.numeric(y))
  plot(cv_fit)
  plot(cv_fit$fit)

  expect_length(cv_fit$cve, length(cv_fit$lambda))
  expect_equal(cv_fit$lambda, cv_fit$fit$lambda, tolerance=1e-1)
  expect_equal(cv_fit$fit$beta[, cv_fit$lambda.min.index], as.numeric(b), tolerance=1e-1)
  expect_equal(cv_fit$fit, fit, tolerance=1e-2)
})

test_that("missing case", {
  X <- matrix(c(14, 25, 4, 7, 1, 2, 23, 11, 30, 18, 19, 26, 10, 22, 21, 27, 9, 5, 13, 17, 24, 12, 16, 29, 15, 20, 28, 3, 8, 6),
              nrow=10, ncol=3)
  b <- matrix(c(-1, 1, 0), nrow=3, ncol=1)
  e <- matrix(0.01 * (-1)^(1:10), nrow=10, ncol=1)
  y <- X %*% b + e
  X_miss <- X
  X_miss[1, 1] <- NA
  X_miss[2, 2] <- NA
  cv_fit <- cv.hmlasso(X_miss, as.numeric(y), nfolds = 5)
  fit <- hmlasso(X_miss, as.numeric(y))
  plot(cv_fit)
  plot(cv_fit$fit)

  expect_length(cv_fit$cve, length(cv_fit$lambda))
  expect_length(unique(cv_fit$foldid), 5)
  expect_equal(cv_fit$lambda, cv_fit$fit$lambda, tolerance=1e-1)
  expect_equal(cv_fit$fit$beta[, cv_fit$lambda.min.index], as.numeric(b), tolerance=5e-1)
  expect_equal(cv_fit$fit, fit, tolerance=1e-2)
})

# test_that("large sample", {
#   # setting
#   n <- 100 # number of samples
#   p <- 10 # number of dimensions
#   r <- 0.5 # correlation levels
#   eps <- 1 # noise level
#   mu <- rep(0, length=p)
#   Sigma <- matrix(r, nrow=p, ncol=p) # correlation among variables is 0
#   diag(Sigma) <- 1 # variance of each variable is 1
#   set.seed(0)
#   X <- mvrnorm(n, mu, Sigma)
#   b <- matrix(0, nrow = p, ncol = 1)
#   # b[1:10,] <- (-1)^(1:10) * (1:10)
#   b[1:4, ] <- c(1, -1, 2, -2)
#   y <- X %*% b + eps * rnorm(n)
#   X_miss <- X
#   X_miss[sample(1:length(X), round(length(X) * 0.1))] <- NA # for random missing pattern
#   cv_fit <- cv.hmlasso(X_miss, as.numeric(y))
#   fit <- hmlasso(X_miss, as.numeric(y))
#   plot(cv_fit)
#   plot(cv_fit$fit)
#
#   expect_length(cv_fit$cve, length(cv_fit$lambda))
#   expect_length(unique(cv_fit$foldid), 10)
#   expect_equal(cv_fit$lambda, cv_fit$fit$lambda, tolerance=1e-1)
#   expect_equal(cv_fit$fit$beta[, cv_fit$lambda.min.index], as.numeric(b), tolerance=5e-1)
#   expect_equal(cv_fit$fit, fit, tolerance=1e-2)
# })
#
# test_that("large sample cocolasso verbose", {
#   # setting
#   n <- 100 # number of samples
#   p <- 10 # number of dimensions
#   r <- 0.5 # correlation levels
#   eps <- 1 # noise level
#   mu <- rep(0, length=p)
#   Sigma <- matrix(r, nrow=p, ncol=p) # correlation among variables is 0
#   diag(Sigma) <- 1 # variance of each variable is 1
#   set.seed(0)
#   X <- mvrnorm(n, mu, Sigma)
#   b <- matrix(0, nrow = p, ncol = 1)
#   # b[1:10,] <- (-1)^(1:10) * (1:10)
#   b[1:4, ] <- c(1, -1, 2, -2)
#   y <- X %*% b + eps * rnorm(n)
#   X_miss <- X
#   X_miss[sample(1:length(X), round(length(X) * 0.5))] <- NA # for random missing pattern
#   cv_fit <- cv.hmlasso(X_miss, as.numeric(y), positify="admm_max", weight_power = 0, verbose=TRUE)
#   fit <- hmlasso(X_miss, as.numeric(y), positify="admm_max", weight_power = 0, verbose=TRUE)
#   plot(cv_fit, ylim=c(0,100))
#   plot(cv_fit$fit, ylim=c(-10, 10))
#
#   expect_length(cv_fit$cve, length(cv_fit$lambda))
#   expect_length(unique(cv_fit$foldid), 10)
#   expect_equal(cv_fit$lambda, cv_fit$fit$lambda, tolerance=1e-1)
#   expect_equal(cv_fit$fit$beta[, cv_fit$lambda.min.index], as.numeric(b), tolerance=2e0)
#   expect_equal(cv_fit$fit, fit, tolerance=1e-2)
# })
#
# test_that("large sample high missing alpha=0.5 frobenius formulation", {
#   # setting
#   n <- 100 # number of samples
#   p <- 10 # number of dimensions
#   r <- 0.5 # correlation levels
#   eps <- 1 # noise level
#   mu <- rep(0, length=p)
#   Sigma <- matrix(r, nrow=p, ncol=p) # correlation among variables is 0
#   diag(Sigma) <- 1 # variance of each variable is 1
#   set.seed(0)
#   X <- mvrnorm(n, mu, Sigma)
#   b <- matrix(0, nrow = p, ncol = 1)
#   # b[1:10,] <- (-1)^(1:10) * (1:10)
#   b[1:4, ] <- c(1, -1, 2, -2)
#   y <- X %*% b + eps * rnorm(n)
#   X_miss <- X
#   X_miss[sample(1:length(X), round(length(X) * 0.7))] <- NA # for random missing pattern
#   cv_fit <- cv.hmlasso(X_miss, as.numeric(y), positify="admm_frob", weight_power = 0.5, verbose=TRUE)
#   fit <- hmlasso(X_miss, as.numeric(y), positify="admm_frob", weight_power = 0.5, verbose=TRUE)
#   plot(cv_fit, ylim=c(0,100))
#   plot(cv_fit$fit, ylim=c(-10, 10))
#
#   expect_length(cv_fit$cve, length(cv_fit$lambda))
#   expect_length(unique(cv_fit$foldid), 10)
#   expect_equal(cv_fit$lambda, cv_fit$fit$lambda, tolerance=1e-1)
#   expect_equal(cv_fit$fit$beta[, cv_fit$lambda.min.index], as.numeric(b), tolerance=2e0)
#   expect_equal(cv_fit$fit, fit, tolerance=1e-2)
# })

# test_that("comparison with ncvreg", {
#   # setting
#   n <- 100 # number of samples
#   p <- 10 # number of dimensions
#   r <- 0.5 # correlation levels
#   eps <- 1 # noise level
#   mu <- rep(0, length=p)
#   Sigma <- matrix(r, nrow=p, ncol=p) # correlation among variables is 0
#   diag(Sigma) <- 1 # variance of each variable is 1
#   set.seed(0)
#   X <- mvrnorm(n, mu, Sigma)
#   b <- matrix(0, nrow = p, ncol = 1)
#   # b[1:10,] <- (-1)^(1:10) * (1:10)
#   b[1:4, ] <- c(1, -1, 2, -2)
#   y <- X %*% b + eps * rnorm(n)
#   X_miss <- X
#   X_miss[sample(1:length(X), round(length(X) * 0.1))] <- NA # for random missing pattern
#   set.seed(0)
#   foldid <- sample(1:10, n, replace=TRUE)
#
#   cv_fit <- cv.hmlasso(X, as.numeric(y), nlambda = 100, foldid = foldid, lambda.min.ratio = 1e-4)
#
#   cv_fit_ncvreg <- cv.ncvreg(X, as.numeric(y), nlambda = 100, fold = foldid, lambda.min = 1e-4, penalty = "lasso")
#
#   expect_equal(cv_fit$lambda, cv_fit_ncvreg$lambda, tolerance=1e-4)
#   expect_equal(cv_fit$cve, cv_fit_ncvreg$cve, tolerance=1e-4)
#   expect_equivalent(cv_fit$fit$beta[, cv_fit$lambda.min.index],
#                     cv_fit_ncvreg$fit$beta[, cv_fit_ncvreg$lambda==cv_fit_ncvreg$lambda.min][-1], tolerance=1e-2)
# })
