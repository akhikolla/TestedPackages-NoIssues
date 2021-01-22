context("iilasso")

test_that("simple iilasso", {
  X <- matrix(c(1,2,3,5,4,7,6,8,9,10),
              nrow=5, ncol=2)
  b <- matrix(c(-1,1), nrow=2, ncol=1)
  e <- matrix(c(0,-0.1,0.1,-0.1,0.1), nrow=5, ncol=1)
  y <- as.numeric(X %*% b + e)
  fit <- lasso(X, y)
  expect_that(fit$beta[,1], equals(c(0,0)))
  expect_that(fit$a0[,1], equals(mean(y)))

  pr <- predict_lasso(fit, X)
  expect_that(pr[,1], equals(rep(mean(y), length=5)))
  expect_that(1 - sum((y - pr[,100])^2) / var(y), equals(1, tolerance=1e-1))
  expect_that(sum((b - fit$beta[,100])^2), equals(0, tolerance=1e-1))
  
  plot_lasso(fit)
})

test_that("iilasso with cross validation", {
  X <- matrix(c(1,2,3,5,4,7,6,8,9,10),
              nrow=5, ncol=2)
  b <- matrix(c(-1,1), nrow=2, ncol=1)
  e <- matrix(c(0,-0.1,0.1,-0.1,0.1), nrow=5, ncol=1)
  y <- as.numeric(X %*% b + e)
  cv_fit <- cv_lasso(X, y, nfolds=5)
  fit <- cv_fit$fit
  expect_that(fit$beta[,1], equals(c(0,0)))
  expect_that(fit$a0[,1], equals(mean(y)))
  expect_that(cv_fit$lambda.min.index, equals(100))

  pr <- predict_lasso(fit, X, cv_fit$lambda.min)
  expect_that(1 - sum((y - pr)^2) / var(y), equals(1, tolerance=1e-1))
  expect_that(sum((b - fit$beta[,100])^2), equals(0, tolerance=1e-1))
  
  plot_cv_lasso(cv_fit)
})

test_that("simple logistic iilasso", {
  X <- matrix(c(1,2,3,5,4,7,6,8,9,10),
              nrow=5, ncol=2)
  b <- matrix(c(-1,1), nrow=2, ncol=1)
  e <- matrix(c(0,-0.1,0.1,-0.1,0.1), nrow=5, ncol=1)
  y <- as.numeric(X %*% b + e)
  y <- ifelse(y>mean(y), 1, 0)
  fit <- lasso(X, y, family="binomial")

  pr <- predict_lasso(fit, X)
  expect_that(pr[,1], equals(rep(mean(y), length=5)))
  expect_that(1 - sum((y - pr[,100])^2) / var(y), equals(1, tolerance=1e-1))
  # expect_that(sum((b - fit$beta[,100])^2), equals(0, tolerance=1e-2))

  plot_lasso(fit)
})

test_that("logistic iilasso with cross validation", {
  X <- matrix(c(1,2,3,5,4,7,6,8,9,10),
              nrow=5, ncol=2)
  b <- matrix(c(-1,1), nrow=2, ncol=1)
  e <- matrix(c(0,-0.1,0.1,-0.1,0.1), nrow=5, ncol=1)
  y <- as.numeric(X %*% b + e)
  y <- ifelse(y>mean(y), 1, 0)
  cv_fit <- cv_lasso(X, y, nfolds=5, family="binomial")
  fit <- cv_fit$fit
  expect_that(fit$beta[,1], equals(c(0,0)))

  pr <- predict_lasso(fit, X, cv_fit$lambda.min)
  pr_cl <- predict_lasso(fit, X, cv_fit$lambda.min, type="class")
  expect_that(1 - sum((y - pr)^2) / var(y), equals(1, tolerance=1e-1))
  expect_that(all(y==pr_cl), equals(TRUE))

  plot_cv_lasso(cv_fit)
})
