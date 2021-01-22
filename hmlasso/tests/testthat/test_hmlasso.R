library(MASS)

context("hmlasso")

test_that("complete, small sample", {
  X <- matrix(c(14, 25, 4, 7, 1, 2, 23, 11, 30, 18, 19, 26, 10, 22, 21, 27, 9, 5, 13, 17, 24, 12, 16, 29, 15, 20, 28, 3, 8, 6),
              nrow=10, ncol=3)
  b <- matrix(c(-1, 1, 0), nrow=3, ncol=1)
  e <- matrix(0.01 * (-1)^(1:10), nrow=10, ncol=1)
  y <- X %*% b + e
  fit <- hmlasso(X, as.numeric(y))
  plot(fit)

  # expect_equal(fit$beta[, 1], rep(0, length=ncol(X)), tolerance=1e-1)
  # expect_equal(fit$a0[1], mean(y), tolerance=1e0)
  expect_equal(fit$beta[, length(fit$lambda)], as.numeric(b), tolerance=1e-1)
  expect_equal(fit$a0[length(fit$lambda)], 0, tolerance=1e-1)

  pr <- predict(fit, X)
  # expect_equal(pr[, 1], rep(mean(y), 10), tolerance=1e0)
  expect_equal(pr[, length(fit$lambda)], as.numeric(y), tolerance=1e-1)
})

test_that("missing, small sample", {
  X <- matrix(c(14, 25, 4, 7, 1, 2, 23, 11, 30, 18, 19, 26, 10, 22, 21, 27, 9, 5, 13, 17, 24, 12, 16, 29, 15, 20, 28, 3, 8, 6),
              nrow=10, ncol=3)
  b <- matrix(c(-1, 1, 0), nrow=3, ncol=1)
  e <- matrix(0.01 * (-1)^(1:10), nrow=10, ncol=1)
  y <- X %*% b + e
  X_miss <- X
  X_miss[1, 1] <- NA
  X_miss[2, 2] <- NA
  fit <- hmlasso(X_miss, as.numeric(y))
  plot(fit)

  # expect_equal(fit$beta[, 1], rep(0, length=ncol(X)), tolerance=1e-1)
  # expect_equal(fit$a0[1], mean(y), tolerance=1e0)
  expect_equal(fit$beta[, length(fit$lambda)], as.numeric(b), tolerance=2e-1)
  expect_equal(fit$a0[length(fit$lambda)], 0, tolerance=2e0)

  pr <- predict(fit, X)
  # expect_equal(pr[, 1], rep(mean(y), 10), tolerance=1e0)
  expect_equal(pr[, length(fit$lambda)], as.numeric(y), tolerance=5e0)
})

test_that("large sample hmlasso", {
  # setting
  n <- 100 # number of samples
  p <- 10 # number of dimensions
  r <- 0.5 # correlation levels
  eps <- 1 # noise level
  mu <- rep(0, length=p)
  Sigma <- matrix(r, nrow=p, ncol=p) # correlation among variables is 0
  diag(Sigma) <- 1 # variance of each variable is 1
  set.seed(0)
  X <- mvrnorm(n, mu, Sigma)
  b <- matrix(0, nrow = p, ncol = 1)
  # b[1:10,] <- (-1)^(1:10) * (1:10)
  b[1:4, ] <- c(1, -1, 2, -2)
  y <- X %*% b + eps * rnorm(n)
  X_miss <- X
  X_miss[sample(1:length(X), round(length(X) * 0.1))] <- NA # for random missing pattern
  fit <- hmlasso(X_miss, as.numeric(y))
  plot(fit)

  # expect_equal(fit$beta[, 1], rep(0, length=ncol(X)), tolerance=1e-1)
  # expect_equal(fit$a0[1], mean(y), tolerance=1e0)
  expect_equal(fit$beta[, length(fit$lambda)], as.numeric(b), tolerance=5e-1)
  expect_equal(fit$a0[length(fit$lambda)], 0, tolerance=1)

  pr <- predict(fit, X)
  expect_equal(pr[, 1], rep(mean(y), n), tolerance=1)
  expect_equal(pr[, length(fit$lambda)], as.numeric(y), tolerance=5e0)
})
