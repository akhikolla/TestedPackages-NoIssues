## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>",
  message = FALSE, warning = FALSE
)

## ---- eval = F----------------------------------------------------------------
#  install.packages("ADtools")  # stable version
#  devtools::install_github("kcf-jackson/ADtools")  # development version

## -----------------------------------------------------------------------------
library(ADtools)
f <- function(X, y) X %*% y
X <- randn(2, 2)
y <- matrix(c(1, 1))
print(list(X = X, y = y, f = f(X, y)))

## -----------------------------------------------------------------------------
# Full Jacobian matrix
f_AD <- auto_diff(f, at = list(X = X, y = y))
f_AD@dx   # returns a Jacobian matrix

## -----------------------------------------------------------------------------
f_AD <- auto_diff(f, at = list(X = X, y = y), wrt = "y")
f_AD@dx   # returns a partial Jacobian matrix

## -----------------------------------------------------------------------------
f_FD <- finite_diff(f, at = list(X = X, y = y))
f_FD

## -----------------------------------------------------------------------------
set.seed(123)
n <- 1000
p <- 3
X <- randn(n, p)
beta <- randn(p, 1)
y <- X %*% beta + rnorm(n)

## -----------------------------------------------------------------------------
gradient_descent <- function(f, vary, fix, learning_rate = 0.01, tol = 1e-6, show = F) {
  repeat {
    df <- auto_diff(f, at = append(vary, fix), wrt = names(vary))
    if (show) print(df@x)
    delta <- learning_rate * as.numeric(df@dx)
    vary <- relist(unlist(vary) - delta, vary)
    if (max(abs(delta)) < tol) break
  }
  vary
}

## -----------------------------------------------------------------------------
lm_loss <- function(y, X, beta) sum((y - X %*% beta)^2)

# Estimate
gradient_descent(
  f = lm_loss, vary = list(beta = rnorm(p, 1)), fix = list(y = y, X = X),  learning_rate = 1e-4
) 
# Truth
t(beta)

## -----------------------------------------------------------------------------
set.seed(123)
n <- 30  # small data
p <- 10
X <- randn(n, p)
beta <- randn(p, 1)
y <- X %*% beta + rnorm(n)

## ---- eval = F----------------------------------------------------------------
#  gibbs_gaussian <- function(X, y, b_0, B_0, alpha_0, delta_0, num_steps = 1e4, burn_ins = ceiling(num_steps / 10)) {
#    # Initialisation
#    init_sigma <- 1 / sqrt(rgamma0(1, alpha_0 / 2, scale = 2 / delta_0))
#  
#    n <- length(y)
#    alpha_1 <- alpha_0 + n
#    sigma_g <- init_sigma
#    inv_B_0 <- solve(B_0)
#    inv_B_0_times_b_0 <- inv_B_0 %*% b_0
#    XTX <- crossprod(X)
#    XTy <- crossprod(X, y)
#    beta_res <- vector("list", num_steps)
#    sigma_res <- vector("list", num_steps)
#  
#    pb <- txtProgressBar(1, num_steps, style = 3)
#    for (i in 1:num_steps) {
#      # Update beta
#      B_g <- solve(sigma_g^(-2) * XTX + inv_B_0)
#      b_g <- B_g %*% (sigma_g^(-2) * XTy + inv_B_0_times_b_0)
#      beta_g <- t(rmvnorm0(1, b_g, B_g))
#  
#      # Update sigma
#      delta_g <- delta_0 + sum((y - X %*% beta_g)^2)
#      sigma_g <- 1 / sqrt(rgamma0(1, alpha_1 / 2, scale = 2 / delta_g))
#  
#      # Keep track
#      beta_res[[i]] <- beta_g
#      sigma_res[[i]] <- sigma_g
#      setTxtProgressBar(pb, i)
#    }
#  
#    # Compute and return the posterior mean
#    sample_ids <- (burn_ins + 1):num_steps
#    beta_pmean <- Reduce(`+`, beta_res[sample_ids]) / length(sample_ids)
#    sigma_pmean <- Reduce(`+`, sigma_res[sample_ids]) / length(sample_ids)
#    list(sigma = sigma_pmean, beta = beta_pmean)
#  }

## ---- eval = F----------------------------------------------------------------
#  gibbs_deriv <- auto_diff(
#    gibbs_gaussian,
#    at = list(
#      b_0 = numeric(p), B_0 = diag(p), alpha_0 = 4, delta_0 = 4,
#      X = X, y = y, num_steps = 50, burn_ins = 5 # Numbers are reduced for CRAN
#    ),
#    wrt = c("b_0", "B_0", "alpha_0", "delta_0")
#  )

