#  borrowr: estimate population average treatment effects with borrowing between data sources.
#  Copyright (C) 2019  Jeffrey A. Verdoliva Boatman
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

#Fit a Bayesian linear model with a normal-inverse gamma prior.
#
#@details Does not include a formula argument because the formula is passed indirectly
#via the design matrix \code{X}.
#@param Y Outcome variable. Must be a column vector.
#@param X Design matrix. Must have the same number of rows as \code{Y}.
#@param X0 A design matrix used for the causal estimator. It is a modified version
#          of \code{X} representing the counterfactual where all observations are
#          assigned to the treatment group 0 and are compliant.
#@param X1 A design matrix used for the causal estimator. It is a modified version
#          of \code{X} representing the counterfactual where all observations are
#          assigned to the treatment group 1 and are compliant.
#@param ndpost The desired number of draws from the posterior distribution of
#              \out{E(Y<sub>1</sub> - Y<sub>0</sub>)}.
# bayes_lm <- function(Y, X, X0, X1, ndpost) {
#   # production version of bayes_lm
#
#   X <- as.matrix(X)
#   X0 <- as.matrix(X0)
#   X1 <- as.matrix(X1)
#   p <- ncol(X)
#   n <- nrow(X)
#
#   # flm <- lm(Y ~ X + 0)
#
#   # uninformative Normal-Inverse Gamma prior specifications
#   # muBeta <- unname(coef(flm))
#   muBeta <- rep(0, p)
#   # rho    <- 0.2
#   rho <- 0
#   Vbeta  <- 100 * (rho * matrix(1, p, p) + (1 - rho) * diag(p))
#
#   # specification for IG(a, b) prior on sigma ^ 2. This
#   # prior is specified using the frequentist point estimate
#   # and asymptotic variance of sigma ^ 2.
#   # ss_mu <- 10 ^ 2    # old: mean for sigma ^ 2
#   # ss_sd <- 2         # old: sd for sigma ^ 2
#   #ss_mu <- summary(flm)$sigma ^ 2
#   #if (is.na(ss_mu))
#     #stop("insufficient number of observations to fit model.")
#   # ss_sd <- sqrt(2 * ss_mu ^ 4 / n)
#   # mean and sd re-paramaterized for inverse gamma dist'n
#   # b <- ss_mu ^ 3 / ss_sd ^ 2 + ss_mu
#   # a <- (b + ss_mu) / ss_mu
#   a <- 0.1
#   b <- 0.1
#
#   # b / (a - 1)
#   # b ^ 2 / ((a - 1) ^ 2 * (a - 2))
#
#   Vbeta_inv <- solve(Vbeta)
#   tX <- t(X)
#   tXX <- tX %*% X
#   tryCatch(solve(tXX),
#     error = function(e) warning("Multicollinearity in design matrix X."))
#
#
#   # posterior paramaters
#   V_star <- solve(Vbeta_inv + tXX)
#   mu_star <- V_star %*% (Vbeta_inv %*% muBeta + tX %*% Y)
#   a_star <- a + n / 2
#   b_star <- b + 1 / 2 * c((t(muBeta) %*% Vbeta_inv %*% muBeta +
#       t(Y) %*% Y - t(mu_star) %*% solve(V_star) %*% mu_star))
#
#   # marginal likelihood ---
#   mu    <- X %*% muBeta # should be == 0 unless muBeta updated
#   shape <- b / a * (diag(n) + X %*% Vbeta %*% tX)
#   log_marg_like <- dmvt(c(Y), c(mu), shape, df = 2 * a)
#
#   # Monte Carlo verification of log_marg_like
#   # lml <- numeric(1e3)
#   # for (jj in seq_len(1e3)) {
#   #   rss <- 1 / rgamma(1, a, b)
#   #   rBeta <- mvtnorm::rmvnorm(1, muBeta, sigma =  rss * Vbeta)
#   #   rmu <- X %*% t(rBeta)
#   #   lml[jj] <- sum(dnorm(c(Y), c(rmu), sd = sqrt(rss), log = TRUE))
#   # }
#
#   # don't name dimensions here. Will be done in the calling function.
#   beta_post <- array(dim = c(ndpost, p))
#
#   colnames(beta_post) <- colnames(X)
#
#   Y0 <- array(dim = c(ndpost, nrow(X0)))
#   Y1 <- array(dim = c(ndpost, nrow(X1)))
#
#   # draw from posterior
#   tryCatch(
#     ss_post <- 1 / rgamma(ndpost, a_star, b_star),
#     error = function(cnd) stop("problem drawing from sigma ^ 2 posterior.")
#   )
#
#
#   for(ii in seq_len(ndpost)) {
#     beta_post[ii, ] <- mu_star + t(chol(ss_post[ii] * V_star)) %*% rnorm(p)
#     mu_post  <- c(X %*% beta_post[ii, ])
#     Y0[ii, ] <- t(X0 %*% beta_post[ii, ])
#     Y1[ii, ] <- t(X1 %*% beta_post[ii, ])
#   }
#
#   # bs <- cov(beta_post)
#   # fr <- vcov(summary(lm(Y ~ X + 0)))
#   # round(bs / fr, 2)
#   # colMeans(beta_post) / coef(lm(Y ~ X + 0))
#
#   pate_post <- apply(Y1 - Y0, 1, bayes_boot_mean)
#
#   out <- list(
#     log_marg_like = log_marg_like,
#     pate_post     = pate_post)
#   out$EY0 <- rowMeans(Y0)
#   out$EY1 <- rowMeans(Y1)
#   out$beta_post <- beta_post
#
#   out
# }

bayes_lm <- function(Y, X, X0, X1, ndpost) {
  # development version of bayes_lm
  multicollinearity_flag <- FALSE
  X <- as.matrix(X)
  tX <- t(X)
  tXX <- tX %*% X
  X0 <- as.matrix(X0)
  X1 <- as.matrix(X1)
  p <- ncol(X)
  n <- nrow(X)

  # tryCatch(tXXinv <-  solve(tXX),
  #  error = function(e) warning("Multicollinearity in design matrix X."))
  tr <- try(tXXinv <- solve(tXX), silent = TRUE)
  if (class(tr)[1L] == "try-error") {
    multicollinearity_flag <- TRUE
    # warning("Multicollinearity in design matrix X.")
    Xm <- colMeans(X)
    Xs <- apply(X, 2, sd)
    #Xm <- t(array(Xm, dim = rev(dim(X))))
    #Xs <- t(array(Xs, dim = rev(dim(X))))
    #dimnames(Xm) <- dimnames(X)
    #dimnames(Xs) <- dimnames(X)
    if ("(Intercept)" %in% colnames(X)) {
      Xm["(Intercept)"] <- 0
      Xs["(Intercept)"] <- 1
      lambda <- diag(c(0, rep(0.1, p - 1)))
    } else {
      lambda <- 0.1 * diag(p)
    }
    X  <- stdize_matrix(X,  Xm, Xs)
    X0 <- stdize_matrix(X0, Xm, Xs) # (X0 - Xm) / Xs
    X1 <- stdize_matrix(X1, Xm, Xs) # (X1 - Xm) / Xs
    tX <- t(X)
    tXX <- tX %*% X + lambda
    tryCatch(tXXinv <- solve(tXX),
      error = function(e) stop("regularization failed."))
  }


  flm <- lm(Y ~ X + 0)

  # uninformative Normal-Inverse Gamma prior specifications
  # muBeta <- unname(coef(flm))
  # muBeta <- rep(0, p)
  # muBeta <- c(coef(flm)[1], rep(0, p - 1))
  muBeta <- c(mean(Y), rep(0, p - 1))
  # rho    <- 0.2
  # rho <- 0
  # Vbeta  <- 100 * (rho * matrix(1, p, p) + (1 - rho) * diag(p))
  # Vbeta <- solve(1 / n * tXX)
  Vbeta <- n * tXXinv # solve(1 / n * tXX)

  # specification for IG(a, b) prior on sigma ^ 2. This
  # prior is specified using the frequentist point estimate
  # and asymptotic variance of sigma ^ 2.
  # ss_mu <- 10 ^ 2    # old: mean for sigma ^ 2
  # ss_sd <- 2         # old: sd for sigma ^ 2
  ss_mu <- sigma(flm) ^ 2
  #if (is.na(ss_mu))
    #stop("insufficient number of observations to fit model.")
  # ss_sd <- sqrt(2 * ss_mu ^ 4 / n)
  ss_sd <- sqrt(2 * ss_mu ^ 4 / 1)
  # mean and sd re-paramaterized for inverse gamma dist'n
  b <- ss_mu ^ 3 / ss_sd ^ 2 + ss_mu
  a <- (b + ss_mu) / ss_mu
  # a <- 0.1
  # b <- 0.1
  # a <- 2.58 / 2
  # b <- 0.28


  # b / (a - 1)
  # b ^ 2 / ((a - 1) ^ 2 * (a - 2))

  # Vbeta_inv <- solve(Vbeta) # could be simplified to '1 / n * tXX'
  Vbeta_inv <- 1 / n * tXX
  # tX <- t(X)
  # tXX <- tX %*% X


  # posterior paramaters
  V_star <- solve(Vbeta_inv + tXX)
  mu_star <- V_star %*% (Vbeta_inv %*% muBeta + tX %*% Y)
  a_star <- a + n / 2
  b_star <- b + 1 / 2 * c((t(muBeta) %*% Vbeta_inv %*% muBeta +
      t(Y) %*% Y - t(mu_star) %*% solve(V_star) %*% mu_star))

  # marginal likelihood ---
  mu    <- X %*% muBeta # should be == 0 unless muBeta updated
  shape <- b / a * (diag(n) + X %*% Vbeta %*% tX)
  log_marg_like <- dmvt(c(Y), c(mu), shape, df = 2 * a)

  # Monte Carlo verification of log_marg_like
  # lml <- numeric(1e3)
  # for (jj in seq_len(1e3)) {
  #   rss <- 1 / rgamma(1, a, b)
  #   rBeta <- mvtnorm::rmvnorm(1, muBeta, sigma =  rss * Vbeta)
  #   rmu <- X %*% t(rBeta)
  #   lml[jj] <- sum(dnorm(c(Y), c(rmu), sd = sqrt(rss), log = TRUE))
  # }

  # don't name dimensions here. Will be done in the calling function.
  beta_post <- array(dim = c(ndpost, p))

  colnames(beta_post) <- colnames(X)

  Y0 <- array(dim = c(ndpost, nrow(X0)))
  Y1 <- array(dim = c(ndpost, nrow(X1)))

  # draw from posterior
  tryCatch(
    ss_post <- 1 / rgamma(ndpost, a_star, b_star),
    error = function(cnd) stop("problem drawing from sigma ^ 2 posterior.")
  )


  for(ii in seq_len(ndpost)) {
    beta_post[ii, ] <- mu_star + t(chol(ss_post[ii] * V_star)) %*% rnorm(p)
    mu_post  <- c(X %*% beta_post[ii, ])
    Y0[ii, ] <- t(X0 %*% beta_post[ii, ])
    Y1[ii, ] <- t(X1 %*% beta_post[ii, ])
  }

  # bs <- cov(beta_post)
  # fr <- vcov(summary(lm(Y ~ X + 0)))
  # round(bs / fr, 2)
  # colMeans(beta_post) / coef(lm(Y ~ X + 0))

  pate_post <- apply(Y1 - Y0, 1, bayes_boot_mean)

  out <- list(
    log_marg_like = log_marg_like,
    pate_post     = pate_post)
  out$EY0 <- rowMeans(Y0)
  out$EY1 <- rowMeans(Y1)
  out$beta_post <- beta_post
  out$multicollinearity_flag <- multicollinearity_flag

  out
}
