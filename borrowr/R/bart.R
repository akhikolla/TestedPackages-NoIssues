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
bart <- function(Y, X, estimate, nprior = 100, ntree = 200, ndpost,...) {

  # cl <- match.call(expand.dots = TRUE)
  #X  <- as.data.frame(X)
  #X0 <- as.data.frame(X0)
  #X1 <- as.data.frame(X1)

  src_var <- attr(X, "src_var")
  trt_var <- attr(X, "trt_var")
  prm_src <- attr(X, "primary_source")
  com_var <- attr(X, "com_var")

  nc <- !is.na(com_var)

  xo <- X[, src_var] == prm_src
  X[, src_var] <- NULL
  # in_prim <- attr(X, "in_prim")

  # think this is unnecessary:
  #X[,  src_var] <- droplevels(X[,  src_var])
  #X0[, src_var] <- droplevels(X0[, src_var])
  #X1[, src_var] <- droplevels(X1[, src_var])

  #X0 <- X[X[, src_var] == primary_source, ]
  #X1 <- X[X[, src_var] == primary_source, ]
  # X0[, trt_var] <- 0
  # X1[, trt_var] <- 1

#  sv <- unique(X[, src_var])
#  nsv <- length(sv)
#  nv <- ncol(X)
#  if (nsv == 1) {
#    X[, src_var]  <- NULL
#    #X0[, src_var] <- NULL
#    #X1[, src_var] <- NULL
#  }

  #X0[, src_var] <- droplevels(X0[, src_var])
  #X1[, src_var] <- droplevels(X1[, src_var])


  # build call to bartModelMatrix, used for predict
  # !!! possible source of bugs?
  # bc <- cl
  # m <- match(names(formals(BART::bartModelMatrix)), names(bc), 0L)
  # bc <- bc[c(1L, m)]
  # bc[[1L]] <- quote(BART::bartModelMatrix)
  # bc$rm.const <- TRUE
  #
  # bmm <- bc
  # bm0 <- bc
  # bm1 <- bc
  #
  # bmm$X <- quote(X)
  # bm0$X <- quote(X0)
  # bm1$X <- quote(X1)
  #
  # bmm <- eval(bmm)
  # bm0 <- eval(bm0)
  # bm1 <- eval(bm1)

  # bmm <- BART::bartModelMatrix(X, cont = FALSE, rm.const = TRUE)
  # bm0 <- eval(bm0)
  # bm1 <- eval(bm1)




  # bc <- cl
  # bc$x.train <- bc$X
  # bc$y.train <- bc$Y
  # m <- match(names(formals(BART::wbart)), names(bc), 0L)
  # bc <- bc[c(1L, m)]
  # bc[[1L]] <- quote(BART::wbart)

  n <- length(Y)

  # scale Y
  sc <- function(x) {
    x_orig <- x

    x <- x - mean(x)
    x <- x / diff(range(x))

    fit <- lm(x ~ x_orig)
    a <- unname(fit$coef[1])
    b <- unname(fit$coef[2])
    list(a = a, b = b)
  }

  aa <- sc(Y)$a
  bb <- sc(Y)$b
  y_star <- c(aa + bb * Y)

  # ger prior params for sigma
  naive_sigma <- summary(lm(y_star ~ ., data = X))$sigma
  nu     <- 3
  lambda <- qchisq(0.1, nu) * naive_sigma ^ 2 / nu
  gamma  <- 1 / (16 * ntree * naive_sigma ^ 2)


  # bfit <- BART::wbart(x.train = X, y.train = y_star, ndpost = ndpost,
  #  rm.const = TRUE, ...)
  # bfit <- quiet(BART::wbart)(x.train = X, y.train = y_star, ndpost = ndpost, ...)
  bfit <- quiet(wbart)(x.train = X, y.train = y_star, ndpost = ndpost, gamma = gamma, ...)
  # bfit <- eval(bc)


  # marginal likelihood  ----
  bart_for_prior  <- bfit

  # draw from prior of gamma * sigma ^ 2
  sigma_mu_prior <- sqrt(gamma * nu * lambda / rgamma(nprior, nu / 2, 1 / 2))

  # create prior tree
  cut_lens  <- sapply(bfit$treedraws$cutpoints, length)
  p         <- length(cut_lens)

  # possible bugs?
  # if(nsv > 1) {
  #   var_probs <- rep(c(1 / (nv * nsv), 1 / nv), c(nsv, nv - 1))
  #   if(length(var_probs) != p)
  #     stop("Check 'var_probs' computation.")
  # } else {
  #   var_probs <- NULL
  # }
  var_probs <- NULL

  # -- find number of shared nodes between observations.
  # -- to do this, treat each tree as a separate draw from the posterior.
  # -- then find the number of shared nodes between observations.

  tree_struct <- make_prior_trees(nprior = nprior,
    ntree     = ntree,
    cut_lens  = cut_lens,
    var_probs = var_probs,
    #depth     = depth,
    sigma_mu  = sigma_mu_prior)
  old_first_line <- paste(nprior, ntree, p, collapse = " ")
  new_first_line <- paste(nprior * ntree, "1", p, collapse = " ")

  bart_for_prior$treedraws$trees  <- gsub(old_first_line, new_first_line, tree_struct)

  # tc    <- textConnection(bart_for_prior$treedraws$tree)
  # trees <- read.table(file = tc, fill = TRUE, row.names = NULL, header = FALSE, col.names = c('node', 'var', 'cut', 'leaf'))
  # close(tc)
  # nds <- trees[!is.na(trees$leaf) & trees$leaf == 0, ]
  # prop.table(table(nds$var))
  # table(nds$cut, nds$var)

  # temp_fitted <- quiet(predict)(bart_for_prior, bmm)
  # xp <- BART::bartModelMatrix(X, numcut = 100, usequants = FALSE,
  #   cont = FALSE, xinfo = matrix(0, 0, 0), rm.const = TRUE)$X
  xp <- bartModelMatrix(X, numcut = 100, usequants = FALSE,
    cont = FALSE, xinfo = matrix(0, 0, 0), rm.const = TRUE)$X
  temp_fitted <- quiet(predict)(bart_for_prior, xp)

  pnum <- rep(seq_len(nprior), each = ntree)

  # iterate over the draws from the prior. For each,
  # compute correlation matrix (a function of the
  # number of shared leaves), then compute the
  # likelihood using that correlation matrix.
  log_marg_likes <- numeric(nprior)
  m <- ntree
  for(ip in seq_len(nprior)) {
    i_temp_fitted <- temp_fitted[pnum == ip, ]
    R <- matchesToCor(i_temp_fitted)$R
    V <- nu * lambda * (m * gamma * R + diag(n))
    log_marg_likes[ip] <- dmvt(y_star, rep(0, n), V / nu, df = nu, log = TRUE)
  }

  log_marg_like <- max(log_marg_likes) + log(mean(exp(log_marg_likes - max(log_marg_likes))))
  log_marg_like <- n * log(bb) + log_marg_like # re-scale the log likelihood


  # m <- match(names(formals(BART::bartModelMatrix)), names(formals(BART::wbart)), 0L)
  # m <- formals(BART::wbart)[m]
  # bmm <- c(quote(BART::bartModelMatrix), m)
  # bmm$X <- quote(X)
  # bmm <- as.call(bmm)

  pate_post <- NA
  EY0 <- NA
  EY1 <- NA
  if(estimate) {
    # Xt <- BART::bartModelMatrix(X, numcut = 100, usequants = FALSE,
    #   cont = FALSE, xinfo = matrix(0, 0, 0), rm.const = TRUE)$X
    # X0 <- X1 <- BART::bartModelMatrix(X, cont = FALSE, rm.const = TRUE)
    # X0 <- X1 <- eval(bmm)$X
    #X0 <- X0[in_prim, ]
    #X1 <- X1[in_prim, ]
    X0 <- X1 <- xp[xo, , drop = FALSE]
    X0[, trt_var] <- 0
    X1[, trt_var] <- 1
    if (nc) {
      X0[, com_var] <- 1
      X1[, com_var] <- 1
    }
    Y0 <- quiet(predict)(bfit, X0)
    Y1 <- quiet(predict)(bfit, X1)
    # Y0 <- quiet(predict)(bfit, bm0)
    # Y1 <- quiet(predict)(bfit, bm1)
    Y0 <- (Y0 - aa) / bb
    Y1 <- (Y1 - aa) / bb
    pate_post <- apply(Y1 - Y0, 1, bayes_boot_mean)
    EY0 <- rowMeans(Y0)
    EY1 <- rowMeans(Y1)
  }

  out <- list(
    log_marg_like = log_marg_like,
    pate_post     = pate_post)
  out$EY0 <- EY0
  out$EY1 <- EY1

  out
}

# set.seed(100)
# n <- 1e2
# dat <- data.frame(
#   src = sample(c("a", "b", "c"), n, replace = TRUE),
#   trt = sample(c(1, 0), n, replace = TRUE),
#   x1  = rnorm(n),
#   x2  = rnorm(n),
#   y   = rnorm(n)
#  )
#
# X <- model.matrix(y ~ trt + x1 + x2 + 0, dat)
# X0 <- X1 <- model.matrix(y ~ trt + x1 + x2 + 0, subset(dat, src == "a"))
# X0[, 'trt'] <- 0
# X1[, 'trt'] <- 1
# Y <- matrix(dat$y, ncol = 1)
#
# debug(bart)
# system.time(foo <- bart(Y, X, X0, X1, nprior = 10, estimate = TRUE, ntree = 200, nskip = 5))
