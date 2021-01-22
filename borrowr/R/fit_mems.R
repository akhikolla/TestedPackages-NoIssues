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


fit_mems <- function(mf, estimator, ndpost, exch_prob, ...) {
  # if(estimator != "bayesian_lm")
  #   stop("'bayesian_lm' is currently the only estimator implemented.")

  formula <- attr(mf, "formula")
  src_var <- attr(mf, "src_var")
  trt_var <- attr(mf, "trt_var")
  com_var <- attr(mf, "com_var")

  nc <- !is.na(com_var)

  ps <- attr(mf, "primary_source")
  sl <- levels(mf[, src_var])
  nl <- length(sl)


  # bad idea:
  # nm <- match(c(src_var, trt_var), names(mf))
  # names(mf)[nm] <- c("src", "trt")

  # sl <- levels(mf[, src_var])
  # ps <- sl[1L]
  # nl <- length(sl)

  il <- replicate(nl - 1, c(TRUE, FALSE), simplify = FALSE)
  il <- c(TRUE, il)
  im <- do.call(expand.grid, il)
  im <- t(as.matrix(im))
  dimnames(im) <- list(
    source = sl,
    MEM    = seq_len(ncol(im))
  )

  prior_mem <- numeric(ncol(im))
  names(prior_mem) <- seq_len(ncol(im))
  names(exch_prob)   <- rownames(im)[-1]

  mem_pate_post <- array(dim = c(ndpost, ncol(im)))
  mem_EY0 <- array(dim = c(ndpost, ncol(im)))
  mem_EY1 <- array(dim = c(ndpost, ncol(im)))
  dimnames(mem_pate_post) <- list(
    draw = seq_len(ndpost),
    MEM  = seq_len(ncol(im))
  )

  dimnames(mem_EY0) <- dimnames(mem_pate_post)
  dimnames(mem_EY1) <- dimnames(mem_pate_post)

  # Y <- matrix(model.response(mf, "numeric"), ncol = 1)
  # Y <- mf[, all.vars(formula[[2]])]
  # Y <- matrix(Y, ncol = 1)

  log_marg_like <- numeric(ncol(im))
  names(log_marg_like) <- paste0("MEM_", seq_len(ncol(im)))

  # in the event of bugs, check sorting and splitting
  Xf <- mf
  Xf <- Xf[order(Xf[, src_var]), ]
  Y <- Xf[, all.vars(formula[[2]])]
  # Y <- matrix(Y, ncol = 1)
  # for bayesian lm:
  if(estimator == "bayesian_lm") {
    Xm <- as.data.frame(model.matrix(formula, Xf))
    X0 <- Xf[Xf[, src_var] == ps, ]
    X1 <- Xf[Xf[, src_var] == ps, ]
    X0[, trt_var] <- 0
    X1[, trt_var] <- 1
    if (nc) {
      X0[, com_var] <- 1
      X1[, com_var] <- 1
    }
    X0 <- as.data.frame(model.matrix(formula, X0))
    X1 <- as.data.frame(model.matrix(formula, X1))

    beta_post_mean <- array(dim = c(ncol(Xm), ncol(im)))
    beta_post_var  <- array(dim = c(ncol(Xm), ncol(im)))
    dimnames(beta_post_mean) <- list(
      coef = colnames(Xm),
      MEM  = seq_len(ncol(im))
    )
    dimnames(beta_post_var) <- list(
      var  = colnames(Xm),
      MEM  = seq_len(ncol(im))
    )

  } else if (estimator == "BART") {
    Xm <- Xf[, c(all.vars(formula[[3]]), src_var), drop = FALSE]
    beta_post_mean <- NA
    beta_post_var  <- NA
  }

  attr(Xm, "src_var") <- src_var
  attr(Xm, "trt_var") <- trt_var
  attr(Xm, "primary_source") <- ps
  attr(Xm, "com_var") <- com_var

  sro <- Xf[, src_var]
  Xs <- split(Xm, sro)
  #Xfs <- lapply(Xfs, as.matrix)
  Ys <- split(Y, sro)
  # X0 <- Xs[[1]]
  # X1 <- Xs[[1]]
  # X0[, trt_var] <- 0 # bayes_lm only?
  # X1[, trt_var] <- 1 # bayes_lm only?
  #sepfits <- list()
  # ...
  # fit the BART model for each level of the source variable.

  # if(bart) {
  #
  # }

  if(estimator == "bayesian_lm") {
    sepfits <- mapply(bayes_lm,
      Y        = Ys,
      X        = Xs,
      X0       = list(X0),
      X1       = list(X1),
      ndpost   = ndpost,
      SIMPLIFY = FALSE)
    multicollinearity_flag <- any(sapply(sepfits, "[[", "multicollinearity_flag"))
    if (multicollinearity_flag)
      warning("Multicollinearity in design matrix X.")
  } else if (estimator == "BART") {
    #for BART only
    message("Fitting BART model to each data source...")
    sepfits <- mapply(bart,
      Y        = Ys,
      X        = Xs,
      #X0       = list(X0),
      #X1       = list(X1),
      estimate = as.list(rep(c(TRUE, FALSE), c(1, nl - 1))),
      ndpost   = ndpost,
      # theta    = 10,
      sparse   = TRUE,
      cont     = FALSE,
      SIMPLIFY = FALSE,
      ...)
  }
  for(mem in seq_len(ncol(im))) {
    if (estimator == "BART")
      message(sprintf("Fitting MEM %d of %d ...", mem, ncol(im)))
    cm <- im[, mem]
    # cn <- names(cm)[which(cm)]
    exch <- 1 * cm[-1]
    prior_mem[mem] <- prod(exch * exch_prob + (1 - exch) * (1 - exch_prob))
    tfits <- sepfits
    if (any(cm[-1L])) {
      tempX <- do.call(rbind, Xs[cm])
      tempY <- do.call(c, Ys[cm])
      if (estimator == "bayesian_lm") {
        tfits[[1]] <- bayes_lm(Y = tempY,
          X        = tempX,
          X0       = X0,
          X1       = X1,
          ndpost   = ndpost)
      } else if (estimator == "BART") {
        tfits[[1]] <- bart(Y = tempY,
          X        = tempX,
          #X0       = X0,
          #X1       = X1,
          estimate = TRUE,
          ndpost   = ndpost,
          # theta    = 10,
          sparse   = TRUE,
          cont     = FALSE,
          ...)
      }
    }
    ll <- sapply(tfits, "[[", "log_marg_like")[c(TRUE, !cm[-1L])]
    log_marg_like[mem] <- sum(ll)
    mem_pate_post[, mem] <- tfits[[1]]$pate_post
    mem_EY0[, mem] <- tfits[[1]]$EY0
    mem_EY1[, mem] <- tfits[[1]]$EY1
    if (estimator == "bayesian_lm") {
      beta_post_mean[, mem] <- colMeans(tfits[[1]]$beta_post)
      beta_post_var[, mem]  <- apply(tfits[[1]]$beta_post, 2, var)
    }
  }

  # if(estimator == "BART") {
  #   Xf <- mf
  #   Xf <- Xf[, all.vars(attr(mf, "terms")[[3]])]
  #   attr(Xf, "src_var") <- src_var
  #   attr(Xf, "trt_var") <- trt_var
  #   attr(Xf, "primary_source") <- attr(mf, "primary_source")
  #   # attr(Xf, "in_prim") <- attr(mf, "in_prim")
  #   Xfo <- Xf[order(Xf[, src_var]), ]
  #   sro <- Xfo[, src_var]
  #   # Xfo[, src_var] <- NULL
  #   Xfs <- split(Xfo, sro)
  #   #Xfs <- lapply(Xfs, as.matrix)
  #   Ys <- split(Y, sro)
  #   #X0 <- Xfs[[1]]
  #   #X1 <- Xfs[[1]]
  #   # X0[, trt_var] <- 0 # now done within the bart function
  #   # X1[, trt_var] <- 1 # now done within the bart function
  #   #sepfits <- list()
  #   # ...
  #   # fit the BART model for each level of the source variable.
  #   message("Fitting BART model to each data source...")
  #   sepfits <- mapply(bart,
  #     Y        = Ys,
  #     X        = Xfs,
  #     #X0       = list(X0),
  #     #X1       = list(X1),
  #     estimate = as.list(rep(c(TRUE, FALSE), c(1, nl - 1))),
  #     ndpost   = ndpost,
  #     beta     = beta,
  #     gf       = gf,
  #     eb       = eb,
  #     # theta    = 10,
  #     sparse   = TRUE,
  #     cont     = FALSE,
  #     SIMPLIFY = FALSE,
  #     ...)
  #   for(mem in seq_len(ncol(im))) {
  #     message(sprintf("Computing MEM %d of %d ...", mem, ncol(im)))
  #     cm <- im[, mem]
  #     # cn <- names(cm)[which(cm)]
  #     tfits <- sepfits
  #     if (any(cm[-1L])) {
  #       tempX <- do.call(rbind, Xfs[cm])
  #       tempY <- do.call(c, Ys[cm])
  #       tfits[[1]] <- bart(Y = tempY,
  #         X        = tempX,
  #         #X0       = X0,
  #         #X1       = X1,
  #         estimate = TRUE,
  #         ndpost   = ndpost,
  #         beta     = beta,
  #         gf       = gf,
  #         eb       = eb,
  #         # theta    = 10,
  #         sparse   = TRUE,
  #         cont     = FALSE,
  #         ...)
  #     }
  #     ll <- sapply(tfits, "[[", "log_marg_like")[c(TRUE, !cm[-1L])]
  #     log_marg_like[mem] <- sum(ll)
  #     mem_pate_post[, mem] <- tfits[[1]]$pate_post
  #   }
  # } else if (estimator == "bayesian_lm") {
  #   for(mem in seq_len(ncol(im))) {
  #     message(sprintf("Computing MEM %d of %d ...", mem, ncol(im)))
  #     cm <- im[, mem]
  #     cn <- names(cm)[which(cm)]
  #     temp_mf <- mf
  #     temp_mf[, src_var][temp_mf[, src_var] %in% cn] <- ps
  #     X <- model.matrix(terms(mf), temp_mf)
  #     X <- apply(X, 2, dropzeros)
  #     if(is.list(X))
  #       X <- do.call(cbind, X)
  #     X0 <- X1 <- X[mf[, src_var] == ps, ]
  #     X0[, trt_var] <- 0
  #     X1[, trt_var] <- 1
  #     fit <- bayes_lm(Y, X, X0, X1, ndpost)
  #     mem_pate_post[, mem] <- fit$pate_post
  #     log_marg_like[mem] <- fit$log_marg_like
  #   }
  # }

  lml <- log_marg_like
  cons  <- -max(lml)
  # post_probs <- exp(cons + log_marg_like) / sum(exp(cons + log_marg_like))
  post_probs <- exp(cons + lml) * prior_mem / sum(exp(cons + lml) * prior_mem)
  # post_probs[which(is.na(post_probs))] <- 0
  pate_post <- sample_posterior(mem_pate_post, post_probs)

  # better to return the matrices rather than the
  # weighted averages.
  EY0 <- sample_posterior(mem_EY0, post_probs)
  EY1 <- sample_posterior(mem_EY1, post_probs)

  out <- list(
    estimator     = estimator,
    pate_post     = pate_post,
    log_marg_like = log_marg_like,
    post_probs    = post_probs,
    MEMs          = im,
    exch_prob    = exch_prob
  )

  out$mem_pate_post <- mem_pate_post
  out$mem_EY0 <- mem_EY0
  out$mem_EY1 <- mem_EY1
  out$EY0 <- EY0
  out$EY1 <- EY1
  out$beta_post_mean <- beta_post_mean
  out$beta_post_var  <- beta_post_var

  out

}


# set.seed(100)
# n <- 100
# dat <- data.frame(
#   src = sample(letters[1:3], n, replace = TRUE),
#   trt = sample(c(1, 0), n, replace = TRUE),
#   x1  = rnorm(n),
#   x2  = rnorm(n),
#   y   = rnorm(n)
# )
#
# dat$src <- relevel(dat$src, ref = "b")
# str(dat)
#
# mf <- model.frame(y ~ src + trt + x1 * x2, dat)
# attr(mf, "src_var") <- "src"
# attr(mf, "trt_var") <- "trt"
# attr(mf, "primary_source") <- "b"
# # debug(fit_mems)
# # undebug(bart)
# # undebug(BART::wbart)
# foo <- fit_mems(mf, "BART", 10, 2)


# bar <- fit_mems(mf, "bayesian_lm")
# foo$log_marg_like
# bar$log_marg_like
# foo$post_probs
# bar$post_probs
