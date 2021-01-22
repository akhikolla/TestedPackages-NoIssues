#' Fit a model using a design matrix
#'
#' @param X matrix of explanatory variables
#' @param y vector of objective variable
#' @param family family of regression: "gaussian" (default) or "binomial"
#' @param impl implementation language of optimization (only "cpp" is supported)
#' @param lambda.min.ratio ratio of max lambda and min lambda (ignored if lambda is specified)
#' @param nlambda the number of lambda (ignored if lambda is specified)
#' @param lambda lambda sequence
#' @param min_eig_th threshold of the minimum eigenvalue in the PSD matrix problem.
#' @param use method to calculate correlation matrix from missing data (default "pairwise.complete.obs")
#' @param positify method for solving PSD matrix
#' @param weight_power weighting power (default 0 meaning no-weighting)
#' @param eig_tol tol parameter in eigs_sym function
#' @param eig_maxitr maxitr parameter in eigs_sym
#' @param mu augmented Lagrangian parameter
#' @param verbose whether output verbose warnings and messages (default FALSE)
#' @param ... parameters for optimization
#'
#' @return lasso model
#' \item{beta}{coefficients}
#' \item{beta_standard}{standardized coefficients}
#' \item{a0}{intercepts}
#' \item{lambda}{regularization parameters}
#' \item{family}{family}
#'
#' @examples
#' X_incompl <- as.matrix(iris[, 1:3])
#' X_incompl[1:5,1] <- NA
#' X_incompl[6:10,2] <- NA
#' y <- iris[, 4]
#' fit <- hmlasso(X_incompl, y, nlambda=50, lambda.min.ratio=1e-2)
#' plot(fit)
#'
#' @export
hmlasso <- function(X, y, family="gaussian", impl="cpp",
                       lambda.min.ratio =  1e-2,
                       nlambda = 100,
                       lambda = NULL,
                       min_eig_th = 1e-6,
                       use = "pairwise.complete.obs",
                       positify="diag",
                       weight_power = 1,
                       eig_tol = 1e-8,
                       eig_maxitr = 1e+8,
                       mu = 1,
                       verbose = FALSE,
                       ...) {

  n <- nrow(X)
  eigen_Gamma <- NULL

  t1 <- Sys.time()
  if (sum(is.na(X)) + sum(is.na(y)) > 0) {
    # missing_lasso if X or y include NA
    # return(missing_lasso(X, y, ...))
    # constract standardized X and y
    X_mean <- apply(X, 2, function(v){mean(v, na.rm=TRUE)})
    X_sd <- apply(X, 2, function(v){
      sqrt((sum(!is.na(v))-1) / sum(!is.na(v)) * var(v, na.rm=TRUE))
    })
    X_tilde <- sweep(sweep(X, 2, X_mean, "-"), 2, X_sd, "/")
    y_tilde <- y - mean(y, na.rm=TRUE)

    # construct covariance matrices
    # Gamma <- (n-1) / n * cov(X_tilde, use=use) #sample covariance
    Gamma <- cor(X_tilde, use=use) #sample covariance # add 20180209
    Gamma <- apply(Gamma, c(1,2), function(v) {
      if (is.na(v) | is.infinite(v) | is.nan(v)) { 0 } else { v }
    })
    # Gamma[is.na(Gamma)] <- 0
    ill_obs <- which((t(!is.na(X_tilde)) %*% (!is.na(X_tilde))) <= 2)
    Gamma[ill_obs] <- 0
    diag(Gamma) <- 1

    min_eig <- eigs_sym(Gamma, k=1, which="SA", opts=list(tol=eig_tol, maxitr=eig_maxitr, retvec=FALSE))$value
    if (min_eig <= min_eig_th) {
      if (verbose) {
        message(paste0("The minimum eigenvalue is smaller than the specified threshold: ", min_eig, "<=", min_eig_th))
      }
      if (positify=="diag") {
        Gamma <- diag_positify(Gamma, min_eig_th, min_eig)
      } else if (positify=="projection") {
        eigen_Gamma <- eigen(Gamma)
        Gamma <- proj_positify(Gamma, min_eig_th, eigen_Gamma)
      } else if (positify=="mean") {
        Gamma <- mean_positify(X_tilde)
      } else if (positify=="admm_max") {
        Gamma <- admm_positify(Gamma, X, epsilon=min_eig_th, mu=mu, norm="max", cor_A=FALSE, cor_B=FALSE,
                               weight_power=weight_power, verbose=verbose)
      } else if (positify=="admm_frob") {
        Gamma <- admm_positify(Gamma, X, epsilon=min_eig_th, mu=mu, norm="frobenius", cor_A=FALSE, cor_B=FALSE,
                               weight_power=weight_power, verbose=verbose)
      }
    }
    # check positive semi-definiteness
    # min_eig2 <- eigs_sym(Gamma, k=1, which="SA")$value
    min_eig2 <- eigs_sym(Gamma, k=1, which="SA", opts=list(tol=eig_tol, maxitr=eig_maxitr, retvec=FALSE))$value
    if (min_eig2 < min_eig_th * 0.99) {
      if (verbose) {
        warning(paste0("Computed min eig = ", min_eig2, ", but min eig threshold = ", min_eig_th))
        warning("Fail to positify twice")
      }
      eigen_Gamma <- eigen(Gamma)
      Gamma <- proj_positify(Gamma, min_eig_th, eigen_Gamma)
      min_eig3 <- eigs_sym(Gamma, k=1, which="SA", opts=list(tol=eig_tol, maxitr=eig_maxitr, retvec=FALSE))$value
      if (min_eig3 < min_eig_th * 0.99) {
        for (kk in 1:10) {
          if (verbose) {
            warning(paste0("Computed min eig = ", min_eig3, ", but min eig threshold = ", min_eig_th))
            warning("Fail to positify third times")
          }
          eigen_Gamma <- eigen(Gamma)
          Gamma <- proj_positify(Gamma, min_eig_th*10^kk, eigen_Gamma)
          min_eig3 <- eigs_sym(Gamma, k=1, which="SA", opts=list(tol=eig_tol, maxitr=eig_maxitr, retvec=FALSE))$value
          if (min_eig3 >= min_eig_th) {
            break
          }
        }
      }
      if (min_eig3 < min_eig_th * 0.99) {
        stop(paste0("Computed min eig = ", min_eig3, ", but min eig threshold = ", min_eig_th, "."))
      }
    }
    # Gamma <- as.matrix(nearPD(Gamma, corr=TRUE)$mat)
    # min_eig <- eigs_sym(Gamma, k=1, which="SA")$value
    # if (min_eig <= min_eig_th) {
    #   Gamma <- Gamma + diag(- min_eig + 1e-2, nrow=nrow(Gamma), ncol=ncol(Gamma))
    # }
    # Gamma <- t(X_tilde) %*% X_tilde / n
    # gamma <- t(X_tilde) %*% y_tilde / n #sample covariance
    gamma <- apply(X_tilde, 2, function(v){
      fill <- which(!is.na(v) & !is.na(as.vector(y_tilde)))
      g <- sum(v[fill] * y_tilde[fill]) / length(fill)
      if (is.nan(g) | is.na(g)) {g <- 0}
      return(g)
    })

    # remove
    # gamma <- apply(X_tilde, 2, function(v){
    #   fill <- which(!is.na(v) & !is.na(as.vector(y_tilde)))
    #   sum(v[fill] * y_tilde[fill]) / length(fill)
    # })
    # gamma[is.na(gamma)] <- 0

  } else {
    # constract standardized X and y
    X_mean <- apply(X, 2, function(v){mean(v, na.rm=TRUE)})
    X_sd <- apply(X, 2, function(v){
      sqrt((sum(!is.na(v))-1) / sum(!is.na(v)) * var(v, na.rm=TRUE))
    })
    X_tilde <- sweep(sweep(X, 2, X_mean, "-"), 2, X_sd, "/")
    y_tilde <- y - mean(y, na.rm=TRUE)

    # construct covariance matrices
    # Gamma <- (n-1) / n * cov(X_tilde, use=use) #sample covariance
    Gamma <- t(X_tilde) %*% X_tilde / n
    gamma <- t(X_tilde) %*% y_tilde / n #sample covariance
  }
  t2 <- Sys.time()
  if (verbose) {
    message(paste0("Covariance matrix estimation takes ", format(t2-t1, digits = 3), " ", units(t2-t1)))
  }

  # if (is.null(lambda)) {
  if (missing(lambda)) {
    lambda_max <- max(gamma)
    lambda <- lambda_max * exp(seq(from=0,
                                   to=log(lambda.min.ratio),
                                   by=log(lambda.min.ratio)/(nlambda-1)))
  }

  # fit
  if (family=="gaussian") {
    fit <- cov_lasso(Gamma, gamma, lambda=lambda, impl=impl, ...)
    if (is.null(dim(fit$beta_standard))) {
      fit$beta <- fit$beta_standard / X_sd
      if (is.null(names(fit$beta))) {
        if (is.null(colnames(X))) {
          names(fit$beta) <- paste0("V",1:length(fit$beta))
        } else {
          names(fit$beta) <- colnames(X)
        }
      }
    } else {
      fit$beta <- sweep(fit$beta_standard, 1, X_sd, FUN="/")
      if (is.null(rownames(fit$beta))) {
        if (is.null(colnames(X))) {
          # rownames(fit$beta) <- paste0("V",1:length(fit$beta))
        } else {
          rownames(fit$beta) <- colnames(X)
        }
      }
    }
    fit$a0 <- mean(y) - t(X_mean) %*% fit$beta

  } else if (family=="binomial") {
    stop("Not yet implemented.")
  } else {
    stop("Specified family is not supported.")
  }
  t3 <- Sys.time()
  if (verbose) {
    message(paste0("Regression takes ", format(t3-t2, digits = 3), " ", units(t3-t2)))
  }


  fit$family <- family
  fit$X_mean <- X_mean
  fit$min_eig <- ifelse(exists("min_eig"), min_eig, 0)
  if (exists("Gamma")) {
    fit$Gamma <- Gamma
  } else {
    fit$Gamma <- 0
  }
  if (exists("gamma")) {
    fit$gamma <- gamma
  } else {
    fit$gamma <- 0
  }
  # fit$Gamma <- ifelse(exists("Gamma"), Gamma, 0)
  # fit$gamma <- ifelse(exists("gamma"), gamma, 0)
  fit$min_eig_th <- min_eig_th
  fit$eigen_Gamma <- eigen_Gamma
  class(fit) <- "hmlasso"
  return(fit)
}


cov_lasso <- function(Gamma, gamma,
                      lambda.min.ratio = 0.0001,
                      nlambda = 100,
                      lambda = NULL,
                      R = NULL,
                      funcR = function(G){abs(G)^2},
                      maxit = 1e+4,
                      eps = 1e-04,
                      warm = "lambda",
                      init.beta = NULL,
                      strong = TRUE,
                      sparse = FALSE,
                      impl = "cpp") {

  # initialize lambda
  if (is.null(lambda)) {
    lambda_max <- max(gamma)
    lambda <- lambda_max * exp(seq(from=0,
                                   to=log(lambda.min.ratio),
                                   by=log(lambda.min.ratio)/(nlambda-1)))
  } else {
    nlambda <- length(lambda)
  }

  # initialize R
  if (is.null(R)) {
    R <- funcR(Gamma)
  }

  # initialize warm
  init.beta <- matrix(rep(0, length=length(gamma)*length(lambda)),
                      nrow=length(gamma), ncol=length(lambda))

  if (impl=="CPP" | impl=="cpp") {
    beta <- covCdaC(Gamma, as.numeric(gamma), lambda, R, init.beta, 0,
                    maxit, eps, warm, strong)
  } else {
    stop("Not implemented")
  }

  # set result
  res <- list()
  res$beta_standard <- beta
  res$lambda <- lambda

  return(res)
}
