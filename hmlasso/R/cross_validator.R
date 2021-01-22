#' Fit a model using a design matrix with cross validation
#'
#' @param X matrix of explanatory variables
#' @param y vector of objective variable
#' @param nfolds the number of folds (ignored if foldid is specified)
#' @param foldid vector indicating id of fold for each sample
#' @param lambda.min.ratio ratio of max lambda and min lambda (ignored if lambda is specified)
#' @param nlambda the number of lambda (ignored if lambda is specified)
#' @param lambda lambda sequence
#' @param unit unit for cross validation error: "sample" (default) or "fold"
#' @param seed random seed of cross validation
#' @param min_eig_th minimum eigenvalue
#' @param use method to calculate correlation matrix from missing data (default "pairwise.complete.obs")
#' @param impute_method imputation method for predictions
#' @param direct_prediction either corrected cross validation is used or not
#' @param adjust_by_tr whether mean (or median) of training data for prediction is used or not
#' @param positify method for solving PSD matrix
#' @param weight_power weighting power (default 0 meaning no-weighting)
#' @param eig_tol tol parameter in eigs_sym function
#' @param eig_maxitr maxitr parameter in eigs_sym
#' @param mu augmented Lagrangian parameter
#' @param verbose whether output verbose warnings and messages (default FALSE)
#' @param ... parameters of hmlasso function
#'
#' @return lasso model
#' \item{fit}{lasso model with hole data}
#' \item{lambda.min}{lambda with minimum cross validation error}
#' \item{lambda.min.index}{index of lambda.min}
#' \item{lambda.1se}{largest lambda such that error is within 1 standard error of the minimum}
#' \item{lambda.1se.index}{index of lambda.1se}
#' \item{foldid}{fold id}
#' \item{cve}{cross validation error}
#' \item{cvse}{cross validation standard error}
#' \item{cvup}{cross validation error + standard error}
#' \item{cvlo}{cross validation error - standard error}
#' \item{pe}{prediction error (for family="binomial")}
#'
#' @examples
#' X_incompl <- as.matrix(iris[, 1:3])
#' X_incompl[1:5,1] <- NA
#' X_incompl[6:10,2] <- NA
#' y <- iris[, 4]
#' cv_fit <- cv.hmlasso(X_incompl, y, nlambda=50, lambda.min.ratio=1e-2)
#' plot(cv_fit)
#' plot(cv_fit$fit)
#'
#' @export
cv.hmlasso <- function(X, y, nfolds = 10,
                     lambda.min.ratio = 1e-2,
                     nlambda = 100,
                     lambda = NULL,
                     foldid = NULL,
                     unit = "sample",
                     seed = 0,
                     min_eig_th = 1e-6,
                     use = "pairwise.complete.obs",
                     impute_method = "mean",
                     direct_prediction = FALSE,
                     adjust_by_tr = FALSE,
                     positify="diag",
                     weight_power = 1,
                     mu = 1,
                     eig_tol = 1e-8,
                     eig_maxitr = 1e+8,
                     verbose = FALSE,
                     ...) {

  if (is.null(foldid)) {
    set.seed(seed)
    foldid <- sample(rep(1:nfolds, floor(nrow(X) / nfolds)), replace=FALSE)
    # foldid <- split(1:nrow(X), folds)
  } else {
    nfolds <- length(unique(foldid))
  }

  # constract standardized X and y
  n <- nrow(X)
  X_mean <- apply(X, 2, function(v){mean(v, na.rm=TRUE)})
  X_sd <- apply(X, 2, function(v){
    sqrt((sum(!is.na(v))-1) / sum(!is.na(v)) * var(v, na.rm=TRUE))
  })
  X_tilde <- sweep(sweep(X, 2, X_mean, "-"), 2, X_sd, "/")

  # replace
  # X_tilde <- apply(X, 2, function(v){(v - mean(v, na.rm=TRUE)) / sqrt((n-1) / n * var(v, na.rm=TRUE))})

  y_tilde <- y - mean(y, na.rm=TRUE)
  # construct covariance matrices
  gamma <- apply(X_tilde, 2, function(v){
    fill <- which(!is.na(v) & !is.na(as.vector(y_tilde)))
    g <- sum(v[fill] * y_tilde[fill]) / length(fill)
    if (is.nan(g) | is.na(g)) {g <- 0}
    return(g)
  })
  # gamma <- t(X_tilde) %*% y_tilde / n #sample covariance
  if (is.null(lambda)) {
    lambda_max <- max(gamma)
    lambda <- lambda_max * exp(seq(from=0,
                                   to=log(lambda.min.ratio),
                                   by=log(lambda.min.ratio)/(nlambda-1)))
  }

  # cross validation
  cv_fit <- lapply(1:nfolds, function(i) {
    if (verbose) {
      message(paste0("fold #", i))
    }

    train_id <- which(foldid != i)
    test_id <- which(foldid == i)

    # fit the model
    # fit <- hmlasso(X[train_id, ], y[train_id], lambda=lambda, positify=positify, weight_power=weight_power, ...)
    fit <- hmlasso(X[train_id, ], y[train_id], lambda=lambda, positify=positify, min_eig_th=min_eig_th,
                   weight_power=weight_power, mu=mu,
                   eig_tol=eig_tol, eig_maxitr=eig_maxitr, use=use, verbose=verbose, ...)

    # evaluate the model
    if (direct_prediction) {
      # direct prediction case
      gamma_new <- apply(X[test_id, ], 2, function(v){
        fill <- which(!is.na(v) & !is.na(as.vector(y[test_id])))
        sum(v[fill] * y[test_id][fill]) / length(fill)
      })
      gamma_new[is.na(gamma_new)] <- 0

      # Gamma_new <- (n-1) / n * cov(X_tilde, use="pairwise.complete.obs") #sample covariance
      # Gamma_new <- (n-1) / n * cov(X_tilde[test_id, ], use="pairwise.complete.obs") #sample covariance
      Gamma_new <- cor(X_tilde[test_id, ], use="pairwise.complete.obs") #sample correlation # add 20180209
      Gamma_new <- apply(Gamma_new, c(1,2), function(v) {
        if (is.na(v) | is.infinite(v) | is.nan(v)) { 0 } else { v }
      })
      # Gamma_new[is.na(Gamma_new)] <- 0
      ill_obs <- which((t(!is.na(X_tilde[test_id, ])) %*% (!is.na(X_tilde[test_id, ]))) <= 2)
      Gamma_new[ill_obs] <- 0
      diag(Gamma_new) <- 1

      # min_eig <- eigs_sym(Gamma_new, k=1, which="SA")$value
      min_eig <- eigs_sym(Gamma_new, k=1, which="SA", opts=list(tol=eig_tol, maxitr=eig_maxitr, retvec=FALSE))$value
      if (min_eig <= min_eig_th) {
        if (adjust_by_tr) {
          stop("incomplete implementation")
          # if (fit$min_eig <= fit$min_eig_th) {
          if (positify=="diag") {
            if (fit$min_eig <= min_eig_th) {
              # a <- - fit$min_eig + fit$min_eig_th
              a <- - fit$min_eig + min_eig_th
              Gamma_new <- (Gamma_new + diag(a, nrow=nrow(Gamma_new), ncol=ncol(Gamma_new))) / (1+a)
            }
          } else if (positify=="projection") {
            if (!is.null(fit$eigen_Gamma)) {
              D_add <- diag(sapply(fit$eigen_Gamma$values, function(v){ifelse(v>min_eig_th, 0, -v+min_eig_th)}))
              eigen_Gamma_new <- eigen(Gamma_new)
              D_new <- diag(eigen_Gamma_new$values) + D_add
              Gamma_new <- eigen_Gamma_new$vectors %*% D_new %*% t(eigen_Gamma_new$vectors)
            }
          } else if (positify=="mean") {
            X_tilde_impute <- sapply(1:ncol(X_tilde), function(j){
              v <- X_tilde[test_id, j]
              v[is.na(v)] <- fit$X_mean[j]
              return(v)
            })
            # X_tilde_impute <- apply(X_tilde[test_id, ], 2, function(v){
            #   v2 <- v
            #   v2[is.na(v2)] <- fit$X_mean
            #   return(v2)
            # })
            Gamma_new <- cor(X_tilde_impute) #sample covariance
          }
        } else {
          if (positify=="diag") {
            Gamma_new <- diag_positify(Gamma_new, min_eig_th, min_eig)
          } else if (positify=="projection") {
            eigen_Gamma <- eigen(Gamma_new)
            Gamma_new <- proj_positify(Gamma_new, min_eig_th, eigen_Gamma)
          } else if (positify=="mean") {
            Gamma_new <- mean_positify(X_tilde[test_id, ])
          } else if (positify=="admm_max") {
            Gamma_new <- admm_positify(Gamma_new, X[test_id, ], epsilon=min_eig_th, mu=mu, norm="max", cor_A=FALSE, cor_B=FALSE,
                                       weight_power=weight_power, verbose=verbose)
          } else if (positify=="admm_frob") {
            Gamma_new <- admm_positify(Gamma_new, X[test_id, ], epsilon=min_eig_th, mu=mu, norm="frobenius", cor_A=FALSE, cor_B=FALSE,
                                       weight_power=weight_power, verbose=verbose)
          }
        }
      }

      mse <- var(y[test_id], na.rm=TRUE) -
        2 * t(gamma_new) %*% fit$beta +
        diag(t(fit$beta) %*% Gamma_new %*% fit$beta)
      pred <- NULL
      actual <- y[test_id]
      resid=NULL
    } else {
      # impute prediction case
      pred <- predict.hmlasso(fit, X[test_id, ], impute_method=impute_method, adjust_by_tr=adjust_by_tr)
      actual <- y[test_id]
      resid <- actual - pred
      if (fit$family == "gaussian") {
        mse <- apply(pred, 2, function(v){
          mean((v - y[test_id])^2)
        })
      } else if (fit$family == "binomial") {
        pred[pred < 0.00001] <- 0.00001
        pred[pred > 0.99999] <- 0.99999
        mse <- apply(pred, 2, function(v) {
          mean(-2 * log(ifelse(actual==1, pred, 1-pred)))
        })
      }
    }
    list(fit=fit, pred=pred, actual=actual, resid=resid, mse=mse,
         train_id=train_id, test_id=test_id)
  })

  # fit for the all data
  if (verbose) {
    message("Whole data fitting")
  }
  # fit <- hmlasso(X, y, lambda=lambda, positify=positify, weight_power=weight_power, ...)
  fit <- hmlasso(X, y, lambda=lambda, positify=positify, min_eig_th=min_eig_th,
                    weight_power=weight_power, mu=mu,
                    eig_tol=eig_tol, eig_maxitr=eig_maxitr, use=use, verbose=verbose, ...)
  # pred <- predict.hmlasso(fit, X)

  res <- list()
  res$lambda <- lambda
  pred_all <- Reduce(rbind, lapply(cv_fit, "[[", "pred"))
  actual_all <- Reduce("c", lapply(cv_fit, "[[", "actual"))
  if (fit$family == "gaussian") {
    if (direct_prediction) {
      resid_all <- NULL
    } else {
      resid_all <- apply(pred_all, 2, function(v){(actual_all - v)^2})
    }
  } else if (fit$family == "binomial") {
    resid_all <- apply(pred_all, 2, function(v) {
      -2 * log(ifelse(actual_all==1, v, 1 - v))
    })
  }
  if (unit=="sample") {
    if (direct_prediction) {
      res$cve <- apply(sapply(cv_fit, "[[", "mse"), 1, mean)
      res$cvse <- 0
    } else {
      res$cve <- apply(resid_all, 2, mean)
      res$cvse <- apply(resid_all, 2, sd) / sqrt(nrow(resid_all))
    }
  } else if (unit=="fold") {
    # res$cve <- apply(sapply(cv_fit, "[[", "resid"), 1, mean)
    # res$cvse <- apply(sapply(cv_fit, "[[", "mse"), 1, sd) / sqrt(nfolds)
    res$cve <- apply(sapply(cv_fit, "[[", "mse"), 1, mean)
    res$cvse <- apply(sapply(cv_fit, "[[", "mse"), 1, sd) / sqrt(nfolds)
  } else {
    stop(paste0(unit, "is not defined"))
  }
  res$cvup <- res$cve + res$cvse
  res$cvlo <- res$cve - res$cvse
  if (fit$family=="binomial") {
    res$pe <- apply(pred_all, 2, function(v){mean(as.numeric(v>0.5)!=actual_all)})
  } else {
    res$pe <- NULL
  }

  res$fit <- fit
  res$lambda.min.index <- which.min(res$cve)
  res$lambda.min <- lambda[res$lambda.min.index]
  res$lambda.1se.index <- suppressWarnings(
    max(Filter(function(i){
      i < res$lambda.min.index
    }, which(res$cve > res$cve[res$lambda.min.index] + res$cvse[res$lambda.min.index])))
  ) # -Inf if there are no lambda.1se.index
  # res$lambda.1se.index <- max(Filter(function(i){
  #   i < res$lambda.min.index
  # }, which(res$cve > res$cve[res$lambda.min.index] + res$cvse[res$lambda.min.index])))
  res$lambda.1se <- res$lambda[res$lambda.1se.index]
  res$foldid <- foldid

  class(res) <- "cv.hmlasso"

  return(res)
}


