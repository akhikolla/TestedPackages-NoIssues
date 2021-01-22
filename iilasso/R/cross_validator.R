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
#' @param cl (not yet implemented)
#' @param ... parameters of lasso function
#' 
#' @return lasso model
#' \item{fit}{lasso model with hole data}
#' \item{lambda.min}{lambda with minimum cross validation error}
#' \item{lambda.min.index}{index of lambda.min}
#' \item{lambda.1se}{largest lambda such that error is within 1 standard error of the minimum}
#' \item{lambda.1se.index}{index of lambda.1se}
#' \item{delta}{delta defined above}
#' \item{foldid}{fold id}
#' \item{cve}{cross validation error}
#' \item{cvse}{cross validation standard error}
#' \item{cvup}{cross validation error + standard error}
#' \item{cvlo}{cross validation error - standard error}
#' \item{pe}{prediction error (for family="binomial")}
#' 
#' @export
#' 
#' @examples 
#' X <- matrix(c(1,2,3,5,4,7,6,8,9,10), nrow=5, ncol=2)
#' b <- matrix(c(-1,1), nrow=2, ncol=1)
#' e <- matrix(c(0,-0.1,0.1,-0.1,0.1), nrow=5, ncol=1)
#' y <- as.numeric(X %*% b + e)
#' cv_fit <- cv_lasso(X, y, nfolds=5)
#' fit <- cv_fit$fit
#' pr <- predict_lasso(fit, X, cv_fit$lambda.min)
#' plot_cv_lasso(cv_fit)
cv_lasso <- function(X, y, nfolds = 10, 
                     lambda.min.ratio = 0.0001,
                     nlambda = 100,
                     lambda = NULL,
                     foldid = NULL,
                     unit = "sample",
                     seed,
                     cl, 
                     ...) {
  
  if (!missing(seed)) {
    set.seed(seed)
  }
  
  if (is.null(foldid)) {
    # set.seed(seed)
    foldid <- sample(rep(1:nfolds, floor(nrow(X) / nfolds)), replace=FALSE)
    # foldid <- split(1:nrow(X), folds)
  }
  
  # constract standardized X and y
  n <- nrow(X)
  X_tilde <- apply(X, 2, function(v){(v - mean(v)) / sqrt((n-1) / n * var(v))})
  y_tilde <- y - mean(y)
  # construct covariance matrices
  gamma <- t(X_tilde) %*% y_tilde / n #sample covariance
  if (is.null(lambda)) {
    lambda_max <- max(gamma)
    lambda <- lambda_max * exp(seq(from=0,
                                   to=log(lambda.min.ratio),
                                   by=log(lambda.min.ratio)/(nlambda-1)))
  }
  
  if (missing(cl)) {
    cv_fit <- lapply(1:nfolds, function(i) {
      train_id <- which(foldid != i)
      test_id <- which(foldid == i)
      fit <- lasso(X[train_id, ], y[train_id], lambda=lambda, ...)
      pred <- predict_lasso(fit, X[test_id, ])
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
      list(fit=fit, pred=pred, actual=actual, resid=resid, mse=mse,
           train_id=train_id, test_id=test_id)
    })
  } else {
    stop("not yet implemented")
    # params <- list(...)
    # print(params)
    # print(ls(all=TRUE))
    # clusterExport(cl, varlist = ls(all=TRUE), envir = environment())
    # # clusterExport(cl, varlist = c("lasso", "predict_lasso"))
    # cv_fit <- parLapply(cl, 1:nfolds, function(i) {
    #   train_id <- sort(unlist(foldid[-i]))
    #   test_id <- foldid[[i]]
    #   fit <- do.call(lasso,
    #                  c(list(X=X[train_id, ]), list(y=y[train_id]), list(lambda=lambda), params))
    #   pred <- predict_lasso(fit, X[test_id, ])
    #   actual <- y[test_id]
    #   resid <- actual - pred
    #   if (fit$family == "gaussian") {
    #     mse <- apply(pred, 2, function(v){
    #       mean((v - y[test_id])^2)
    #     })
    #   } else if (fit$family == "binomial") {
    #     pred[pred < 0.00001] <- 0.00001
    #     pred[pred > 0.99999] <- 0.99999
    #     mse <- apply(pred, 2, function(v) {
    #       mean(-2 * log(ifelse(actual==1, pred, 1-pred)))
    #     })
    #   }
    #   list(fit=fit, pred=pred, actual=actual, resid=resid, mse=mse,
    #        train_id=train_id, test_id=test_id)
    # })
  }
  fit <- lasso(X, y, lambda=lambda, ...)
  # pred <- predict_lasso(fit, X)

  # browser()
    
  res <- list()
  res$lambda <- lambda
  pred_all <- Reduce(rbind, lapply(cv_fit, "[[", "pred"))
  actual_all <- Reduce("c", lapply(cv_fit, "[[", "actual"))
  if (fit$family == "gaussian") {
    resid_all <- apply(pred_all, 2, function(v){(actual_all - v)^2})
  } else if (fit$family == "binomial") {
    resid_all <- apply(pred_all, 2, function(v) {
      -2 * log(ifelse(actual_all==1, v, 1 - v))
    })
  }
  if (unit=="sample") {
    res$cve <- apply(resid_all, 2, mean)
    res$cvse <- apply(resid_all, 2, sd) / sqrt(nrow(resid_all))
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
  res$lambda.min.index <- as.numeric(which.min(res$cve))
  res$lambda.min <- lambda[res$lambda.min.index]
  res$lambda.1se.index <- as.numeric(max(Filter(function(i){
    i < res$lambda.min.index
  }, which(res$cve > res$cve[res$lambda.min.index] + res$cvse[res$lambda.min.index]))))
  res$lambda.1se <- res$lambda[res$lambda.1se.index]
  res$foldid <- foldid
  
  return(res)
}


