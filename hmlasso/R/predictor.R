#' Predict responses
#'
#' @param object hmlasso model
#' @param newx matrix of explanatory variables
#' @param s selected lambda (default: all)
#' @param impute_method imputation method for predictions (default: "mean")
#' @param adjust_by_tr whether mean (or median) of training data for prediction is used or not
#' @param ... parameters of predict function
#'
#' @examples
#' X_incompl <- as.matrix(iris[, 1:3])
#' X_incompl[1:5,1] <- NA
#' X_incompl[6:10,2] <- NA
#' y <- iris[, 4]
#' cv_fit <- cv.hmlasso(X_incompl, y, nlambda=50, lambda.min.ratio=1e-2)
#' predict(cv_fit$fit, X_incompl)
#'
#' @export
predict.hmlasso <- function(object, newx, s=NULL, impute_method="mean", adjust_by_tr=FALSE, ...) {

  if (object$family=="gaussian") {
    if (any(is.na(newx))) {
      if (is.null(s)) {
        active_var <- names(which(apply(object$beta, 1, function(v){any(v!=0)})))
        if (impute_method == "mean" | impute_method == "median") {
          if (adjust_by_tr) {
            compl <- impute_data(newx[, active_var], method=NULL, values=object$X_mean)
          } else {
            compl <- impute_data(newx[, active_var], method=impute_method, values=NULL)
          }
        } else {
          stop("Specified imputation method is not supported.")
        }
        # newx[, active_var] <- compl
        for (var in active_var) {
          newx[, var] <- compl[, var]
        }
        newx[is.na(newx)] <- 0
        t(apply(newx %*% object$beta, 1, function(v){v + as.numeric(object$a0)}))
      } else {
        active_var <- names(which(object$beta[, object$lambda==s]!=0))
        if (impute_method == "mean" | impute_method == "median") {
          if (adjust_by_tr) {
            compl <- impute_data(newx[, active_var], values=object$X_mean)
          } else {
            compl <- impute_data(newx[, active_var], method=impute_method)
          }
        } else {
          stop("Specified imputation method is not supported.")
        }
        # newx[, active_var] <- compl
        for (var in active_var) {
          newx[, var] <- compl[, var]
        }
        newx[is.na(newx)] <- 0
        newx %*% object$beta[, object$lambda==s, drop=FALSE] + object$a0[object$lambda==s]
      }
    } else {
      if (is.null(s)) {
        t(apply(newx %*% object$beta, 1, function(v){v + as.numeric(object$a0)}))
      } else {
        newx %*% object$beta[, object$lambda==s, drop=FALSE] + object$a0[object$lambda==s]
      }
    }

  } else {
    stop("specified family is not supported")
  }
}

impute_data <- function(data, method="mean", values=NULL) {
  if (is.null(values)) {
    # cat("impute using ts\n")
    imputed_data <- apply(data, 2, function(vec){impute_vec(vec, method=method)})
  } else {
    # cat("impute using tr\n")
    imputed_data <- sapply(1:ncol(data), function(j){
      as.matrix(impute_vec(data[, j], value=values[j]),
                nrow=nrow(data), ncol=1)
    })
  }
  colnames(imputed_data) <- colnames(data)
  rownames(imputed_data) <- rownames(data)
  return(imputed_data)
}

impute_vec <- function(vec, method="mean", value=NULL) {
  v <- vec
  if (!is.null(value)) {
    v[is.na(vec)] <- value
  } else if (method=="mean") {
    v[is.na(vec)] <- mean(vec, na.rm=TRUE)
  } else if (method=="median") {
    v[is.na(vec)] <- median(vec, na.rm=TRUE)
  } else {
    stop(paste0(method, " is not implemented"))
  }
  return(v)
}
