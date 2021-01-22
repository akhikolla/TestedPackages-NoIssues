#' Predict responses
#' 
#' @param fit IILasso model
#' @param newx matrix of explanatory variables
#' @param s selected lambda (default: all)
#' @param type prediction type for logistic lasso: "response" (default) or "class"
#' 
#' @return prediction matrix (if s is NULL) or vector (if s is specified)
#' 
#' @export
#' 
#' @examples 
#' X <- matrix(c(1,2,3,5,4,7,6,8,9,10), nrow=5, ncol=2)
#' b <- matrix(c(-1,1), nrow=2, ncol=1)
#' e <- matrix(c(0,-0.1,0.1,-0.1,0.1), nrow=5, ncol=1)
#' y <- as.numeric(X %*% b + e)
#' fit <- lasso(X, y)
#' pr <- predict_lasso(fit, X)
#' plot_lasso(fit)
predict_lasso <- function(fit, newx, s=NULL, type="response") {
  
  # colnames
  varnames <- colnames(newx)
  samplenames <- rownames(newx)
  
  if (fit$family=="gaussian") {
    # if (sum(is.na(X)) > 0) {
    #   return(missing_predict_lasso(fit, newx, s))
    # }
    if (is.null(s)) {
      pr <- t(apply(newx %*% fit$beta, 1, function(v){v + as.numeric(fit$a0)}))
      rownames(pr) <- samplenames
      colnames(pr) <- fit$lambda
    } else {
      pr <- as.numeric(newx %*% fit$beta[, fit$lambda==s, drop=FALSE] + fit$a0[fit$lambda==s])
      names(pr) <- samplenames
    }
    
  } else if (fit$family=="binomial") {
    if (type=="response") {
      if (is.null(s)) {
        pr <- 1 / (1 + exp(- t(apply(newx %*% fit$beta, 1, function(v){v + as.numeric(fit$a0)}))))
        rownames(pr) <- samplenames
        colnames(pr) <- fit$lambda
      } else {
        pr <- as.numeric(1 / (1 + exp(- (newx %*% fit$beta[, fit$lambda==s, drop=FALSE] + fit$a0[fit$lambda==s]))))
        names(pr) <- samplenames
      }
    } else if (type=="class") {
      if (is.null(s)) {
        res <- 1 / (1 + exp(- t(apply(newx %*% fit$beta, 1, function(v){v + as.numeric(fit$a0)}))))
        pr <- apply(res, 2, function(v){as.numeric(v > 0.5)})
        rownames(pr) <- samplenames
        colnames(pr) <- fit$lambda
      } else {
        res <- 1 / (1 + exp(- (newx %*% fit$beta[, fit$lambda==s, drop=FALSE] + fit$a0[fit$lambda==s])))
        pr <- as.numeric(res > 0.5)
        names(pr) <- samplenames
      }
    }
    
  } else {
    stop("specified family is not supported")
  }
  
  return(pr)
}

