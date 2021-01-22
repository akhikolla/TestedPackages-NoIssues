#' Fit a model using a design matrix
#' 
#' @param X matrix of explanatory variables
#' @param y vector of objective variable
#' @param family family of regression: "gaussian" (default) or "binomial"
#' @param impl implementation language of optimization: "cpp" (default) or "r"
#' @param lambda.min.ratio ratio of max lambda and min lambda (ignored if lambda is specified)
#' @param nlambda the number of lambda (ignored if lambda is specified)
#' @param lambda lambda sequence
#' @param warm warm start direction: "lambda" (default) or "delta"
#' @param ... parameters for optimization
#' 
#' @return lasso model
#' \item{beta}{coefficients}
#' \item{beta_standard}{standardized coefficients}
#' \item{a0}{intercepts}
#' \item{lambda}{regularization parameters}
#' \item{alpha}{alpha defined above}
#' \item{delta}{delta defined above}
#' \item{family}{family}
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
lasso <- function(X, y, family="gaussian", impl="cpp", 
                  lambda.min.ratio = 0.0001,
                  nlambda = 100,
                  lambda = NULL,
                  warm = "lambda", ...) {
  # colnames
  varnames <- colnames(X)

  # stop if X or y include NA
  if (sum(is.na(X)) + sum(is.na(y)) > 0) {
    stop("There are missing values")
  }
  
  # constract standardized X and y
  n <- nrow(X)
  X_mean <- apply(X, 2, mean)
  X_sd <- apply(X, 2, function(v){
    sqrt((n-1) / n * var(v))
  })
  X_tilde <- sweep(sweep(X, 2, X_mean, "-"), 2, X_sd, "/")
  y_tilde <- y - mean(y)
  
  # construct covariance matrices
  Gamma <- (n-1) / n * cov(X_tilde) #sample covariance
  # Gamma <- t(X_tilde) %*% X_tilde / n
  gamma <- t(X_tilde) %*% y_tilde / n #sample covariance
  
  if (is.null(lambda)) {
    lambda_max <- max(gamma)
    lambda <- lambda_max * exp(seq(from=0,
                                   to=log(lambda.min.ratio),
                                   by=log(lambda.min.ratio)/(nlambda-1)))
  }
  
  # fit
  if (family=="gaussian") {
    if (warm == "delta") {
      init.beta <- cov_lasso(Gamma, gamma, delta=0, 
                             lambda=lambda, impl=impl)$beta_standard
      fit <- cov_lasso(Gamma, gamma, lambda=lambda, warm=warm, init.beta=init.beta, impl=impl, ...)
    } else if (warm == "lambda"){
      fit <- cov_lasso(Gamma, gamma, lambda=lambda, warm=warm, impl=impl, ...)
    }
    if (is.null(dim(fit$beta_standard))) {
      fit$beta <- fit$beta_standard / X_sd
      # if (is.null(names(fit$beta))) {
      #   names(fit$beta) <- paste0("V",1:length(fit$beta))
      # }
    } else {
      fit$beta <- sweep(fit$beta_standard, 1, X_sd, FUN="/")
      # if (is.null(rownames(fit$beta))) {
      #   rownames(fit$beta) <- paste0("V",1:nrow(fit$beta))
      # }
    }
    fit$a0 <- mean(y) - t(X_mean) %*% fit$beta
    
  } else if (family=="binomial") {
    if (warm == "delta") {
      init.beta <- logit_lasso(X_tilde, y, delta=0,
                               lambda=lambda, impl=impl)$beta_standard
      fit <- logit_lasso(X_tilde, y, lambda=lambda, warm=warm, init.beta=init.beta, impl=impl, ...)
    } else if (warm == "lambda") {
      fit <- logit_lasso(X_tilde, y, lambda=lambda, warm=warm, impl=impl, ...)
    }
    if (is.null(dim(fit$beta_standard))) {
      fit$beta <- fit$beta_standard[-1, ] / X_sd
      # if (is.null(names(fit$beta))) {
      #   names(fit$beta) <- paste0("V",1:length(fit$beta))
      # }
    } else {
      fit$beta <- sweep(fit$beta_standard[-1, , drop=FALSE], 1, X_sd, FUN="/")
      # if (is.null(rownames(fit$beta))) {
      #   rownames(fit$beta) <- paste0("V",1:nrow(fit$beta))
      # }
    }
    fit$a0 <- fit$beta_standard[1, ] -
      apply(fit$beta_standard[-1, ], 2, function(v){sum(v * X_mean / X_sd)})
    
  } else {
    stop("specified family is not supported")
  }
  
  fit$family <- family
  colnames(fit$beta) <- lambda
  rownames(fit$beta) <- varnames
  names(fit$a0) <- lambda

  return(fit)
}

#' Set up a lambda sequence
#' 
#' @param X matrix of explanatory variables
#' @param y vector of objective variable
#' @param family family of regression: "gaussian" (default) or "binomial"
#' @param lambda.min.ratio ratio of max lambda and min lambda
#' @param nlambda the number of lambda (ignored if lambda is specified)
#' 
#' @return lambda
#' 
#' @export
#' 
#' @examples
#' X <- matrix(c(1,2,3,5,4,7,6,8,9,10), nrow=5, ncol=2)
#' b <- matrix(c(-1,1), nrow=2, ncol=1)
#' e <- matrix(c(0,-0.1,0.1,-0.1,0.1), nrow=5, ncol=1)
#' y <- as.numeric(X %*% b + e)
#' setup_lambda(X, y)
setup_lambda <- function(X, y, family="gaussian", lambda.min.ratio=1e-4, nlambda=100) {
  
  n <- nrow(X)
  p <- ncol(X)
  fit <- glm(as.numeric(y)~1, family=family)

  if (family=="gaussian") {
    lambda_max <- max(abs(t(X) %*% fit$residuals))  / n
  } else {
    lambda_max <- max(abs(t(X) %*% residuals(fit, "working"))) / n
  }
  lambda <- lambda_max * exp(seq(from=0,
                                 to=log(lambda.min.ratio),
                                 by=log(lambda.min.ratio)/(nlambda-1)))
  return(lambda)
}

#' Fit a linear regression model using a covariance matrix
#' 
#' @param Gamma covariance matrix of explanatory variables
#' @param gamma covariance vector of explanatory and objective variables
#' @param impl implementation language of optimization: "cpp" (default) or "r"
#' @param lambda.min.ratio ratio of max lambda and min lambda
#' @param nlambda the number of lambda (ignored if lambda is specified)
#' @param lambda lambda sequence
#' @param warm warm start direction: "lambda" (default) or "delta"
#' @param delta ratio of regularization (exclusive penalty / l1 penalty) (default: 0)
#' @param alpha mixing parameter of regularization of l1 and exclusive penalty terms (delta = (1 - alpha) / alpha)
#' @param R matrix using exclusive penalty term
#' @param funcR function of R (input: X, output: R)
#' @param maxit max iteration (default: 1e+4)
#' @param eps convergence threshold for optimization (default: 1e-4)
#' @param init.beta initial values of beta
#' @param strong whether use strong screening (default) or not
#' @param sparse whether use sparse matrix or not (default)
#' @param abs (experimental) whether use absolute value of beta (default) or not
#' 
#' @return lasso model
#' \item{beta_standard}{standardized coefficients}
#' \item{lambda}{regularization parameters}
#' \item{alpha}{alpha defined above}
#' \item{delta}{delta defined above}
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
cov_lasso <- function(Gamma, gamma,
                      lambda.min.ratio = 0.0001,
                      nlambda = 100,
                      lambda = NULL,
                      delta = 0,
                      alpha = NULL,
                      R = NULL,
                      funcR = function(G){abs(G)^2},
                      maxit = 1e+4,
                      eps = 1e-04,
                      warm = "lambda",
                      init.beta = NULL,
                      strong = TRUE,
                      sparse = FALSE,
                      impl = "cpp",
                      abs = TRUE) {
  
  # initialize lambda
  if (is.null(lambda)) {
    lambda_max <- max(gamma)
    lambda <- lambda_max * exp(seq(from=0,
                                   to=log(lambda.min.ratio),
                                   by=log(lambda.min.ratio)/(nlambda-1)))
  } else {
    nlambda <- length(lambda)
  }
  
  # initialize alpha & delta
  if (!is.null(alpha)) {
    if (alpha <= 0 | alpha>1) {
      stop("0 < alpha <= 1")
    } else {
      delta <- (1 - alpha) / alpha
    }
  }
  
  # initialize R
  if (is.null(R)) {
    R <- funcR(Gamma)
  }
  
  # initialize warm
  if (warm == "lambda") {
    init.beta <- matrix(rep(0, length=length(gamma)*length(lambda)), 
                        nrow=length(gamma), ncol=length(lambda))
  } else if (warm == "delta") {
    if (is.null(init.beta)) {
      stop("need init.beta if warm is delta")
    }
  }
  
  if (impl=="CPP" | impl=="cpp") {
    if (abs) {
      beta <- covCdaC(Gamma, as.numeric(gamma), lambda, R, init.beta, delta, 
                      maxit, eps, warm, strong)
    } else {
      beta <- covCdaC2(Gamma, as.numeric(gamma), lambda, R, init.beta, delta, 
                      maxit, eps, warm, strong)
    }
  } else if (impl=="R" | impl=="r") {
    if (abs) {
      beta <- cov_cda_r(Gamma, gamma, lambda, R, init.beta, delta,
                        maxit, eps, warm, strong, sparse)
    } else {
      beta <- cov_cda_r2(Gamma, gamma, lambda, R, init.beta, delta,
                        maxit, eps, warm, strong, sparse)
    }
  }

  # set result
  res <- list()
  res$beta_standard <- beta
  res$lambda <- lambda
  res$alpha <- 1 / (1 + delta)
  res$delta <- delta
  
  return(res)
}

#' Fit a logistic regression model using a design matrix
#' 
#' @param X_tilde standardized matrix of explanatory variables
#' @param y vector of objective variable
#' @param impl implementation language of optimization: "cpp" (default) or "r"
#' @param lambda.min.ratio ratio of max lambda and min lambda
#' @param nlambda the number of lambda (ignored if lambda is specified)
#' @param lambda lambda sequence
#' @param warm warm start direction: "lambda" (default) or "delta"
#' @param delta ratio of regularization (exclusive penalty / l1 penalty) (default: 0)
#' @param alpha mixing parameter of regularization of l1 and exclusive penalty terms (delta = (1 - alpha) / alpha)
#' @param R matrix using exclusive penalty term
#' @param funcR function of R (input: X, output: R)
#' @param maxit max iteration (default: 1e+4)
#' @param eps convergence threshold for optimization (default: 1e-4)
#' @param init.beta initial values of beta
#' @param strong whether use strong screening (default) or not
#' @param sparse whether use sparse matrix or not (default)
#' @param abs (experimental) whether use absolute value of beta (default) or not
#' 
#' @return lasso model
#' \item{beta_standard}{standardized coefficients}
#' \item{lambda}{regularization parameters}
#' \item{alpha}{alpha defined above}
#' \item{delta}{delta defined above}
#' 
#' @export
#' 
#' @examples
#' X <- matrix(c(1,2,3,5,4,7,6,8,9,10), nrow=5, ncol=2)
#' b <- matrix(c(-1,1), nrow=2, ncol=1)
#' e <- matrix(c(0,-0.1,0.1,-0.1,0.1), nrow=5, ncol=1)
#' y <- as.numeric(X %*% b + e)
#' y <- ifelse(y>mean(y), 1, 0)
#' fit <- lasso(X, y, family="binomial")
#' pr <- predict_lasso(fit, X)
#' plot_lasso(fit)
logit_lasso <- function(X_tilde, y,
                        lambda.min.ratio = 0.0001,
                        nlambda = 100,
                        lambda = NULL,
                        delta = 0,
                        alpha = NULL,
                        R = NULL,
                        funcR = function(G){abs(G)^2},
                        maxit = 1e+4,
                        eps = 1e-04,
                        warm = "lambda",
                        init.beta = NULL,
                        strong = FALSE,
                        sparse = FALSE, 
                        impl="cpp",
                        abs=TRUE) {
  
  # initialize lambda
  if (is.null(lambda)) {
    lambda_max <- max(gamma)
    lambda <- lambda_max * exp(seq(from=0,
                                   to=log(lambda.min.ratio),
                                   by=log(lambda.min.ratio)/(nlambda-1)))
  } else {
    nlambda <- length(lambda)
  }
  
  # initialize alpha & delta
  if (!is.null(alpha)) {
    if (alpha <= 0 | alpha>1) {
      stop("0 < alpha <= 1")
    } else {
      delta <- (1 - alpha) / alpha
    }
  }
  
  # initialize R
  if (is.null(R)) {
    n <- nrow(X_tilde)
    R <- funcR((n-1)/n*cov(X_tilde))
  }
  
  # initialize warm
  if (warm == "lambda") {
    init.beta <- matrix(rep(0, length=length(gamma)*length(lambda)), 
                        nrow=length(gamma), ncol=length(lambda))
  } else if (warm == "delta") {
    if (is.null(init.beta)) {
      stop("need init.beta if warm is delta")
    }
  }
  
  if (impl=="CPP" | impl=="cpp") {
    if (abs) {
      beta <- logitCdaC(X_tilde, as.numeric(y), lambda, R, init.beta, delta, 
                        maxit, eps, warm, strong)
    } else {
      beta <- logitCdaC2(X_tilde, as.numeric(y), lambda, R, init.beta, delta, 
                        maxit, eps, warm, strong)
    }
  } else if (impl=="R" | impl=="r") {
    if (abs) {
      beta <- logit_cda_r(X_tilde, y, lambda, R, init.beta, delta,
                          maxit, eps, warm, strong, sparse)
    } else {
      beta <- logit_cda_r2(X_tilde, y, lambda, R, init.beta, delta,
                          maxit, eps, warm, strong, sparse)
    }
  }
  
  # set result
  res <- list()
  res$beta_standard <- beta
  res$lambda <- lambda
  res$alpha <- 1 / (1 + delta)
  res$delta <- delta
  
  return(res)
}
