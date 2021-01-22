#' soft thresholding function
#' 
#' @param z z
#' @param g gamma
#' 
#' @return value
#' 
#' @keywords internal
soft_threshold <- function(z, g) {
  if (g < abs(z)) {
    if(z > 0) {
      z - g
    } else {
      z + g
    }
  } else {
    0
  }
}

#' update rule function
#' 
#' @param u u
#' @param l1 l1
#' @param l2 l2
#' @param v v
#' 
#' @return value
#' 
#' @keywords internal
update_lasso <- function(u, l1, l2, v) {
  if (abs(u) <= l1) {
    return(0)
  } else {
    return(sign(u) * (abs(u) - l1) / (v * (1 + l2)))
  }
}

#' Optimize a linear regression model by coordinate descent algorithm using a covariance matrix with R
#' 
#' @param Gamma covariance matrix of explanatory variables
#' @param gamma covariance vector of explanatory and objective variables
#' @param lambda lambda sequence
#' @param warm warm start direction: "lambda" (default) or "delta"
#' @param delta ratio of regularization between l1 and exclusive penalty terms
#' @param R matrix using exclusive penalty term
#' @param maxit max iteration
#' @param eps convergence threshold for optimization
#' @param init.beta initial values of beta
#' @param strong whether use strong screening or not
#' @param sparse whether use sparse matrix or not
#' 
#' @return standardized beta
#' 
#' @keywords internal
cov_cda_r <- function(Gamma, gamma, lambda, R, init.beta, delta, 
                      maxit, eps, warm, strong, sparse) {
  
  # initialize
  p <- nrow(Gamma)
  nlambda <- length(lambda)
  if (sparse) {
    beta <- sparseMatrix(1, 1, x=0, dims=c(p, nlambda))
    beta_prev <- sparseMatrix(1, 1, x=0, dims=c(p, 1))
  } else {
    beta <- matrix(rep(0, length=p*nlambda), nrow=p, ncol=nlambda)
    beta_prev <- matrix(rep(0, length=p), nrow=p, ncol=1)
  }
  r <- gamma
  
  # coordinate descent algorithm
  for (k in 1:nlambda) {
    cat(paste0("lambda: ", k, "\n"))
    if (warm == "delta") {
      # beta-initialize
      beta[, k] <- init.beta[, k]
    } else if (warm == "lambda") {
      # warm-start
      beta[, k] <- beta_prev
    }
    
    # strong screening
    if (strong & (k!=1)) {
      candid <- abs(gamma - Gamma %*% beta[, k-1]) >= 2 * lambda[k] - lambda[k-1]
    } else {
      candid <- rep(TRUE, length=p)
    }
    
    if (delta==0) { # normal lasso

      # cda for active set
      for (i in 1:maxit) {
        candid2 <- (beta[, k]!=0) & candid

        for (j in which(candid2)) {
          z <- gamma[j, ] - Gamma[j, -j] %*% beta[-j, k]
          g <- lambda[k]
          beta[j, k] <- soft_threshold(z, g)
        }
        residual <- max(abs(beta[, k] - beta_prev))
        beta_prev <- beta[, k]
        if (residual < eps) {
          cat(paste0("1 iteration: ", i, "\n"))
          break
        }
      }

      # cda for strong set
      for (i in 1:maxit) {
        for (j in which(candid)) {
          z <- gamma[j, ] - Gamma[j, -j] %*% beta[-j, k]
          g <- lambda[k]
          beta[j, k] <- soft_threshold(z, g)
        }
        residual <- max(abs(beta[, k] - beta_prev))
        beta_prev <- beta[, k]
        if (residual < eps) {
          cat(paste0("2 iteration: ", i, "\n"))
          break
        }
      }
      
      # cda for all variables
      for (i in 1:maxit) {
        for (j in 1:p) {
          z <- gamma[j, ] - Gamma[j, -j] %*% beta[-j, k]
          # z <- r[j, ]
          g <- lambda[k]
          beta[j, k] <- soft_threshold(z, g)
          # shift <- beta[j, k] - beta_prev[j]
          # if (shift != 0) {
          #   r[-j] <- r[-j] - shift * Gamma[-j, j, drop=FALSE]
          # }
        }
        residual <- max(abs(beta[, k] - beta_prev))
        beta_prev <- beta[, k]
        if (residual < eps) {
          cat(paste0("3 iteration: ", i, "\n"))
          break
        }
      }
    } else { # exclusive lasso
      # cda for active set
      for (i in 1:maxit) {
        candid2 <- (beta[, k]!=0) & candid
        for (j in which(candid2)) {
          z <- gamma[j, ] - Gamma[j, -j] %*% beta[-j, k]
          g <- (1 + delta * R[j, -j] %*% abs(beta[-j, k])) * lambda[k]
          beta[j, k] <- soft_threshold(z, g) / (Gamma[j, j] + delta * lambda[k] * R[j, j])
        }
        residual <- max(abs(beta[, k] - beta_prev))
        beta_prev <- beta[, k]
        if (residual < eps) {
          break
        }
        if (sum(is.infinite(beta_prev) | is.na(beta_prev)) > 0) {
          break
        }
      }
      
      # cda for strong set
      for (i in 1:maxit) {
        for (j in which(candid)) {
          z <- gamma[j, ] - Gamma[j, -j] %*% beta[-j, k]
          g <- (1 + delta * R[j, -j] %*% abs(beta[-j, k])) * lambda[k]
          beta[j, k] <- soft_threshold(z, g) / (Gamma[j, j] + delta * lambda[k] * R[j, j])
        }
        residual <- max(abs(beta[, k] - beta_prev))
        beta_prev <- beta[, k]
        if (residual < eps) {
          break
        }
        if (sum(is.infinite(beta_prev) | is.na(beta_prev)) > 0) {
          break
        }
      }
      
      # cda for all variables
      for (i in 1:maxit) {
        for (j in 1:p) {
          z <- gamma[j, ] - Gamma[j, -j] %*% beta[-j, k]
          g <- (1 + delta * R[j, -j] %*% abs(beta[-j, k])) * lambda[k]
          beta[j, k] <- soft_threshold(z, g) / (Gamma[j, j] + delta * lambda[k] * R[j, j])
        }
        residual <- max(abs(beta[, k] - beta_prev))
        beta_prev <- beta[, k]
        if (residual < eps) {
          break
        }
        if (sum(is.infinite(beta_prev) | is.na(beta_prev)) > 0) {
          break
        }
      }
    }
    
    # break if beta is inf
    if (sum(is.infinite(beta_prev) | is.na(beta_prev)) > 0) {
      break
    }
  }
  
  return(beta)
  
}

#' Optimize a logistic regression model by coordinate descent algorithm using a design matrix with R
#' 
#' @param X_tilde standardized matrix of explanatory variables
#' @param y vector of objective variable
#' @param lambda lambda sequence
#' @param warm warm start direction: "lambda" (default) or "delta"
#' @param delta ratio of regularization between l1 and exclusive penalty terms
#' @param R matrix using exclusive penalty term
#' @param maxit max iteration
#' @param eps convergence threshold for optimization
#' @param init.beta initial values of beta
#' @param strong whether use strong screening or not
#' @param sparse whether use sparse matrix or not
#' 
#' @return standardized beta
#' 
#' @keywords internal
logit_cda_r <- function(X_tilde, y, lambda, R, init.beta, delta, 
                        maxit, eps, warm, strong, sparse) {

  # add a column of 1s for intercept
  n <- nrow(X_tilde)
  X_tilde <- cbind(rep(1, length=n), X_tilde)
  colnames(X_tilde) <- c("(Intercept)", colnames(X_tilde)[2:ncol(X_tilde)])

  # initialize
  p <- ncol(X_tilde)
  nlambda <- length(lambda)
  if (sparse) {
    beta <- sparseMatrix(1, 1, x=0, dims=c(p, nlambda))
    # beta[1, ] <- log(mean(y) / (1 - mean(y)))
    rownames(beta) <- colnames(X_tilde)
    beta_tilde <- sparseMatrix(1, 1, x=0, dims=c(p, 1))
    R <- rbind(matrix(0, nrow=1, ncol=p),
               cbind(matrix(0, nrow=(p-1), ncol=1), R))
  } else {
    beta <- matrix(0, nrow=p, ncol=nlambda)
    # beta[1, ] <- log(mean(y) / (1 - mean(y)))
    rownames(beta) <- colnames(X_tilde)
    beta_tilde <- matrix(0, nrow=p, ncol=1)
    R <- rbind(matrix(0, nrow=1, ncol=p),
               cbind(matrix(0, nrow=(p-1), ncol=1), R))
  }
  
  candid <- c(TRUE, rep(FALSE, length=(p-1)))
  
  # initialization
  eta <- rep(0, length=nrow(X_tilde))

  # coordinate descent algorithm
  for (k in 1:nlambda) {
    cat(paste0("lambda: ", k, "\n"))
    
    if (warm == "delta") {
      # beta-initialize
      beta[, k] <- init.beta[, k]
    } else if (warm == "lambda") {
      # warm-start
      if (k != 1) {
        beta[, k] <- beta[, (k - 1)]
      }
      # beta[, k] <- beta_tilde
    }
    
    iter <- 0
    while(iter < maxit) {
      while(iter < maxit) {
        iter <- iter + 1
        
        # quadratic approximation for logistic log-likelihood
        # cat(paste0("approx: ", t, "\n"))
        beta_tilde <- beta[, k]
        eta <- X_tilde %*% beta_tilde
        pi <- sapply(eta, function(e) {
          if (e > 10) {
            return(1)
          } else if (e < -10) {
            return(0)
          } else {
            return(1 / (1 + exp(-e)))
          }
        })
        w <- pi * (1 - pi)
        w <- ifelse(w > 1e-4,  w, 1e-4)
        s <- y - pi
        r = s / w
        
        max_change <- 0
        
        for (j in 1:p) {
          if (candid[j]) {
            # xwr <- t(X_tilde[, j, drop=FALSE]) %*% diag(as.numeric(w)) %*% r
            xwr <- sum(X_tilde[, j] * as.numeric(w) * r)
            # xwx <- t(X_tilde[, j, drop=FALSE]) %*% diag(as.numeric(w)) %*% X_tilde[, j, drop=FALSE]
            xwx <- sum(as.numeric(w) * X_tilde[, j]^2)
            
            if (delta == 0) { # normal lasso
              u <- as.numeric(xwr / n + (xwx/n) * beta_tilde[j])
              # u <- as.numeric(xwr / n + (xwx/n) * beta[j, k])
              v <- as.numeric(xwx / n)
              l1 <- lambda[k]
              l2 <- 0
            } else { # exclusive lasso
              u <- as.numeric(xwr / n + (xwx/n) * beta_tilde[j])
              # u <- as.numeric(xwr / n + (xwx/n) * beta[j, k])
              v <- as.numeric(xwx / n) + delta * lambda[k] * R[j, j]
              l1 <- lambda[k] + delta * lambda[k] * (R[j, -j] %*% abs(beta[-j, k]))
              l2 <- 0
            }
            
            if (j==1) { # intercept without penalty
              beta[j, k] <- update_lasso(u, 0, l2, 1)
            } else {
              beta[j, k] <- update_lasso(u, l1, l2, v)
            }
            
            shift <- beta[j, k] - beta_tilde[j]
            if (shift != 0) {
              si = shift * X_tilde[, j]
              r = r - si
              eta = eta + si
            }
            
            if (max_change < abs(shift) * sqrt(v)) {
              max_change <- abs(shift) * sqrt(v)
            }
          }
        }
        
        # check for convergence
        # max_change <- max(abs(beta[, k] - beta_tilde) * sqrt(v))
        beta_tilde <- beta[, k]
        if (max_change < eps) {
          cat(paste0("inner iteration: ", iter, "\n"))
          break
        }
      }
      
      # check for violation
      violations <- 0
      for (j in 1:p) {
        if (!candid[j]) {
          z <- sum(X_tilde[, j] * as.numeric(s)) / n
          l1 <- lambda[k]
          if (abs(z) > l1) {
            candid[j] <- TRUE
            violations <- violations + 1
          }
        }
      }
      if (violations==0) {
        cat(paste0("outer iteration: ", iter, "\n"))
        break
      }
    }

    # break if beta is inf
    if (sum(is.infinite(beta_tilde) | is.na(beta_tilde)) > 0) {
      return(beta)
    }
    # break if maxit reached
    if (iter == maxit) {
      return(beta)
    }
  }
  
  return(beta)
}


#' (Experimental) Optimize a ULasso linear regression model by coordinate descent algorithm using a covariance matrix with R
#' 
#' @param Gamma covariance matrix of explanatory variables
#' @param gamma covariance vector of explanatory and objective variables
#' @param lambda lambda sequence
#' @param warm warm start direction: "lambda" (default) or "delta"
#' @param delta ratio of regularization between l1 and exclusive penalty terms
#' @param R matrix using exclusive penalty term
#' @param maxit max iteration
#' @param eps convergence threshold for optimization
#' @param init.beta initial values of beta
#' @param strong whether use strong screening or not
#' @param sparse whether use sparse matrix or not
#' 
#' @return standardized beta
#' 
#' @keywords internal
cov_cda_r2 <- function(Gamma, gamma, lambda, R, init.beta, delta, 
                      maxit, eps, warm, strong, sparse) {
  
  # initialize
  p <- nrow(Gamma)
  nlambda <- length(lambda)
  if (sparse) {
    beta <- sparseMatrix(1, 1, x=0, dims=c(p, nlambda))
    beta_prev <- sparseMatrix(1, 1, x=0, dims=c(p, 1))
  } else {
    beta <- matrix(rep(0, length=p*nlambda), nrow=p, ncol=nlambda)
    beta_prev <- matrix(rep(0, length=p), nrow=p, ncol=1)
  }
  r <- gamma
  
  # coordinate descent algorithm
  for (k in 1:nlambda) {
    cat(paste0("lambda: ", k, "\n"))
    if (warm == "delta") {
      # beta-initialize
      beta[, k] <- init.beta[, k]
    } else if (warm == "lambda") {
      # warm-start
      beta[, k] <- beta_prev
    }
    
    # strong screening
    if (strong & (k!=1)) {
      candid <- abs(gamma - Gamma %*% beta[, k-1]) >= 2 * lambda[k] - lambda[k-1]
    } else {
      candid <- rep(TRUE, length=p)
    }
    
    if (delta==0) { # normal lasso
      
      # cda for active set
      for (i in 1:maxit) {
        candid2 <- (beta[, k]!=0) & candid
        
        for (j in which(candid2)) {
          z <- gamma[j, ] - Gamma[j, -j] %*% beta[-j, k]
          g <- lambda[k]
          beta[j, k] <- soft_threshold(z, g)
        }
        residual <- max(abs(beta[, k] - beta_prev))
        beta_prev <- beta[, k]
        if (residual < eps) {
          cat(paste0("1 iteration: ", i, "\n"))
          break
        }
      }
      
      # cda for strong set
      for (i in 1:maxit) {
        for (j in which(candid)) {
          z <- gamma[j, ] - Gamma[j, -j] %*% beta[-j, k]
          g <- lambda[k]
          beta[j, k] <- soft_threshold(z, g)
        }
        residual <- max(abs(beta[, k] - beta_prev))
        beta_prev <- beta[, k]
        if (residual < eps) {
          cat(paste0("2 iteration: ", i, "\n"))
          break
        }
      }
      
      # cda for all variables
      for (i in 1:maxit) {
        for (j in 1:p) {
          z <- gamma[j, ] - Gamma[j, -j] %*% beta[-j, k]
          # z <- r[j, ]
          g <- lambda[k]
          beta[j, k] <- soft_threshold(z, g)
          # shift <- beta[j, k] - beta_prev[j]
          # if (shift != 0) {
          #   r[-j] <- r[-j] - shift * Gamma[-j, j, drop=FALSE]
          # }
        }
        residual <- max(abs(beta[, k] - beta_prev))
        beta_prev <- beta[, k]
        if (residual < eps) {
          cat(paste0("3 iteration: ", i, "\n"))
          break
        }
      }
    } else { # exclusive lasso
      # cda for active set
      for (i in 1:maxit) {
        candid2 <- (beta[, k]!=0) & candid
        for (j in which(candid2)) {
          z <- gamma[j, ] - Gamma[j, -j] %*% beta[-j, k] - delta * R[j, -j] %*% beta[-j, k] * lambda[k]
          g <- lambda[k]
          beta[j, k] <- soft_threshold(z, g) / (Gamma[j, j] + delta * lambda[k] * R[j, j])
        }
        residual <- max(abs(beta[, k] - beta_prev))
        beta_prev <- beta[, k]
        if (residual < eps) {
          break
        }
        if (sum(is.infinite(beta_prev) | is.na(beta_prev)) > 0) {
          break
        }
      }
      
      # cda for strong set
      for (i in 1:maxit) {
        for (j in which(candid)) {
          z <- gamma[j, ] - Gamma[j, -j] %*% beta[-j, k] - delta * R[j, -j] %*% beta[-j, k] * lambda[k]
          g <- lambda[k]
          beta[j, k] <- soft_threshold(z, g) / (Gamma[j, j] + delta * lambda[k] * R[j, j])
        }
        residual <- max(abs(beta[, k] - beta_prev))
        beta_prev <- beta[, k]
        if (residual < eps) {
          break
        }
        if (sum(is.infinite(beta_prev) | is.na(beta_prev)) > 0) {
          break
        }
      }
      
      # cda for all variables
      for (i in 1:maxit) {
        for (j in 1:p) {
          z <- gamma[j, ] - Gamma[j, -j] %*% beta[-j, k] - delta * R[j, -j] %*% beta[-j, k] * lambda[k]
          g <- lambda[k]
          beta[j, k] <- soft_threshold(z, g) / (Gamma[j, j] + delta * lambda[k] * R[j, j])
        }
        residual <- max(abs(beta[, k] - beta_prev))
        beta_prev <- beta[, k]
        if (residual < eps) {
          break
        }
        if (sum(is.infinite(beta_prev) | is.na(beta_prev)) > 0) {
          break
        }
      }
    }
    
    # break if beta is inf
    if (sum(is.infinite(beta_prev) | is.na(beta_prev)) > 0) {
      break
    }
  }
  
  return(beta)
  
}


#' (Experimental) Optimize a ULasso logistic regression model by coordinate descent algorithm using a design matrix with R
#' 
#' @param X_tilde standardized matrix of explanatory variables
#' @param y vector of objective variable
#' @param lambda lambda sequence
#' @param warm warm start direction: "lambda" (default) or "delta"
#' @param delta ratio of regularization between l1 and exclusive penalty terms
#' @param R matrix using exclusive penalty term
#' @param maxit max iteration
#' @param eps convergence threshold for optimization
#' @param init.beta initial values of beta
#' @param strong whether use strong screening or not
#' @param sparse whether use sparse matrix or not
#' 
#' @return standardized beta
#' 
#' @keywords internal
logit_cda_r2 <- function(X_tilde, y, lambda, R, init.beta, delta, 
                        maxit, eps, warm, strong, sparse) {
  
  # add a column of 1s for intercept
  n <- nrow(X_tilde)
  X_tilde <- cbind(rep(1, length=n), X_tilde)
  colnames(X_tilde) <- c("(Intercept)", colnames(X_tilde)[2:ncol(X_tilde)])
  
  # initialize
  p <- ncol(X_tilde)
  nlambda <- length(lambda)
  if (sparse) {
    beta <- sparseMatrix(1, 1, x=0, dims=c(p, nlambda))
    # beta[1, ] <- log(mean(y) / (1 - mean(y)))
    rownames(beta) <- colnames(X_tilde)
    beta_tilde <- sparseMatrix(1, 1, x=0, dims=c(p, 1))
    R <- rbind(matrix(0, nrow=1, ncol=p),
               cbind(matrix(0, nrow=(p-1), ncol=1), R))
  } else {
    beta <- matrix(0, nrow=p, ncol=nlambda)
    # beta[1, ] <- log(mean(y) / (1 - mean(y)))
    rownames(beta) <- colnames(X_tilde)
    beta_tilde <- matrix(0, nrow=p, ncol=1)
    R <- rbind(matrix(0, nrow=1, ncol=p),
               cbind(matrix(0, nrow=(p-1), ncol=1), R))
  }
  
  candid <- c(TRUE, rep(FALSE, length=(p-1)))
  
  # initialization
  eta <- rep(0, length=nrow(X_tilde))
  
  # coordinate descent algorithm
  for (k in 1:nlambda) {
    cat(paste0("lambda: ", k, "\n"))
    
    if (warm == "delta") {
      # beta-initialize
      beta[, k] <- init.beta[, k]
    } else if (warm == "lambda") {
      # warm-start
      if (k != 1) {
        beta[, k] <- beta[, (k - 1)]
      }
      # beta[, k] <- beta_tilde
    }
    
    iter <- 0
    while(iter < maxit) {
      while(iter < maxit) {
        iter <- iter + 1
        
        # quadratic approximation for logistic log-likelihood
        # cat(paste0("approx: ", t, "\n"))
        beta_tilde <- beta[, k]
        eta <- X_tilde %*% beta_tilde
        pi <- sapply(eta, function(e) {
          if (e > 10) {
            return(1)
          } else if (e < -10) {
            return(0)
          } else {
            return(1 / (1 + exp(-e)))
          }
        })
        w <- pi * (1 - pi)
        w <- ifelse(w > 1e-4,  w, 1e-4)
        s <- y - pi
        r = s / w
        
        max_change <- 0
        
        for (j in 1:p) {
          if (candid[j]) {
            # xwr <- t(X_tilde[, j, drop=FALSE]) %*% diag(as.numeric(w)) %*% r
            xwr <- sum(X_tilde[, j] * as.numeric(w) * r)
            # xwx <- t(X_tilde[, j, drop=FALSE]) %*% diag(as.numeric(w)) %*% X_tilde[, j, drop=FALSE]
            xwx <- sum(as.numeric(w) * X_tilde[, j]^2)
            
            if (delta == 0) { # normal lasso
              u <- as.numeric(xwr / n + (xwx/n) * beta_tilde[j])
              # u <- as.numeric(xwr / n + (xwx/n) * beta[j, k])
              v <- as.numeric(xwx / n)
              l1 <- lambda[k]
              l2 <- 0
            } else { # exclusive lasso
              u <- as.numeric(xwr / n + (xwx/n) * beta_tilde[j] - delta * lambda[k] * (R[j, -j] %*% beta[-j, k]))
              # u <- as.numeric(xwr / n + (xwx/n) * beta[j, k])
              v <- as.numeric(xwx / n) + delta * lambda[k] * R[j, j]
              l1 <- lambda[k]
              l2 <- 0
            }
            
            if (j==1) { # intercept without penalty
              beta[j, k] <- update_lasso(u, 0, l2, 1)
            } else {
              beta[j, k] <- update_lasso(u, l1, l2, v)
            }
            
            shift <- beta[j, k] - beta_tilde[j]
            if (shift != 0) {
              si = shift * X_tilde[, j]
              r = r - si
              eta = eta + si
            }
            
            if (max_change < abs(shift) * sqrt(v)) {
              max_change <- abs(shift) * sqrt(v)
            }
          }
        }
        
        # check for convergence
        # max_change <- max(abs(beta[, k] - beta_tilde) * sqrt(v))
        beta_tilde <- beta[, k]
        if (max_change < eps) {
          cat(paste0("inner iteration: ", iter, "\n"))
          break
        }
      }
      
      # check for violation
      violations <- 0
      for (j in 1:p) {
        if (!candid[j]) {
          z <- sum(X_tilde[, j] * as.numeric(s)) / n
          l1 <- lambda[k]
          if (abs(z) > l1) {
            candid[j] <- TRUE
            violations <- violations + 1
          }
        }
      }
      if (violations==0) {
        cat(paste0("outer iteration: ", iter, "\n"))
        break
      }
    }
    
    # break if beta is inf
    if (sum(is.infinite(beta_tilde) | is.na(beta_tilde)) > 0) {
      return(beta)
    }
    # break if maxit reached
    if (iter == maxit) {
      return(beta)
    }
  }
  
  return(beta)
}
