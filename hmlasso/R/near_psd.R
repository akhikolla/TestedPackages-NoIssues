max_norm <- function(M1, M2=matrix(0, nrow=nrow(M1), ncol=ncol(M1)), W=matrix(1, nrow=nrow(M1), ncol=ncol(M1))) {
  max(abs(W * (M1 - M2)))
}

f_norm <- function(M1, M2=matrix(0, nrow=nrow(M1), ncol=ncol(M1)), W=matrix(1, nrow=nrow(M1), ncol=ncol(M1))) {
  sqrt(sum((W * (M1 - M2))^2))
}

proj <- function(A, epsilon=1e-6, correlation=TRUE, maxit=10) {
  eig <- eigen(A)
  eig_values <- sapply(eig$values, function(v){max(v, epsilon)})
  Anew <- eig$vectors %*% diag(eig_values) %*% t(eig$vectors)
  if (correlation) {
    Anew <- Anew / mean(diag(Anew))
    diag(Anew) <- 1
  }

  return(Anew)
}

vecl <- function(A) {
  A[upper.tri(A, diag=TRUE)]
}

matl <- function(v) {
  n <- round(sqrt(8*length(v)+1)/2 - 1/2)
  A <- matrix(0, nrow=n, ncol=n)
  A[upper.tri(A, diag=TRUE)] <- v
  A <- A + t(A) - diag(diag(A))
  return(A)
}

l1 <- function(x, mu) {
  if (sum(abs(x)) <= mu) {
    return(x)
  }

  u <- abs(x)
  u_order <- order(u, decreasing = TRUE)
  rho <- 0
  for (j in 1:length(u)) {
    z <- u[u_order[j]] - (1/j) * (sum(u[u_order[1:j]]) - mu)
    if (z > 0) {
      rho <- j
      next
    } else {
      break
    }
  }

  theta <- 1/rho * (sum(u[u_order[1:rho]]) - mu)
  w <- sapply(x, function(v){sign(v) * max(abs(v)-theta, 0)}) # modify 2018/07/02

  return(w)
}

linf <- function(x, mu) {
  u <- abs(x)
  u_order <- order(u, decreasing = TRUE)
  rho <- 0
  for (j in 1:length(u)) {
    z <- u[u_order[j]] - (1/j) * (sum(u[u_order[1:j]]) - mu/2)
    if (z > 0) {
      rho <- j
      next
    } else {
      break
    }
  }

  theta <- 1/rho * (sum(u[u_order[1:rho]]) - mu/2)
  y <- sapply(x, function(v){
    if (abs(v) <= theta) {
      return(v)
    } else {
      return(sign(v) * theta)
    }
  })

  return(y)
}

linfw <- function(x, mu, w=rep(1, times=length(x)), eps=1e-6) {
  w[abs(w) < eps] <- eps

  u <- abs(x)
  u_order <- order(w * u, decreasing = TRUE)
  rho <- 0
  for (j in 1:length(u)) {
    z <- w[u_order[j]] * u[u_order[j]] - (sum(u[u_order[1:j]]) - mu/2) / sum(1/w[u_order[1:j]])
    if (z > 0) {
      rho <- j
      next
    } else {
      break
    }
  }

  theta <- (sum(u[u_order[1:rho]]) - mu/2) / sum(1/w[u_order[1:rho]])
  y <- sapply(1:length(x), function(j) {
    if (w[j] * abs(x[j]) <= theta) {
      return(x[j])
    } else {
      return(sign(x[j]) * theta / w[j])
    }
  })

  return(y)
}

# positify by imputing mean
mean_positify <- function(X_tilde) {
  X_tilde_impute <- apply(X_tilde, 2, function(v){
    v2 <- v
    v2[is.na(v2)] <- mean(v, na.rm=TRUE)
    return(v2)
  })
  Gamma <- cor(X_tilde_impute) #sample covariance
  return(Gamma)
}

# positify by adding a constant to diagonals
diag_positify <- function(Gamma, min_eig_th, min_eig = NULL) {
  if (is.null(min_eig)) {
    min_eig <- eigs_sym(Gamma, k=1, which="SA", opts=list(tol=1e-8, maxitr=1e+6, retvec=FALSE))$value
  }
  a <- - min_eig + min_eig_th
  Gamma <- (Gamma + diag(a, nrow=nrow(Gamma), ncol=ncol(Gamma))) / (1+a)
  return(Gamma)
}

# positify by projection on to the space of matrices with eigen values greater than zero
proj_positify <- function(Gamma, min_eig_th, eigen_Gamma = NULL) {
  if (is.null(eigen_Gamma)) {
    eigen_Gamma <- eigen(Gamma)
  }
  D <- diag(eigen_Gamma$values)
  D[D<min_eig_th] <- min_eig_th
  D <- diag(diag(D))
  Gamma <- eigen_Gamma$vectors %*% D %*% t(eigen_Gamma$vectors)
  return(Gamma)
}

# positify by admm for non-weighting formulation
admm_positify <- function(Sigmahat, X=NULL, weight_power = 1, epsilon=1e-6, mu=1, maxit=1e+4, tol=1e-6,
                          norm="max", cor_A=FALSE, cor_B=FALSE,
                          # H=matrix(1, nrow=nrow(Sigmahat), ncol=ncol(Sigmahat)),
                          tau=2, m=10, verbose=FALSE) {

  # initialize
  Ahat <- Sigmahat
  Bhat <- matrix(0, nrow=nrow(Sigmahat), ncol=ncol(Sigmahat))
  Lambda <- matrix(0, nrow=nrow(Sigmahat), ncol=ncol(Sigmahat))
  H <- (t(!is.na(X)) %*% (!is.na(X)))^weight_power
  H <- (H / max(H, na.rm=TRUE))^weight_power
  H <- apply(H, c(1,2), function(v) {
    if (is.na(v) | is.infinite(v) | is.nan(v)) {epsilon} else { if (v<epsilon) {epsilon} else {v} }
  })
  w <- vecl(H) / max(H, na.rm=TRUE)

  if (min(eigen(Sigmahat)$values) > 0) {
    return(Sigmahat)
  }

  if (norm=="max") {
    # iterate
    for (i in 1:maxit) {
      # A step
      Ahat_before <- Ahat
      Ahat <- proj(Bhat + Sigmahat + mu * Lambda, epsilon, cor_A)
      # B step
      Bhat_before <- Bhat
      cvec <- vecl(Ahat - Sigmahat - mu * Lambda)
      Bhat <- matl(linfw(cvec, mu, w))
      if (cor_B) {
        diag(Bhat) <- 0
      }
      # Lambda step
      Lambda_before <- Lambda
      R <- Ahat - Bhat - Sigmahat
      Lambda <- Lambda - 1/mu * R
      # adjust mu
      S <- (Bhat - Bhat_before) / mu
      if (f_norm(R) > m * f_norm(S)) {
        mu <- mu / tau
      } else if (f_norm(S) > m * f_norm(R)) {
        mu <- mu * tau
      }

      # check terminal condition
      change <- max(abs(Ahat_before - Ahat),
                    abs(Bhat_before - Bhat),
                    abs(Lambda_before - Lambda),
                    abs(Ahat - Sigmahat - Bhat))
      if (change < tol) {
        if (verbose) {
          message(paste0("ADMM total iteration: ", i))
          message(paste0("ADMM converges: ", change, "<", tol))
        }
        break
      }
      if (i == maxit) {
        if (verbose) {
          message(paste0("ADMM total iteration: ", i))
          warning(paste0("ADMM does not converge: ", change, ">=", tol))
          # cat(paste0(max(abs(Ahat_before - Ahat)), "\n"))
          # cat(paste0(max(abs(Bhat_before - Bhat)), "\n"))
          # cat(paste0(max(abs(Lambda_before - Lambda)), "\n"))
          # cat(paste0(max(abs(Ahat - Sigmahat - Bhat)), "\n"))
        }
      }
    }
  } else if (norm=="frobenius") {
    # iterate
    for (i in 1:maxit) {
      # A step
      Ahat_before <- Ahat
      Ahat <- proj(Bhat + Sigmahat + mu * Lambda, epsilon, cor_A)
      # B step
      Bhat_before <- Bhat
      Bhat <- (Ahat - Sigmahat - mu * Lambda) / (mu * (H^weight_power)^2 + 1)
      if (cor_B) {
        diag(Bhat) <- 0
      }
      # Lambda step
      Lambda_before <- Lambda
      R <- Ahat - Bhat - Sigmahat
      Lambda <- Lambda - 1/mu * R
      # adjust mu
      S <- (Bhat - Bhat_before) / mu
      if (f_norm(R) > 10 * f_norm(S)) {
        mu <- mu / 2
      } else if (f_norm(S) > 10 * f_norm(R)) {
        mu <- mu * 2
      }

      # check terminal condition
      change <- max(abs(Ahat_before - Ahat),
                    abs(Bhat_before - Bhat),
                    abs(Lambda_before - Lambda),
                    abs(Ahat - Sigmahat - Bhat))
      if (change < tol) {
        if (verbose) {
          message(paste0("ADMM total iteration: ", i))
          message(paste0("ADMM converges: ", change, "<", tol))
        }
        break
      }
      if (i == maxit) {
        if (verbose) {
          message(paste0("ADMM total iteration: ", i))
          warning(paste0("ADMM does not converge: ", change, ">=", tol))
          # cat(paste0(max(abs(Ahat_before - Ahat)), "\n"))
          # cat(paste0(max(abs(Bhat_before - Bhat)), "\n"))
          # cat(paste0(max(abs(Lambda_before - Lambda)), "\n"))
          # cat(paste0(max(abs(Ahat - Sigmahat - Bhat)), "\n"))
        }
      }
    }
  }

  return(Ahat)
}
