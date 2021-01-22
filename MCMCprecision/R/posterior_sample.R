# use base::eigen to sample from the posterior for the stationary distribution
posterior.sample <- function (i, tab, epsilon = 0, digits = 8){
  M <- ncol(tab)

  # 1.) sample from conjugate posterior with prior: Dirichlet(0,..,0)
  P <- matrix(rgamma(M^2, tab + epsilon, 1),
                    nrow = M, ncol = M)
  sel <- rowSums(P) > 0
  P[sel,] <- P[sel,,drop=FALSE]/rowSums(P[sel,,drop=FALSE])

  # 2.) get estimate for stationary distribution (largest eigenvalue = 1)
  decomp <- eigen(t(P))
  idx <- which.max(Re(decomp$values))
  if (round(decomp$values[idx],digits) != 1){
    return (rep(NA, M))
  } else {
    ev <- Re(decomp$vectors[,idx])
    return (ev/sum(ev))
  }
}
