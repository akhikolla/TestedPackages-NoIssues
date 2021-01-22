rsvd <- function (M, n, q=2, oversampling=2, transpose=FALSE, verbose=FALSE) {
  ## --- randomized SVD according to Halko, Martinsson & Tropp (2009, p. 9) ---  
  
  ## We can apply the rSVD algorithm either to A = M or to A = t(M), depending on the format of M.
  ## Preliminary testing suggested that the original algorithm (A = M) is suitable for a matrix with many columns,
  ## while the transpose algorithm (A = t(M)) works better if the matrix has many rows and a limited number of columns.
  ## With the current implementation, which uses SVD rather than QR decomposition to obtain an orthonormal basis,
  ## there does not seem to be a substantial difference.
  
  dsm.is.canonical(M) # ensure that M is a suitable matrix (we don't need to enforce canonical format)
  nR <- nrow(M)
  nC <- ncol(M)
  if (n < 1 || n > min(nR, nC)) stop(sprintf("number of singular components out of range n = 1 ... %d", min(nR, nC)))
  
  k2 <- min(oversampling * n, nR, nC)   # = 2*k in the paper
  if (verbose) cat(sprintf("Randomized SVD reduction%s to %d => %d dimensions:\n", if (transpose) " (transposed)" else "", k2, n))

  if (!transpose) {

    if (verbose) cat(" - sampling range of A\n") # -- original algorithm applied to A = M
    Omega <- matrix(rnorm(nC*k2), nC, k2)
    Y <- M %*% Omega
    rm(Omega)
    if (q >= 1) for (i in 1:q) {
      if (verbose) cat(sprintf(" - power iteration #%d\n", i))
      Y <- M %*% crossprod(M, Y)
    }
    if (verbose) cat(sprintf(" - orthonormal basis of %d x %d matrix\n", nrow(Y), ncol(Y)))
    Q <- svd(Y, nu=k2, nv=0)$u  # orthonormal basis of rg(Y); SVD is faster than and as accurate as QR decomposition
    rm(Y)
    B <- crossprod(Q, M)
    if (verbose) cat(sprintf(" - SVD decomposition of %d x %d matrix\n", nrow(B), ncol(B)))
    SVD <- svd(B, nu=n, nv=n)   # SVD of B, truncated to n target dimensions
    rm(B)
    if (verbose) cat(" - composing final result\n")
    return(list(
      u = Q %*% SVD$u, # U = Q * \hat{U}
      v = SVD$v,       # V
      d = SVD$d[1:n])) # diag(Sigma)
    
  } else {

    if (verbose) cat(" - sampling range of A\n") # -- transposed algorithm for A = t(M)
    Omega <- matrix(rnorm(k2*nR), k2, nR) # = t(Omega)
    Y <- Omega %*% M                      # = t(A * Omega)
    rm(Omega)
    if (q >= 1) for (i in 1:q) {
      if (verbose) cat(sprintf(" - power iteration #%d\n", i))
      Y <- tcrossprod(Y, M) %*% M         # = t( (A * t(A))^i * A * Omega) ) = t(Y)
    }
    if (verbose) cat(sprintf(" - orthonormal basis of %d x %d matrix\n", ncol(Y), nrow(Y)))
    Q <- svd(Y, nu=0, nv=k2)$v            # orthonormal basis of rg(Y) = column space of t(Y); Q is _not_ transposed
    rm(Y)
    B <- M %*% Q                          # = t( t(Q) * A ) = t(B)
    if (verbose) cat(sprintf(" - SVD decomposition of %d x %d matrix\n", nrow(B), ncol(B)))
    SVD <- svd(B, nu=n, nv=n)             # t(B) = V * Sigma * t(\hat{U}), truncated to n target dimensions
    rm(B)
    if (verbose) cat(" - composing final result\n")
    ## we now have A = U * Sigma * t(V) with U = Q * \hat{U}, where \hat{U} = SVD$v (_not_ transposed!) and V = SVD$u;
    ## the approximate SVD of M = t(A) is therefore M = V * Sigma * t(U)
    return(list(
      u = SVD$u,       # V
      v = Q %*% SVD$v, # U = Q * \hat{U}
      d = SVD$d[1:n])) # diag(Sigma)

  }
}
