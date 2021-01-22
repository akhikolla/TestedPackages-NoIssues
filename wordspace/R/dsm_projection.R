dsm.projection <- function (model, n, method=c("svd", "rsvd", "asvd", "ri", "ri+svd"), oversampling=NA, q=2, rate=.01, power=1, with.basis=FALSE, verbose=FALSE) {
  method <- match.arg(method)
  M <- find.canonical.matrix(model)
  info <- dsm.is.canonical(M) # to find out whether M is sparse or dense
  
  nR <- nrow(M)
  nC <- ncol(M)
  
  if (is.na(n)) n <- min(nR, nC)
  if (method == "ri") {
    if (2*n > nC) stop(sprintf("random indexing from %d to %d dimensions makes no sense", nC, n))
  } else {
    if (n > min(nR, nC)) stop("number of target dimensions exceeds dimensionality of DSM matrix")
  }

  if (is.na(oversampling)) {
    oversampling <- switch(method, svd=1, rsvd=2, asvd=10, ri=1, "ri+svd"=20)
  }
  if (with.basis && method %in% c("ri", "ri+svd")) stop("with.basis=TRUE is not supported for RI-based models")

  B <- NULL
  R2 <- NULL
  sigma <- NULL
  if (method == "svd") {
    ## --- standard SVD algorithm (dense or sparse) ---

    if (verbose) cat(sprintf("SVD reduction to %d dimensions:\n", n))
    if (verbose) cat(" - SVD decomposition\n")
    if (info$sparse) {
      SVD <- sparsesvd(M, rank=n) # use SVDLIBC algorithm for truncated SVD of sparse matrix
      if (ncol(SVD$u) < n) {      # in case poorly conditioned singular components have been removed
        if (verbose) cat(sprintf(" - %d poorly conditioned singular components have been dropped\n", n - ncol(SVD$u)))
        n <- ncol(SVD$u)
      }
    } else {
      SVD <- svd(M, nu=n, nv=(if (with.basis) n else 0)) # we don't need right singular vectors for the dimensionality reduction
    }
    if (verbose) cat(" - composing final matrix\n")
    sigma <- SVD$d[1:n]
    if (power == 1) {
      S <- scaleMargins(SVD$u, cols=sigma) # dimensionality-reduced matrix
    } else {
      S <- scaleMargins(SVD$u, cols=sigma^power)
    }
    if (with.basis) B <- SVD$v
    R2 <- sigma^2 / norm(M, "F")^2 # proportion of "variance" captured by SVD dimensions
    rm(SVD)

  } else if (method == "asvd") {
    ## --- approximated SVD based on random sample of rows (DEPRECATED) ---

    sub.size <- min(n * oversampling, nR)
    if (verbose) cat(sprintf("Approximate SVD reduction to %d dimensions, based on %d rows:\n", n, sub.size))
    sub.idx <- sort(sample(1:nR, sub.size, replace=FALSE))
    M.sub <- M[sub.idx, ]
    if (verbose) cat(" - SVD decomposition\n")
    SVD <- svd(M.sub, nu=0, nv=n) # here we only need the right singular vectors
    if (verbose) cat(" - composing final matrix\n")
    S <- M %*% SVD$v  # V projects columns to first n latent dimensions
    if (with.basis) B <- SVD$v
    rm(SVD)
    S <- as.matrix(S) # make sure result is an ordinary dense matrix (for efficient further processing)
    R2 <- colNorms(M, "euclidean")^2 / norm(M, "F")^2 # this should be the proportion of explained "variance"
    sigma <- SVD$d[1:n]
    if (power != 1) S <- scaleMargins(S, cols=sigma^(power - 1))
    
  } else if (method == "rsvd") {
    ## --- randomized SVD according to Halko, Martinsson & Tropp (2009, p. 9) ---

    ## preliminary testing suggests there is no substantial difference between the original and transposed rSVD algorith, so we currently always use the original version
    SVD <- rsvd(M, n=n, q=q, oversampling=oversampling, transpose=FALSE, verbose=verbose)

    sigma <- SVD$d
    if (power == 1) {
      S <- scaleMargins(SVD$u, cols=sigma)
    } else {
      S <- scaleMargins(SVD$u, cols=sigma^power)
    }
    if (with.basis) B <- SVD$v
    R2 <- sigma^2 / norm(M, "F")^2
    rm(SVD)

  } else if (method %in% c("ri", "ri+svd")) {
    ## --- straightforward random indexing (with specified fill rate), optionally followed by rSVD ---
    if (rate < 0 || rate > 1) stop("RI rate= must be between 0 and 1")
    
    ## TODO: -- check references on statistical guarantees, appropriate fill rates, etc.--
    if (method == "ri+svd") {
      nRI <- n * oversampling       # number of intermediate random dimensions
      nRI <- max(2*n, min(nRI, floor(nC / 2)))
    } else {
      nRI <- n                      # number of random dimensions
    }
    if (verbose) cat(sprintf("Random Indexing in %d dimensions:\n", nRI))
    
    if (!info$sparse) {
      ## if original matrix can be stored in dense representation, RI should not pose any memory problems
      n.fill <- as.integer(nC * rate) # number of nonzero entries in each random vector
      if (n.fill < 1) n.fill <- 1
      scale <- 1 / sqrt(n.fill)       # scale random vectors so they are normalised

      if (verbose) cat(sprintf(" - generating %d random vectors with %d of %d nonzero elements\n", nRI, n.fill, nC))
      Q <- matrix(0, nRI, nC)
      for (.d in 1:nRI) {
        .idx <- sort(sample.int(nC, n.fill))
        Q[.d, .idx] <- scale * (1 - 2*rbinom(n.fill, 1, .5))
      }

      if (verbose) cat(" - projecting into random subspace\n")
      S <- tcrossprod(M, Q)
      rm(Q)
    } else {
      ## efficient C implementation of RI for sparse matrix (which is guaranteed to be in canonical format)
      S <- CPP_random_indexing_sparse(nR, nC, M@p, M@i, M@x, nRI, rate, verbose)
    }
    
    if (method == "ri+svd") {
      S <- dsm.projection(S, "svd", n, verbose=verbose, with.basis=FALSE) # use plain SVD since matrix is already dense
    }
    attr(S, "R2") <- R2 <- NULL # partial R2 not accurate for non-orthogonal projection
  } else {
    stop("dimensionality reduction method '", method, "' has not been implemented yet")
  }

  dimnames(S) <- list(rownames(M), paste(method, 1:n, sep=""))
  if (with.basis) {
    dimnames(B) <- list(colnames(M), colnames(S))
    attr(S, "basis") <- B
  }
  if (!is.null(R2)) attr(S, "R2") <- R2
  if (!is.null(sigma)) attr(S, "sigma") <- sigma
  return(S)
}
