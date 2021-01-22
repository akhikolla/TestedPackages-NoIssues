pair.distances <- function (w1, w2, M, ..., transform=NULL, rank=c("none", "fwd", "bwd", "avg"), avg.method=c("arithmetic", "geometric", "harmonic"), batchsize=10e6, verbose=FALSE) {
  rank <- match.arg(rank)
  avg.method <- match.arg(avg.method)
  if (!is.null(transform) && !is.function(transform)) stop("transform= must be a vectorized function expecting a single argument")
  w1 <- as.character(w1)
  w2 <- as.character(w2)
  stopifnot(length(w1) == length(w2))
  is.dist <- inherits(M, "dist.matrix") || isTRUE(attr(M, "dist.matrix"))
  if (!is.dist) M <- find.canonical.matrix(M) # ensure DSM matrix is in canonical format, or extract from DSM object
  
  if (rank != "none") {
    ## case 1: dispatch to specialised function for computing neighbour ranks
    if (rank == "fwd") {
      r <- pair.ranks(w1, w2, M, ..., rank="fwd", is.dist=is.dist, batchsize=batchsize, verbose=verbose)
    } else if (rank == "bwd") {
      r <- pair.ranks(w1, w2, M, ..., rank="bwd", is.dist=is.dist, batchsize=batchsize, verbose=verbose)
    } else {
      r1 <- pair.ranks(w1, w2, M, ..., rank="fwd", is.dist=is.dist, batchsize=batchsize, verbose=verbose)
      r2 <- pair.ranks(w1, w2, M, ..., rank="bwd", is.dist=is.dist, batchsize=batchsize, verbose=verbose)
      r <- switch(avg.method,
                  arithmetic = (r1 + r2) / 2,
                  geometric = sqrt(r1 * r2),
                  harmonic = ifelse(is.finite(pmax(r1, r2)), 2 * r1 * r2 / (r1 + r2), Inf))
    }
    if (!is.null(transform)) transform(r) else r
  } else {
    ## case 2: compute regular distances or similarities (rank == "none")
    n <- length(w1)
    types1 <- unique(w1)
    types2 <- unique(w2)
    n.types1 <- as.double(length(types1))
    n.types2 <- as.double(length(types2))
    
    ## if there are too many distinct types (such that intermediate distance matrix would have > chunksize elements), 
    ## partition input recursively and combine result vectors (unless M is a pre-computed distance matrix)
    if (!is.dist) {
      n.elements <- n.types1 * n.types2 # size of distance matrix to be computed
      split.batch <- n.elements > batchsize && n >= 4
      if (verbose) cat(sprintf("%s- pair.distances(): %d pairs, %d x %d types = %.1fM elements %s\n", paste(rep(" ", verbose), collapse=""), n, length(types1), length(types2), n.elements/1e6, if (split.batch) "" else "***"))
      if (split.batch) {
        pivot <- floor(n/2)
        verbose.val <- if (verbose) verbose + 2 else FALSE
        res1 <- pair.distances(w1[1:pivot], w2[1:pivot], M, ..., batchsize=batchsize, verbose=verbose.val)
        res2 <- pair.distances(w1[(pivot+1):n], w2[(pivot+1):n], M, ..., batchsize=batchsize, verbose=verbose.val)
        is.similarity <- isTRUE(attr(res1, "similarity")) # pass through similarity marker
        res <- structure(c(res1, res2), similarity=is.similarity)
        if (!is.null(transform)) return(transform(res)) else return(res)
      }
    }
    
    if (is.dist) {
      ## case 2a: look up distances or similarities directly in pre-computed matrix M
      distances <- M  # NB: is M is sparse, it cannot be indexed with a matrix of labels below
    }
    else {
      ## case 2b: compute distance matrix as a superset of the required distances between rows of M
      distances <- dist.matrix(M, byrow=TRUE, terms=types1, terms2=types2, skip.missing=TRUE, ...)
    }
    is.similarity <- isTRUE(attr(distances, "similarity"))
    miss.val <- if (is.similarity) -Inf else Inf

    res <- rep(miss.val, n)
    w1.row <- match(w1, rownames(distances)) # row index of w1
    w2.col <- match(w2, colnames(distances)) # column index of w1
    is.known <- !is.na(w1.row) & !is.na(w2.col)
    res[is.known] <- distances[cbind(w1.row[is.known], w2.col[is.known])]
    res <- structure(res, names=paste(w1, w2, sep="/"), similarity=is.similarity)
    if (!is.null(transform)) return(transform(res)) else return(res)
  }
}

## need specialised implemenation for neighbour ranks, which loops over w1 types instead of word pairs,
## in order to avoid repeated expensive computation of all-neighbour rankings
pair.ranks <- function (w1, w2, M, ..., rank=c("fwd", "bwd"), is.dist=FALSE, batchsize=10e6, verbose=FALSE) {
  rank <- match.arg(rank)
  ## already ensured by pair.distances()
  # w1 <- as.character(w1)
  # w2 <- as.character(w2)
  # stopifnot(length(w1) == length(w2))
  res <- rep(Inf, length(w1))

  if (is.dist) {
    ## case 1: pre-computed distance or similarity matrix
    
    is.similarity <- isTRUE(attr(M, "similarity"))
    is.symmetric <- isTRUE(attr(M, "symmetric")) # whether to adjust for word as its own neighbour
    is.sparse <- dsm.is.canonical(M)$sparse
    if (is.sparse && !is.similarity) stop("only non-negative similarity matrix supported in sparse format")
    is.known <- w1 %in% rownames(M) & w2 %in% colnames(M) # word pairs found in pre-computed distance matrix
    
    u1 <- if (rank == "fwd") w1 else w2 # compute rank of u2 among neighbours of u1
    u2 <- if (rank == "fwd") w2 else w1 # but we still have to decide between rows and columns of M below
    u1.types <- sort(unique(u1[is.known])) # list of u1 types we need to process
    n.types <- length(u1.types)
    items.per.batch <- ceiling(batchsize / (if (rank == "fwd") ncol(M) else nrow(M))) # number of u1 types we can process per batch

    for (i.start in seq(1, n.types, items.per.batch)) {
      i.end <- min(i.start + items.per.batch - 1, n.types)
      if (verbose) cat(sprintf(" - pair.ranks(as.dist, '%s'): types #%d .. #%d of %d (%s .. %s)\n", 
                               rank, i.start, i.end, n.types, u1.types[i.start], u1.types[i.end]))
      batch.types <- u1.types[i.start:i.end] # extract and rank row ("fwd") or column ("bwd") vectors for these types from M
      distances <- if (rank == "fwd") t(as.matrix(M[batch.types, , drop=FALSE])) else as.matrix(M[, batch.types, drop=FALSE])
      if (is.similarity) distances <- -distances # to rank by decreasing similarity
      ranks <- apply(distances, 2, function (x) {
        res <- rank(x, ties.method="min")
        if (is.sparse) res[x == 0] <- Inf # empty cells in sparse matrix cannot be neighbours
        res
      }) # ranks is column-major: look up ranks[u2, u1]
      idx <- (u1 %in% batch.types) & is.known
      adj.rank <- if (is.symmetric) 1 else 0         # adjust ranks for u1 as its own first neighbour if M is symmetric
      res[idx] <- ranks[cbind(u2[idx], u1[idx])] - adj.rank
    }
  } else {
    ## case 2: compute distances between rows of DSM matrix
    
    u1 <- if (rank == "fwd") w1 else w2 # compute rank of u2 among neighbours of u1
    u2 <- if (rank == "fwd") w2 else w1
    known.words <- rownames(M)
    is.known <- u1 %in% known.words & u2 %in% known.words # word pairs found in DSM
    u1.types <- sort(unique(u1[is.known])) # list of u1 types found in DSM
    n.types <- length(u1.types)
    items.per.batch <- ceiling(batchsize / nrow(M)) # number of u1 types we can process in a single batch
    
    for (i.start in seq(1, n.types, items.per.batch)) {
      i.end <- min(i.start + items.per.batch - 1, n.types)
      if (verbose) cat(sprintf(" - pair.ranks(): types #%d .. #%d of %d (%s .. %s)\n", 
                               i.start, i.end, n.types, u1.types[i.start], u1.types[i.end]))
      batch.types <- u1.types[i.start:i.end]
      distances <- dist.matrix(M, byrow=TRUE, terms=batch.types, terms2=NULL, skip.missing=FALSE, ...)
      is.similarity <- isTRUE(attr(distances, "similarity"))
      if (is.similarity) distances <- -distances
      ranks <- apply(distances, 1, function (x) rank(x, ties.method="min")) # ranks is column-major: ranks[u2, u1]
      idx <- (u1 %in% batch.types) & is.known
      res[idx] <- ranks[cbind(u2[idx], u1[idx])] - 1 # adjust for u1 as its own neighbour (cannot be used for cross-distances!)
    }
  }
    
  structure(res, names=paste(w1, w2, sep="/"), similarity=FALSE)
}
