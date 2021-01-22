nearest.neighbours <- function (M, term, n=10, M2=NULL, byrow=TRUE, drop=TRUE, skip.missing=FALSE, dist.matrix=FALSE, ..., batchsize=50e6, verbose=FALSE) {
  is.dist <- inherits(M, "dist.matrix") || isTRUE(attr(M, "dist.matrix"))
  if (is.dist) {
    ## case 1: M is a pre-computed distance matrix
    ##  - byrow determines if target term is looked up in rows or columns of M
    ##  - if M is not marked symmetric, we need to treat this as a cross-distance computation
    cross.distance <- !isTRUE(attr(M, "symmetric"))
    if (!is.null(M2)) stop("M2 cannot be specified if M is a pre-computed distance matrix")
    if (dist.matrix && cross.distance) stop("pre-computed distance matrix must be symmetric for dist.matrix=TRUE")
  } else {
    ## case 2: M is a matrix of row or column vectors
    ##   - compute cross-distances if M2 is specified, otherwise (usually symmetric) distances in M
    ##   - byrow determines whether row or column vectors of M and M2 are used
    cross.distance <- !is.null(M2)
    M <- find.canonical.matrix(M) # ensure that M is a suitable matrix, or extract from DSM object
    if (cross.distance) {
      M2 <- find.canonical.matrix(M2) # ensure that M2 is a suitable matrix, or extract from DSM object
      if (byrow && ncol(M) != ncol(M2)) stop("M and M2 are not conformable (must have the same number of columns)")
      if (!byrow && nrow(M) != nrow(M2)) stop("M and M2 are not conformable (must have the same number of rows)")
    } else {
      M2 <- M                           # so we can always compute distances against M2
    }
  }

  nn.of.vector <- is.numeric(term)
  if (nn.of.vector) {
    ## case 1: targets are given as row vectors
    if (is.dist) stop("cannot find nearest neighbours of vector if M is a pre-computed distance matrix")
    if (!is.matrix(term)) term <- matrix(term, nrow=1) # convert plain vector into row vector
    if (byrow && ncol(term) != ncol(M2)) stop("target vectors in term= do not have the right number of dimensions")
    if (!byrow && ncol(term) != nrow(M2)) stop("target vectors in term= do not have the right number of dimensions")
    n.terms <- nrow(term)
    if (is.null(rownames(term))) rownames(term) <- paste0("VEC", 1:n.terms)
  } else {
    ## case 2: targets are given as strings (to be looked up in rows or columns of M)
    term.labels <- if (byrow) rownames(M) else colnames(M)
    found <- term %in% term.labels
    if (any(!found) && !skip.missing) stop("target term(s) not found in M: ", paste(term[!found], collapse=", "))
    term <- term[found]
    n.terms <- length(term)
  }
  if (n.terms == 0) return(NULL)
  
  ## unless we're working on a pre-computed dist.matrix, process vector of lookup terms in moderately sized batches
  if (!is.dist) {
    n.cand <- if (byrow) nrow(M2) else ncol(M2) # neighbour candidates
    if (n.cand > 1 && as.double(n.terms) * n.cand > batchsize) {
      items.per.batch <- ceiling(batchsize / n.cand)
      res.list <- lapply(seq(1, n.terms, items.per.batch), function (i.start) {
        i.end <- min(i.start + items.per.batch - 1, n.terms)
        if (verbose) cat(sprintf(" - nearest.neighbours(): terms #%d .. #%d of %d (size = %.1fM)\n", i.start, i.end, n.terms, (i.end-i.start+1) * n.cand / 1e6))
        term.batch <- if (nn.of.vector) term[i.start:i.end, , drop=FALSE] else term[i.start:i.end]
        nearest.neighbours(M, term.batch, n=n, M2=if (cross.distance) M2 else NULL, drop=FALSE, skip.missing=FALSE, byrow=byrow, dist.matrix=dist.matrix, ..., batchsize=Inf, verbose=verbose)
      })
      return(do.call(c, res.list))
    }
  }

  ## items to look up
  items <- if (nn.of.vector) rownames(term) else term 
  
  ## prepare distance matrix DM between all targets and candidates,
  ## arranging so that it is always accessed by column (which is more efficient)
  if (is.dist) {
    ## pre-computed distance matrix: extract relevant rows or columns
    if (byrow) {
      items.ok <- items[items %in% rownames(M)]
      DM <- t(M[items.ok, , drop=FALSE]) # now accessed by column
    } else {
      items.ok <- items[items %in% colnames(M)]
      DM <- M[, items.ok, drop=FALSE]
    }
    ## there should be methods [.dist.matrix and [<-.dist.matrix so that we don't need to reconstruct a dist.matrix object here, 
    ## but this would probably make row and column access considerably slower; and a sparse distance matrix works differently anyway
    DM <- as.distmat(DM, similarity=isTRUE(attr(M, "similarity"))) # won't be symmetric, even if M is
  } else {
    ## compute distance matrix between specified items and all other targets
    if (nn.of.vector) {
      M.term <- if (byrow) term else t(term) # M.term = matrix of target vectors (rows or columns)
    } else {
      M.term <- if (byrow) M[term, , drop=FALSE] else M[, term, drop=FALSE]
    }
    ## it's more efficient to have the smaller matrix M.term first because it will be the inner loop
    ## when computing a general distance matrix (resulting in better cache coherence)
    DM <- t(dist.matrix(M=M.term, M2=M2, byrow=byrow, ...)) # items correspond to columns, regardless of <byrow>
  }
  similarity <- isTRUE(attr(DM, "similarity"))
  sparse <- dsm.is.canonical(DM)$sparse # may only happen for pre-computed similarity matrix
  if (sparse && !similarity) stop("only non-negative similarity matrix supported in sparse format")
  
  result <- lapply(items, function (.t) {
    if (sparse) {
      ## sparse pre-computed similarity matrix
      neighbours <- DM[, .t]
      neighbours <- neighbours[neighbours > 0] # only non-zero cells are candidates for neighbours
    } else {
      ## in all other cases
      neighbours <- DM[, .t]
    }
    neighbours <- sort(neighbours, decreasing=similarity)
    if (!nn.of.vector && !cross.distance) {
      neighbours <- head(neighbours, n + 1) # remove target from list of nearest neighbours
      neighbours <- neighbours[names(neighbours) != .t] # this should remove at most 1 element
    }
    neighbours <- head(neighbours, n)
    if (dist.matrix) {
      ## case 1: compute distance matrix between nearest neighbours (including target term)
      nn.terms <- names(neighbours)
      if (is.dist) {
        ## pre-computed distance matrix M must be symmetric (i.e. !cross.distance)
        nn.terms <- c(.t, nn.terms) # prepend target term to list of neighbours (has been removed above)
        nn.dist <- M[nn.terms, nn.terms] 
        nn.dist <- as.distmat(nn.dist, symmetric=TRUE, similarity=similarity)
      } else {
        ## dense (cross-)distance matrix computed on the fly
        nn.matrix <- if (byrow) M2[nn.terms, , drop=FALSE] else t(M2[, nn.terms, drop=FALSE]) # matrix of row vectors for all neighbours
        nn.terms <- c(.t, nn.terms) # add target term
        if (nn.of.vector) {
          nn.matrix <- rbind2(term[.t, , drop=FALSE], nn.matrix) # add specified target vector
          ## NB: rbind() doesn't dispatch to Matrix package; need rbind2() or rBind()
        } else {
          nn.matrix <- rbind2(if (byrow) M[.t, , drop=FALSE] else t(M[, .t, drop=FALSE]), nn.matrix) # add target vector from M
        }
        rownames(nn.matrix) <- nn.terms
        nn.dist <- dist.matrix(nn.matrix, byrow=TRUE, ...)
      }
      attr(nn.dist, "selected") <- nn.terms == .t # mark target as selected
      nn.dist
    } else {
      ## case 2: return specified number of nearest neighbours
      neighbours
    }
  })
  names(result) <- items
  
  if (drop && length(result) == 1) result[[1]] else result
}
