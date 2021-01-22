rbind.dsm <- function (..., term.suffix=NULL, deparse.level=1) {
  .Deprecated(msg="The 'rbind' method for DSM objects is experimental. It may be removed or modified in a future release of the 'wordspace' package.")
  
  models <- list(...) # should be one or more objects of class "dsm" if rbind() dispatches here
  if (!is.null(term.suffix) && length(term.suffix) != length(models)) stop("term.suffix must provide as many strings as there are DSMs")

  models.info <- lapply(models, check.dsm, validate=TRUE) # validate models and extract dimensions etc.

  have.M.vec <- sapply(models.info, function (i) i$M$ok) # are raw frequencies available?
  have.M <- all(have.M.vec)
  if (any(have.M.vec) && !have.M) stop("either all DSM objects must contain raw frequency data (matrix M), or none of them")

  have.S.vec <- sapply(models.info, function (i) i$S$ok) # are scored matrices available?
  have.S <- all(have.S.vec)
  if (any(have.S.vec) && !have.S) warning("some but not all DSM objects contain score matrix S, dropped from result")

  if (!(have.M || have.S)) stop("neither raw frequencies M nor score matrices S found consistently across all DSM objects")

  any.locked <- any(sapply(models.info, function (i) i$locked)) # are any of the input DSMs locked already?
  N.vec <- sapply(models.info, function (i) i$N) # extract sample sizes

  cols.merged <- .combine.marginals(lapply(models, function (m) m$cols), margin="column", mode="same") # check feature dimensions, then combine column marginals
  marginals.inconsistent <- attr(cols.merged, "adjusted") # if marginals or sample sizes differ between DSMs, they're adjusted to the maximum value to ensure consistency
  if (have.M) {
    if (diff(range(N.vec)) >= 1) marginals.inconsistent <- TRUE
    N <- max(N.vec)
  } else {
    N <- NA
  }
  if (marginals.inconsistent && !have.S) warning("DSM objects have inconsistent column marginals / sample sizes, should calculate scores before combining them")

  rows.merged <- .bind.marginals(lapply(models, function (m) m$rows), margin="row", term.suffix=term.suffix)
  
  ## TODO: rBind is memory-inefficient because it recursively combines two matrices at a time
  ##  - replace by custom implementation (using triplet representation?), which need not preserve row/col names
  ##  - standard rbind() for matrices seems compact and fast
  res <- list(
    rows = rows.merged,
    cols = cols.merged,
    globals = models[[1]]$globals,
    locked = marginals.inconsistent || any.locked
  )
  if (have.M) {
    res$globals$N <- N
    any.sparse <- any(sapply(models.info, function (i) i$M$sparse)) # are there any sparse matrices?
    res$M <- do.call(if (any.sparse) rBind else rbind, lapply(models, function (m) m$M))
    dimnames(res$M) <- list(res$rows$term, res$cols$term)
  }
  if (have.S) {
    any.sparse <- any(sapply(models.info, function (i) i$S$sparse)) # are there any sparse matrices?
    res$S <- do.call(if (any.sparse) rBind else rbind, lapply(models, function (m) m$S))
    dimnames(res$S) <- list(res$rows$term, res$cols$term)
  }

  class(res) <- c("dsm", "list")
  return(res)
}

cbind.dsm <- function (..., term.suffix=NULL, deparse.level=1) {
  .Deprecated(msg="The 'cbind' method for DSM objects is experimental. It may be removed or modified in a future release of the 'wordspace' package.")
  stop("not yet implemented")
}

merge.dsm <- function (x, y, ..., rows=TRUE, all=FALSE, term.suffix=NULL) {
  .Deprecated(msg="The 'merge' method for DSM objects is deprecated. It will be removed in the next release of the 'wordspace' package and may be re-introduced later with different semantics.")
  models <- list(x, y, ...)
  n.models <- length(models)
  if (!is.null(term.suffix) && length(term.suffix) != n.models) stop("term.suffix must provide as many strings as there are DSMs")
  models.info <- lapply(models, check.dsm, validate=TRUE) # validate models and extract dimensions etc.

  if (all) stop("all=TRUE is not yet implemented")
  if (!rows) stop("rows=FALSE is not yet implemented")

  have.M.vec <- sapply(models.info, function (i) i$have.M) # are raw frequencies available?
  have.M <- all(have.M.vec)
  if (any(have.M.vec) && !have.M) stop("either all DSM objects must contain raw frequency data (matrix M), or none of them")

  have.S.vec <- sapply(models.info, function (i) i$have.S) # are scored matrices available
  have.S <- all(have.S.vec)
  if (any(have.S.vec) && !have.S) warning("some but not all DSM objects contain score matrix S, dropped from result")

  if (!(have.M || have.S)) stop("neither raw frequencies M nor score matrices S found consistently across all DSM objects")

  any.locked <- any(sapply(models.info, function (i) i$locked)) # are any of the input DSMs locked already?
  any.sparse <- any(sapply(models.info, function (i) i$sparse)) # are there any sparse matrices?
  N.vec <- sapply(models, function (i) i$N) # extract sample sizes

  # bind rows of DSM objects, preserving only terms that are shared by all DSM
  cols.merged <- .combine.marginals(lapply(models, function (m) m$cols), margin="column", mode="intersect")
  marginals.inconsistent <- attr(cols.merged, "adjusted") # if marginals differ between DSMs, they're adjusted to the maximum value to ensure consistency
  if (have.M) {
    if (diff(range(N.vec)) >= 1) marginals.inconsistent <- TRUE
    N <- max(N.vec)
  } else {
    N <- NA
  }
  if (marginals.inconsistent && !have.S) warning("DSM objects have inconsistent column marginals / sample sizes, should calculate scores before combining them")
  
  if (any.sparse) {
    ## for sparse matrices, extract/reorder columns, then call rbind() to merge the models
    ## TODO: this wastes huge amounts of memory; needs to be reimplemented in a more sophisticated way!
    adjusted.models <- lapply(models, function (.m) subset(.m, select=na.omit( match(cols.merged$term, .m$cols$term) )))
    return( do.call(rbind, c(adjusted.models, list(term.suffix=term.suffix))) )
  }
  else {
    ## for dense matrices, build merged DSM directly, filling in a pre-allocated matrix (more memory-efficient)
    rows.merged <- .bind.marginals(lapply(models, function (m) m$rows), margin="row", term.suffix=term.suffix)
    n.rows <- sapply(models.info, function (i) i$nrow)
    first.row <- cumsum(c(1, n.rows)) # rows offsets of individual DSMs in combined matrix M
  
    if (have.M) M <- matrix(nrow=nrow(rows.merged), ncol=nrow(cols.merged), dimnames=list(rows.merged$term, cols.merged$term))
    if (have.S) S <- matrix(nrow=nrow(rows.merged), ncol=nrow(cols.merged), dimnames=list(rows.merged$term, cols.merged$term))
    for (i in 1:n.models) {
      model <- models[[i]]
      col.idx <- na.omit( match(cols.merged$term, model$cols$term) ) # extract columns of i-th DSM matrix that match the common terms
      if (have.M) M[ (first.row[i]):(first.row[i]+n.rows[i]-1), ] <- model$M[ , col.idx]
      if (have.S) S[ (first.row[i]):(first.row[i]+n.rows[i]-1), ] <- model$S[ , col.idx]
    }
    
    ## construct and return merged DSM
    res <- list(rows = rows.merged,
                cols = cols.merged,
                globals = models[[1]]$globals,
                locked = marginals.inconsistent || any.locked)
    if (have.M) {
      res$M <- M
      res$globals$N <- res$N <- N
    }
    if (have.S) res$S <- S

    class(res) <- c("dsm", "list")
    return(res)
  }
}  
