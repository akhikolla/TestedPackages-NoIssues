subset.dsm <- function (x, subset=NULL, select=NULL,
                        recursive=FALSE, drop.zeroes=FALSE, matrix.only=FALSE,
                        envir=parent.frame(), run.gc = FALSE,
                        ...) {
  info <- check.dsm(x, validate=TRUE) # make sure that rows/columns are consistent
  if (recursive && matrix.only) stop("matrix.only=TRUE cannot be combined with recursive=TRUE")

  while (recursive) {
    force(envir) # fix parent environment to context of initial function call
    y <- do.call(subset.dsm, list(x=x, subset=substitute(subset), select=substitute(select), recursive=FALSE, drop.zeroes=drop.zeroes, matrix.only=FALSE, envir=envir))
    y.info <- check.dsm(y, validate=TRUE)
    if (y.info$nrow == info$nrow && y.info$ncol == info$ncol) return(y)
    x <- y
    info <- y.info
    if (run.gc) gc(verbose=FALSE) # avoid more than one duplicate during recursion
  }
  
  condition <- substitute(subset)
  row.idx <- eval(condition, c(x$rows, x$globals), envir)
  if (is.null(row.idx)) {
    row.idx <- 1:info$nrow
    update.nz.cols <- FALSE
  } else {
    ## todo: check validity (either Boolean of correct length or numeric vector with indexes in range)
    update.nz.cols <- TRUE # rows may have been deleted, so nnzero counts for columns need to be updated
  }
  
  condition <- substitute(select)
  col.idx <- eval(condition, c(x$cols, x$globals), envir)
  if (is.null(col.idx)) {
    col.idx <- 1:info$ncol
    update.nz.rows <- FALSE
  } else {
    ## todo: check validity (either Boolean of correct length or numeric vector with indexes in range)
    update.nz.rows <- TRUE # columns may have been deleted, so nnzero counts for rows need to be updated
  }

  if (is.logical(row.idx)) row.idx <- which(row.idx) # make sure we have numeric indices for additional subsetting
  if (is.logical(col.idx)) col.idx <- which(col.idx)

  if (!matrix.only && identical(row.idx, seq_len(info$nrow)) && identical(col.idx, seq_len(info$ncol))) {
    ## no change -> return unmodified DSM object (saves unnecessary copy in last iteration of recursive=TRUE)
    return(x)
  }
    
  if (drop.zeroes) {
    M <- if (info$S$ok) x$S else x$M  # primary data matrix (use scores if available, which may be sparser than frequencies)
    M.sub <- M[row.idx, col.idx, drop=FALSE] # proposed subset matrix
    nnzero.rows <- rowNorms(M.sub, method="minkowski", p=0)
    nnzero.cols <- colNorms(M.sub, method="minkowski", p=0)
    rm(M.sub) # will take the real subset matrix in the code below (it would be better to just keep M.sub if no zeroes need to be deleted)
    keep.rows <- nnzero.rows > 0
    row.idx <- row.idx[keep.rows] # drop rows without nonzero entries
    nnzero.rows <- nnzero.rows[keep.rows] # updated nnzero counts for final subset
    keep.cols <- nnzero.cols > 0
    col.idx <- col.idx[keep.cols] # drop columns without nonzero entries
    nnzero.cols <- nnzero.cols[keep.cols]
    rm(keep.rows, keep.cols)
  }

  if (matrix.only) {
    M <- if (info$S$ok) x$S else x$M  # matrix.only=TRUE: just return subset of the appropriate matrix
    return(M[row.idx, col.idx, drop=FALSE])
  }

  ## for small result sets, it is more memory-efficient to construct a new DSM object from scratch
  y <- list(rows=x$rows[row.idx, , drop=FALSE],   # mandatory components
            cols=x$cols[col.idx, , drop=FALSE],
            globals=x$globals)
  if (info$M$ok) {
    y$M <- x$M[row.idx, col.idx, drop=FALSE]
    if (!is.na(info$M$nonneg)) attr(y$M, "nonneg") <- info$M$nonneg
  }
  if (info$S$ok) {
    y$S <- x$S[row.idx, col.idx, drop=FALSE]
    if (!is.na(info$S$nonneg)) attr(y$S, "nonneg") <- info$S$nonneg
  }

  ## if rows and/or columns may have been deleted, update the relevant nonzero counts
  if (drop.zeroes) {
    y$rows$nnzero <- nnzero.rows  # we've already computed these above
    y$cols$nnzero <- nnzero.cols
  } else {
    if (update.nz.rows || update.nz.cols) {
      M <- if (info$S$ok) y$S else y$M # nonzero counts are based on the primary data matrix
      if (update.nz.rows) y$rows$nnzero <- rowNorms(M, method="minkowski", p=0)
      if (update.nz.cols) y$cols$nnzero <- colNorms(M, method="minkowski", p=0)
    }
  }

  structure(y, class=c("dsm", "list"))
}
