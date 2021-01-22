dimnames.dsm <- function (x) {
  res <- check.dsm(x, validate=TRUE) # row/column names must be consistent
  if (res$S$ok) {
    dimnames(x$S)
  } else if (res$M$ok) {
    dimnames(x$M)
  } else stop("no co-occurrence matrix available in DSM object")
}

`dimnames<-.dsm` <- function (x, value) {
  info <- check.dsm(x) # don't need to check row/column names since they'll be overwritten
  if (length(value) != 2) stop("value must be a list of length two with row and column names, respectively")

  valR <- value[[1]] # validate new row and column names
  if (!is.character(valR)) stop("rownames must be a character vector")
  nR <- length(valR)
  if (nR != info$nrow) stop(sprintf("rownames should be a vector of length %d, but has %d elements", info$nrow, nR))

  valC <- value[[2]]
  if (!is.character(valC)) stop("colnames must be a character vector")
  nC <- length(valC)
  if (nC != info$ncol) stop(sprintf("rownames should be a vector of length %d, but has %d elements", info$ncol, nC))
  
  x$rows$term <- valR  # set in row/column info tables
  x$cols$term <- valC
  if (info$M$ok) dimnames(x$M) <- value # set both dimnames in single call to avoid overhead
  if (info$S$ok) dimnames(x$S) <- value
  
  x
}
