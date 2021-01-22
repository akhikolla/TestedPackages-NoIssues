as.matrix.dsm <- function (x, what=c("auto", "M", "S"), ...) {
  what <- match.arg(what)
  info <- check.dsm(x, validate=TRUE) # ensure that row/column names are correct
  if (what == "auto") {
    if (info$S$ok) what <- "S" 
    else if (info$M$ok) what <- "M"
    else stop("neither M nor S are available")
  }

  if (what == "M") {
    if (!info$M$ok) stop("co-occurrence matrix M is not available")
    dsm.canonical.matrix(x$M)
  }
  else if (what == "S") {
    if (!info$S$ok) stop("score matrix S is not available")
    dsm.canonical.matrix(x$S)
  }
  else stop(sprintf("internal error -- unsupported what = '%s'", what))
}
