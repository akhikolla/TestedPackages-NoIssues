scaleMargins <- function (M, rows=NULL, cols=NULL, duplicate=TRUE) {
  info <- dsm.is.canonical(M)
  if (!info$canonical) M <- dsm.canonical.matrix(M)
  nr <- nrow(M)
  nc <- ncol(M)
  
  if (is.null(rows)) {
    rows <- rep(1, nr)
  } else {
    if (length(rows) == 1) rows <- rep(rows, nr)
    if (length(rows) != nr) stop("rows= must either be a scalar or conformable with the rows of M=")
  }

  if (is.null(cols)) {
    cols <- rep(1, nc)
  } else {
    if (length(cols) == 1) cols <- rep(cols, nc)
    if (length(cols) != nc) stop("cols= must either be a scalar or conformable with the columns of M=")
  }
  
  if (info$sparse) {
    CPP_scale_margins_sparse(M, rows=rows, cols=cols, duplicate=duplicate)
  } else {
    CPP_scale_margins_dense(M, rows=rows, cols=cols, duplicate=duplicate)
  }
}
