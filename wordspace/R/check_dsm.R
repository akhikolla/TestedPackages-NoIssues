check.dsm <- function (model, validate=FALSE, nonneg.check=FALSE) {
  stopifnot(inherits(model, "dsm"))
  slots <- names(model)

  if ("M" %in% slots) {
    M <- dsm.is.canonical(model$M, nonneg.check=nonneg.check)
    M$ok <- TRUE
  } else {
    M <- data.frame(ok=FALSE)
  }
  rownames(M) <- "M"
  
  if ("S" %in% slots) {
    S <- dsm.is.canonical(model$S, nonneg.check=nonneg.check)
    S$ok <- TRUE
  } else {
    S <- data.frame(ok=FALSE)
  }
  rownames(S) <- "S"

  stopifnot(M$ok || S$ok) # need to have either frequency matrix or score matrix (or both)
  stopifnot(all(c("rows", "cols", "globals") %in% slots))
  required <- if (M$ok) c("term", "f") else c("term") # required columns in $rows and $cols
  stopifnot(all(required %in% colnames(model$rows)))
  stopifnot(all(required %in% colnames(model$cols)))

  N <- if ("N" %in% names(model$globals)) model$globals$N else NA
  if (M$ok && !is.finite(N)) stop("missing information on sample size N (but frequency matrix M is present)")

  is.locked <- if ("locked" %in% names(model$globals)) model$globals$locked else FALSE
  
  n.rows <- nrow(model$rows)
  n.cols <- nrow(model$cols)

  if (M$ok) {
    stopifnot(nrow(model$M) == n.rows)
    stopifnot(ncol(model$M) == n.cols)
    if (validate) {
      stopifnot(all(rownames(model$M) == model$rows$term))
      stopifnot(all(colnames(model$M) == model$cols$term))
    }
  }

  if (S$ok) {
    stopifnot(nrow(model$S) == n.rows)
    stopifnot(ncol(model$S) == n.cols)
    if (validate) {
      stopifnot(all(rownames(model$S) == model$rows$term))
      stopifnot(all(colnames(model$S) == model$cols$term))
    }
  }
  
  list(nrow=n.rows, ncol=n.cols, N=N, M=M, S=S, locked=is.locked)
}
