dsm.is.canonical <- function (x, nonneg.check = FALSE) {
  if (is.matrix(x) && is.numeric(x)) {
    sparse <- FALSE  # regular dense numeric matrix
    canonical <- TRUE
  } else {
    if (!is(x, "dMatrix")) stop("first argument must be a dense or sparse numeric matrix")
    sparse <- if (is(x, "sparseMatrix")) TRUE else FALSE
    canonical <- if (sparse && is(x, "dgCMatrix")) TRUE else FALSE
  }
  if (nonneg.check) {
    if (!sparse || is(x, "dgeMatrix") || is(x, "dgCMatrix") || is(x, "dgRMatrix")) {
      nonneg <- signcount(x, "nonneg") # efficient non-negativity check is supported for these formats
    } else {
      nonneg <- !any(x < 0) # caution: all(x >= 0) would result in dense matrix!
    }
  } else {
    nonneg <- attr(x, "nonneg")
    if (is.null(nonneg)) nonneg <- NA
  }
  data.frame(sparse=sparse, canonical=canonical, nonneg=nonneg)
}

dsm.canonical.matrix <- function (x, triplet = FALSE, annotate = FALSE, nonneg.check = FALSE) {
  if (nonneg.check && !annotate) {
    warning("nonneg.check=TRUE ignored without annotate=TRUE")
    nonneg.check <- FALSE
  }
  flags <- dsm.is.canonical(x, nonneg.check=FALSE) # nonneg check deferred to when we have a canonical matrix
  if (flags$sparse) {
    if (triplet) {
      if (nonneg.check && flags$canonical) {
        flags$nonneg <- signcount(x, "nonneg") # efficient non-negativity check possible before conversion to triplet form
        nonneg.check <- FALSE
      }
      if (!is(x, "dgTMatrix")) x <- as(x, "dgTMatrix")
    } else {
      if (!flags$canonical) x <- as(x, "dgCMatrix")
    }
  } else {
    if (!flags$canonical) x <- as.matrix(x)
  }
  if (annotate) {
    attr(x, "sparse") <- flags$sparse
    attr(x, "nonneg") <- if (nonneg.check) dsm.is.canonical(x, nonneg.check=TRUE)$nonneg else flags$nonneg
  }
  x
}
