.complain.about.words <- function (w, what="words") {
  n <- length(w)
  if (n >= 10) {
    sprintf("%d missing %s: %s, ...", n, what, paste(w[1:6], collapse=", "))
  }
  else {
    sprintf("missing %s: %s", what, paste(w, collapse=", "))
  }
}

dsm <- function (M = NULL, target = NULL, feature = NULL, score = NULL, rowinfo = NULL, colinfo = NULL, N = NA, globals = list(), raw.freq = FALSE, sort = FALSE, verbose = FALSE) {
  if (!is.null(M)) {

    if (!is.null(score)) stop("you must specify either M= or score=, target= and feature=")
    flags <- dsm.is.canonical(M)
    if (verbose) cat(sprintf("Compiling DSM from %d x %d co-occurrence matrix ...\n", nrow(M), ncol(M)))

    if (!flags$canonical) {
      if (verbose) cat(sprintf(" - converting %s matrix into canonical format\n", if (flags$sparse) "sparse" else "dense"))
      M <- dsm.canonical.matrix(M)
    }

    if (!is.null(target)) {
      if (length(target) != nrow(M)) stop("target= must provide one label for each row of M")
      rownames(M) <- target
    } else {
      if (is.null(rownames(M))) stop("target= must be specified if M doesn't have row labels")
    }

    if (!is.null(feature)) {
      if (length(feature) != ncol(M)) stop("feature= must provide one label for each column of M")
      colnames(M) <- feature
    } else {
      if (is.null(colnames(M))) stop("feature= must be specified if M doesn't have column labels")
    }
    
  } else {

    if (is.null(score) || is.null(target) || is.null(feature)) stop("you must specify either M= or score=, target= and feature=")
    n.items <- length(target)
    ok <- (length(score) == n.items || length(score) == 1) && length(feature) == n.items
    if (!ok) stop("score=, target= and feature= must be vectors of the same length")
    stopifnot(is.numeric(score))
    if (is.factor(target)) target <- as.character(target)
    if (is.factor(feature)) feature <- as.character(feature)
    stopifnot(is.character(target) && is.character(feature))
    if (verbose) cat(sprintf("Compiling DSM from %.2fM items in triplet representation ...\n", n.items / 1e6))
    
    if (verbose) cat(" - target & feature terms\n")
    t.dict <- unique(target)  # compile lists of target and feature types
    f.dict <- unique(feature)
  
    if (sort) {
      t.dict <- sort(t.dict)
      f.dict <- sort(f.dict)
    }
  
    row.idx <- match(target, t.dict) # row/column indices of triplets
    col.idx <- match(feature, f.dict)
  
    if (verbose) cat(" - building sparse matrix\n")
    M <- sparseMatrix(i=row.idx, j=col.idx, x=score, giveCsparse=TRUE) # may either be M or S in the final object
    rm(row.idx, col.idx) # free memory
    rownames(M) <- t.dict
    colnames(M) <- f.dict

    flags <- list(sparse=TRUE)
  }

  if (raw.freq) {
    if (verbose) cat(" - checking non-negative frequency counts\n")
    if (!signcount(M, "nonneg")) stop("raw frequency counts must be non-negative")
    attr(M, "nonneg") <- TRUE
  }

  if (verbose) cat(" - collecting target and feature information\n")
  if (is.null(rowinfo)) {
    rowinfo <- data.frame(term=as.character(rownames(M)), stringsAsFactors=FALSE)
  } else {
    if (!("term" %in% colnames(rowinfo))) stop("rowinfo= must specify target types in column 'term'")
    rowinfo$term <- as.character(rowinfo$term)
    idx <- match(rownames(M), rowinfo$term)
    if (any(is.na(idx))) stop(.complain.about.words(rownames(M)[is.na(idx)], "target types in rowinfo"))
    rowinfo <- rowinfo[idx, , drop=FALSE]
  }

  if (is.null(colinfo)) {
    colinfo <- data.frame(term=as.character(colnames(M)), stringsAsFactors=FALSE)
  } else {
    if (!("term" %in% colnames(colinfo))) stop("colinfo= must specify feature types in column 'term'")
    colinfo$term <- as.character(colinfo$term)
    idx <- match(colnames(M), colinfo$term)
    if (any(is.na(idx))) stop(.complain.about.words(colnames(M)[is.na(idx)], "feature types in colinfo"))
    colinfo <- colinfo[idx, , drop=FALSE]
  }

  if (verbose) cat(" - computing marginal statistics\n")
  if (is.na(N) && !is.null(globals$N)) N <- globals$N  # use sample size from globals unless specified with N=
  if (!("nnzero" %in% colnames(rowinfo))) rowinfo$nnzero <- rowNorms(M, method="minkowski", p=0) # we now have efficient nonzero counts with "Hamming length"
  if (!("nnzero" %in% colnames(colinfo))) colinfo$nnzero <- colNorms(M, method="minkowski", p=0)
  if (raw.freq) {
    if (!("f" %in% colnames(rowinfo))) rowinfo$f <- rowSums(M) # M must be non-negative at this point
    if (!("f" %in% colnames(colinfo))) colinfo$f <- colSums(M)
    if (is.na(N)) N <- sum(M)
  }

  if (flags$sparse) {
    n.nzero <- signcount(M, "nnzero")
    if (verbose) cat(sprintf("%d x %d matrix with %d nonzero entries (= %.1f%%)\n", nrow(M), ncol(M), n.nzero, 100 * n.nzero / prod(dim(M))))
  } else {
    if (verbose) cat(sprintf("%d x %d matrix with %.2fM cells\n", nrow(M), ncol(M), prod(dim(M)) / 1e6))
  }

  globals$N <- N
  globals$locked <- FALSE
  if (raw.freq) {
    dsm.obj <- list(M=M, rows=rowinfo, cols=colinfo, globals=globals)
  } else {
    dsm.obj <- list(S=M, rows=rowinfo, cols=colinfo, globals=globals)
  }
  
  structure(dsm.obj, class=c("dsm", "list"))
}
