context.vectors <- function (M, contexts, split="\\s+", drop.missing=TRUE, row.names=NULL) {
  M <- find.canonical.matrix(M) # ensure that M is a suitable matrix, or extract matrix from DSM
  known.terms <- rownames(M)
  nR <- nrow(M)
  nC <- ncol(M)
  if (is.null(row.names)) {
    row.names <- if (is.null(names(contexts))) 1:length(contexts) else names(contexts)
  } else {
    if (length(row.names) != length(contexts)) stop("row.names= must have same length as contexts=")
  }
  if (is.character(contexts)) {
    tokens.list <- strsplit(contexts, split, perl=TRUE)
  } else {
    if (!is.list(contexts)) stop("contexts= must be either a character vector or a list")
    tokens.list <- contexts
  }

  CM <- t(vapply(tokens.list, function (tokens) {
    weights <- NULL
    if (is.character(tokens)) {
      ## context = vector of tokens
      idx <- na.omit(match(tokens, known.terms)) # row numbers of known terms in M (possibly repeated))
    } else if (is.logical(tokens)) {
      ## context = logical index vector into M (deprecated)
      if (length(tokens) != nR) stop("invalid logical index vector in contexts= (wrong length)")
      idx <- which(tokens)
    } else if (is.numeric(tokens)) {
      terms <- names(tokens)
      if (is.character(terms)) {
        ## context = weighted bag of words = vector of weights labelled with terms
        idx <- match(terms, known.terms)
        ok <- !is.na(idx)
        idx <- idx[ok]
        weights <- tokens[ok]
      }
      else if (is.integer(tokens)) {
        ## context = index of row numbers into M (deprecated)
        idx <- tokens
        if (length(idx) > 0 && (min(idx) < 1 || max(idx) > nR)) stop("invalid integer index vector in contexts= (row number out of range)")
      }
      else stop("invalid numeric vector without labels in contexts=")
    } else stop("invalid specification in contexts= (must be character, numeric or logical vector)")
    if (length(idx) > 0) {
      if (is.null(weights)) {
        colMeans(M[idx, , drop=FALSE]) # unweighted centroid vector
      }
      else {
        colSums(scaleMargins(M[idx, , drop=FALSE], rows=weights)) / sum(weights) # weighted centroid vector
      }
    } else {
      if (drop.missing) rep(NA, nC) else rep(0, nC) # return null vector for context without known tokens
    }
  }, FUN.VALUE=numeric(nC), USE.NAMES=FALSE))
  rownames(CM) <- row.names
  if (!is.null(colnames(M))) colnames(CM) <- colnames(M)
  if (drop.missing) {
    idx.miss <- is.na(CM[, 1]) # assuming there were no NAs or NaNs in M
    CM[!idx.miss, , drop=FALSE]
  } else {
    CM
  }
}
