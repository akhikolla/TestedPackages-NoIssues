.check.norm <- function (method, p) {
  if (method == "minkowski" && p == Inf) {
    method <- "maximum"
    p <- 2
  }
  if (method == "minkowski" && (p < 0 || !is.finite(p))) stop("Minkowski p-norm con only be computed for 0 <= p < Inf")
  
  ## internal codes for selected norm (must match C code in <row_norms.c>)
  code <- switch(method, euclidean=0, maximum=1, manhattan=2, minkowski=3)
  if (is.null(code)) stop("unknown norm selected (internal error)")
  
  list(code=as.integer(code), p=as.double(p))
}

rowNorms <- function (M, method = "euclidean", p = 2) {
  method <- match.arg(method, c("euclidean", "maximum", "manhattan", "minkowski"))
  norm <- .check.norm(method, p)

  info <- dsm.is.canonical(M)
  if (!info$canonical) M <- dsm.canonical.matrix(M)

  if (info$sparse) {
    result <- CPP_row_norms_sparse(nrow(M), ncol(M), M@p, M@i, M@x, norm$code, norm$p)
  } else {
    result <- CPP_row_norms_dense(M, norm$code, norm$p)
  }

  names(result) <- rownames(M)
  result
}

colNorms <- function (M, method = "euclidean", p = 2) {
  method <- match.arg(method, c("euclidean", "maximum", "manhattan", "minkowski"))
  norm <- .check.norm(method, p)
  
  info <- dsm.is.canonical(M)
  if (!info$canonical) M <- dsm.canonical.matrix(M)

  if (info$sparse) {  
    result <- CPP_col_norms_sparse(nrow(M), ncol(M), M@p, M@i, M@x, norm$code, norm$p)
  } else {
    result <- CPP_col_norms_dense(M, norm$code, norm$p)
  }

  names(result) <- colnames(M)
  result
}
