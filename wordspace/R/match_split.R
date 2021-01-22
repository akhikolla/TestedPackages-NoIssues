match.split <- function (x, f, values=NULL, groups=NULL, nomatch=NA_integer_) {
  if (length(x) != length(f)) stop("arguments x and f must be vectors of the same length")
  if (missing(groups)) {
    f <- as.factor(f)
    groups <- levels(f)
  } else {
    f <- factor(f, levels=groups)
  }
  n <- length(x)
  
  idx.split <- split(1:n, f)
  if (is.null(values)) {
    x.split <- lapply(idx.split, function (.idx) unique(x[.idx]))
    values <- Reduce(intersect, x.split)
    if (length(values) < 1) stop("no values are shared between all groups")
  }

  ## find positions of first matches in each group, then obtain original positions with .idx[...]
  res <- do.call(cbind, lapply(idx.split, function (.idx) .idx[ match(values, x[.idx]) ]))

  if (!missing(nomatch)) res[is.na(res)] <- as.integer(nomatch) # replace NA's for "no match" by user-specified value
  rownames(res) <- values
  return(res)
}
