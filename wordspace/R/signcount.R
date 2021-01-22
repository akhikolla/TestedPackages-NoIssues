signcount <- function (x, what=c("counts", "nonneg", "nnzero")) {
  what <- match.arg(what)
  if (is.double(x)) {
    res <- CPP_signcount(x)
  }
  else if (is.integer(x)) {
    res <- CPP_signcount_int(x)
  }
  else if (is(x, "dgeMatrix")) {
    res <- CPP_signcount(x@x)
  }
  else if (is(x, "dgCMatrix") || is(x, "dgRMatrix")) {
    n.val <- prod(dim(x))
    res <- CPP_signcount(x@x)
    res[2] <- res[2] + n.val - sum(res) # add structural zeroes to count
  }
  else {
    stop("'x' must be a numeric vector or matrix, a dense Matrix or a sparseMatrix in compressed representation")
  }
  ## <res> now contains c(pos, zero, neg)
  switch(
    what,
    counts = structure(res, names=c("pos", "zero", "neg")),
    nonneg = (res[3] == 0),
    nnzero = res[1] + res[3],
    stop("internal error -- unsupported value for what= argument"))
}
