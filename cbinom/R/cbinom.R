dcbinom <- function(x, size, prob, log = FALSE){
  if (!is.numeric(x) || !is.numeric(size) || !is.numeric(prob))
    stop("Non-numeric argument to mathematical function")
  fmat <- dcblp(x, size, prob)
  # columns of fmat are:
  #  log(pcbinom(x + h)), log(pcbinom(x - h)), 2*h if h < x < size + 1 - h
  #  log(pcbinom(x + h)), log(pcbinom(x)), h if x <= h
  #  log(pcbinom(x), log(pcbinom(x - h), h if x > size + 1 - h
  #   where h is small
  suppressWarnings({
    if (!log) {
      ans <- (exp(fmat[, 1]) - exp(fmat[, 2]))/fmat[,3]
      ans[ans < 0] <- 0
    } else {
      ans <- log((fmat[, 1] - fmat[, 2])/fmat[, 3]) +
        pcbinom(x, size, prob, log.p = T)
      ans[is.na(ans)] <- -Inf
    }
  })
  return(ans)
}
pcbinom <- function(q, size, prob, lower.tail = TRUE, log.p = FALSE){
  if (
    suppressWarnings(sum(is.na(q)) != sum(is.na(as.numeric(q)))) ||
    suppressWarnings(sum(is.na(size)) != sum(is.na(as.numeric(size)))) ||
    suppressWarnings(sum(is.na(prob)) != sum(is.na(as.numeric(prob))))
  ) stop("Non-numeric argument to mathematical function")
  if (is.null(lower.tail) || is.na(lower.tail)){
    lower.tail = TRUE
  } else {
    lower.tail <- as.logical(lower.tail)
    if (is.na(lower.tail[1])){
      lower.tail <- TRUE
      warning("NAs introduced by coercion")
    }
  }
  lower.tail <- lower.tail[1]
  if (is.null(log.p)){
    log.p = TRUE
  } else if (is.na(log.p)){
    log.p = TRUE
  } else {
    log.p <- as.logical(log.p)
    if (is.na(log.p[1])){
      log.p <- TRUE
      warning("NAs introduced by coercion")
    }
  }
  log.p <- log.p[1]
  ans <- pcbinomC(q, size, prob, log.p)
  if (!lower.tail){
    if (!log.p){
      ans <- 1 - ans
    } else {
      ans <- log(1 - exp(ans))
    }
  }
  if (any(is.nan(ans))){
    warning("NaNs produced")
  }
  return(ans)
}
qcbinom <- function(p, size, prob, lower.tail = TRUE, log.p = FALSE){
  if (
    suppressWarnings(sum(is.na(p)) != sum(is.na(as.numeric(p)))) ||
    suppressWarnings(sum(is.na(size)) != sum(is.na(as.numeric(size)))) ||
    suppressWarnings(sum(is.na(prob)) != sum(is.na(as.numeric(prob))))
  ) stop("Non-numeric argument to mathematical function")
  if (is.null(lower.tail) || is.na(lower.tail)){
    lower.tail = TRUE
  } else {
    lower.tail <- as.logical(lower.tail)
    if (is.na(lower.tail[1])){
      lower.tail <- TRUE
      warning("NAs introduced by coercion")
    }
  }
  lower.tail <- lower.tail[1]
  if (is.null(log.p)){
    log.p = TRUE
  } else if (is.na(log.p)){
    log.p = TRUE
  } else {
    log.p <- as.logical(log.p)
    if (is.na(log.p[1])){
      log.p <- TRUE
      warning("NAs introduced by coercion")
    }
  }
  log.p <- log.p[1]
  suppressWarnings({
    if (!log.p){
      if (lower.tail) p <- log(p) else p <- log(1 - p)
    } else {
      if (!lower.tail) p <- log(1 - exp(p))
    }
  })
  ans <- qcbinomC(p, size, prob, rcb = F)
  if (any(is.nan(ans))) warning("NaNs produced")
  return(ans)
}
rcbinom <- function(n, size, prob){
  if (suppressWarnings(sum(is.na(size)) != sum(is.na(as.numeric(size)))) ||
      suppressWarnings(sum(is.na(prob)) != sum(is.na(as.numeric(prob)))))
    stop("invalid arguments")
  if (length(n) == 1){
    if(!is.finite(n) || n < 0) stop("invalid arguments")
  } else {
    n <- length(n)
  }
  u <- runif(n, 0, 1)
  ans <- qcbinomC(log(u), size, prob, rcb = T)
  if (any(is.na(ans))) warning("NAs produced")
  return(ans)
}
