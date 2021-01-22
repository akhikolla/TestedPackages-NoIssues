#'

#' Moore-Penrose generalized inverse
#'
#' @param X matrix
#' @param tol tolerance for determining bad entries
#' @examples
#' # create a positive definite 5x5 matrix
#' x <- crossprod(matrix(rnorm(25),5))
#' # make it singular
#' x[,2] <- x[,3]+x[,5]
#' geninv(x)
#' @return A matrix of the same dimension as \code{X} is returned, the Moore-Penrose generalized inverse.
#' @export
geninv <- function(X, tol=.Machine$double.eps^(2/3)) {
  stopifnot(is.numeric(X), length(dim(X)) == 2, is.matrix(X))
  nm <- colnames(X)
  s <- svd(X)
  p <- s$d > max(tol * s$d[1L], 0)
  inv <- if (all(p)) 
    s$v %*% (1/s$d * t(s$u))
  else if (!any(p)) 
    array(0, dim(X)[2L:1L])
  else {
    s$v[, p, drop = FALSE] %*% ((1/s$d[p]) * 
                                t(s$u[, p, drop = FALSE]))
  }
#  badvars <- if(any(!p)) names(collin(X))
  structure(inv, dimnames=list(nm,nm))
}

nazero <- function(x) ifelse(is.na(x),0,x)


#' Extract the mixed proportional hazard distribution
#'
#' @description
#' Various functions for extracting the proportional hazard distribution.
#'
#' \code{mphdist} extracts the hazard distribution.
#' @param pset
#' a parameter set of class \code{"mphcrm.pset"}, typically
#' \code{opt[[1]]$par}, where \code{opt} is returned from \code{\link{mphcrm}}.
#' If given a list of results, extracts the
#' first in the list.
#' @return
#' A matrix.
#' @examples
#' # load a dataset and a precomputed fitted model
#' data(durdata)
#' best <- fit[[1]]
#' mphdist(best)
#' mphmoments(best)
#' mphcov.log(best)
#' @export
mphdist <- function(pset) {
  if(inherits(pset,'mphcrm.list')) pset <- pset[[1]]$par
  if(inherits(pset,'mphcrm.opt')) pset <- pset$par
  mus <- exp(sapply(pset$parset, function(pp) pp$mu))
  if(is.null(colnames(mus))) {
    mus <- t(mus)
    colnames(mus) <- names(pset$parset)
  }
  rownames(mus) <- sprintf('point %2d',seq_len(nrow(mus)))
  prob <- a2p(pset$pargs)
  cbind(prob,mus)
}

#' Extract the mixed proportional log hazard distribution
#' @rdname mphdist
#' @description \code{mphdist.log} extracts the log hazard distribution.
#' @export
mphdist.log <- function(pset) {
  if(inherits(pset,'mphcrm.list')) pset <- pset[[1]]$par
  if(inherits(pset,'mphcrm.opt')) pset <- pset$par
  mus <- sapply(pset$parset, function(pp) pp$mu)
  if(is.null(colnames(mus))) {
    mus <- t(mus)
    colnames(mus) <- names(pset$parset)
  }
  rownames(mus) <- sprintf('point %2d',seq_len(nrow(mus)))
  prob <- a2p(pset$parg)
  cbind(prob,mus)
}

#' Extract moments of the mixed proportional hazard distribution
#' @rdname mphdist
#' @description
#' \code{mphmoments} returns the first and second moments of the hazard distribution.
#' @export
mphmoments <- function(pset) {
  dist <- mphdist(pset)
  mean <- colSums(dist[,1]*dist[,-1,drop=FALSE])
  if(length(mean) == 1) 
    variance <- sum(dist[,1]*(dist[,-1]-mean)^2) 
  else 
    variance <- rowSums(apply(dist, 1, function(x) x[1]*(x[-1]-mean)^2))
  sd <- sqrt(variance)
  cbind(mean,variance,sd)
}
#' Extract moments of the mixed proportional log hazard distribution
#' @rdname mphdist
#' @description
#' \code{mphmoments.log} returns the first and second moments of the log hazard distribution.
#' @export
mphmoments.log <- function(pset) {
  dist <- mphdist.log(pset)
  mean <- colSums(dist[,1]*dist[,-1,drop=FALSE])
  if(length(mean) == 1) 
    variance <- sum(dist[,1]*(dist[,-1]-mean)^2) 
  else
    variance <- rowSums(as.matrix(apply(dist, 1, function(x) x[1]*(x[-1]-mean)^2)))
  sd <- sqrt(variance)
  cbind(mean,variance,sd)
}

#' Extract covariance matrix of the proportional hazard distribution
#' @rdname mphdist
#' @description
#' \code{mphcov} returns the variance/covariance matrix of the hazard distribution.
#' @export
mphcov <- function(pset) {
  dist <- mphdist(pset)
  mean <- colSums(dist[,1]*dist[,-1,drop=FALSE])
  crossprod(sqrt(dist[,1])*(dist[,-1,drop=FALSE]-rep(mean,each=nrow(dist))))
}

#' Extract medians of the proportional hazard distribution
#' @rdname mphdist
#' @description
#' \code{mphmedian} returns the medians of the hazard distribution.
#' @export
mphmedian <- function(pset) {
  dist <- mphdist(pset)
  sapply(colnames(dist)[-1], function(tr) {
    if(nrow(dist)==1) return(dist[1,tr])
    oo <- order(dist[,tr])
    do <- dist[oo,,drop=FALSE]
    cp <- cumsum(do[,'prob'])
    above <- which(cp > 0.5)[1]
    if(above == 1) return(do[1,tr])
    p <- (cp[above]-0.5)/(cp[above]-cp[above-1])
    as.numeric((1-p)*do[above,tr] + p*do[above-1,tr])
  })
}



#' Extract covariance matrix of the proportional hazard distribution
#' @rdname mphdist
#' @description
#' \code{mphcov.log} returns the variance/covariance matrix of the log hazard distribution.
#' @export
mphcov.log <- function(pset) {
  dist <- mphdist.log(pset)
  mean <- colSums(dist[,1]*dist[,-1,drop=FALSE])
  crossprod(sqrt(dist[,1])*(dist[,-1,drop=FALSE]-rep(mean,each=nrow(dist))))
}

#' Extract standard errors of the estimated parameters
#' @param x
#' The Fisher matrix, typically from \code{opt[[1]]$fisher}, where \code{opt} is returned
#' from \code{\link{mphcrm}}.
#' @param tol tolerance for \link{geninv}
#' @return A named vector of standard errors is returned.
#' @export
se <- function(x,tol=.Machine$double.eps) {
  if(is.matrix(x)) return(sqrt(diag(geninv(x,tol))))
  if(!is.null(x$fisher)) return(sqrt(diag(geninv(x$fisher,tol))))
  stop("Can't find a matrix to invert")
}

#' Calculate various pseudo R2's 
#'
#' @description
#' There are several variants of pseudo \eqn{R^2} that can be computed for a likelihood
#' estimation. They all relate the log likelihood of the estimated model to the log likelihood
#' of the null model.
#'
#' The ones included here are McFadden's, Adjusted McFadden's, Cox & Snell's, and Nagelkerke, Cragg, and Uhler's.
#' 
#' @param opt returned value from \code{\link{mphcrm}}.
#' @return
#' A matrix is returned, with one row for each iteration containing the various pseudo \eqn{R^2}s.
#' @export
pseudoR2 <- function(opt) {
  nullmod <- opt[['nullmodel']]
  t(sapply(opt[-length(opt)], function(op) {
    c(mcfadden=1-op$value/nullmod$value,
      adjmcfadden=1-(op$value - length(unlist(op$par)))/nullmod$value,
      coxsnell=1-exp((nullmod$value-op$value)*2/op$nspells),
      ncu=(1-exp((nullmod$value-op$value)*2/op$nspells))/(1-exp(nullmod$value*2/op$nspells)))
  }))
}

#' Prettyprint a time interval
#' @description
#' Converts a time in seconds to a short string e.g. \code{"3m4s"}.
#' @param t
#' numeric. time in seconds.
#' @return A character string is returned.
#' @examples
#' timestr(1.3)
#' timestr(73)
#' timestr(4684)
#' @export
timestr <- function(t) {
  dec <- t - as.integer(t)
  t <- as.integer(t-dec)
  s <- t %% 60
  t <- t %/% 60
  m <- t %% 60
  h <- t %/% 60
  if(h > 0) {
    str <- sprintf('%dh%dm',h,m)
  } else if(m > 0) {
    str <- sprintf('%dm%.0fs',m,s)
  } else {
    str <- sprintf('%.1fs',s+dec)
  }
  str
}

#' Collapse levels of a factor
#' @description
#' Combines levels of a factor into new levels
#' @param f
#' factor.
#' @param newlevels
#' list. The names of \code{newlevels} are the new levels. Each list element is a list
#' of old levels in the factor \code{f} which should be combined into the new level
#' @examples
#' # create a factor with levels 30:60
#' age <- factor(sample(30:60, 200, replace=TRUE))
#' # combine 35-40 into a single level, 41-50 into a single level, and 51-60 into a single level
#  # levels 30-34 are left as is. Note the backticks, necessary because these must be parsed as names.
#' g <- smashlevels(age, list(`35-40` = 35:40, `41-50` = 41:50, `51-60` = 51:60))
#' table(g)
#' # If the syntax permits, the backticks can be avoided.
#' h <- smashlevels(age, list(young=30:34, pushing40 = 35:40, pushing50 = 41:50, fossilized = 51:120))
#' table(h)
#' @export
smashlevels <- function(f, newlevels) {
  f <- as.factor(f)
  olev <- levels(f)
  for(nl in names(newlevels)) olev[match(newlevels[[nl]], olev)] <- nl
  levels(f) <- olev
  f
}
