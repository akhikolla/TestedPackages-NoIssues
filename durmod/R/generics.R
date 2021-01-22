#' @method coef mphcrm.pset
#' @export
coef.mphcrm.pset <- function(object, ...) {
  pars <- unlist(lapply(object$parset, function(pp) {
    c(pp$pars,unlist(pp$facs))
  }))
  mus <- unlist(lapply(object$parset, function(pp) {
    pp$mu
  }))
  probs <- a2p(object$pargs)
  names(probs) <- paste('P',seq_along(probs),sep='')
  c(pars,mus,probs)
}



# There are four classes:
# mphcrm.list, a list as returned from mphcrm
# each element is of class mphcrm.opt
# The mphcrm.opt contains a par-element of class mphcrm.pset
#


#' @method print mphcrm.list
#' @export
print.mphcrm.list <- function(x,...) {
  for(i in seq_along(x)) {
    oo <- x[[i]]
    N <- length(oo$par$pargs)+1
    ll <- oo$value
    nm <- names(x)[i]
    cat(sprintf('%s: estimate with %d points, log-likelihood: %.4f\n', nm, N,ll))
    if(i == 1) {cat('\n');print(oo[[1]]); cat('\n')}
  }
}

#' @method print mphcrm.opt
#' @export 
print.mphcrm.opt <- function(x,...) {
  N <- length(x$par$pargs)+1
  ll <- x$value
  cat(sprintf('Estimate with %d points, log-likelihood: %.4f\n', N, ll))
  print(x$par)
}

#' @method summary mphcrm.opt
#' @export
summary.mphcrm.opt <- function(object, ...,tol=.Machine$double.eps) {
  val <- flatten(object$par)
  dist <- grep('(\\.mu[0-9]+|pargs[0-9]*)$',names(val))
  se <- if(!is.null(object$fisher)) se(object,tol)[-dist] else NA
  tval <- val[-dist]/se
  rdf <- object$nobs
  list(loglik=object$value,
       coefs=cbind(value=val[-dist], 
                   se=se,
                   t=tval,
                   `Pr(>|t|)`=2*pt(abs(tval), rdf,lower.tail=FALSE)),
       moments=mphmoments(object$par))
}

#' @method summary mphcrm.pset
#' @export
summary.mphcrm.pset <- function(object,...) {
  val <- flatten(object)
  dist <- grep('(\\.mu[0-9]+|pargs[0-9]*)$',names(val))
  list(value = val[-dist], 
       mphdist=mphdist(object))
}

#' @method print mphcrm.pset
#' @export
print.mphcrm.pset <- function(x, ...) {
  val <- flatten(x)
  dist <- grep('(\\.mu[0-9]+|pargs[0-9]*)$',names(val))
  print(val[-dist])
  cat('\nProportional hazard distribution\n')
  print(round(mphdist(x),8))
}

#' Convert a structured coefficient set to a vector
#' @description \code{\link{mphcrm}} stores coefficients in a list,
#'   not in a vector. This is because they should be treated
#'   differently according to whether they are probabilities,
#'   proportional hazards, or coefficients for factor levels or
#'   ordinary covariates.  \code{flatten} extracts them as a named
#'   vector. \code{unflatten} puts them back in structured form.
#' @details \code{flatten}/\code{unflatten} is just a thinly disguised
#'   \code{\link[base]{unlist}}/\code{\link[utils]{relist}}, but uses
#'   slightly more readable names.
#' @param x parameter set as typically found in \code{opt[[1]]\$par},
#'   where \code{opt} is returned from \code{mphcrm}.
#' @param exclude For internal use
#' @export
flatten <- function(x, exclude=attr(x,'exclude')) {
  class(x) <- setdiff(class(x),'mphcrm.pset')
  vec <- unlist(as.relistable(x))
  # fix the names so more readable
  # parset.t1.{pars,mu,facs}.x <- t1.x
  newnames <- gsub('^parset\\.(.*)\\.(pars|mu|facs)\\.(.*)','\\1.\\3',names(vec))
  names(vec) <- newnames

  skel <- attr(vec,'skeleton')
  class(skel) <- c('mphcrm.pset',class(skel))
  attr(vec,'skeleton') <- skel
  if(length(exclude) == 0) return(vec)
  # remove the excluded parameters
  exc <- sapply(exclude, function(ee) ee[[1]])
  if(is.character(exc)) exc <- which(names(vec) %in% exc)
  exc <- intersect(exc,seq_along(vec))
  if(length(exc)==0) return(structure(vec,exclude=NULL))
  vals <- vec[exc]
  vec <- vec[-exc]
  attr(vec,'skeleton') <- skel
  exc <- structure(mapply(function(a,b) list(a,b), exc, vals,SIMPLIFY=FALSE), names=names(vals))
  structure(vec,exclude=exc)
}


#' @rdname flatten
#' @param flesh
#' vector of class \code{"relistable"}, as returned from \code{\link{flatten}}.
#' @param skeleton
#' For internal use
#' @export
unflatten <- function(flesh, skeleton=attr(flesh, 'skeleton'), exclude=attr(flesh,'exclude')) {
  class(skeleton) <- setdiff(class(skeleton),'mphcrm.pset')
  ne <- length(exclude)
  if(ne == 0) {
    pset <- relist(flesh,skeleton)
    class(pset) <- c('mphcrm.pset',class(pset))
    return(pset)
  }

  # insert the excluded parameters before relisting
  exc <- sapply(exclude, function(ee) ee[[1]])
  vals <- sapply(exclude, function(ee) ee[[2]])
  vec <- numeric(length(flesh)+ne)
  idx <- seq_along(vec)
  ok <- setdiff(idx,exc)
  vec[ok] <- flesh
  vec[exc] <- vals
  names(vec)[ok] <- names(flesh)
  names(vec)[exc] <- names(vals)
  pset <- relist(vec,skeleton)
  class(pset) <- c('mphcrm.pset',class(pset))
  structure(pset,exclude=as.list(exc))
}

exclude <- function(x,exclude,...) {
  # the generic doesn't have a ... argument, so call it directly
  relist(getS3method('unlist','mphcrm.pset')(x,exclude=exclude))
}

#' @method logLik mphcrm.opt
#' @param useobs Use number of observations for computing degrees of freedom, not number of spells.
##' @export
logLik.mphcrm.opt <- function(object, ..., useobs=FALSE) {
  val <- object$value
  attr(val, 'nall') <- object$nobs
  attr(val, 'nobs') <- object$nspells 
  if(useobs) N <- object$nobs else N <- object$nspells
  attr(val, 'df') <- length(flatten(object$par))
  structure(val, class='logLik')
}

#' @method logLik mphcrm.list
##' @export
logLik.mphcrm.list <- function(object,..., useobs=FALSE) {
  logLik(object[[1L]],...,useobs)
}

#' @method coef mphcrm.pset
##' @export
coef.mphcrm.pset <- function(object,...)  structure(flatten(object),skeleton=NULL)

#' @method coef mphcrm.opt
##' @export
coef.mphcrm.opt <- function(object,...) coef(object$par)

#' @method vcov mphcrm.opt
##' @export
vcov.mphcrm.opt <- function(object,...)  geninv(object$fisher)
