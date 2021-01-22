
# binomial coefficients

# not for the user
binombase = function(n, k){
  n_k = subC(n, k)
  if ( gtC(k, n_k) ){
    k = n_k
  }
  rm(n_k)
  binomC(n, k)
}

#' @title Binomial Coefficients for vli Objects
#' @author Javier Leiva Cuadrado
#' @param n object of class vli or 32 bits integer
#' @param k object of class vli or 32 bits integer
#' @return object of class vli
#' @description \code{binom} computes binomial coefficients of vli (Very Large Integer) objects. That is, given two positive integers \code{n} and \code{k} with \code{n >= k}, the function \code{binom(n, k)} returns the number of ways to choose a subset of \code{k} elements, disregarding their order, from a set of \code{n} elements.
#' @examples
#' x <- as.vli("100")
#' binom(x, 20)
#' @name 06. Binomial coefficients
#' @rdname binom
#' @export binom
#'
binom <- function(n, k) UseMethod("binom")

#' @rdname binom
#' @method binom default
#' @export binom default
#'
binom.default <- function(n, k) stop("The x object passed as argument is neither a vli object nor a 32 bits integer")

#' @rdname binom
#' @method binom numeric
#' @export binom numeric
#'
binom.numeric <- function(n, k){
  if( (abs(n) < 2147483648) & (n > 0) ){
    n = vliC(toString(n))
    if ( ltC(n, .pkgenv$one) ) stop("binom is only defined for positive integer values")
    if ( !is.vli(k) ){
      if( (abs(k) < 2147483648) & (k > 0) ) k = vliC(toString(k))
      else stop("binom is only defined for positive integer values")
    }
    else if ( ltC(k, .pkgenv$one) ) stop("loge is only defined for positive integer values")
    if ( gtC(k, n) ) stop("n has to be greater or equal to k")
    binombase(n, k)
  }
  else stop("binom is only defined for positive integer numbers")
}

#' @rdname binom
#' @method binom vli
#' @export binom vli
#'
binom.vli <- function(n, k){
  if ( ltC(n, .pkgenv$one) ) stop("binom is only defined for positive integer values")
  if ( !is.vli(k) ){
    if( (abs(k) < 2147483648) & (k > 0) ) k = vliC(toString(k))
    else stop("binom is only defined for positive integer values")
  }
  else if ( ltC(k, .pkgenv$one) ) stop("binom is only defined for positive integer values")
  if ( gtC(k, n) ) stop("n has to be greater or equal to k")
  binombase(n, k)
}
