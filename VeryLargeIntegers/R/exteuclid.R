

# extended euclidean algorithm

# not for the user
exteuclidbase <- function(x, y){
  s = .pkgenv$zero
  t = .pkgenv$one
  r = y
  old_s = .pkgenv$one
  old_t = 0
  old_r = x
  while ( neqC(r, .pkgenv$zero) ){
    qu = divbaseC(old_r, r)[[1]]
    z = r
    r = old_r - (qu * r)
    old_r = z
    rm(z)
    z = s
    s = old_s - (qu * s)
    old_s = z
    rm(z)
    z = t
    t = old_t - (qu * t)
    old_t = z
    rm(z)
  }
  out = vli(3)
  out[[1]] = old_r # gcd
  out[[2]] = old_s # first Bezout coefficent
  out[[3]] = old_t # second Bezout coefficent
  return(out)
}

#' @title Extended Euclidean Algorithm for vli Objects
#' @author Javier Leiva Cuadrado
#' @param x object of class vli or 32 bits integer
#' @param y object of class vli or 32 bits integer
#' @return list of 3 objects of class vli: the first is the greatest common divisor of \code{x} and \code{y}, and the other two are the Bezout's coefficients
#' @description Computation of the Extended Euclidean algorithm for vli (Very Large Integers) objects. Given two positive integers, \code{x} and \code{y}, the Extended Euclidean algorithm looks for two integers \code{a} and \code{b} (called Bezout's coefficients) such that \code{(a * x) + (b * y) = 1}. To do this, the algorithm needs to compute the greatest common divisor of \code{x} and \code{y}, so it is also returned by the function.
#' @details The returned object is a list of 3 elements. To access the numbers, it is necessary to use the list operator \code{[[i]]}, where "\code{i}" has to be 1 for the greatest common divisor, 2 for the first Bezout coefficient and 3 for the second Bezout coefficient (see the example).
#' @examples x <- as.vli("232636113097")
#' y <- as.vli("52442092785616")
#' result <- exteuclid(x, y)
#' ( result[[2]] * x ) + ( result[[3]] * y )
#' @name 11. Extended Euclidean algorithm
#' @rdname exteuclid
#' @export exteuclid
#'
exteuclid <- function(x, y) UseMethod("exteuclid")

#' @rdname exteuclid
#' @method exteuclid default
#' @export exteuclid default
#'
exteuclid.default <- function(x, y) stop("x and y have to be specified as vli objects or 32 bits integers")

#' @rdname exteuclid
#' @method exteuclid numeric
#' @export exteuclid numeric
#'
exteuclid.numeric <- function(x, y){
  if ( abs(x) < 2147483648 ){
    if ( x >= 0 ){
      x = vliC(toString(x))
    }
    else stop("exteuclid is only defined for positive integer numbers")
  }
  else stop("The first object passed as argument is neither a vli object nor a 32 bits integer")
  if ( !is.vli(y) ){
    if ( is.numeric(y) & (abs(y) < 2147483648) ){
      if ( y >= 0 ){
        y = vliC(toString(y))
      }
      else stop("exteuclid is only defined for positive integer numbers")
    }
    else stop("The second object passed as argument is neither a vli object nor a 32 bits integer")
  }
  else if ( y$sign == -1 ) stop("exteuclid is only defined for positive integer numbers")
  return(exteuclidbase(x,y))
}

#' @rdname exteuclid
#' @method exteuclid vli
#' @export exteuclid vli
#'
exteuclid.vli <- function(x, y){
  if ( x$sign == -1 ) stop("exteuclid is only defined for positive integer numbers")
  if ( !is.vli(y) ){
    if ( is.numeric(y) & (abs(y) < 2147483648) ){
      if ( y >= 0 ){
        y = vliC(toString(y))
      }
      else stop("exteuclid is only defined for positive integer numbers")
    }
    else stop("The second object passed as argument is neither a vli object nor a 32 bits integer")
  }
  else if ( y$sign == -1 ) stop("exteuclid is only defined for positive integer numbers")
  return(exteuclidbase(x,y))
}
