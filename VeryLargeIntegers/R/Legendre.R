

# Legendre's formula (Given an integer n and a prime p, find the largest x such that p^x divides n!)

# not for the user
Legendrebase <- function(n, p){
  out = .pkgenv$zero
  while ( neqC(n, .pkgenv$zero) ){
    n = divbaseC(n, p)[[1]]
    out = sumC(out, n)
  }
  out
}

#' @title Legrendre's Formula for vli Objects
#' @author Javier Leiva Cuadrado
#' @param n a positive integer; object of class vli or 32 bits integer
#' @param p a prime number; object of class vli or 32 bits integer
#' @return object of class vli
#' @description Given a positive integer \code{n} and a prime \code{p}, the Legendre's Formula finds the largest integer \code{x} such that \code{p^x} divides the factorial of \code{n}, \code{n!}.
#' @examples p <- as.vli(577)
#' is.prime(p)
#' Legendre(12222, p)
#' @name 13. Legrendre's Formula
#' @rdname Legendre
#' @export Legendre
#'
Legendre <- function(n, p) UseMethod("Legendre")

#' @rdname Legendre
#' @method Legendre default
#' @export Legendre default
#'
Legendre.default <- function(n, p) stop("n and p have to be passed as vli objects or 32 bits integers")

#' @rdname Legendre
#' @method Legendre numeric
#' @export Legendre numeric
#'
Legendre.numeric <- function(n, p){
  if ( abs(n) < 2147483648 ){
    n = vliC(toString(n))
  }
  else stop("The n object passed as argument is neither a vli object nor a 32 bits integer")
  if ( !is.vli(p) ){
    if ( is.numeric(p) & (abs(p) < 2147483648) ){
      p = vliC(toString(p))
    }
    else stop("The p object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if ( !is.primeMRbase(p, iter=100) ) stop("The p object passed as argument is not a prime number")
  Legendrebase(n, p)
}

#' @rdname Legendre
#' @method Legendre vli
#' @export Legendre vli
#'
Legendre.vli <- function(n, p){
  if ( !is.vli(p) ){
    if ( is.numeric(p) & (abs(p) < 2147483648) ){
      p = vliC(toString(p))
    }
    else stop("The p object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if ( !is.primeMRbase(p, iter=100) ) stop("The p object passed as argument is not a prime number")
  Legendrebase(n, p)
}
