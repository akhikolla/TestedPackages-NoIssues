
# Jacobi Symbol

# not for the user
jacobbase <- function(a, n){
  out = .pkgenv$one
  if ( eqC(a, .pkgenv$zero) ) return(.pkgenv$zero)
  if ( ltC(a, .pkgenv$zero) ){
    a$sign = ( -1 ) * a$sign
    if ( eqC(divbaseC(n, .pkgenv$four)[[2]], .pkgenv$three) ) out$sign = (-1) * out$sign
  }
  if ( eqC(a, .pkgenv$one) ) return(out)
  while ( neqC(a, .pkgenv$zero) ){
    if ( ltC(a, .pkgenv$zero) ){
      a$sign = ( -1 ) * a$sign
      if ( eqC(divbaseC(n, .pkgenv$four)[[2]], .pkgenv$three) ) out$sign = (-1) * out$sign
    }
    while ( eqC(divbaseC(a, .pkgenv$two)[[2]], .pkgenv$zero) ){
      a = divbaseC(a, .pkgenv$two)[[1]]
      if ( eqC(divbaseC(n, .pkgenv$eight)[[2]], .pkgenv$three) | eqC(divbaseC(n, .pkgenv$eight)[[2]], .pkgenv$five) ) out$sign = (-1) * out$sign
    }
    t = a
    a = n
    n = t
    rm(t)
    if ( eqC(divbaseC(a, .pkgenv$four)[[2]], .pkgenv$three) & eqC(divbaseC(n, .pkgenv$four)[[2]], .pkgenv$three) ) out$sign = (-1) * out$sign
    a = divbaseC(a, n)[[2]]
    if ( gtC(a, divbaseC(n, .pkgenv$two)[[1]]) ) a = subC(a, n)
  }
  if ( eqC(n, .pkgenv$one) ) return(out)
  return(.pkgenv$zero)
}

#' @title Computation of the Jacobi Symbol for vli Objects
#' @author Javier Leiva Cuadrado
#' @param a object of class vli or 32 bits integer
#' @param n positive odd integer; object of class vli or 32 bits integer
#' @return object of class vli with value \code{-1}, \code{0} or \code{1}.
#' @description Computation of the Jacobi Symbol for vli (Very Large Integers) objects. The Jacobi Symbol is a generalization of the Legendre Symbol, not being necessary that \code{n} be a prime number.
#'
#' It is needed in many algorithms of modular arithmetic, computational number theory and cryptography. For example, it is used by the present package in the Solovay-Strassen probabilistic primality test.
#' @examples x <- as.vli("342635653456")
#' y <- as.vli("3210591001")
#' Jacobi(x, y)
#' @name 16. Jacobi Symbol
#' @rdname Jacobi
#' @export Jacobi
#'
Jacobi <- function(a, n) UseMethod("Jacobi")

#' @rdname Jacobi
#' @method Jacobi default
#' @export Jacobi default
#'
Jacobi.default <- function(a, n) stop("Both arguments have to be passed as vli objects or 32 bits integers")

#' @rdname Jacobi
#' @method Jacobi numeric
#' @export Jacobi numeric
#'
Jacobi.numeric <- function(a, n){
  if ( abs(a) < 2147483648 ){
    a = vliC(toString(a))
  }
  else stop("The first object passed as argument is neither a vli object nor a 32 bits integer")
  if ( !is.vli(n) ){
    if ( is.numeric(n) & (abs(n) < 2147483648) ){
      n = vliC(toString(n))
      if ( is.even(n) | leqC(n, .pkgenv$zero)) stop("The Jacobi Symbol is only defined for positive odd values of n")
    }
    else stop("The second object passed as argument is neither a vli object nor a 32 bits integer")
  }
  jacobbase(a, n)
}

#' @rdname Jacobi
#' @method Jacobi vli
#' @export Jacobi vli
#'
Jacobi.vli <- function(a, n){
  if ( !is.vli(n) ){
    if ( is.numeric(n) & (abs(n) < 2147483648) ){
      n = vliC(toString(n))
      if ( is.even(n) | leqC(n, .pkgenv$zero) ) stop("The Jacobi Symbol is only defined for positive odd values of n")
    }
    else stop("The second object passed as argument is neither a vli object nor a 32 bits integer")
  }
  jacobbase(a, n)
}
