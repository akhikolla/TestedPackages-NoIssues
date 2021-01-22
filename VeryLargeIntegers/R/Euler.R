
# Euler's phi function

# not for the user
phibase <- function(x){
  out = x
  n = x
  p = .pkgenv$two
  while ( leqC( mulbaseC(p, p) , x ) ){
    if ( eqC( (divbaseC(n, p)[[2]]), .pkgenv$zero ) ){
      while ( eqC( (divbaseC(n, p)[[2]]), .pkgenv$zero ) ){
        n = divbaseC(n, p)[[1]]
      }
      out = subC( out, (divbaseC(out, p)[[1]]) )
    }
    p = sumC(p, .pkgenv$one)
  }
  if ( gtC(n, .pkgenv$one) ){
    out = subC(out, (divbaseC(out, n)[[1]]) )
  }
  return(out)
}

#' @title Euler's Phi Function for vli Objects
#' @author Javier Leiva Cuadrado
#' @param x positive integer; object of class vli or 32 bits integer
#' @return object of class vli
#' @description Euler's Phi Function for vli (Very Large Integers) objects. Given a positive integer \code{x}, the Euler's Phi Function returns the number of positive integers up to \code{x} that are relatively prime to \code{x}.
#' @details The returned value by the \code{phi} function is equal to the order of the group of units of the ring \code{Z/Zn} (the multiplicative group of integers modulo \code{n}). It is also called Euler's Totient Function, and plays a major part in Number Theory and in the RSA Cryptosystem.
#' @examples
#' \dontrun{
#' x <- as.vli("24352")
#' phi(x)
#' }
#' @name 17. Euler's phi function
#' @rdname Euler
#' @export phi
#'
phi <- function(x) UseMethod("phi")

#' @rdname Euler
#' @method phi default
#' @export phi default
#'
phi.default <- function(x) stop("The object passed as argument is neither a vli object nor a 32 bits integer")

#' @rdname Euler
#' @method phi numeric
#' @export phi numeric
#'
phi.numeric <- function(x){
  if ( abs(x) < 2147483648 ){
    if ( x >= 0 ){
      x = vliC(toString(x))
    }
    else stop("x has to be a positive integer")
  }
  else stop("The object passed as argument is neither a vli object nor a 32 bits integer")
  phibase(x)
}

#' @rdname Euler
#' @method phi vli
#' @export phi vli
#'
phi.vli <- function(x){
  if ( x$sign == -1 ) stop("x has to be a positive integer")
  phibase(x)
}
