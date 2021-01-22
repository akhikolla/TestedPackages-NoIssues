

# greatest common divisor

# not for the user
gcdbase <- function(x, y){
  while ( neqC(x, .pkgenv$zero) ){
    z = x
    x = divbaseC(y, x)[[2]]
    y = z
    rm(z)
  }
  y
}

#' @title Greatest Common Divisor for vli Objects
#' @author Javier Leiva Cuadrado
#' @param x object of class vli or 32 bits integer
#' @param y object of class vli or 32 bits integer
#' @return object of class vli
#' @description \code{gcd} computes and returns the greatest common divisor of two vli (Very Large Integers) objects.
#' @examples x <- as.vli("1225312091263347514461245")
#' y <- as.vli("357590484262521")
#' gcd(x, y)
#' @name 09. Greatest common divisor
#' @rdname gcd
#' @export gcd
#'
gcd <- function(x, y) UseMethod("gcd")

#' @rdname gcd
#' @method gcd default
#' @export gcd default
#'
gcd.default <- function(x, y) stop("x and y have to be specified as vli objects or 32 bits integers")

#' @rdname gcd
#' @method gcd numeric
#' @export gcd numeric
#'
gcd.numeric <- function(x, y){
  if ( abs(x) < 2147483648 ){
    if ( x >= 0 ){
      x = vliC(toString(x))
    }
    else stop("gcd is only defined for positive integer numbers")
  }
  else stop("The first object passed as argument is neither a vli object nor a 32 bits integer")
  if ( !is.vli(y) ){
    if ( is.numeric(y) & (abs(y) < 2147483648) ){
      if ( y >= 0 ){
        y = vliC(toString(y))
      }
      else stop("gcd is only defined for positive integer numbers")
    }
    else stop("The second object passed as argument is neither a vli object nor a 32 bits integer")
  }
  else if ( y$sign == -1 ) stop("gcd is only defined for positive integer numbers")
  return(gcdbase(x,y))
}

#' @rdname gcd
#' @method gcd vli
#' @export gcd vli
#'
gcd.vli <- function(x, y){
  if ( x$sign == -1 ) stop("gcd is only defined for positive integer numbers")
  if ( !is.vli(y) ){
    if ( is.numeric(y) & (abs(y) < 2147483648) ){
      if ( y >= 0 ){
        y = vliC(toString(y))
      }
      else stop("gcd is only defined for positive integer numbers")
    }
    else stop("The second object passed as argument is neither a vli object nor a 32 bits integer")
  }
  else if ( y$sign == -1 ) stop("gcd is only defined for positive integer numbers")
  return(gcdbase(x,y))
}
