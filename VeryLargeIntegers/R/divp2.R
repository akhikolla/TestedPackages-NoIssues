

# x vli or integer, k positive integer, return vli x/(2^k)

#' @title Efficient Division by a Power of 2
#' @author Javier Leiva Cuadrado
#' @param x dividend; object of class vli or 32 bits integer
#' @param k exponent of the divisor (the divisor will be \code{2^k}); 32 bits integer
#' @return object of class vli
#' @description \code{divp2} efficiently divides an object of class vli by a power of 2.
#' @details Given two integers \code{x} (vli or 32 bits integer) and \code{k} (32 bits integer), the function \code{divp2(x, k)} computes and returns \code{x/(2^k)} as an object of class vli.
#' @examples # Dividing a random 500 digits integer by 2^10 = 1024
#' x <- rvlidigits(500)
#' x
#' divp2(x, 10)
#' @name 05. Efficent division by a power of 2
#' @rdname divp2
#' @method divp2 vli
#' @export divp2 vli
#'
divp2 <- function(x, k) UseMethod("divp2")

#' @rdname divp2
#' @method divp2 default
#' @export divp2 default
#'
divp2.default = function(x, k) stop("x must be a 32 bits integer or a vli class object and k must be a positive integer")

#' @rdname divp2
#' @method divp2 numeric
#' @export divp2 numeric
#'
divp2.numeric = function(x, k){
  if ( abs(x) < 2147483648 ){
    x = vliC(toString(x))
  }
  else stop("The x object passed as argument is neither a vli object nor a 32 bits integer")
  if ( !is.numeric(k) ) stop("The exponent k has to be a 32 bits integer number")
  else if ( !(abs(k) < 2147483648) ) stop("The exponent k has to be a 32 bits integer number")
  else{
    if ( k < 0 ) stop("divp2 is only defined for positive integer exponents k")
  }
  x$value = divp2C(x$value, as.integer(abs(k)))
  if ( eqC(x, .pkgenv$zero) ) x$sign = 1
  else if( x$sign == -1 ) x = subC(x, .pkgenv$one)
  x
}

#' @rdname divp2
#' @method divp2 vli
#' @export divp2 vli
#'
divp2.vli = function(x,k){
  if ( !is.numeric(k) ) stop("The exponent k has to be a 32 bits integer number")
  else if ( !(abs(k) < 2147483648) ) stop("The exponent k has to be a 32 bits integer number")
  else{
    if ( k < 0 ) stop("divp2 is only defined for positive integer exponents k")
  }
  x$value = divp2C(x$value, as.integer(abs(k)))
  if ( eqC(x, .pkgenv$zero) ) x$sign = 1
  else if ( x$sign == -1 ) x = subC(x, .pkgenv$one)
  x
}
