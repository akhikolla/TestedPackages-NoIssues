
# Integer roots


# square root

#' @title Integer roots for vli Objects
#' @author Javier Leiva Cuadrado
#' @param x base of the root; object of class vli or 32 bits integer
#' @return object of class vli
#' @description Computation of integer roots and their remainders of vli (Very Large Integers) objects. Functions \code{sqrt} and \code{rootk} returns respectively the integer square root and the integer k-th root of the given value. Functions \code{sqrtrem} and \code{rootkrem} returns the corresponding remainder.
#' @examples x <- as.vli("4124135")
#' sqrt(x)
#' sqrtrem(x)
#' sqrt(x)^2 + sqrtrem(x) == x
#' \dontrun{
#' rootk(as.vli("1492346293864978561249785"), 5)
#' }
#' @name 03. Roots
#' @rdname roots
#' @method sqrt vli
#' @export sqrt vli
#'
sqrt.vli = function(x){
  if ( x$sign == -1 ) stop("sqrt is only defined for positive integer numbers")
  rootC(x)[[1]]
}

# square root remainder

#' @rdname roots
#' @export sqrtrem
#'
sqrtrem <- function(x) UseMethod("sqrtrem")

#' @rdname roots
#' @method sqrtrem default
#' @export sqrtrem default
#'
sqrtrem.default = function(x) stop("The x object passed as argument is neither a vli object nor a 32 bits integer")

#' @rdname roots
#' @method sqrtrem numeric
#' @export sqrtrem numeric
#'
sqrtrem.numeric = function(x){
  if ( (abs(x) < 2147483648) & (x >= 0) ){
    x = vliC(toString(x))
    rootC(x)[[2]]
  }
  else  stop("sqrtrem is only defined for positive integer numbers")
}

#' @rdname roots
#' @method sqrtrem vli
#' @export sqrtrem vli
#'
sqrtrem.vli = function(x){
  if( x$sign == -1 ) stop("sqrtrem is only defined for positive integer numbers")
  rootC(x)[[2]]
}

# k-th root

#' @rdname roots
#' @export rootk
#'
rootk <- function(x, k) UseMethod("rootk")


#' @rdname roots
#' @method rootk default
#' @export rootk default
#'
rootk.default = function(x, k) stop("rootk is only defined for vli objects or 32 bits integers")


#' @rdname roots
#' @param k index of the root; object of class vli or 32 bits integer
#' @method rootk numeric
#' @export rootk numeric
#'
rootk.numeric = function(x, k){
  if ( !is.vli(k) ){
    if ( is.numeric(k) & (abs(k) < 2147483648) ){
      if ( k > 0 ){
        k = vliC(toString(k))
      }
      else stop("rootk is only defined for positive integer values of k")
    }
    else stop("The k object passed as argument is neither a vli object nor a 32 bits integer")
  }
  else{
    if( (k$sign == -1) | eqC(k, .pkgenv$zero) ) stop("rootk is only defined for positive integer values of k")
  }
  s = 1
  if ( abs(x) < 2147483648 ){
    if ( x >= 0 ){
      x = vliC(toString(x))
    }
    else {
      if ( is.even(k) ) stop("rootk is only defined for positive integer numbers when k is even")
      else {s = -1; x = vliC(toString(x)); x$sign=1}
    }
  }
  else stop("x argument is neither a vli object nor a 32 bits integer")
  if ( x$sign == -1 ){
    if ( is.even(k) ) stop("rootk is only defined for positive integer numbers when k is even")
    else {s = -1; x$sign = 1}
  }
  return(vlivC(s, rootkC(x, k)[[1]]$value))
}

#' @rdname roots
#' @method rootk vli
#' @export rootk vli
#'
rootk.vli = function(x, k){
  if ( !is.vli(k) ){
    if ( is.numeric(k) & (abs(k) < 2147483648) ){
      if ( k > 0 ){
        k = vliC(toString(k))
      }
      else stop("rootk is only defined for positive integer values of k")
    }
    else stop("The k object passed as argument is neither a vli object nor a 32 bits integer")
  }
  else{
    if( (k$sign == -1) | eqC(k, .pkgenv$zero) ) stop("rootk is only defined for positive integer values of k")
  }
  s = 1
  if ( x$sign == -1 ){
    if ( is.even(k) ) stop("rootk is only defined for positive integer numbers when k is even")
    else {s = -1; x$sign = 1}
  }
  return( vlivC(s, rootkC(x, k)[[1]]$value) )
}


# k-th root remainder

#' @rdname roots
#' @export rootkrem
#'
rootkrem <- function(x, k) UseMethod("rootkrem")

#' @rdname roots
#' @method rootkrem default
#' @export rootkrem default
#'
rootkrem.default = function(x, k) stop("rootk is only defined for vli objects or 32 bits integers")

#' @rdname roots
#' @method rootkrem numeric
#' @export rootkrem numeric
#'
rootkrem.numeric = function(x, k){
  if ( !is.vli(k) ){
    if ( is.numeric(k) & (abs(k) < 2147483648) ){
      if ( k > 0 ){
        k = vliC(toString(k))
      }
      else stop("rootkrem is only defined for positive integer values of k")
    }
    else stop("The k object passed as argument is neither a vli object nor a 32 bits integer")
  }
  else{
    if( (k$sign == -1) | eqC(k, .pkgenv$zero) ) stop("rootk is only defined for positive integer values of k")
  }
  s = 1
  if ( abs(x) < 2147483648 ){
    if ( x >= 0 ){
      x = vliC(toString(x))
    }
    else {
      if ( is.even(k) ) stop("rootkrem is only defined for positive integer numbers when k is even")
      else {s = -1; x = vliC(toString(x)); x$sign=1}
    }
  }
  else stop("x argument is neither a vli object nor a 32 bits integer")
  if ( x$sign == -1 ){
    if ( is.even(k) ) stop("rootkrem is only defined for positive integer numbers when k is even")
    else {s = -1; x$sign = 1}
  }
  return(vlivC(s, rootkC(x, k)[[2]]$value))
}

#' @rdname roots
#' @method rootkrem vli
#' @export rootkrem vli
#'
rootkrem.vli = function(x, k){
  if ( !is.vli(k) ){
    if ( is.numeric(k) & (abs(k) < 2147483648) ){
      if ( k > 0 ){
        k = vliC(toString(k))
      }
      else stop("rootkrem is only defined for positive integer values of k")
    }
    else stop("The k object passed as argument is neither a vli object nor a 32 bits integer")
  }
  else{
    if( (k$sign == -1) | eqC(k, .pkgenv$zero) ) stop("rootkrem is only defined for positive integer values of k")
  }
  s = 1
  if ( x$sign == -1 ){
    if ( is.even(k) ) stop("rootkrem is only defined for positive integer numbers when k is even")
    else {s = -1; x$sign = 1}
  }
  return(vlivC(s, rootkC(x, k)[[2]]$value))
}
