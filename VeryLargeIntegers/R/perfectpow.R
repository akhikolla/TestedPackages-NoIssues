

# perfect power test

# not for the user
perfectpowbase <- function(x){
  s = x$sign
  x = abs(x)
  b = .pkgenv$two
  while ( leqC(expC(.pkgenv$two, b), x) ){
    a = .pkgenv$one
    c = x
    while ( geqC( subC(c, a), .pkgenv$two) ){
      m = divbaseC( sumC(a, c), .pkgenv$two)[[1]]
      p = expC(m, b)
      x1 = sumC(x, .pkgenv$one)
      if ( gtC(p, x1) ) p = x1
      if ( eqC(p, x) ){
        if ( s == -1 ){
          if ( !is.even(b) ){m$sign = -1; return(list(m, b))}
          else return(list(.pkgenv$zero, .pkgenv$zero))
        }
        else return(list(m, b))
      }
      if ( ltC(p, x) ) a = m
      else c = m
    }
    b = sumC(b, .pkgenv$one)
  }
  return(list(.pkgenv$zero, .pkgenv$zero))
}

#' @title Perfect Power Tools for vli Objects
#' @author Javier Leiva Cuadrado
#' @param x object of class vli or 32 bits integer
#' @return \code{is.perfectpow(x)} returns a Boolean
#'
#' \code{perfectpow(x)} returns a list of two objects of class vli
#' @description A positive integer is a perfect power if it can be expressed as an integer power of another positive integer. That is, a positive integer \code{x} is a perfect power if there exist two positive integers \code{a} and \code{b} such that \code{x = a^b} (note that \code{a} and \code{b} might not be unique).
#' @details The function \code{is.perfectpow(x)} returns \code{TRUE} if there exist two positive integers \code{a} and \code{b} such that \code{x = a^b}, and returns \code{FALSE} if there not exist.
#'
#' The function \code{perfectpow(x)} returns a list of two vli objects, \code{a} and \code{b}, such that \code{x = a^b}. If there not exist such numbers, the two vli objects will be equal to zero. Although the concept is usually defined only for positive integers, the function has been also programmed to work with negative integers.
#' @examples x <- as.vli("234925792")
#' is.perfectpow(x)
#'
#' x <- as.vli("77808066022325383192121677734375")
#' is.perfectpow(x)
#' res <- perfectpow(x)
#' res
#' res[[1]]^res[[2]]
#' @name 12. Perfect power
#' @rdname perfectpow
#' @export perfectpow
#'
perfectpow <- function(x) UseMethod("perfectpow")

#' @rdname perfectpow
#' @method perfectpow default
#' @export perfectpow default
#'
perfectpow.default <- function(x) stop("x has to be passed as vli object or 32 bits integer")

#' @rdname perfectpow
#' @method perfectpow numeric
#' @export perfectpow numeric
#'
perfectpow.numeric <- function(x){
  if ( abs(x) < 2147483648 ){
    x = vliC(toString(x))
  }
  else stop("The x object passed as argument is neither a vli object nor a 32 bits integer")
  perfectpowbase(x)
}

#' @rdname perfectpow
#' @method perfectpow vli
#' @export perfectpow vli
#'
perfectpow.vli <- function(x) perfectpowbase(x)


# is perfectpower

#' @rdname perfectpow
#' @export is.perfectpow
#'
is.perfectpow <- function(x) UseMethod("is.perfectpow")

#' @rdname perfectpow
#' @method is.perfectpow default
#' @export is.perfectpow default
#'
is.perfectpow.default <- function(x) stop("x has to be passed as vli object or 32 bits integer")

#' @rdname perfectpow
#' @method is.perfectpow numeric
#' @export is.perfectpow numeric
#'
is.perfectpow.numeric <- function(x){
  if ( abs(x) < 2147483648 ){
    x = vliC(toString(x))
  }
  else stop("The x object passed as argument is neither a vli object nor a 32 bits integer")
  neqC( perfectpowbase(x)[[1]], .pkgenv$zero)
}

#' @rdname perfectpow
#' @method is.perfectpow vli
#' @export is.perfectpow vli
#'
is.perfectpow.vli <- function(x) neqC( perfectpowbase(x)[[1]], .pkgenv$zero)

