
# vli Logarithms



# base10-logarithm

#' @title Integer Logarithms for vli Objects
#' @author Javier Leiva Cuadrado
#' @param x object of class vli or 32 bits integer
#' @return object of class vli
#' @description Computation of integer logarithms and their remainders for objects of class vli.
#'
#' Functions \code{log}, \code{log10} and \code{loge} return respectively the integer generalized logarithm, the integer base-10 logarithm and the integer natural logarithm of the given values. Functions \code{logrem} and \code{log10rem} returns the corresponding remainder.
#' @examples x <- as.vli("3873899469432")
#' log(x, base = 5)
#' logrem(x, base = 5)
#' ( 5^log(x, base = 5) ) + logrem(x, base = 5) == x
#' x <- as.vli("149234629386497858748773210293261249785")
#' log10(x)
#' @name 04. Logarithms
#' @rdname logarithms
#' @method log10 vli
#' @export log10 vli
#'
log10.vli = function(x){
  if ( ltC(x, .pkgenv$one) ) stop("log10 is only defined for positive integer values")
  if ( eqC(x, .pkgenv$one) ) return(.pkgenv$zero)
  vliC(toString(x$length - 1))
}


# base10-logarithm remainder

#' @rdname logarithms
#' @export log10rem
#'
log10rem <- function(x) UseMethod("log10rem")

#' @rdname logarithms
#' @method log10rem default
#' @export log10rem default
#'
log10rem.default <- function(x) stop("The x object passed as argument is neither a vli object nor a 32 bits integer")

#' @rdname logarithms
#' @method log10rem numeric
#' @export log10rem numeric
#'
log10rem.numeric <- function(x){
  if( (abs(x) < 2147483648) & (x > 0) ){
    x = vliC(toString(x))
    if ( ltC(x, .pkgenv$one) ) stop("log10rem is only defined for positive integer numbers")
    subC(x, expC((vlivC(1, c(10))), log10(x)) )
  }
  else stop("log10rem is only defined for positive integer numbers")
}

#' @rdname logarithms
#' @method log10rem vli
#' @export log10rem vli
#'
log10rem.vli <- function(x){
  if (ltC(x, .pkgenv$one)) stop("log10rem is only defined for positive integer numbers")
  subC(x, expC((vlivC(1, c(10))), log10(x)) )
}


# generalized logarithm

#' @rdname logarithms
#' @param base base of the logarithm; object of class vli or 32 bits integer
#' @method log vli
#' @export log vli
#'
log.vli = function(x, base){
  if ( !is.vli(base) ){
    if ( is.numeric(base) & (abs(base) < 2147483648) ){
      if ( base > 0 ){
        base = vliC(toString(base))
      }
      else stop("log.vli is only defined for positive integer base")
    }
    else stop("The base object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if ( ltC(x, .pkgenv$one) ) stop("log.vli is only defined for positive integer numbers")
  if ( eqC(x, .pkgenv$one) ) return(.pkgenv$zero)
  if ( eqC(base, vlivC(1,c(10))) ) return(vliC(toString(x$length - 1)))
  n = .pkgenv$zero
  t = .pkgenv$one
  while ( gtC(x, t) ){
    t = mulbaseC(t, base)
    n = sumC(n, .pkgenv$one)
  }
  subC(n, .pkgenv$one)
}


# generalized logarithm remainder

#' @rdname logarithms
#' @export logrem
#'
logrem <- function(x, base) UseMethod("logrem")

#' @rdname logarithms
#' @method logrem default
#' @export logrem default
#'
logrem.default <- function(x, base) stop("The x object passed as argument is neither a vli object nor a 32 bits integer")

#' @rdname logarithms
#' @method logrem numeric
#' @export logrem numeric
#'
logrem.numeric <- function(x, base){
  if ( !is.vli(base) ){
    if ( is.numeric(base) & (abs(base) < 2147483648) ){
      if ( base > 0 ){
        base = vliC( toString(base) )
      }
      else stop("logrem is only defined for positive integer base")
    }
    else stop("The base object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if( (abs(x) < 2147483648) & (x > 0) ){
    x = vliC(toString(x))
    if ( ltC(x, .pkgenv$one) ) stop("logrem is only defined for positive integer numbers")
    subC(x, expC((base), log(x, base)) )
  }
  else stop("logrem is only defined for positive integer numbers")
}

#' @rdname logarithms
#' @method logrem vli
#' @export logrem vli
#'
logrem.vli <- function(x, base){
  if ( !is.vli(base) ){
    if ( is.numeric(base) & (abs(base) < 2147483648) ){
      if ( base > 0 ){
        base = vliC( toString(base) )
      }
      else stop("logrem is only defined for positive integer base")
    }
    else stop("The base object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if ( ltC(x, .pkgenv$one) ) stop("logrem is only defined for positive integer numbers")
  subC(x, expC((base), log(x, base)) )
}


# base-e logarithm aproximation

#' @rdname logarithms
#' @export loge
#'
loge <- function(x) UseMethod("loge")

#' @rdname logarithms
#' @method loge default
#' @export loge default
#'
loge.default <- function(x) stop("The x object passed as argument is neither a vli object nor a 32 bits integer")

#' @rdname logarithms
#' @method loge numeric
#' @export loge numeric
#'
loge.numeric <- function(x){
  if( (abs(x) < 2147483648) & (x > 0) ){
    x = vliC(toString(trunc(x)))
    if ( ltC(x, .pkgenv$one) ) stop("loge is only defined for positive integer values")
    if ( eqC(x, .pkgenv$one) ) return(.pkgenv$zero)
    vliC(toString(trunc(log(as.numeric(collapseC(x$value)),base=exp(1)))))
  }
  else stop("loge is only defined for positive integer numbers")
}

#' @rdname logarithms
#' @method loge vli
#' @export loge vli
#'
loge.vli <- function(x){
  if ( ltC(x, .pkgenv$one) ) stop("loge is only defined for positive integer values")
  if ( eqC(x, .pkgenv$one) ) return(.pkgenv$zero)
  vliC(toString(trunc(log(as.numeric(collapseC(x$value)),base=exp(1)))))
}
