

# least common multiple

# not for the user
lcmulbase <- function(x, y){
  divbaseC(mulbaseC(x, y), gcdbase(x, y))[[1]]
}

#' @title Least Common Multiple for vli Objects
#' @author Javier Leiva Cuadrado
#' @param x object of class vli or 32 bits integer
#' @param y object of class vli or 32 bits integer
#' @return object of class vli
#' @description Computation of the least common multiple of two vli (Very Large Integers) objects.
#' @examples x <- as.vli("125634750214756")
#' y <- as.vli("761048412524216246")
#' lcmul(x, y)
#' @name 10. Least common multiple
#' @rdname lcmul
#' @export lcmul
#'
lcmul <- function(x, y) UseMethod("lcmul")

#' @rdname lcmul
#' @method lcmul default
#' @export lcmul default
#'
lcmul.default <- function(x, y) stop("x and y have to be passed as vli objects or 32 bits integers")

#' @rdname lcmul
#' @method lcmul numeric
#' @export lcmul numeric
#'
lcmul.numeric <- function(x, y){
  if ( abs(x) < 2147483648 ){
    if ( x >= 0 ){
      x = vliC(toString(x))
    }
    else stop("lcmul is only defined for positive integer numbers")
  }
  else stop("The first object passed as argument is neither a vli object nor a 32 bits integer")
  if ( !is.vli(y) ){
    if ( is.numeric(y) & (abs(y) < 2147483648) ){
      if ( y >= 0 ){
        y = vliC(toString(y))
      }
      else stop("lcmul is only defined for positive integer numbers")
    }
    else stop("The second object passed as argument is neither a vli object nor a 32 bits integer")
  }
  else if ( y$sign == -1 ) stop("lcmul is only defined for positive integer numbers")
  return(lcmulbase(x,y))
}

#' @rdname lcmul
#' @method lcmul vli
#' @export lcmul vli
#'
lcmul.vli <- function(x, y){
  if ( x$sign == -1 ) stop("lcmul is only defined for positive integer numbers")
  if ( !is.vli(y) ){
    if ( is.numeric(y) & (abs(y) < 2147483648) ){
      if ( y >= 0 ){
        y = vliC(toString(y))
      }
      else stop("lcmul is only defined for positive integer numbers")
    }
    else stop("The second object passed as argument is neither a vli object nor a 32 bits integer")
  }
  else if ( y$sign == -1 ) stop("lcmul is only defined for positive integer numbers")
  return(lcmulbase(x,y))
}
