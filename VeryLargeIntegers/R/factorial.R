

# factorial of an integer n

#' @title Factorial of a vli Object
#' @author Javier Leiva Cuadrado
#' @param n object of class vli or 32 bits integer
#' @return object of class vli
#' @description \code{factvli} computes and returns the factorial of a vli (Very Large Integers) object. Given a positive integer \code{n}, the factorial of \code{n}, \code{n!}, is defined as the product of all the positive integers from \code{1} to \code{n}.
#' @examples
#' \dontrun{
#' n <- as.vli("420")
#' factvli(n)
#' }
#' @name 07. Factorial
#' @rdname factorial
#' @export factvli
#'
factvli <- function(n) UseMethod("factvli")

#' @rdname factorial
#' @method factvli default
#' @export factvli default
#'
factvli.default <- function(n) stop("The n object passed as argument is neither a vli object nor a 32 bits integer")

#' @rdname factorial
#' @method factvli numeric
#' @export factvli numeric
#'
factvli.numeric <- function(n){
  if( (abs(n) < 2147483648) ){
    if ( n > 0 ) n = vliC(toString(n))
    else stop("factvli is only defined for positive integer numbers")
  }
  else stop("The n object passed as argument is neither a vli object nor a 32 bits integer")
  factbaseC(n)
}

#' @rdname factorial
#' @method factvli vli
#' @export factvli vli
#'
factvli.vli <- function(n){
  if ( ltC(n, .pkgenv$one) ) stop("factvli is only defined for positive integer values")
  factbaseC(n)
}
