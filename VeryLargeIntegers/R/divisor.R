

# find a random divisor

# not for the user
divisorbase <- function(n, iter){
  if ( is.primeMRbase(n, iter) ){
    return(n)
  }
  if ( is.even(n) ) return(.pkgenv$two)
  if ( eqC(divbaseC(n, .pkgenv$three)[[2]], .pkgenv$zero) ) return(.pkgenv$three)
  if ( eqC(divbaseC(n, .pkgenv$five)[[2]], .pkgenv$zero) ) return(.pkgenv$five)
  found = FALSE
  while ( !found ){
    x = rvliunifbase(.pkgenv$zero, n)
    m = gcdbase(x, n)
    if ( (gtC(m, .pkgenv$one)) && (neqC(m, n)) ) found = TRUE
  }
  m
}

#' @title Finding a Random Divisor of a vli Object
#' @author Javier Leiva Cuadrado
#' @param n object of class vli or 32 bits integer
#' @param iter number of iterations for testing if the given number is prime; numeric
#' @return object of class vli
#' @description \code{divisor} returns a randomly chosen divisor of a given number.
#' @details The algorithm determines if the given number is prime or composite by usign the Miller-Rabin Probabilistic Primality Test. If it is prime, it returns the number itself. If it is composite, it returns a randomly chosen divisor. The number of iterations is configurable to set the desired accuracy. A too low number of iterations could cause an infinite loop because of being looking for a divisor of a prime number.
#' @examples r <- rvliprime(100)
#' r
#' x <- r * 51
#' x
#' divisor(x, iter = 100)
#' @name 14. Finding a random divisor
#' @rdname divisor
#' @export divisor
#'
divisor <- function(n, iter=100) UseMethod("divisor")

#' @rdname divisor
#' @method divisor default
#' @export divisor default
#'
divisor.default <- function(n, iter=100) stop("n has to be specified as a vli object or a 32 bits integer")

#' @rdname divisor
#' @method divisor numeric
#' @export divisor numeric
#'
divisor.numeric <- function(n, iter=100){
  if ( abs(n) < 2147483648 ){
    n = vliC(toString(n))
  }
  else stop("The n object passed as argument is neither a vli object nor a 32 bits integer")
  if ( is.numeric(iter) ){
    if ( iter < 1 ) stop("The number of iterations, iter, has to be a numeric object greater than zero")
  }
  else stop("The number of iterations, iter, has to be a numeric object greater than zero")
  divisorbase(n, iter)
}

#' @rdname divisor
#' @method divisor vli
#' @export divisor vli
#'
divisor.vli <- function(n, iter=100){
  if ( is.numeric(iter) ){
    if ( iter < 1 ) stop("The number of iterations, iter, has to be a numeric object greater than zero")
  }
  else stop("The number of iterations, iter, has to be a numeric object greater than zero")
  divisorbase(n, iter)
}
