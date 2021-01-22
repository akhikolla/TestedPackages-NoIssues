

# count primes up to n

# not for the user
countbase <- function(n, iter, bar){
  out = .pkgenv$one
  i = .pkgenv$three
  j = 2
  prog = 3
  numn = as.numeric(collapseC(n$value))
  if ( bar ) pb = txtProgressBar(min = 0, max = numn, style = 3)
  while ( ltC(i, n) ){
    if ( is.primeMRbase(i, iter) ){
      out = sumC(out, .pkgenv$one)
      j = j + 1
    }
    i = sumC(i, .pkgenv$one)
    if ( bar ){
      prog = prog + 1
      setTxtProgressBar(pb, prog)
    }
  }
  if ( bar ) cat("\n")
  out
}

#' @title Counting the Number of Primes Up to a Given Bound
#' @author Javier Leiva Cuadrado
#' @param n upper bound of the interval in which we want to count the number of primes; object of class vli or 32 bits integer
#' @param iter number of iterations for each number being tested; numeric
#' @param bar to choose if display or not a progress bar; boolean
#' @description The function \code{primescount} returns the number of primes found up to a given bound. The implemented algorithm uses the Miller-Rabin Primality Test to determine whether a number is prime or not. The number of iterations is configurable, to set the desired accuracy.
#' @examples
#' \dontrun{
#' ## Counting primes up to 200
#' primescount(n = 200, iter = 10, bar = TRUE)
#'
#' ## Computing the approximation of pi(x)
#' pi(200)
#'
#' ## Showing the numbers by using the Solovay-Strassen test
#' primes(n = 200, iter = 10, test = "SS", bar = TRUE)
#' }
#' @name 22. Counting the number of primes
#' @rdname primescount
#' @export primescount
#'
primescount <- function(n, iter=10, bar = TRUE) UseMethod("primescount")

#' @rdname primescount
#' @method primescount default
#' @export primescount default
#'
primescount.default <- function(n, iter=10, bar = TRUE) stop("n has to be specified as a vli object or a 32 bits integer")

#' @rdname primescount
#' @method primescount numeric
#' @export primescount numeric
#'
primescount.numeric <- function(n, iter=10, bar = TRUE){
  if ( abs(n) < 2147483648 ){
    if ( n > 1 ){
      n = vliC(toString(n))
    }
    else stop("n has to be greater than one")
  }
  else stop("The n argument is neither a vli object nor a 32 bits integer")
  if ( ltC(n, .pkgenv$two) ) stop("n has to be greater than one")
  if ( is.numeric(iter) ){
    if ( iter < 1 ) stop("The number of iterations, iter, has to be a numeric object greater than zero")
  }
  else stop("The number of iterations, iter, has to be a numeric object greater than zero")
  if ( !is.logical(bar) ) stop("The bar argument has to be an object of type 'logical'")
  countbase(n, iter, bar)
}

#' @rdname primescount
#' @method primescount vli
#' @export primescount vli
#'
primescount.vli <- function(n, iter=10, bar = TRUE){
  if ( ltC(n, .pkgenv$two) ) stop("n has to be greater than one")
  if ( is.numeric(iter) ){
    if ( iter < 1 ) stop("The number of iterations, iter, has to be a numeric object greater than zero")
  }
  else stop("The number of iterations, iter, has to be a numeric object greater than zero")
  countbase(n, iter, bar)
}
