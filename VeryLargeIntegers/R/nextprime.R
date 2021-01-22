

# Next prime

# not for the user
nextprimebase <- function(n, iter, test){
  t = n
  found = FALSE
  if ( test == 'MR' ){
    while ( !found ){
      t = sumC(t, .pkgenv$one)
      found = is.primeMRbase(t, iter)
    }
    return(t)
  }
  if ( test == 'F' ){
    while ( !found ){
      t = sumC(t, .pkgenv$one)
      found = is.primeFbase(t, iter)
    }
    return(t)
  }
  if ( test == 'SS' ){
    while ( !found ){
      t = sumC(t, .pkgenv$one)
      found = is.primeSSbase(t, iter)
    }
    return(t)
  }
}

#' @title Next Prime Number
#' @author Javier Leiva Cuadrado
#' @param n object of class vli or 32 bits integer
#' @param iter number of iterations for testing whether or not each number is prime; numeric
#' @param test chosen test: "F" for the Fermat Test, "SS" for the Solovay-Strassen Test or "MR" (by default) for the Miller-Rabin Test; character
#' @return object of class vli
#' @description The function \code{nextprime} computes and returns the smallest prime number that is greater than the given number.
#' @details The number of iterations is configurable to set the desired accuracy. A small number of iterations might cause not finding a prime number.
#' @examples n <- as.vli("982234568923564")
#' x <- nextprime(n)
#' x
#' is.prime(x)
#' @name 20. Next prime number
#' @rdname nextprime
#' @export nextprime
#'
nextprime <- function(n, iter=10, test='MR') UseMethod("nextprime")


#' @rdname nextprime
#' @method nextprime default
#' @export nextprime default
#'
nextprime.default <- function(n, iter=10, test='MR') stop("The object n passed as argument is neither a vli object nor a 32 bits integer")


#' @rdname nextprime
#' @method nextprime numeric
#' @export nextprime numeric
#'
nextprime.numeric <- function(n, iter=10, test='MR'){
  if ( abs(n) < 2147483648 ){
    if ( n > 1 ){
      n = vliC(toString(n))
    }
    else stop("n has to be greater than one")
  }
  else stop("The n argument is neither a vli object nor a 32 bits integer")
  if ( ltC(n, .pkgenv$two) ) stop("n has to be greater than one")
  if ( test != 'F' & test != 'MR' & test != 'SS' ) stop("The test argument has to be 'F' for Fermat test, 'SS' for Solovay-Strassen test or 'MR' (by default) for Miller-Rabin test")
  if ( is.numeric(iter) ){
    if ( iter < 1 ) stop("The number of iterations, iter, has to be a numeric object greater than zero")
  }
  else stop("The number of iterations, iter, has to be a numeric object greater than zero")
  nextprimebase(n, iter, test)
}


#' @rdname nextprime
#' @method nextprime vli
#' @export nextprime vli
#'
nextprime.vli <- function(n, iter=10, test='MR'){
  if ( ltC(n, .pkgenv$two) ) stop("n has to be greater than one" )
  if ( test != 'F' & test != 'MR' & test != 'SS' ) stop("The test argument has to be 'F' for Fermat test, 'SS' for Solovay-Strassen test or 'MR' (by default) for Miller-Rabin test")
  if ( is.numeric(iter) ){
    if ( iter < 1 ) stop("The number of iterations, iter, has to be a numeric object greater than zero")
  }
  else stop("The number of iterations, iter, has to be a numeric object greater than zero")
  nextprimebase(n, iter, test)
}
