

# primes up to n

# not for the user
primesbase <- function(n, test, iter, bar){
  out = c(noquote("2"))
  i = .pkgenv$three
  j = 2
  prog = 3
  if ( bar ){
    numn = as.numeric(collapseC(n$value))
    pb = txtProgressBar(min = 0, max = numn, style = 3)
  }
  if ( test == 'MR' ){
    while ( ltC(i, n) ){
      if ( is.primeMRbase(i, iter) ){
        out[j] = noquote(collapseC(i$value))
        j = j+1
      }
      i = sumC(i, .pkgenv$one)
      if ( bar ){
        prog = prog + 1
        setTxtProgressBar(pb, prog)
      }
    }
  }
  else if ( test == 'F' ){
    while ( ltC(i, n) ){
      if ( is.primeFbase(i, iter) ){
        out[j] = noquote(collapseC(i$value))
        j = j+1
      }
      i = sumC(i, .pkgenv$one)
      if ( bar ){
        prog = prog + 1
        setTxtProgressBar(pb, prog)
      }
    }
  }
  if ( test == 'SS' ){
    while ( ltC(i, n) ){
      if ( is.primeSSbase(i, iter) ){
        out[j] = noquote(collapseC(i$value))
        j = j+1
      }
      i = sumC(i, .pkgenv$one)
      if ( bar ){
        prog = prog + 1
        setTxtProgressBar(pb, prog)
      }
    }
  }
  if ( bar ) cat("\n")
  out
}

#' @title Finding All Primes Up to a Given Bound
#' @author Javier Leiva Cuadrado
#' @param n upper bound of the interval in which look for primes; object of class vli or 32 bits
#' @param test chosen test for each number: "F" for the Fermat Test, "SS" for the Solovay-Strassen Test or "MR" (by default) for the Miller-Rabin Test; character
#' @param iter number of iterations for each number being tested; numeric
#' @param bar to choose if display or not a progress bar; boolean
#' @return vector of objects of class "noquote"
#' @description The function \code{primes} displays a vector with all prime numbers up to a given bound. Computation can be made by using different probabilistic primality tests at the user's choice (Fermat Test, Miller-Rabin Test or Solovay-Strassen Test). The number of iterations is also configurable, to set the desired accuracy.
#' @examples primes(n = 600, iter = 10, test = "MR", bar = TRUE)
#' @name 19. Finding all primes
#' @rdname primes
#' @export
#'
primes <- function(n, test='MR', iter=10, bar = TRUE) UseMethod("primes")

#' @rdname primes
#' @method primes default
#' @export primes default
#'
primes.default <- function(n, test='MR', iter=10, bar = TRUE) stop("n has to be specified as a vli object or a 32 bits integer")

#' @rdname primes
#' @method primes numeric
#' @export primes numeric
#'
primes.numeric <- function(n, test='MR', iter=10, bar = TRUE){
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
  if ( !is.logical(bar) ) stop("The bar argument has to be an object of type 'logical'")
  primesbase(n, test, iter, bar)
}

#' @rdname primes
#' @method primes vli
#' @export primes vli
#'
primes.vli <- function(n, test='MR', iter=10, bar = TRUE){
  if ( ltC(n, .pkgenv$two) ) stop("n has to be greater than one")
  if ( test != 'F' & test != 'MR' & test != 'SS' ) stop("The test argument has to be 'F' for Fermat test, 'SS' for Solovay-Strassen test or 'MR' (by default) for Miller-Rabin test")
  if ( is.numeric(iter) ){
    if ( iter < 1 ) stop("The number of iterations, iter, has to be a numeric object greater than zero")
  }
  else stop("The number of iterations, iter, has to be a numeric object greater than zero")
  primesbase(n, test, iter, bar)
}
