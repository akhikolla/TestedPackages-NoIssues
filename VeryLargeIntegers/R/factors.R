

# find all prime factors

# not for the user
factorsbase <- function(n, iter, output){
  if ( is.primeMRbase(n, iter) ) return(n)
  out = c(noquote("0"))
  tot = .pkgenv$one
  x = n
  i=1
  while ( neqC(tot, n) ){
    m = divisor(x)
    tot = mulbaseC(tot, m)
    x = divbaseC(x, m)[[1]]
    out[i] = noquote(collapseC(m$value))
    i = i+1
  }
  if ( output == "print" ) return(out)
  else if ( output == "list" ){
    str = strsplit(as.character(out), split=" ")
    l = length(str)
    outl = vli(l)
    for ( j in 1:l ){outl[[j]] = vliC(str[[j]])}
    return(outl)
  }
}

#' @title Factorization of vli Objects
#' @author Javier Leiva Cuadrado
#' @param n integer to be factorized; vli class object or 32 bits integer
#' @param iter number of iterations for testing if the given number is prime; numeric
#' @param output chosen way for objects being returned: \code{'list'} to return the result as a list of vli objects or \code{'print'} (by default) to simply display the result on the screen; character
#' @return list of objects of class vli or the result displayed on the screen, depending on the \code{output} argument
#' @description \code{factors} returns all the prime factors of a given number.
#' @details The implemented algorithm is based in a Monte Carlo method for integer factorization called Pollard's Rho Algorithm.
#'
#' It determines if the given number is prime or composite by usign the Miller-Rabin Probabilistic Primality Test. If it is prime, it returns the number itself. If it is composite, it calls iteratively the \code{divisor} function until all the prime factors of the given number are found.
#'
#' It is a Monte Carlo method, therefore it is not deterministic. The number of iterations is configurable, to set the desired accuracy. A too low number of iterations could cause an infinite loop because of being looking for a divisor of a prime number.
#' @examples
#' x <- as.vli("584843")
#' factors(x, iter = 100)
#' @name 15. Factorization
#' @rdname factors
#' @export factors
#'
factors <- function(n, iter=10, output="print") UseMethod("factors")

#' @rdname factors
#' @method factors default
#' @export factors default
#'
factors.default <- function(n, iter=10, output="print") stop("n has to be specified as a vli object or a 32 bits integer")

#' @rdname factors
#' @method factors numeric
#' @export factors numeric
#'
factors.numeric <- function(n, iter=10, output="print"){
  if ( abs(n) < 2147483648 ){
    n = vliC(toString(n))
  }
  else stop("The n object passed as argument is neither a vli object nor a 32 bits integer")
  if ( is.numeric(iter) ){
    if ( iter < 1 ) stop("The number of iterations, iter, has to be a numeric object greater than zero")
  }
  else stop("The number of iterations, iter, has to be a numeric object greater than zero")
  factorsbase(n, iter, output)
}

#' @rdname factors
#' @method factors vli
#' @export factors vli
#'
factors.vli <- function(n, iter=10, output="print"){
  if ( is.numeric(iter) ){
    if ( iter < 1 ) stop("The number of iterations, iter, has to be a numeric object greater than zero")
  }
  else stop("The number of iterations, iter, has to be a numeric object greater than zero")
  if ( output!="print" & output!="list" ) stop("output argument has to be 'list' to return the result as a list of vli objects or 'print' (by default) to simply display the result on the screen")
  factorsbase(n, iter, output)
}
