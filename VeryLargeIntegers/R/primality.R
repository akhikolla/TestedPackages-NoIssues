

# Probabilistic Primality Tests


# Fermat test

# not for the user
is.primeFbase <- function(x, iter){
  x[1] = 1 # abs(x)
  if ( leqC(x, .pkgenv$one) | eqC(x, .pkgenv$four) ) return(FALSE)
  if ( leqC(x, .pkgenv$three ) | eqC(x, .pkgenv$five ) ) return(TRUE)
  if ( is.even(x) ) return(FALSE)
  if ( eqC(divbaseC(x, .pkgenv$three)[[2]], .pkgenv$zero) ) return(FALSE)
  if ( eqC(divbaseC(x, .pkgenv$five)[[2]], .pkgenv$zero) ) return(FALSE)
  for ( i in 1:iter ){
    a = rvliunifbase(.pkgenv$two, subC(x, .pkgenv$one))
    if ( neqC(powmodbase(a, subC(x, .pkgenv$one), x), .pkgenv$one) ) return(FALSE)
  }
  return(TRUE)
}

#' @title Probabilistic Primality Tests for vli Objects
#' @author Javier Leiva Cuadrado
#' @param x number to be tested; object of class vli or 32 bits integer
#' @param iter number of iterations; numeric
#' @return boolean:
#' if the test determines that the given number is prime it returns \code{TRUE}
#' if the test determines that the given number is composite it returns \code{FALSE}
#' @description Functions to compute different probabilistic primality tests for vli (Very Large Integer) objects.
#'
#'The function \code{is.primeF} computes the Fermat Primality Test.
#' @details Probabilistic primality tests are algorithms that determine if an integer is prime or composite. They are not deterministic tests so there is a probability of error (it is never reported a prime number as composite, but it is possible for a composite number to be reported as prime). This probability of error can be calculated and reduced as much as we want by increasing the number of iterations of the test.
#'
#' Each test is different, therefore they have different computational efficiency and one could be better than other for testing some numbers. However, the Miller-Rabin test is the most accurated of all three and, because of that, it is the test chosen by default in every function that needs primality testing in the present package.
#'
#' The Fermat Primality Test detects composite numbers by using the Fermat's Little Theorem, which says that, if \code{p} is prime, for any integer \code{a} satisfaying \code{gcd(a, p) = 1} we have that \code{a^(p-1) = 1 }(\code{mod p}).
#' Each iteration randomly pick an integer \code{a}. The more iterations are computed, the greater probability to find an \code{a} that does not verify such conditions and, therefore, it reveals that \code{p} is composite. However, there are some composite numbers \code{p} that have the property that \code{a^(p-1) = 1 }(\code{mod p}) for every \code{a} coprime to \code{p}. These numbers are called Carmichael numbers or Fermat pseudoprimes, and it is not possible for the Fermat Test to detect that they are composite numbers. But there are only 105212 such numbers up to \code{10^15} (approximately 1 Carmichael number per each 10.000.000.000 integer numbers). The first five are: 561, 1105, 1729, 2465 and 2821.
#'
#' As a conclusion, we can say that if the chosen \code{x} number is prime, the Fermat test returns \code{TRUE}. If it is an odd composite (but not a Carmichael number), it returns \code{FALSE} with probability at least \code{1/2^k}, where \code{k} is the number of computed iterations.
#' @examples
#' \dontrun{
#' ## Testing a 32 bits integer using the Miller-Rabin Test
#' is.primeMR(2845127, iter = 10)
#'
#' ## Testing an object of class vli using the Fermat Test
#' x <- as.vli("2801401243675128975602569907852141")
#' is.primeF(x, iter = 100)
#'
#' ## Testing the same object of class vli using the general
#' ##     is.prime function and the Solovay-Strassen Test
#' is.prime(x, iter = 100, test = "SS")
#' }
#' @name 18. Probabilistic primality tests
#' @rdname primality
#' @export is.primeF
#'
is.primeF <- function(x, iter=10) UseMethod("is.primeF")

#' @rdname primality
#' @method is.primeF default
#' @export is.primeF default
#'
is.primeF.default <- function(x, iter=10) stop("x has to be specified as a vli object or a 32 bits integer")

#' @rdname primality
#' @method is.primeF numeric
#' @export is.primeF numeric
#'
is.primeF.numeric <- function(x, iter=10){
  if ( abs(x) < 2147483648 ){
    x = vliC(toString(x))
  }
  else stop("The x object passed as argument is neither a vli object nor a 32 bits integer")
  if ( is.numeric(iter) ){
    if ( iter < 1 ) stop("The number of iterations, iter, has to be a numeric object greater than zero")
  }
  else stop("The number of iterations, iter, has to be a numeric object greater than zero")
  is.primeFbase(x, iter)
}

#' @rdname primality
#' @method is.primeF vli
#' @export is.primeF vli
#'
is.primeF.vli <- function(x, iter=10){
  if ( is.numeric(iter) ){
    if ( iter < 1 ) stop("The number of iterations, iter, has to be a numeric object greater than zero")
  }
  else stop("The number of iterations, iter, has to be a numeric object greater than zero")
  is.primeFbase(x, iter)
}


# Miller-Rabin test

# not for the user
MRtest <- function(x, d){
  a = rvliunifbase(.pkgenv$two, subC(x, .pkgenv$one))
  n = powmodbase(a, d, x)
  if ( eqC(n, .pkgenv$one) | eqC(n, subC(x, .pkgenv$one)) ) return(TRUE)
  while ( neqC(d, subC(x, .pkgenv$one)) ){
    n = mulmodbase(n, n, x)
    d = d * .pkgenv$two
    if ( eqC(n, .pkgenv$one) ) return(FALSE)
    if ( eqC(n, subC(x, .pkgenv$one)) ) return(TRUE)
  }
  return(FALSE)
}

# not for the user
is.primeMRbase <- function(x, iter){
  x[1] = 1  # abs(x)
  if ( leqC(x, .pkgenv$one) | eqC(x, .pkgenv$four) ) return(FALSE)
  if ( leqC( x, .pkgenv$three) | eqC(x, .pkgenv$five) ) return(TRUE)
  if ( is.even(x) ) return(FALSE)
  if ( eqC(divbaseC(x, .pkgenv$three)[[2]], .pkgenv$zero) ) return(FALSE)
  if ( eqC(divbaseC(x, .pkgenv$five)[[2]], .pkgenv$zero) ) return(FALSE)
  d = subC(x, .pkgenv$one)
  while ( eqC( (divbaseC(d, .pkgenv$two)[[2]]), .pkgenv$zero) ){
    d = divbaseC(d, .pkgenv$two)[[1]]
  }
  for ( i in 1:iter ){
    if ( MRtest(x, d) == FALSE ) return(FALSE)
  }
  return(TRUE)
}


#' @description The function \code{is.primeMR} computes the Miller-Rabin Primality Test.
#' @details The Miller-Rabin Primality Test is a more sophisticated version of the Fermat test. If the chosen \code{x} number is prime, the test returns \code{TRUE}. If \code{x} is an odd composite the algorithm returns \code{TRUE} (that is, it fails) with probability less than \code{1/4^k}, where \code{k} is the number of computed iterations. In cases of very big numbers, the probability is even smaller.
#' @rdname primality
#' @export is.primeMR
#'
is.primeMR <- function(x, iter=10) UseMethod("is.primeMR")

#' @rdname primality
#' @method is.primeMR default
#' @export is.primeMR default
#'
is.primeMR.default <- function(x, iter=10) stop("x has to be specified as a vli object or a 32 bits integer")

#' @rdname primality
#' @method is.primeMR numeric
#' @export is.primeMR numeric
#'
is.primeMR.numeric <- function(x, iter=10){
  if ( abs(x) < 2147483648 ){
    x = vliC(toString(x))
  }
  else stop("The x object passed as argument is neither a vli object nor a 32 bits integer")
  if ( is.numeric(iter) ){
    if ( iter < 1 ) stop("The number of iterations, iter, has to be a numeric object greater than zero")
  }
  else stop("The number of iterations, iter, has to be a numeric object greater than zero")
  is.primeMRbase(x, iter)
}

#' @rdname primality
#' @method is.primeMR vli
#' @export is.primeMR vli
#'
is.primeMR.vli <- function(x, iter=10){
  if ( is.numeric(iter) ){
    if ( iter < 1 ) stop("The number of iterations, iter, has to be a numeric object greater than zero")
  }
  else stop("The number of iterations, iter, has to be a numeric object greater than zero")
  is.primeMRbase(x, iter)
}


# Solovay-Strassen test

# not for the user
is.primeSSbase <- function(x, iter){
  if ( ltC(x, .pkgenv$two) ) return(FALSE)
  if ( neqC(x, .pkgenv$two) & eqC(divbaseC(x, .pkgenv$two)[[2]], .pkgenv$zero) ) return(FALSE)
  for ( i in 1:iter ){
    a = rvliunifbase(.pkgenv$one, x)
    j = divbaseC(sumC(x, jacobbase(a, x)), x)[[2]]
    m = powmod(a, divbaseC(subC(x, .pkgenv$one), .pkgenv$two)[[1]], x)
    if ( eqC(j, .pkgenv$zero) | neqC(m, j) ) return(FALSE)
  }
  return(TRUE)
}

#' @description The function \code{is.primeSS} computes the Solovay-Strassen Primality Test.
#' @details The Solovay-Strassen test is based in a known algebraic property of the Jacobi symbol. The probabily of failure is also less than \code{1/2^k}, where \code{k} is the number of computed iterations. However, unlike it happens with the Fermat test, there are not odd composite numbers that can not be detected with enough iterations of the Solovay-Strassen test.
#' @rdname primality
#' @export is.primeSS
#'
is.primeSS <- function(x, iter=10) UseMethod("is.primeSS")

#' @rdname primality
#' @method is.primeSS default
#' @export is.primeSS default
#'
is.primeSS.default <- function(x, iter=10) stop ("x has to be specified as a vli object or a 32 bits integer")

#' @rdname primality
#' @method is.primeSS numeric
#' @export is.primeSS numeric
#'
is.primeSS.numeric <- function(x, iter=10){
  if ( abs(x) < 2147483648 ){
    x = vliC(toString(x))
  }
  else stop("The x object passed as argument is neither a vli object nor a 32 bits integer")
  if ( is.numeric(iter) ){
    if ( iter < 1 ) stop("The number of iterations, iter, has to be a numeric object greater than zero")
  }
  else stop("The number of iterations, iter, has to be a numeric object greater than zero")
  is.primeSSbase(x, iter)
}

#' @rdname primality
#' @method is.primeSS vli
#' @export is.primeSS vli
#'
is.primeSS.vli <- function(x, iter=10){
  if ( is.numeric(iter) ){
    if ( iter < 1 ) stop("The number of iterations, iter, has to be a numeric object greater than zero")
  }
  else stop("The number of iterations, iter, has to be a numeric object greater than zero")
  is.primeSSbase(x, iter)
}


# is prime, general function

#' @param test chosen test: "F" for the Fermat Test, "SS" for the Solovay-Strassen Test or "MR" (by default) for the Miller-Rabin Test; character
#' @description The function \code{is.prime} is a general function that computes the test specified in the \code{test} argument.
#' @rdname primality
#' @export is.prime
#'
is.prime <- function(x, iter=10, test='MR'){
  if ( test == 'MR' ){is.primeMR(x, iter)}
  else if ( test == 'F' ){is.primeF(x, iter)}
  else if ( test == 'SS' ){is.primeSS(x, iter)}
  else stop("The test argument has to be 'F' for Fermat test, 'SS' for Solovay-Strassen test or 'MR' (by default) for Miller-Rabin test")
}
