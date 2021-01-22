

# Random vli numbers generators



# by number of digits

#' @title Random Generators of vli Objects
#' @author Javier Leiva Cuadrado
#' @param d number of digits of the vli class object being generated; numeric
#' @return objects of class vli in all cases:
#'
#' \code{rvlidigits(d)} returns a object of class vli belonging to the interval \code{[0, 10^d)}
#' @description Random generators of vli (Very Large Integer) objects following different probability distributions.
#' @details The function \code{rvlidigits(d)} returns a vli object of \code{d} digits randomly generated following the uniform distribution. It is the most efficient way of generating random vli objects.
#' @examples rvlidigits(2000)
#' rvliunif(3425, as.vli("2061341345304562604342"))
#' rvlibin(100, 0.6)
#' rvlinegbin(as.vli("1000000"), 0.5)
#' rvliprime(as.vli("100000"), iter = 10, test = "MR")
#' @name 24. Random generators
#' @rdname random
#' @export rvlidigits
#'
rvlidigits <- function(d){
  if ( !is.numeric(d) ) stop("Argument d has to be an integer number greater than 0")
  if ( d < 1 ) stop("Argument d has to be an integer number greater than 0")
  randomvliC(d)
}


# from uniform distribution

# not for the user
rvliunifbase <- function(x, y){
  sumC(divbaseC( randomvliC(y$length), subC(y,x) )[[2]], x)
}

#' @param x lower bound for the object of class vli being generated; object of class vli or 32 bits integer
#' @param y upper bound for the object of class vli being generated; object of class vli or 32 bits integer
#' @return \code{rvliunif(x, y)} returns a object of class vli belonging to the interval \code{[x, y)}
#' @details The function \code{rvliunif(x, y)} returns a vli object randomly generated following the Uniform distribution with parameters \code{x} and \code{y}.
#' @rdname random
#' @export rvliunif
#'
rvliunif <- function(x, y) UseMethod("rvliunif")

#' @rdname random
#' @method rvliunif default
#' @export rvliunif default
#'
rvliunif.default <- function(x, y) stop("The x object passed as argument is neither a vli object nor a 32 bits integer")

#' @rdname random
#' @method rvliunif numeric
#' @export rvliunif numeric
#'
rvliunif.numeric <- function(x, y){
  if ( abs(x) < 2147483648 ){
    x = vliC(toString(x))
  }
  else stop("The x object passed as argument is neither a vli object nor a 32 bits integer")
  if ( !is.vli(y) ){
    if ( is.numeric(y) & (abs(y) < 2147483648) ){
      y = vliC(toString(y))
    }
    else stop("The y object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if ( geqC(x, y) ) stop("y has to be greater than x")
  return(rvliunifbase(x, y))
}

#' @rdname random
#' @method rvliunif vli
#' @export rvliunif vli
#'
rvliunif.vli <- function(x, y){
  if ( !is.vli(y) ){
    if ( is.numeric(y) & (abs(y) < 2147483648) ){
      y = vliC(toString(y))
    }
    else stop("The y object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if ( geqC(x, y) ) stop("y has to be greater than x")
  return(rvliunifbase(x, y))
}


# from binomial distribution

# not for the user
rvlibinbase <- function(n, p){
  out = .pkgenv$zero
  sum = 0
  lim = vlivC(1, c(1, 0000))
  t = divbaseC(n, lim)
  q = as.integer(collapseC(t[[1]]$value))
  r = as.integer(collapseC(t[[2]]$value))
  rm(t); rm(lim);
  if ( q > 0 ){
    for ( i in 1:q ){
      for ( j in 1:10000 ){
        if ( runif(1) < p ) sum = as.integer(sum + 1)
      }
      if ( sum > 2147470000 ){
        out = sumC(out, vliC(toString(sum)))
        sum = 0
      }
    }
  }
  if ( r > 0 ){
    for ( k in 1:r ){
      if ( runif(1) < p ) sum = sum + 1
    }
  }
  sumC(out, vliC(toString(sum)))
}

#' @param n number of independent Bernoulli trials; object of class vli 32 bits integer
#' @param p probability of success; numeric
#' @return \code{rvlibin(n, p)} returns a object of class vli belonging to the interval \code{[0, n]}
#' @details The function \code{rvlibin(n, p)} returns a vli object randomly generated following the Binomial distribution with parameters \code{n} and \code{p}, where \code{n} is the number of Bernoulli trials and \code{p} the probability of success.
#' @rdname random
#' @export rvlibin
#'
rvlibin <- function(n, p) UseMethod("rvlibin")

#' @rdname random
#' @method rvlibin default
#' @export rvlibin default
#'
rvlibin.default <- function(n, p) stop("The number of independent Bernoulli trials, n, has to be a vli object or a 32 bits integer")

#' @rdname random
#' @method rvlibin numeric
#' @export rvlibin numeric
#'
rvlibin.numeric <- function(n, p){
  if ( (n < 2147483648) & (n > 0) ){
    n = vliC(toString(n))
  }
  else stop("The number of independent Bernoulli trials, n, has to be a vli object or a 32 bits integer with value greater than zero")

  if ( !is.numeric(p) ){
    stop("The probability of success, p, has to be a numeric objetc with value between 0 and 1")
  }
  else if ( (p > 1) | (p < 0) ) stop("The probability of success, p, has to be a numeric object with value between 0 and 1")
  rvlibinbase(n, p)
}

#' @rdname random
#' @method rvlibin vli
#' @export rvlibin vli
#'
rvlibin.vli <- function(n, p){
  if ( !is.numeric(p) ){
    stop("The probability of success, p, has to be a numeric object with value between 0 and 1")
  }
  else if ( (p > 1) | (p < 0) ) stop("The probability of success, p, has to be a numeric object with value between 0 and 1")
  rvlibinbase(n, p)
}


# from negative binomial distribution

# not for the user
rvlinegbinbase <- function(s, p){
  lim = vlivC(1, c(100, 0000))
  d = divbaseC(s, lim)
  parts = as.integer(d[[1]])
  rem = as.integer(d[[2]])
  succ = 0
  if ( parts > 0 ){
    for ( i in 1:parts ){
      succ = succ + negbinC(1000000, p)
    }
  }
  succ = succ + negbinC(rem, p)
  succ
}

#' @param s number of successes; vli class object or 32 bits integer
#' @return \code{rvlinegbin(x, y)} returns a object of class vli belonging to the interval \code{[n, Inf)}
#' @details The function \code{rvlinegbin(x, y)} returns a vli object randomly generated following the Negative Binomial distribution with parameters \code{s} and \code{p}, where \code{s} is the number of successes and \code{p} the probability of success.
#' @rdname random
#' @export rvlinegbin
#'
rvlinegbin <- function(s, p) UseMethod("rvlinegbin")

#' @rdname random
#' @method rvlinegbin default
#' @export rvlinegbin default
#'
rvlinegbin.default <- function(s, p) stop("The number of successes, s, has to be a vli object or a 32 bits integer")

#' @rdname random
#' @method rvlinegbin numeric
#' @export rvlinegbin numeric
#'
rvlinegbin.numeric <- function(s, p){
  if ( (s < 2147483648) & (s > 0) ){
    s = vliC(toString(s))
  }
  else stop("The number of successes, s, has to be a vli object or a 32 bits integer with value greater than zero")

  if ( !is.numeric(p) ){
    stop("The probability of success, p, has to be a numeric objetc with value between 0 and 1")
  }
  else if ( (p > 1) | (p < 0) ) stop("The probability of success, p, has to be a numeric object with value between 0 and 1")
  rvlinegbinbase(s, p)
}

#' @rdname random
#' @method rvlinegbin vli
#' @export rvlinegbin vli
#'
rvlinegbin.vli <- function(s, p){
  if ( !is.numeric(p) ){
    stop("The probability of success, p, has to be a numeric object with value between 0 and 1")
  }
  else if ( (p > 1) | (p < 0) ) stop("The probability of success, p, has to be a numeric object with value between 0 and 1")
  rvlinegbinbase(s, p)
}


# random prime number in [2,y)

# not for the user
rvliprimebase <- function(y, iter, test){
  found = FALSE
  if ( test == 'MR' ){
    while ( !found ){
      t = rvliunifbase(.pkgenv$two, y)
      found = is.primeMRbase(t, iter)
    }
    return(t)
  }
  if ( test == 'F' ){
    while ( !found ){
      t = rvliunifbase(.pkgenv$two, y)
      found = is.primeFbase(t, iter)
    }
    return(t)
  }
  if ( test == 'SS' ){
    while ( !found ){
      t = rvliunifbase(.pkgenv$two, y)
      found = is.primeSSbase(t, iter)
    }
    return(t)
  }
}

#' @param iter number of iterations for each number to be tested; numeric
#' @param test chosen primality test: "F" for the Fermat Test, "SS" for the Solovay-Strassen Test or "MR" (by default) for the Miller-Rabin Test; character
#' @return \code{rvliprime(y, iter, test)} returns a object of class vli with the value of a prime number belonging to the interval \code{[2, y)}
#' @details The function \code{rvliprime(y, iter, test)} returns a vli object randomly chosen from the set of primes up to \code{y}.
#' @rdname random
#' @export rvliprime
#'
rvliprime <- function(y, iter=10, test='MR') UseMethod("rvliprime")

#' @rdname random
#' @method rvliprime default
#' @export rvliprime default
#'
rvliprime.default <- function(y, iter=10, test='MR') stop("The object y passed as argument is neither a vli object nor a 32 bits integer")

#' @rdname random
#' @method rvliprime numeric
#' @export rvliprime numeric
#'
rvliprime.numeric <- function(y, iter=10, test='MR'){
  if ( abs(y) < 2147483648 ){
    if ( y > 1 ){
      y = vliC(toString(y))
    }
    else stop("y has to be greater than one")
  }
  else stop("The y argument is neither a vli object nor a 32 bits integer")
  if ( ltC(y, .pkgenv$two) ) stop("y has to be greater than one")
  if ( test != 'F' & test != 'MR' & test != 'SS' ) stop("The test argument has to be 'F' for Fermat test, 'SS' for Solovay-Strassen test or 'MR' (by default) for Miller-Rabin test")
  if ( is.numeric(iter) ){
    if ( iter < 1 ) stop("The number of iterations, iter, has to be a numeric object greater than zero")
  }
  else stop("The number of iterations, iter, has to be a numeric object greater than zero")
  rvliprimebase(y, iter, test)
}

#' @rdname random
#' @method rvliprime vli
#' @export rvliprime vli
#'
rvliprime.vli <- function(y, iter=10, test='MR'){
  if ( ltC(y, .pkgenv$two) ) stop("y has to be greater than one")
  if ( test != 'F' & test != 'MR' & test != 'SS' ) stop("The test argument has to be 'F' for Fermat test, 'SS' for Solovay-Strassen test or 'MR' (by default) for Miller-Rabin test")
  if ( is.numeric(iter) ){
    if ( iter < 1 ) stop("The number of iterations, iter, has to be a numeric object greater than zero")
  }
  else stop("The number of iterations, iter, has to be a numeric object greater than zero")
  rvliprimebase(y, iter, test)
}




