

# Fibonacci numbers tools


# First m Fibonacci numbers

# not for the user
Fibobase <- function(m, output){
  out = vli(m)
  out[[1]] = .pkgenv$zero
  out[[2]] = .pkgenv$one
  for ( i in 3:m ){
    out[[i]] = sumC(out[[i-1]], out[[i-2]])
  }
  if ( output == "list" ) out
  else if ( output == "print" ){
    outp = c(noquote("0"))
    for ( j in 1:m ){outp[[j]] = noquote(collapseC(out[[j]]$value))}
    outp
  }
}

#' @title Fibonacci Numbers Tools for vli Objects
#' @author Javier Leiva Cuadrado
#' @param m object of class vli or 32 bits integer
#' @param output chosen way for objects being returned: \code{'list'} to return the result as a list of vli objects or \code{'print'} (by default) to simply display the result on the screen; character
#' @return The function \code{Fibonacci(m, output)} returns a list of objects of class vli or the result displayed on the screen, depending on the \code{output} argument.
#'
#' The function \code{nthFibonacci(n)} returns a object of class vli.
#'
#' The function \code{is.Fibonacci(x)} returns a boolean.
#' @description The Fibonacci Sequence is defined as follows:
#'
#' \code{x[1] = 0},
#'
#' \code{x[2] = 1},
#'
#' \code{...}
#'
#' \code{x[n] = x[n-1] + x[n-2]}.
#'
#' A positive integer is said to be a Fibonacci Number if it is an element of the Fibonacci Sequence.
#'
#' The function \code{Fibonacci(m, output)} computes and displays the first \code{m} elements of the Fibonacci Sequence.
#'
#' The function \code{nthFibonacci(n)} computes and displays the \code{n}-th element of the Fibonacci Sequence.
#'
#' The function \code{is.Fibonacci(x)} says whether or not \code{x} is a Fibonacci Number.
#' @examples Fibonacci(200)
#'
#' n <- as.vli("50000")
#' nthFibonacci(n)
#'
#' x <- as.vli("5358359254990966640871840")
#' is.Fibonacci(x)
#'
#' y <- x + 1
#' is.Fibonacci(y)
#' @name 23. Fibonacci numbers
#' @rdname Fibonacci
#' @export Fibonacci
#'
Fibonacci <- function(m, output="print") UseMethod("Fibonacci")

#' @rdname Fibonacci
#' @method Fibonacci default
#' @export Fibonacci default
#'
Fibonacci.default <- function(m, output="print") stop("m has to be a numeric object")

#' @rdname Fibonacci
#' @method Fibonacci numeric
#' @export Fibonacci numeric
#'
Fibonacci.numeric <- function(m, output="print"){
  if ( m < 1 ) stop("m has to be a numeric object greater than zero")
  if ( output!="print" & output!="list" ) stop("output argument has to be 'list' to return the result as a list of vli objects or 'print' (by default) to simply display the result on the screen")
  Fibobase(m, output)
}


# n-th Fibonacci number

# not for the user
nthFibobase <- function(n){
  a = .pkgenv$zero
  b = .pkgenv$one
  i = .pkgenv$three
  if ( eqC(n, .pkgenv$one) ) return(a)
  if ( eqC(n, .pkgenv$two) ) return(b)
  while ( leqC(i, n) ){
    c = sumC(a, b)
    a = b
    b = c
    i = sumC(i, .pkgenv$one)
  }
  b
}

#' @param n vli class object or 32 bits integer
#' @rdname Fibonacci
#' @export nthFibonacci
#'
nthFibonacci <- function(n) UseMethod("nthFibonacci")

#' @rdname Fibonacci
#' @method nthFibonacci default
#' @export nthFibonacci default
#'
nthFibonacci.default <- function(n) stop("n has to be a vli object or a 32 bits integer")

#' @rdname Fibonacci
#' @method nthFibonacci numeric
#' @export nthFibonacci numeric
#'
nthFibonacci.numeric <- function(n){
  if ( abs(n) < 2147483648){
    if ( n < 1 ) stop("n has to be greater than zero")
    n = vliC(toString(n))
  }
  else stop("The object passed as argument is neither a vli object nor a 32 bits integer")
  nthFibobase(n)
}

#' @rdname Fibonacci
#' @method nthFibonacci vli
#' @export nthFibonacci vli
#'
nthFibonacci.vli <- function(n){
  if ( n < 1 ) stop("n has to be greater than zero")
  nthFibobase(n)
}


# is Fibonacci

# not for the user
is.Fibobase <- function(x){
  if ( eqC(x, .pkgenv$zero) ) return(TRUE)
  p = mulbaseC( mulbaseC(.pkgenv$five, x), x )
  p1 = sumC(p, .pkgenv$four)
  s1 = sqrt(p1)
  p2 = subC(p, .pkgenv$four)
  s2 = sqrt(p2)
  ( eqC((s1*s1), p1)  |  eqC((s2*s2), p2) )
}

#' @param x vli class object or 32 bits integer
#' @rdname Fibonacci
#' @export is.Fibonacci
#'
is.Fibonacci <- function(x) UseMethod("is.Fibonacci")

#' @rdname Fibonacci
#' @method is.Fibonacci default
#' @export is.Fibonacci default
#'
is.Fibonacci.default <- function(x) stop("x has to be a vli object or a 32 bits integer")

#' @rdname Fibonacci
#' @method is.Fibonacci numeric
#' @export is.Fibonacci numeric
#'
is.Fibonacci.numeric <- function(x){
  if ( abs(x) < 2147483648 ){
    if ( x < 0 ) stop("x has to be a positive integer")
    x = vliC(toString(x))
  }
  else stop("The object passed as argument is neither a vli object nor a 32 bits integer")
  is.Fibobase(x)
}

#' @rdname Fibonacci
#' @method is.Fibonacci vli
#' @export is.Fibonacci vli
#'
is.Fibonacci.vli <- function(x){
  if ( ltC(x, .pkgenv$zero) ) stop("x has to be a positive integer")
  is.Fibobase(x)
}
