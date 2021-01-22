

## Basic arithmetic an logical operators




# Arithmetic operators


# managing binary inputs for basic arithmetic operators
# not for the user
bininput = function(x, y){
  if ( is.vli(x) ){
    if( missing(y) ){1}
    else if ( is.vli(y) ){2}
    else if ( is.numeric(y) ){
      if ( abs(y) < 2147483648 ){3}
      else {4}
    }
    else {5}
  }
  else if ( is.numeric(x) ){
    if ( abs(x) < 2147483648 ){6}
    else {7}
  }
  else {8}
}

#' @title Basic Arithmetic and Logical Operators for vli Objects
#' @author Javier Leiva Cuadrado
#' @param x object of class vli or 32 bits integer
#' @param y object of class vli or 32 bits integer
#' @return objects of class vli with the arithmetic operators; booleans with the logical operators
#' @description Basic arithmetic and logical operators for vli (Very Large Integers) objects.
#' @details As in the creation of vli objects (through the function \code{as.vli}), punctuation signs will be ignored (see the last example).
#'
#' The algorithm implemented for the operator "\code{*}" computes the product with a trivial method when imput numbers have less than 40 digits and with the Karatsuba algorithm for fast multiplications when they are larger.
#' @examples x <- as.vli("712376544526091241")
#' x ^ 61
#' x / as.vli("4225234")
#' x > -x
#' x <= 10000000
#' 13.2415 - as.vli(132415)
#' @name 02. Arithmetic and logic
#' @rdname arithmeticlogic
#' @method + vli
#' @export + vli
#'
'+.vli' = function(x, y) {
  input = bininput(x, y)
  if ( input == 1 ){x}
  else if ( input == 2 ) sumC(x, y)
  else if ( input == 3 ) sumC(x, vliC(toString(as.integer(y))))
  else if ( input == 4 | input == 5 ) stop("The second object passed as argument is neither a vli object nor a 32 bits integer")
  else if ( input == 6) sumC(vliC(toString(x)), y)
  else if ( input == 7 | input == 8 ) stop("The first object passed as argument is neither a vli object nor a 32 bits integer")
}

#' @rdname arithmeticlogic
#' @method - vli
#' @export - vli
#'
'-.vli' = function(x, y){
  input = bininput(x, y)
  if ( input == 1 ){
    if ( x$sign == 1 ) x$sign = -1
    else x$sign = 1
    x
  }
  else if ( input == 2) subC(x, y)
  else if ( input == 3) subC(x, vliC(toString(y)))
  else if ( input == 4 | input == 5 ) stop("The second object passed as argument is neither a vli object nor a 32 bits integer")
  else if ( input == 6 ) subC(vliC(toString(x)), y)
  else if ( input == 7 | input == 8 ) stop("The first object passed as argument is neither a vli object nor a 32 bits integer")
}

#' @rdname arithmeticlogic
#' @method * vli
#' @export * vli
#'
'*.vli' = function(x, y){
  input = bininput(x, y)
  if ( input == 2 ){
    return(mulbaseC(x, y))
  }
  else if ( input == 3) return(mulbaseC(x, vliC(toString(y))) )
  else if ( input == 4 | input == 5 ) stop("The second object passed as argument is neither a vli object nor a 32 bits integer")
  else if ( input == 6) return(mulbaseC(vliC(toString(x)), y))
  else if ( input == 7 | input == 8) stop("The first object passed as argument is neither a vli object nor a 32 bits integer")
}

#' @rdname arithmeticlogic
#' @method / vli
#' @export / vli
#'
'/.vli' = function(x, y){
  input = bininput(x, y)
  if ( input == 2 ){
    if ( y$length == 1 ){
      if ( y$value == c(0) ) stop("division-by-zero error")
      else return(divbaseC(x, y)[[1]])
    }
    else return(divbaseC(x, y)[[1]])
  }
  else if ( input == 3 ){
    if ( y == 0 ) stop("division-by-zero error")
    else return(divbaseC(x, vliC(toString(y)))[[1]])
  }
  else if ( input == 4 | input == 5 ) stop("The second object passed as argument is neither a vli object nor a 32 bits integer")
  else if ( input == 6 ){
    if ( y$value == c(0) ) stop("division-by-zero error")
    else return(divbaseC(vliC(toString(x)), y)[[1]])
  }
  else if ( input == 7 | input == 8 ) stop("The first object passed as argument is neither a vli object nor a 32 bits integer")
}

#' @rdname arithmeticlogic
#' @method %% vli
#' @export %% vli
#'
'%%.vli' = function(x, y){
  input = bininput(x, y)
  if (input == 2) return(divbaseC(x, y)[[2]])
  else if ( input == 3 ) return(divbaseC(x,vliC(toString(y)))[[2]])
  else if ( input == 4 | input == 5 ) stop("The second object passed as argument is neither a vli object nor a 32 bits integer")
  else if ( input == 6 ) return(divbaseC(vliC(toString(x)), y)[[2]])
  else if ( input == 7 | input == 8 ) stop("The first object passed as argument is neither a vli object nor a 32 bits integer")
}

#' @rdname arithmeticlogic
#' @method abs vli
#' @export abs vli
#'
'abs.vli' = function(x){
  x$sign = 1
  x
}

#' @rdname arithmeticlogic
#' @method ^ vli
#' @export ^ vli
#'
'^.vli' = function(x, y){
  input = bininput(x, y)
  if ( input == 2 ) expC(x,y)
  else if ( input == 3 ) expC(x, vliC(toString(y)))
  else if ( input == 4 | input == 5 ) stop("The second object passed as argument is neither a vli object nor a 32 bits integer number")
  else if ( input == 6 ) expC(vliC(toString(x)), y)
  else if ( input == 7 | input == 8 ) stop("The first object passed as argument is neither a vli object nor a 32 bits integer number")
}



# logical operators

#' @rdname arithmeticlogic
#' @method > vli
#' @export > vli
#'
'>.vli' <- function(x, y){
  input = bininput(x, y)
  if ( input == 2 ) gtC(x, y)
  else if ( input == 3 ) gtC(x, vliC(toString(y)))
  else if ( input == 4 | input == 5 ) stop("The second object passed as argument is neither a vli object nor a 32 bits integer")
  else if ( input == 6 ) gtC(vliC(toString(x)), y)
  else if ( input == 7 | input == 8 ) stop("The first object passed as argument is neither a vli object nor a 32 bits integer")
}

#' @rdname arithmeticlogic
#' @method < vli
#' @export < vli
#'
'<.vli' <- function(x, y){
  input = bininput(x, y)
  if ( input == 2) ltC(x, y)
  else if ( input == 3 ) ltC(x, vliC(toString(y)))
  else if ( input == 4 | input == 5 ) stop("The second object passed as argument is neither a vli object nor a 32 bits integer")
  else if ( input == 6 ) ltC(vliC(toString(x)), y)
  else if ( input == 7 | input == 8 ) stop("The first object passed as argument is neither a vli object nor a 32 bits integer")
}

#' @rdname arithmeticlogic
#' @method >= vli
#' @export >= vli
#'
'>=.vli' <- function(x, y){
  input = bininput(x, y)
  if ( input == 2) geqC(x,y)
  else if ( input == 3 ) geqC(x, vliC(toString(y)))
  else if ( input == 4 | input == 5 ) stop("The second object passed as argument is neither a vli object nor a 32 bits integer")
  else if ( input == 6 ) geqC(vliC(toString(x)), y)
  else if ( input == 7 | input == 8 ) stop("The first object passed as argument is neither a vli object nor a 32 bits integer")
}

#' @rdname arithmeticlogic
#' @method <= vli
#' @export <= vli
#'
'<=.vli' <- function(x, y){
  input = bininput(x, y)
  if ( input == 2) leqC(x, y)
  else if ( input == 3) leqC(x, vliC(toString(y)))
  else if ( input == 4 | input == 5) stop("The second object passed as argument is neither a vli object nor a 32 bits integer")
  else if ( input == 6) leqC(vliC(toString(x)), y)
  else if ( input == 7 | input == 8) stop("The first object passed as argument is neither a vli object nor a 32 bits integer")
}

#' @rdname arithmeticlogic
#' @method == vli
#' @export == vli
#'
'==.vli' <- function(x, y){
  input = bininput(x,y)
  if ( input == 2 ) eqC(x, y)
  else if ( input == 3 ) eqC(x, vliC(toString(y)))
  else if ( input == 4 | input == 5 ) stop("The second object passed as argument is neither a vli object nor a 32 bits integer")
  else if ( input == 6 ) eqC(vliC(toString(x)), y)
  else if ( input == 7 | input == 8 ) stop("The first object passed as argument is neither a vli object nor a 32 bits integer")
}

#' @rdname arithmeticlogic
#' @method != vli
#' @export != vli
#'
'!=.vli' <- function(x, y){
  input = bininput(x, y)
  if ( input == 2 ) neqC(x, y)
  else if ( input == 3 ) neqC(x, vliC(toString(y)))
  else if ( input == 4 | input == 5 ) stop("The second object passed as argument is neither a vli object nor a 32 bits integer")
  else if ( input == 6 ) neqC(vliC(toString(x)), y)
  else if ( input == 7 | input == 8 ) stop("The first object passed as argument is neither a vli object nor a 32 bits integer")
}
