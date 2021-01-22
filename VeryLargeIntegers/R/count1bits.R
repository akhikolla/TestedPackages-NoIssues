

# number or bits with value 1 in the base 2 representation of x

#' @title Counting the Number of 1-Bits in vli Objects
#' @author Javier Leiva Cuadrado
#' @param x object of class vli
#' @return integer
#' @description Counting the number of 1-bits in the base 2 expression of vli (Very Large Integer) objects.
#' @examples x <- as.vli("69158247560284795612")
#' count1bits(x)
#' @name 25. Counting 1 bits
#' @rdname count1bits
#' @export count1bits
#'
count1bits <- function(x) UseMethod("count1bits")

#' @rdname count1bits
#' @method count1bits default
#' @export count1bits default
#'
count1bits.default <- function(x) stop("The object passed as argument is neither a vli object nor a 32 bits integer")

#' @rdname count1bits
#' @method count1bits numeric
#' @export count1bits numeric
#'
count1bits.numeric <- function(x){
  if ( is.numeric(x) & (abs(x) < 2147483648) ){
    x = vliC(toString(x))
  }
  else stop("The object passed as argument is neither a vli object nor a 32 bits integer")
  onebitsC(x)
}

#' @rdname count1bits
#' @method count1bits vli
#' @export count1bits vli
#'
count1bits.vli <- function(x){
  onebitsC(x)
}
