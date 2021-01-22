
# creating am object of class vli

#' @title Very Large Integers Basics
#' @author Javier Leiva Cuadrado
#' @param n value for the vli object being created; character or numeric
#' @description vli is a S3 class that allows to store and operate with arbitrarily large integers. Each object of class vli has 3 attributes (\code{sign}, \code{length} and \code{value}) that can be accessed as shown in the examples. The (absolute) value of the number is stored in a numeric vector to avoid truncation.
#' @details In \code{as.vli(n)}, if \code{n} is numeric, it must be a 32 bits integer to avoid the loss of precision. The idea is to use numeric objects only for small numbers. In other case, character objects are prefered.
#' The function \code{as.integer(x)}, where \code{x} a vli object, only works when the absolute value of \code{x} is up to 2.147.483.648 (32 bits). In other case it returns an error.
#' The function \code{asnumeric(x)} could cause loss of precision if the value of \code{x} is big.
#' The function \code{vli(m)} initialize a list of \code{m} objects of class vli.
#' Punctuation signs are ignored in the creation of vli objects (see the last example).
#' @examples
#' ## Creating a new vli object
#' x <- as.vli("-89027148538375418689123052")
#'
#' ## Printing a vli object
#' print(x)
#'
#' ## Testing the class
#' is.vli(x)
#'
#' ## Coercing into a character object
#' as.character(x)
#'
#' ## Accessing to the attributes of the vli object
#' x$sign
#' x$value
#' x$length
#'
#' ## Punctuation signs are ignored
#' as.vli("2345.25")
#' @name 01. Basics
#' @rdname vli
#' @export as.vli
#'
as.vli <- function(n) UseMethod("as.vli")

#' @rdname vli
#' @method as.vli default
#' @export as.vli default
#'
as.vli.default <- function(n) stop("n argument must verify to be one of the next: 1) A character object containing only an integer number and (optionally) its sign; 2) A 32 bits integer number")

#' @rdname vli
#' @method as.vli character
#' @export as.vli character
#'
as.vli.character <- function(n){
  if(
    !tryCatch(
      is.numeric(as.numeric(n)),
      warning = function(w){
        FALSE
      },
      error = function(e){
        FALSE
      }
    )
    | (grepl("e", n))
  )
  stop("n argument must be either a character object containing only an integer number and (optionally) its sign or a 32 bits integer number")
  vliC(n)
}

#' @rdname vli
#' @method as.vli numeric
#' @export as.vli numeric
#'
as.vli.numeric <- function(n){
  if ( abs(n) >= 2147483648 )
    stop("n argument must be either a character object containing only an integer number and (optionally) its sign or a 32 bits integer number")
  vliC(toString(n))
}


# vli to numeric

#' @title Very Large Integers
#' @param x object of class vli
#' @rdname vli
#' @export asnumeric
#'
asnumeric <- function(x) UseMethod("asnumeric", x)

#' @rdname vli
#' @method asnumeric default
#' @export asnumeric default
#'
asnumeric.default <- function(x) as.numeric(x)

#' @rdname vli
#' @method asnumeric vli
#' @export asnumeric vli
#'
asnumeric.vli <- function(x) (x$sign * as.numeric(collapseC(x$value)))


# vli to integer

#' @rdname vli
#' @method as.integer vli
#' @export as.integer vli
#'
as.integer.vli <- function(x, ...){
  manint = vlivC(1, c(21, 4748, 3648))
  if ( gtC(abs(x), manint) ) stop("Only vli objets with absolute value up to 2147483648 can be coerced to integers")
  else return( x$sign * as.integer(collapseC(x$value)) )
}


# vli to character

#' @rdname vli
#' @method as.integer vli
#' @export as.character vli
#'
as.character.vli <- function(x, ...){
  out = collapseC(x$value)
  if ( x$sign == -1 ) out = paste0("-", out)
  out
}


# Initializing a list of m objects of class vli

#' @rdname vli
#' @param m number of vli objects being initialized; numeric
#' @export vli
#'
vli <- function(m){
  if ( !is.numeric(m) ) stop("The list size, m, must be a numeric object")
  lapply(character(m), as.vli)
}


# printing

#' @rdname vli
#' @param ... further arguments passed to or from other methods
#' @method print vli
#' @export print vli
#'
print.vli <- function(x, ...){
  cat("Very Large Integer: \n")
  print.noquote(printvliC(x))
}


# class belonging

#' @rdname vli
#' @export is.vli
#'
is.vli <- function(x){
  class(x) == "vli"
}







