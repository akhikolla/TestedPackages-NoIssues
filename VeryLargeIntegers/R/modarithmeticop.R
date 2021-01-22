

# modular arithmetic operators


# modular sum

#' @title Basic Modular-Arithmetic Operators for vli Objects
#' @author Javier Leiva Cuadrado
#' @param x vli class object or 32 bits integer
#' @param y vli class object or 32 bits integer
#' @param mod vli class object or 32 bits integer
#' @param n vli class object or 32 bits integer
#' @return object of class vli
#' @description Basic modular-arithmetic operators for vli (Very Large Integers) objects.
#' @details The functions \code{summod}, \code{submod} and \code{mulmod} compute respectively the sum, the substraction and the multiplication of \code{x} and \code{y} under modulo \code{mod}.
#' @examples x <- as.vli("8925378246957826904701")
#' y <- as.vli("347892325634785693")
#' mod <- as.vli(21341)
#'
#' summod(x, y, mod)
#'
#' mulmod(x, invmod(x, n = 123), mod = 123) == 1
#'
#' z <- divmod(x, y, mod)
#' mulmod(z, y, mod) == x %% mod
#' @name 08. Modular-arithmetic
#' @rdname modarithmeticop
#' @export summod
#'
summod <- function(x, y, mod) UseMethod("summod")

#' @rdname modarithmeticop
#' @method summod default
#' @export summod default
#'
summod.default <- function(x, y, mod) stop ("x, y and mod have to be specified as vli objects or 32 bits integers")

#' @rdname modarithmeticop
#' @method summod numeric
#' @export summod numeric
#'
summod.numeric <- function(x, y, mod){
  if ( abs(x) < 2147483648 ) x = vliC(toString(x))
  else stop("The x object passed as argument is neither a vli object nor a 32 bits integer")
  if ( !is.vli(y) ){
    if ( is.numeric(y) & (abs(y) < 2147483648) ){
      y = vliC(toString(y))
    }
    else stop("The y object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if ( !is.vli(mod) ){
    if ( is.numeric(mod) & (abs(mod) < 2147483648) ){
      mod = vliC(toString(mod))
    }
    else stop("The mod object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if ( eqC(mod, .pkgenv$zero) ) stop("mod argument can not be equal to zero")
  return( divbaseC( sumC(x,y), mod)[[2]] )
}

#' @rdname modarithmeticop
#' @method summod vli
#' @export summod vli
#'
summod.vli <- function(x, y, mod){
  if ( !is.vli(y) ){
    if ( is.numeric(y) & (abs(y) < 2147483648) ){
      y = vliC(toString(y))
    }
    else stop("The y object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if ( !is.vli(mod) ){
    if ( is.numeric(mod) & (abs(mod) < 2147483648) ){
      mod = vliC(toString(mod))
    }
    else stop("The mod object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if ( eqC(mod, .pkgenv$zero) ) stop("mod argument can not be equal to zero")
  return( divbaseC( sumC(x,y), mod)[[2]] )
}


# modular substraction

#' @rdname modarithmeticop
#' @export submod
#'
submod <- function(x, y, mod) UseMethod("submod")

#' @rdname modarithmeticop
#' @method submod default
#' @export submod default
#'
submod.default <- function(x, y, mod) stop ("x, y and mod have to be specified as vli objects or 32 bits integers")

#' @rdname modarithmeticop
#' @method submod numeric
#' @export submod numeric
#'
submod.numeric <- function(x, y, mod){
  if ( abs(x) < 2147483648 ) x = vliC(toString(x))
  else stop("The x object passed as argument is neither a vli object nor a 32 bits integer")
  if ( !is.vli(y) ){
    if ( is.numeric(y) & (abs(y) < 2147483648) ){
      y = vliC(toString(y))
    }
    else stop("The y object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if ( !is.vli(mod) ){
    if ( is.numeric(mod) & (abs(mod) < 2147483648) ){
      mod = vliC(toString(mod))
    }
    else stop("The mod object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if ( eqC(mod, .pkgenv$zero) ) stop("mod argument can not be equal to zero")
  return( divbaseC( subC(x,y), mod)[[2]] )
}

#' @rdname modarithmeticop
#' @method submod vli
#' @export submod vli
#'
submod.vli <- function(x, y, mod){
  if ( !is.vli(y) ){
    if ( is.numeric(y) & (abs(y) < 2147483648) ){
      y = vliC(toString(y))
    }
    else stop("The y object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if ( !is.vli(mod) ){
    if ( is.numeric(mod) & (abs(mod) < 2147483648) ){
      mod = vliC(toString(mod))
    }
    else stop("The mod object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if ( eqC(mod, .pkgenv$zero) ) stop("mod argument can not be equal to zero")
  return( divbaseC( subC(x,y), mod)[[2]] )
}


# modular multiplication

# not for the user
mulmodbase <- function(x, y, mod){
  if ( max(c(x$length, y$length) ) < 40 ){
    return(divbaseC( ( mulC( divbaseC(x,mod)[[2]], divbaseC(y,mod)[[2]] ) ),mod)[[2]] )
  }
  else{
    return(divbaseC( ( karC( divbaseC(x,mod)[[2]], divbaseC(y,mod)[[2]] ) ),mod)[[2]] )
  }
}

#' @rdname modarithmeticop
#' @export mulmod
#'
mulmod <- function(x, y, mod) UseMethod("mulmod")

#' @rdname modarithmeticop
#' @method mulmod default
#' @export mulmod default
#'
mulmod.default <- function(x, y, mod) stop("x, y and mod have to be specified as vli objects or 32 bits integers")

#' @rdname modarithmeticop
#' @method mulmod numeric
#' @export mulmod numeric
#'
mulmod.numeric <- function(x, y, mod){
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
  if ( !is.vli(mod) ){
    if ( is.numeric(mod) & (abs(mod) < 2147483648) ){
      mod = vliC(toString(mod))
    }
    else stop("The mod object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if ( eqC(mod, .pkgenv$zero) ) stop("mod argument can not be equal to zero")
  return(mulmodbase(x, y, mod))
}

#' @rdname modarithmeticop
#' @method mulmod vli
#' @export mulmod vli
#'
mulmod.vli <- function(x, y, mod){
  if ( !is.vli(y) ){
    if ( is.numeric(y) & (abs(y) < 2147483648) ){
      y = vliC(toString(y))
    }
    else stop("The y object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if ( !is.vli(mod) ){
    if ( is.numeric(mod) & (abs(mod) < 2147483648) ){
      mod = vliC(toString(mod))
    }
    else stop("The mod object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if ( eqC(mod, .pkgenv$zero) ) stop("mod argument can not be equal to zero")
  return(mulmodbase(x, y, mod))
}


# not for the user
is.even <- function(x){
  ((tail(x$value, 1)) %% 2) == 0
}


# modular n-th power

#not for the user
powmodbase <- function(x, n, mod){
  res = .pkgenv$one
  x = divbaseC(x, mod)[[2]]
  while( gtC(n, .pkgenv$zero) ){
    if ( !is.even(n) ){
      res = mulmodbase(res, x, mod)
    }
    n = divbaseC(n, .pkgenv$two)[[1]]
    x = mulmod(x, x, mod)
  }
  res
}

#' @details The function \code{powmod} computes the \code{n}-th power of \code{x} under modulo \code{mod}.
#' @rdname modarithmeticop
#' @export powmod
#'
powmod <- function(x, n, mod) UseMethod("powmod")

#' @rdname modarithmeticop
#' @method powmod default
#' @export powmod default
#'
powmod.default <- function(x, n, mod) stop("x, y and mod have to be specified as vli objects or 32 bits integers")

#' @rdname modarithmeticop
#' @method powmod numeric
#' @export powmod numeric
#'
powmod.numeric <- function(x, n, mod){
  if ( abs(x) < 2147483648 ){
    x = vliC(toString(x))
  }
  else stop("The x object passed as argument is neither a vli object nor a 32 bits integer")
  if ( !is.vli(n) ){
    if ( is.numeric(n) & (abs(n) < 2147483648) ){
      n = vliC(toString(n))
    }
    else stop("The y object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if ( !is.vli(mod) ){
    if ( is.numeric(mod) & (abs(mod) < 2147483648) ){
      mod = vliC(toString(mod))
    }
    else stop("The mod object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if ( eqC(mod, .pkgenv$zero) ) stop("mod argument can not be equal to zero")
  return(powmodbase(x, n, mod))
}

#' @rdname modarithmeticop
#' @method powmod vli
#' @export powmod vli
#'
powmod.vli <- function(x, n, mod){
  if ( !is.vli(n) ){
    if ( is.numeric(n) & (abs(n) < 2147483648) ){
      n = vliC(toString(n))
    }
    else stop("The y object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if ( !is.vli(mod) ){
    if ( is.numeric(mod) & (abs(mod) < 2147483648) ){
      mod = vliC(toString(mod))
    }
    else stop("The mod object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if ( eqC(mod, .pkgenv$zero) ) stop("mod argument can not be equal to zero")
  return(powmodbase(x, n, mod))
}


# modular multiplicative inverse of x in Zn:  y=x^(-1) such that x*y = 1 (mod n)

# not for the user
invmodbase <- function(x, n){
  if ( neqC( gcdbase(x, n), .pkgenv$one) ){
    stop("x and n are not coprimes so it does not exists a multiplicative inverse of x in the ring of integer modulo n")
  }
  else  return(exteuclidbase(x, n)[[2]])
}

#' @details The function \code{invmod} returns the modular multiplicative inverse of \code{x} in Z\code{n}; that is,  \code{y = x^(-1)} such that \code{ x * y = 1 (}mod\code{ n)}.
#' @rdname modarithmeticop
#' @export invmod
#'
invmod <- function(x, n) UseMethod("invmod")

#' @rdname modarithmeticop
#' @method invmod default
#' @export invmod default
#'
invmod.default <- function(x, n) stop("x and n have to be specified as vli objects or 32 bits integers")

#' @rdname modarithmeticop
#' @method invmod numeric
#' @export invmod numeric
#'
invmod.numeric <- function(x, n){
  if ( abs(x) < 2147483648 ){
    if ( x >= 0 ){
      x = vliC(toString(x))
    }
    else stop("invmod is only defined for positive integer numbers")
  }
  else stop("The first object passed as argument is neither a vli object nor a 32 bits integer")
  if ( !is.vli(n) ){
    if ( is.numeric(n) & (abs(n) < 2147483648) ){
      if ( n >= 0 ){
        n = vliC(toString(n))
      }
      else stop("invmod is only defined for positive integer numbers")
    }
    else stop("The second object passed as argument is neither a vli object nor a 32 bits integer")
  }
  else if ( n$sign == -1 ) stop("invmod is only defined for positive integer numbers")
  if ( eqC(n, .pkgenv$zero) ) stop("n can not be equal to zero")
  return(invmodbase(x, n))
}

#' @rdname modarithmeticop
#' @method invmod vli
#' @export invmod vli
#'
invmod.vli <- function(x, n){
  if ( x$sign == -1 ) stop("invmod is only defined for positive integer numbers")
  if ( !is.vli(n) ){
    if ( is.numeric(n) & (abs(n) < 2147483648) ){
      if ( n >= 0 ){
        n = vliC(toString(n))
      }
      else stop("invmod is only defined for positive integer numbers")
    }
    else stop("The second object passed as argument is neither a vli object nor a 32 bits integer")
  }
  else if ( n$sign == -1 ) stop("invmod is only defined for positive integer numbers")
  if ( eqC(n, .pkgenv$zero) ) stop("n can not be equal to zero")
  return(invmodbase(x, n))
}


# modular division. divmod(x,y,mod) return z such that (y * z) % mod = x % mod

# not for the user
divmodbase <- function(x, y, mod){
  if ( neqC( gcdbase(y, mod), .pkgenv$one) ){
    stop("y and mod are not coprimes so it does not exists a multiplicative inverse of y in the ring of integer modulo mod and modular division is not defined")
  }
  mulbaseC(invmodbase(y, mod), x) %% mod
}

#' @details The function \code{divmod} returns the modular division of \code{x} over \code{y}; that is, \code{z} such that \code{y * z (}mod \code{mod) = x (}mod \code{mod)}.
#' @rdname modarithmeticop
#' @export divmod
#'
divmod <- function(x, y, mod) UseMethod("divmod")

#' @rdname modarithmeticop
#' @method divmod default
#' @export divmod default
#'
divmod.default <- function(x, y, mod) stop("x, y and mod have to be specified as vli objects or 32 bits integers")

#' @rdname modarithmeticop
#' @method divmod numeric
#' @export divmod numeric
#'
divmod.numeric <- function(x, y, mod){
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
  if ( !is.vli(mod) ){
    if ( is.numeric(mod) & (abs(mod) < 2147483648) ){
      mod = vliC(toString(mod))
    }
    else stop("The mod object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if ( eqC(mod, .pkgenv$zero) ) stop("mod argument can not be equal to zero")
  return(divmodbase(x, y, mod))
}

#' @rdname modarithmeticop
#' @method divmod vli
#' @export divmod vli
#'
divmod.vli <- function(x, y, mod){
  if ( !is.vli(y) ){
    if ( is.numeric(y) & (abs(y) < 2147483648) ){
      y = vliC(toString(y))
    }
    else stop("The y object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if ( !is.vli(mod) ){
    if ( is.numeric(mod) & (abs(mod) < 2147483648) ){
      mod = vliC(toString(mod))
    }
    else stop("The mod object passed as argument is neither a vli object nor a 32 bits integer")
  }
  if ( eqC(mod, .pkgenv$zero) ) stop("mod argument can not be equal to zero")
  return(divmodbase(x, y, mod))
}

