
# function to calc Fibonacci numbers
calc_fib = function(n) {
  if (n < 0L) stop("n should be >= 0")
  if (n == 0L) return(0L)
  if (n == 1L || n == 2L) return(1L)
  x = rep(1L, n)
  for (i in 3L:n) x[[i]] = x[[i - 1]] + x[[i - 2]]
  x[[n]]
}

#* Fibonacci numbers calculation
#* @param n Integer number
#* @get /fib
#* @serializer unboxedJSON
function(n) {
  n = as.integer(n)
  if (is.na(n)) {
    stop("\"n\"must be integer number.")
  }
  calc_fib(n)
}
