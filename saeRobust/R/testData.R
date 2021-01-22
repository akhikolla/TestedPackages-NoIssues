#' Construction of test data
#'
#' @param ... matrices
#'
#' @export
#'
#' @references Weihs / Mersmann / Ligges (2014): Foundations of Statistical
#'   Algorithms: With References to R Packages
#'
#' @rdname constTestData
#' @examples
#' ## Examples from Weihs et. al. (2014) p. 108
#' library("Matrix")
#' testMatX(Matrix(998), Matrix(998))
#' Z <- Matrix(c(998, 0, 0, 0), 2, 2)
#' testMatX(Z, Z)
testMatX <- function(...) {

    Z <- list(...)
    p <- nrow(Z[[1]])
    m <- length(Z) + 1
    n <- m * p

    I <- do.call(rbind, replicate(m, Diagonal(p), simplify = FALSE))

    listOfCols <- lapply(1:m, function(i) {
        if (i == m) {
            ind <- c(rep(0, length.out = n - p), rep(1, length.out = p))
            Reduce(rbind, replicate(m, Z[[1]])) - I * ind
        } else {
            ind <- c(rep(1, length.out = n - i * p), rep(0, length.out = i * p))
            Reduce(rbind, replicate(m, Z[[i]])) + I * ind
        }
    })

    Matrix(do.call(cbind, listOfCols), forceCheck = TRUE)

}

#' @param x a matrix
#' @param beta a vector with parameters
#' @rdname constTestData
#' @export
#' @examples
#' testResponse0(testMatX(Matrix(1)))
testResponse0 <- function(x, beta = rep(1, ncol(x))) {
    as.numeric(x %*% beta)
}

#' @param y0 a response vector (numeric)
#' @param k values in 1 to 4 (integer)
#' @param .sd the standard deviation used for random numbers
#' @export
#' @rdname constTestData
#' @examples
#' library("magrittr")
#' Matrix(1) %>% testMatX %>% testResponse0 %>% testResponse
testResponse <- function(y0, k = 1:4, .sd = sd(y0)) {

    l2norm <- . %>% .^2 %>% sum %>% sqrt

    i <- -21
    r0 <- rnorm(length(y0), sd = .sd)

    for (j in seq(19, -19, -2)) {
        if (l2norm(r0) >= 2^j * l2norm(y0)) { i <- j; break}
    }

    ck <- 2^c(-21 - i, -3 - i, 1 - i, 19 - i)

    y0s <- do.call(cbind, replicate(length(k), Matrix(y0, ncol = 1)))
    r0s <- do.call(cbind, lapply(ck[k], . %>% { . * r0 } %>% Matrix(ncol = 1)))
    y0s + r0s

}

#' @param n dimension
#'
#' @export
#' @rdname constTestData
testRook <- function(n) {
  spdep::nb2mat(spdep::cell2nb(n, 1, "rook"), style = "W")
}



