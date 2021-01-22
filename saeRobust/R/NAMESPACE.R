#' @importFrom assertthat assert_that
#' @importFrom magrittr %>% %<>%
#' @importFrom stats coef coefficients dnorm fitted.values model.frame model.matrix model.response pnorm predict qnorm quantile residuals rnorm sd weights getCall formula
#' @importFrom pbapply pbreplicate
#' @import aoos
#' @import ggplot2
#' @import Matrix
#' @importMethodsFrom Matrix update
#' @import methods
#' @import modules
#' @import Rcpp
#' @useDynLib saeRobust
NULL

globalVariables(c(".", ".self", "psi", "convCrit", "maxIter", "maxIterParam", "k", "K"))

# needed for S4-dispatch:
setOldClass("rfh")
setOldClass("rfhVariance")
setOldClass(c("fitrsfh", "fitrfh"))
setOldClass(c("fitrtfh", "fitrfh"))
setOldClass(c("fitrstfh", "fitrfh"))


# Diagonal <- function(n, x = NULL) {
#   # This wrapper exists to stop the Matrix package to rise warnings when multiplying
#   # diagonals. I don't see a reason for these warnings and am convinced it is a
#   # bug in Matrix.
#   if (missing(n)) as.matrix(Matrix::Diagonal(x = x))
#   else if (!missing(n) & is.null(x)) as.matrix(Matrix::Diagonal(n = n))
#   else as.matrix(Matrix::Diagonal(n = n, x = x))
# }
