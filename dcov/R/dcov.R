#' Distance covariance
#'
#' This method implements the method to compute the value of distance covariance
#' proposed by \cite{Székely et al.(2007)} and \cite{Székely and Rizzo(2013)}
#' by Armadillo library. For distance covariance between two one dimensional
#' variables, the fast algorithm proposed by \cite{Huo and Székely(2016)} is used.
#' @aliases dcov
#' @aliases dcor
#' @param x the matrix of x
#' @param y the matrix of y
#' @param type "V" or "U", for V- or U-statistics of distance covariance or
#' correlation. The default value is "V".
#' @note Note that the result of \code{dcov(x,y,"V")} and \code{dcor(x,y,"V")}
#'  is same with the result of energy::dcov(x,y)^2 and energy::dcor(x,y)^2.
#'  The result of \code{dcov(x,y,'U')} and \code{dcor(x,y,'U')} is same with
#'  the result of \code{energy::dcovU(x,y)} and \code{energy::bcdcor(x,y)}.
#' @seealso dcov2d
#' @references Székely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and testing dependence by correlation of distances. The annals of statistics, 35(6), 2769-2794.
#' @references Székely, G. J., & Rizzo, M. L. (2013). The distance correlation t-test of independence in high dimension. Journal of Multivariate Analysis, 117, 193-213.
#' @references Huo, X., & Székely, G. J. (2016). Fast computing for distance covariance. Technometrics, 58(4), 435-447.
#' @useDynLib dcov
#' @import Rcpp
#' @export
#' @rdname dcov
#' @examples
#' x = matrix(rnorm(200),100,2)
#' y = matrix(rnorm(200),100,2)
#' dcov(x,y)
#' dcor(x,y)
#'
dcov<-function(x,y,type=c('V','U')){

  type <- match.arg(type)
  if(!is.matrix(x)) x = as.matrix(x)
  if(!is.matrix(y)) y = as.matrix(y)
  if(nrow(x)!=nrow(y)) stop("x and y should have same numeber of rows.")
  .Call('_dcov_dcov', PACKAGE = 'dcov', x, y, type)

}

#' @rdname dcov
#' @export

dcor<-function(x,y,type=c('V','U')){

  type <- match.arg(type)
  if(!is.matrix(x)) x = as.matrix(x)
  if(!is.matrix(y)) y = as.matrix(y)
  if(nrow(x)!=nrow(y)) stop("x and y should have same numeber of rows.")
  .Call('_dcov_dcor', PACKAGE = 'dcov', x, y, type)

}
