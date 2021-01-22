#' Fast distance covariance for two bivariate variables
#'
#' This method implements the fast algorithm proposed by \cite{Huo and Székely}. The
#' result of \code{dcov2d} and \code{dcor2d} is same with the result of
#' \code{energy::dcov2d} and \code{energy::dcor2d}
#'
#' @aliases dcov2d
#' @aliases dcor2d
#' @param x the vector of x
#' @param y the vector of y
#' @param type "V" or "U", for V- or U-statistics of distance covariance or
#' correlation. The default value is "V".
#' @references Székely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and testing dependence by correlation of distances. The annals of statistics, 35(6), 2769-2794.
#' @references Székely, G. J., & Rizzo, M. L. (2013). The distance correlation t-test of independence in high dimension. Journal of Multivariate Analysis, 117, 193-213.
#' @references Huo, X., & Székely, G. J. (2016). Fast computing for distance covariance. Technometrics, 58(4), 435-447.
#' @export
#' @rdname dcov2d
#' @examples
#' x = rnorm(200)
#' y = rnorm(200)
#' dcov2d(x,y)
#' dcor2d(x,y)


dcov2d<-function(x,y,type=c("V","U")){
  type <- match.arg(type)
  if (!is.vector(x) || !is.vector(y)) {
    if (NCOL(x) > 1 || NCOL(y) > 1)
      stop("this method is only for univariate x and y")
  }
  if(length(x)!=length(y)){
    stop("x and y should have the same length.")
  }
  dcov1v1(x,y,type)
}

#' @export
#' @rdname dcov2d


dcor2d<-function(x,y,type=c("V","U")){
  type <- match.arg(type)
  if (!is.vector(x) || !is.vector(y)) {
    if (NCOL(x) > 1 || NCOL(y) > 1)
      stop("this method is only for univariate x and y")
  }
  dcor1v1(x,y,type)
}


