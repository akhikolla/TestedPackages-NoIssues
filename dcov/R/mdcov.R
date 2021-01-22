#' Marginal distance covariance
#' This function implements the method of calculating distance covariance between y
#' and each column in x
#' @aliases mdcov
#' @aliases mdcor
#' @param y the matrix of y
#' @param x the matrix of x, distance covariance is calculated for each variable
#' in x with y.
#' @param type "V" or "U", for V- or U-statistics of distance covariance or
#' correlation. The default value is "V".
#' @export
#' @rdname mdcov
#' @examples
#' n = 200; p = 10
#' y = matrix(rnorm(n*2),n,2)
#' x = matrix(rnorm(n*p),n,p)
#' res1 = mdcov(y,x)
#' res2 = numeric(p)
#' for(j in 1:p){res2[j] = dcov::dcov(y,x[,j])}
#' # res1 is same with res2
#' res1 - res2
#' res3 = mdcor(y,x)
#' res4 = numeric(p)
#' for(j in 1:p){res4[j] = dcov::dcor(y,x[,j])}
#' # res3 is same with res4
#' res3-res4

mdcov<-function(y,x,type=c('V','U')){

  type <- match.arg(type)
  if(!is.matrix(x)) x = as.matrix(x)
  if(!is.matrix(y)) y = as.matrix(y)
  if(nrow(x)!=nrow(y)) stop("x and y should have same numeber of rows.")
  .Call('_dcov_mdcov', PACKAGE = 'dcov', y, x, type)

}

#' @export
#' @rdname mdcov

mdcor<-function(y,x,type=c('V','U')){

  type <- match.arg(type)
  if(!is.matrix(x)) x = as.matrix(x)
  if(!is.matrix(y)) y = as.matrix(y)
  if(nrow(x)!=nrow(y)) stop("x and y should have same numeber of rows.")
  .Call('_dcov_mdcor', PACKAGE = 'dcov', y, x, type)

}
