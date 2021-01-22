#' Partial distance covariance
#'
#' This method implements the method to compute the value of partial distance covariance
#' proposed by \cite{Székely and Rizzo, 2014}.
#' @aliases pdcov
#' @aliases pdcor
#' @param x the matrix of x
#' @param y the matrix of y
#' @param z the matrix of z. Given the value of z, pdcov or pdcor between x and y is calcuated.
#' @param type "V" or "U", for V- or U-statistics of partial distance covariance or
#' correlation. The default value is "U".
#' @references Székely, G. J., & Rizzo, M. L. (2014). Partial distance correlation with methods for dissimilarities. The Annals of Statistics, 42(6), 2382-2412.
#' @export
#' @rdname pdcov
#' @examples
#' z = matrix(rnorm(400),200,2)
#' x = matrix(rnorm(400),200,2)*z
#' y = matrix(rnorm(400),200,2)*z
#' pdcov(x,y,z)
#' pdcor(x,y,z)
#'

pdcov<-function(x,y,z,type=c("U","V")){
  type <- match.arg(type)
  if(!is.matrix(x)) x = as.matrix(x)
  if(!is.matrix(y)) y = as.matrix(y)
  if(!is.matrix(z)) z = as.matrix(z)
  if(nrow(x)!=nrow(y)||nrow(x)!=nrow(z)) 
    stop("x, y and z should have same numeber of rows.")
  .Call('_dcov_pdcov', PACKAGE = 'dcov', x, y, z, type)
}

#' @rdname pdcov
#' @export

pdcor<-function(x,y,z,type=c("U","V")){
  type <- match.arg(type)
  if(!is.matrix(x)) x = as.matrix(x)
  if(!is.matrix(y)) y = as.matrix(y)
  if(!is.matrix(z)) z = as.matrix(z)
  if(nrow(x)!=nrow(y)||nrow(x)!=nrow(z)) 
    stop("x, y and z should have same numeber of rows.")
  .Call('_dcov_pdcor', PACKAGE = 'dcov', x, y, z, type)
}


