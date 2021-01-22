#' Projection covariance between two random vectors
#' This function implements the projection correlation in \cite{Zhu et al.(2017)}
#' @param x the matrix of x
#' @param y the matrix of y
#' @references Zhu, L., Xu, K., Li, R., & Zhong, W. (2017). Projection correlation between two random vectors. Biometrika, 104(4), 829-843.
#' @export
#' @examples
#' x = matrix(rnorm(200),100,2)
#' y = matrix(rnorm(200),100,2)
#' pcov(x,y)
pcov<-function(x,y){
  if(!is.matrix(x)) x = as.matrix(x)
  if(!is.matrix(y)) y = as.matrix(y)
  if(nrow(x)!=nrow(y)) stop("x and y should have same numeber of rows.")
  .Call('_dcov_pcovCpp', PACKAGE = 'dcov', x, y)
}
