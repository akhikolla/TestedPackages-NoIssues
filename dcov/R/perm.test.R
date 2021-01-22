#' Permutation test of distance correlation and partial distance correlation
#'
#' Simple independence test based on data permutation using distance correlation
#' and partial distance correlation.
#' @aliases dcor.test
#' @aliases pdcor.test
#' @param x the data of x
#' @param y the data of y
#' @param z the data of controlling variables. Given z, pdcor between x and y is
#' calculated.
#' @param R the number of replicates
#' @param type "U" or "V"
#' @examples
#' n = 200
#' z = rnorm(n)
#' x = rnorm(n)*z
#' y = rnorm(n)*z
#' res1 = dcor.test(x,y,R=500)
#' res2 = pdcor.test(x,y,z,R=500)
#'

#' @export
#' @rdname perm_test
dcor.test<-function(x,y,R=500,type=c('V','U')){

  type <- match.arg(type)
  if(!is.matrix(x)) x = as.matrix(x)
  if(!is.matrix(y)) y = as.matrix(y)
  if(nrow(x)!=nrow(y)) stop("x and y should have same numeber of rows.")
  if(R<1) stop('error : R should be larger than 0.')
  .Call('_dcov_dcor_test', PACKAGE = 'dcov', x, y, R, type)

}

#' @export
#' @rdname perm_test
pdcor.test<-function(x,y,z,R=500,type=c('U','V')){

  type <- match.arg(type)
  if(!is.matrix(x)) x = as.matrix(x)
  if(!is.matrix(y)) y = as.matrix(y)
  if(!is.matrix(z)) z = as.matrix(z)
  if(nrow(x)!=nrow(y)||nrow(x)!=nrow(z)) 
    stop("x, y and z should have same numeber of rows.")
  if(R<1) stop('error : R should be larger than 0.')
  .Call('_dcov_pdcor_test', PACKAGE = 'dcov', x, y, z, R, type)

}

