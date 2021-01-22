#' Centering method
#' This method implements the double centering and U-centering
#' during computing distance covariance.
#' @aliases centering
#' @aliases centering_from_data
#' @param x the matrix of x
#' @param D the pairwise distance matrix
#' @param type "V" or "U". "V" for double centering. "U" for U-centering.
#' @rdname centering
#' @export
#' @examples
#' x = matrix(rnorm(200),100,2)
#' D = as.matrix(dist(x))
#' A = centering(D,'U')
#' A = centering_from_data(x)
#'


centering<-function(D, type = c("V","U")){

  type <- match.arg(type)
  if(!is.matrix(D)||nrow(D)!=ncol(D)){
    stop("D should be a square matrix.")
  }
  Dx = D
  .Call('_dcov_centering', PACKAGE = 'dcov', Dx, type)

  return(Dx)

}


#' @export
#' @rdname centering


centering_from_data<-function(x, type = c("V","U")){

  type <- match.arg(type)
  if(!is.matrix(x)) x = as.matrix(x)
  D = matrix(0,nrow(x),nrow(x))
  .Call('_dcov_centering_from_data', PACKAGE = 'dcov', x, D, type)
  return(D)

}
