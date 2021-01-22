#'Cumulative Distribution Function of the Dependent QFs Ratio
#'
#'The function computes the CDF of the ratio of two dependent and possibly indefinite quadratic forms.
#'
#'@param q vector of quantiles.
#'@param lambdas vector of eigenvalues of the matrix (A-qB).
#'@param A matrix of the numerator QF. If not specified but \code{B} is passed, it is assumed to be the identity.
#'@param B matrix of the numerator QF. If not specified but \code{A} is passed, it is assumed to be the identity.
#'@param eps requested absolute error.
#'@param maxit_comp maximum number of iterations.
#'@param lambdas_tol maximum value admitted for the weight skewness for both the numerator and the denominator. When it is not NULL (default), elements of lambdas such that the ratio max(lambdas)/lambdas is greater than the specified value are removed.
#'
#'@details
#'The distribution function of the following ratio of dependent quadratic forms is computed:
#'\deqn{P\left(\frac{Y^TAY }{Y^TBY}<q\right),}
#'where \eqn{Y\sim N(0, I)}.
#'
#'The transformation to the following indefinite quadratic form is exploited:
#'\deqn{P\left(Y^T(A-qB)Y <0\right).}
#'
#'The following inputs can be provided:
#'
#'\itemize{
#'   \item{vector \code{lambdas} that contains the eigenvalues of the matrix \eqn{(A-qB)}}. Input \code{q} is ignored.
#'   \item{matrix \code{A} and/or matrix \code{B}: in these cases \code{q} is required to be not null and an eventual missing specification of one matrix make it equal to the identity.}
#'}
#'
#'
#'
#'@return
#' The values of the CDF at quantiles \code{q}.
#'
#'
#'@export


pQF_depratio <- function(q=NULL, lambdas=NULL,
                         A=NULL,B=NULL,
                         eps = 1e-6, maxit_comp = 1e5,
                         lambdas_tol=NULL) {
  # no value specified
  if(is.null(lambdas) && is.null(A)&&is.null(B)){
    stop("At least one among 'lambdas', 'A' and 'B' must be specified")
  }
  #only lambdas specified
  if(!is.null(lambdas) && ((!is.null(q))|(!is.null(A))|(!is.null(B)))){
    warning("When weights 'lambdas' are provided 'q', 'A' and 'B' are ignored")
  }
  if(!is.null(lambdas)){
    lambdas_eval <- lambdas
  }

  #lambdas not specified and at least one between A and B provided
  if(is.null(lambdas) && (!is.null(A)||(!is.null(B)))){
    if(!is.null(A)){
      if(ncol(A)!=nrow(A)){
        stop("'A' must be a square matrix")
      }
    }
    if(!is.null(B)){
      if(ncol(B)!=nrow(B)){
        stop("'B' must be a square matrix")
      }
    }
    if(is.null(A)){
      lambdas_B<-eigen(B)$val
      lambdas_eval<-1-q*lambdas_B
    }else if(is.null(B)){
      lambdas_eval<-eigen(A)$val-q
    }else{
      if(!all.equal(dim(A),dim(B))){
        stop("'A' and 'B' must have the same dimensions")
      }
      lambdas_eval<-eigen(A-q*B)$val
    }
  }

  if(sum(c(eps, maxit_comp)<0)!=0){
    stop("All the parameters 'eps' and 'maxit_comp' must be strictly positive")
  }
  if(length(lambdas_eval[lambdas_eval>0]) <= 2){
    h <- .55
  }else{
    if(length(lambdas_eval[lambdas_eval>0]) == 3){
      h <- 1/3
    }else{
      h <- .1
    }
  }


  lambdas_1<-lambdas_eval[lambdas_eval>0];
  lambdas_2<-lambdas_eval[lambdas_eval<0];
  lambdas_2<--lambdas_2;

  if(!is.null(lambdas_tol)){
    message(paste0("Removed ", sum(max(lambdas_1)/lambdas_1 >= lambdas_tol) + sum(max(lambdas_2)/lambdas_2 >= lambdas_tol), " weights using the specified 'lambdas_tol'."))
    lambdas_1<-lambdas_1[max(lambdas_1)/lambdas_1<lambdas_tol];
    lambdas_2<-lambdas_2[max(lambdas_2)/lambdas_2<lambdas_tol];
  }

  eigen_skew <- max(c(max(lambdas_1) / min(lambdas_1), max(lambdas_2) / min(lambdas_2)))
  if(eigen_skew > 1e5){
    warning(paste0("The skewness of the weights, i.e. max(lambdas)/min(lambdas) is ",eigen_skew,".\n
                   Consider to specify an appropriate value for argument 'lambdas_tol'."))
  }

  delta <- .1i
  step_delta <- .05
  maxit_delta<-20

  F_val <- .Call(`_QF_pQF_depratio_c`, lambdas_1, lambdas_2,
                 h, delta, eps, maxit_comp,
                 maxit_delta)
  return(F_val)
}
