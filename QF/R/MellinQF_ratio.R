#'Mellin Transform of the Independent Positive QFs Ratio
#'
#'The function computes the Mellin transform of the ratio of independent and positive definite quadratic forms producing a \code{MellinQF_ratio} object.
#'The output can be used to evaluate the density, cumulative and quantile functions of the target quadratic form.
#'
#'@param lambdas_num vector of positive weights for the numerator.
#'@param lambdas_den vector of positive weights for the denominator.
#'@param etas_num vector of non-centrality parameters for the numerator. Default all zeros (central chi square).
#'@param etas_den vector of non-centrality parameters for the denominator. Default all zeros (central chi square).
#'@param eps required absolute error for density and cumulative functions.
#'@param rho distribution total probability mass for which it is desired to keep the error \code{eps}.
#'@param maxit_comp maximum number of iterations.
#'@param eps_quant required numerical error for quantile computation.
#'@param maxit_quant maximum number of iterations before stopping the quantile computation.
#'@param lambdas_tol maximum value admitted for the weight skewness for both the numerator and the denominator. When it is not NULL (default), elements of lambdas such that the ratio max(lambdas)/lambdas is greater than the specified value are removed.
#'
#'
#'@details
#'
#'The Mellin transform of the ratio of two independent quadratic forms having positive weights \code{lambdas_num} and \code{lambdas_den}
#' and non-centrality parameters \code{etas_num} and \code{etas_den} is computed
#'exploiting the density formulation by Ruben (1962).  The numerical error is controlled in order to provide the requested precision (\code{eps}) for the
#'interval of quantiles that contains the specified total probability \code{rho}.
#'
#'The argument \code{eps_quant} controls the relative precision requested for the computation of quantiles that determine the range in which the error \code{eps} is
#'guaranteed, whereas \code{maxit_quant} sets the maximum number of Newton-Raphson iterations of the algorithm.
#'
#'@return
#'
#'The function returns an object of the class \code{MellinQF_ratio} that contains information on the Mellin transform
#'of the ratio of two linear combinations of positively weighted chi-square random variables. This information can be used in order to
#'evaluate the density, cumulative distribution and quantile functions of the desired quadratic form.
#'
#'An object of the class \code{MellinQF_ratio} has the following components:
#'
#' * \code{range_q}: the range of quantiles that contains the specified mass of probability \code{rho} in which it
#' is possible to compute density and CDF preserving the error level \code{eps}.
#'
#' * \code{Mellin}: a list containing the values of the Mellin transform (\code{Mellin}),
#' the corresponding evaluation points (\code{z}), the integration step \code{delta} and the lowest numerator weight (\code{lambda_min}).
#'
#' * the inputs \code{rho}, \code{lambdas_num}, \code{lambdas_den}, \code{etas_num}, \code{etas_den}, \code{eps} needed for CDF, PDF and quantile function computation.
#'
#'
#'@source
#'
#'Ruben, Harold. "Probability content of regions under spherical normal distributions, IV:
#'The distribution of homogeneous and non-homogeneous quadratic functions of normal variables."
#'The Annals of Mathematical Statistics 33.2 (1962): 542-570.
#'
#'@seealso
#'The function \code{\link{print.MellinQF_ratio}} can be used to summarize the basic information on the Mellin transform.
#'
#'The object can be used in the function \code{\link{dQF_ratio}} to compute the density function of the QFs ratio,
#'\code{\link{pQF_ratio}} for the CDF and \code{\link{qQF_ratio}} for the quantile function.
#'
#'@examples
#'
#'\donttest{
#' library(QF)
#' # Definition of the QFs
#' lambdas_QF_num <- c(rep(7, 6),rep(3, 2))
#' etas_QF_num <- c(rep(6, 6), rep(2, 2))
#' lambdas_QF_den <- c(0.6, 0.3, 0.1)
#' # Computation Mellin transform
#' eps <- 1e-7
#' rho <- 0.999
#' Mellin_ratio <- compute_MellinQF_ratio(lambdas_QF_num, lambdas_QF_den,
#'                                        etas_QF_num, eps = eps, rho = rho)
#' print(Mellin_ratio)
#' }
#'
#'
#'
#'@export


compute_MellinQF_ratio <- function(lambdas_num, lambdas_den,
                                   etas_num = rep(0, length(lambdas_num)),
                                   etas_den = rep(0, length(lambdas_den)),
                                   eps = 1e-6, rho= 1-1e-4,
                                   maxit_comp = 1e5,
                                   eps_quant = 1e-6, maxit_quant = 1e4, lambdas_tol = NULL) {
  if(sum(lambdas_num < 0) != 0 || sum(lambdas_den < 0) != 0 ||
     sum(etas_num < 0) != 0 || sum(etas_den < 0) != 0){
    stop("All elements of 'lambdas_num', 'lambdas_den', 'etas_num' and 'etas_den' must be positive")
  }
  lambdas_num <- lambdas_num[lambdas_num > 0]
  lambdas_den <- lambdas_den[lambdas_den > 0]
  etas_num <- etas_num[lambdas_num > 0]
  etas_den <- etas_den[lambdas_den > 0]

  # select eigenvalues according to maximum skewnell level lambdas_tol
  if(!is.null(lambdas_tol)){
    message(paste0("Removed ", sum(max(lambdas_num) / lambdas_num >= lambdas_tol) + sum(max(lambdas_den) / lambdas_den >= lambdas_tol), " weights using the specified 'lambdas_tol'."))
    lambdas_num <- lambdas_num[max(lambdas_num) / lambdas_num < lambdas_tol]
    etas_num <- etas_num[max(lambdas_num) / lambdas_num < lambdas_tol]
    lambdas_den <- lambdas_den[max(lambdas_den) / lambdas_den < lambdas_tol]
    etas_den <- etas_den[max(lambdas_den) / lambdas_den < lambdas_tol]
  }

  eigen_skew <- max(c(max(lambdas_num) / min(lambdas_num), max(lambdas_den) / min(lambdas_den)))
  if(eigen_skew > 1e5){
    warning(paste0("The skewness of the weights, i.e. max(lambdas)/min(lambdas) is ",eigen_skew,".\n
                   Consider to specify an appropriate value for argument 'lambdas_tol'."))
  }

  n_lambdas_num <- length(lambdas_num)
  n_lambdas_den <- length(lambdas_den)
  if(rho < 0 || rho > 1){
    stop("'rho' must be between 0 and 1")
  }
  if(sum(c(eps, maxit_comp,
           eps_quant, maxit_quant)<0)!=0){
    stop("All the parameters 'delta', 'step_delta', 'eps', 'maxit',
           'eps_quant', 'maxit_quant' must be strictly positive")
  }
  if(n_lambdas_num <= 2){
    h <- .55
  }else{
    if(n_lambdas_num == 3){
      h <- 1/3
    }else{
      h <- .1
    }
  }
  delta <- .1i
  step_delta <- .05
  maxit_delta <- 50

  mellin <-     .Call(`_QF_get_mellin_QF_ratio`, lambdas_num, lambdas_den,
                      etas_num, etas_den, rho, h, delta, eps, eps_quant,
                      maxit_comp, maxit_quant, maxit_delta, step_delta)

  mellin[["lambdas_num"]] <- lambdas_num
  mellin[["lambdas_den"]] <- lambdas_den
  mellin[["etas_num"]] <- etas_num
  mellin[["etas_den"]] <- etas_den
  mellin[["eps"]] <- eps

  attr(mellin, "class") <- "MellinQF_ratio"
  return(mellin)
}

#' Print \code{MellinQF_ratio}
#'
#'It allows to visualize useful information on the \code{MellinQF_ratio} object.
#'
#'@param x \code{MellinQF_ratio} object.
#'@param digits number of digits to display in the output.
#'@param ... potentially more arguments passed to methods.
#'
#'@export

print.MellinQF_ratio <- function(x, digits = 3L, ...){
  cat("-------------\n")
  cat("Mellin transform of a ratio of positive and independent quadratic forms having weights 'lambdas_num':\n")
  cat(paste0("Min: ",round(min(x$lambdas_num), digits = digits),"\n"))
  cat(paste0("Mean: ",round(mean(x$lambdas_num), digits = digits),"\n"))
  cat(paste0("Max: ",round(max(x$lambdas_num), digits = digits),"\n"))
  cat("'lambdas_den':\n")
  cat(paste0("Min: ",round(min(x$lambdas_den), digits = digits),"\n"))
  cat(paste0("Mean: ",round(mean(x$lambdas_den), digits = digits),"\n"))
  cat(paste0("Max: ",round(max(x$lambdas_den), digits = digits),"\n"))
  cat("and non-centrality parameters 'etas_num':\n")
  cat(paste0("Min: ",round(min(x$etas_num), digits = digits),"\n"))
  cat(paste0("Mean: ",round(mean(x$etas_num), digits = digits),"\n"))
  cat(paste0("Max: ",round(max(x$etas_num), digits = digits),"\n"))
  cat("'etas_den':\n")
  cat(paste0("Min: ",round(min(x$etas_den), digits = digits),"\n"))
  cat(paste0("Mean: ",round(mean(x$etas_den), digits = digits),"\n"))
  cat(paste0("Max: ",round(max(x$etas_den), digits = digits),"\n"))
  cat("-------------\n")
  cat("It can be used to compute density and cumulative functions when:\n")
  cat(paste0(round(x$range_q[1], digits = digits),"<q<",round(x$range_q[2], digits = digits),";\n"))
  cat("and the quantile function when:\n")
  cat(paste0((1 - x$rho) / 2, "<p<", 1 - (1 - x$rho) / 2,";\n"))
  cat(paste0("with an error 'eps': ", x$eps,"\n"))
  cat("-------------\n")
}





