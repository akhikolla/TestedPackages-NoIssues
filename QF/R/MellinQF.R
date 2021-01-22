#'Mellin Transform of a Positive QF
#'
#'The function computes the Mellin transform of a positive definite quadratic form producing a \code{MellinQF} object.
#'The output can be used to evaluate the density, cumulative and quantile functions of
#'the target quadratic form.
#'
#'@param lambdas vector of positive weights.
#'@param etas vector of non-centrality parameters. Default all zeros (central chi square).
#'@param eps required absolute error for density and cumulative functions.
#'@param rho distribution total probability mass for which it is desired to keep the error \code{eps}.
#'@param maxit_comp maximum number of iterations.
#'@param eps_quant required numerical error for quantile computation.
#'@param maxit_quant maximum number of iterations before stopping the quantile computation.
#'@param lambdas_tol maximum value admitted for the weight skewness. When it is not NULL (default), elements of lambdas such that the ratio max(lambdas)/lambdas is greater than the specified value are removed.
#'
#'
#'@details
#'The quadratic form having positive weights \code{lambdas} and non-centrality parameters \code{etas} is considered:
#'\deqn{Q=\sum_{i=1}^r \lambda_i\chi^2_{1,\eta_i}. }
#'
#'Its Mellin transform is computed by exploiting the density formulation by Ruben (1962).
#'The numerical error is controlled in order to provide the requested precision (\code{eps}) for the
#'interval of quantiles that contains the specified total probability \code{rho}.
#'
#'The argument \code{eps_quant} controls the relative precision requested for the
#'computation of quantiles that determine the range in which the error \code{eps} is
#'guaranteed, whereas \code{maxit_quant} sets the maximum number of Newton-Raphson iterations of the algorithm.
#'
#'@return
#'
#'The function returns an object of the class \code{MellinQF} that contains information on the Mellin transform
#'of a linear combination of positively weighted chi-square random variables. This information can be used in order to
#'evaluate the density, cumulative distribution and quantile functions.
#'
#' An object of the class \code{MellinQF} has the following components:
#'
#' * \code{range_q}: the range of quantiles that contains the specified mass of probability \code{rho} in which it
#' is possible to compute density and CDF preserving the error level \code{eps}.
#'
#' * \code{Mellin}: a list containing the values of the Mellin transform (\code{Mellin}),
#' the corresponding evaluation points (\code{z}), the integration step \code{delta} and the lowest weight (\code{lambda_min}).
#'
#' * the inputs \code{rho}, \code{lambdas}, \code{etas}, \code{eps} needed for CDF, PDF and quantile function computation.
#'
#'
#'@source
#'
#'Ruben, Harold. "Probability content of regions under spherical normal distributions, IV:
#'The distribution of homogeneous and non-homogeneous quadratic functions of normal variables."
#'The Annals of Mathematical Statistics 33.2 (1962): 542-570.
#'
#'@seealso
#'The function \code{\link{print.MellinQF}} can be used to summarize the basic information on the Mellin transform.
#'
#'The object can be used in the function \code{\link{dQF}} to compute the density function of the QF,
#'\code{\link{pQF}} for the CDF and \code{\link{qQF}} for the quantile function.
#'
#'
#'@examples
#'
#'\donttest{
#' library(QF)
#' # Definition of the QF
#' lambdas_QF <- c(rep(7, 6),rep(3, 2))
#' etas_QF <- c(rep(6, 6), rep(2, 2))
#' # Computation Mellin transform
#' eps <- 1e-7
#' rho <- 0.999
#' Mellin <- compute_MellinQF(lambdas_QF, etas_QF, eps = eps, rho = rho)
#' print(Mellin)
#' }
#'
#'
#'@export


compute_MellinQF <- function(lambdas, etas = rep(0, length(lambdas)),
                             eps = 1e-6, rho= 1-1e-4,
                          maxit_comp = 1e5,
                          eps_quant = 1e-6, maxit_quant = 1e4, lambdas_tol = NULL) {
  if(sum(lambdas<0)!=0 || sum(etas<0)!=0){
    stop("All elements of 'lambdas' and 'etas' must be positive")
  }

  if(rho < 0 || rho > 1){
    stop("'rho' must be between 0 and 1")
  }
  if(sum(c(eps, maxit_comp,
           eps_quant, maxit_quant)<0)!=0){
    stop("All the parameters 'delta', 'step_delta', 'eps', 'maxit',
           'eps_quant', 'maxit_quant' must be strictly positive")
  }

  lambdas <- lambdas[lambdas > 0]
  etas <- etas[lambdas > 0]

  # select eigenvalues according to maximum skewnell level lambdas_tol
  if(!is.null(lambdas_tol)){
    message(paste0("Removed ", sum(max(lambdas) / lambdas >= lambdas_tol), " weights using the specified 'lambdas_tol'."))
    lambdas <- lambdas[max(lambdas) / lambdas < lambdas_tol]
    etas <- etas[max(lambdas) / lambdas < lambdas_tol]
  }

  eigen_skew <- max(lambdas) / min(lambdas)
  if(eigen_skew > 1e5){
    warning(paste0("The skewness of the weights, i.e. max(lambdas)/min(lambdas) is ",eigen_skew,".\n
                   Consider to specify an appropriate value for argument 'lambdas_tol'."))
  }

  n_lambdas <- length(lambdas)
  if(n_lambdas <= 2){
    h <- .55
  }else{
    if(n_lambdas == 3){
      h <- 1/3
    }else{
      h <- .1
    }
  }
  delta <- .1i
  step_delta <- .05
  maxit_delta <- 50


  mellin <- .Call(`_QF_get_mellin_QF`, lambdas, etas,
                  rho, h, delta, eps, eps_quant,
                  maxit_comp, maxit_quant, maxit_delta,
                  step_delta)

  mellin[["lambdas"]] <- lambdas
  mellin[["etas"]] <- etas
  mellin[["eps"]] <- eps

  attr(mellin, "class") <- "MellinQF"
  return(mellin)
}

#' Printing \code{MellinQF}
#'
#'It allows to visualize useful information on the \code{MellinQF} object.
#'
#'@param x \code{MellinQF} object.
#'@param digits number of digits to display in the output.
#'@param ... potentially more arguments passed to methods.
#'
#'@export

print.MellinQF <- function(x, digits = 3L, ...){
  cat("-------------\n")
  cat("Mellin transform of a quadratic form having weights 'lambdas':\n")
  cat(paste0("Min: ",round(min(x$lambdas), digits = digits),"\n"))
  cat(paste0("Mean: ",round(mean(x$lambdas), digits = digits),"\n"))
  cat(paste0("Max: ",round(max(x$lambdas), digits = digits),"\n"))
  cat("and non-centrality parameters 'etas':\n")
  cat(paste0("Min: ",round(min(x$etas), digits = digits),"\n"))
  cat(paste0("Mean: ",round(mean(x$etas), digits = digits),"\n"))
  cat(paste0("Max: ",round(max(x$etas), digits = digits),"\n"))
  cat("-------------\n")
  cat("It can be used to compute density and cumulative functions when:\n")
  cat(paste0(round(x$range_q[1], digits = digits),"<q<",round(x$range_q[2], digits = digits),";\n"))
  cat("and the quantile function when:\n")
  cat(paste0((1 - x$rho) / 2, "<p<", 1 - (1 - x$rho) / 2,";\n"))
  cat(paste0("with an error 'eps': ", x$eps,"\n"))
  cat("-------------\n")
}






