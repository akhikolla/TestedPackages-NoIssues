#' MLE for stationary distribution of discrete MCMC variables
#'
#' Maximum-likelihood estimation of stationary distribution \eqn{\pi} based on (a) a sampled trajectory \eqn{z} of a model-indicator variable or (b) a sampled transition count matrix \eqn{N}.
#'
#' @inheritParams stationary
#'
#' @param method Different types of MLEs:
#' \itemize{
#'   \item \code{"iid"}: Assumes i.i.d. sampling of the model indicator variable \code{z} and estimates \eqn{\pi} as the relative frequencies each model was sampled.
#'   \item \code{"rev"}: Estimate stationary distribution under the constraint that the transition matrix is reversible (i.e., fulfills detailed balance) based on the iterative fixed-point algorithm proposed by Trendelkamp-Schroer et al. (2015)
#'   \item \code{"eigen"}: Computes the first left-eigenvector (normalized to sum to 1) of the sampled transition matrix
#'  }
#' @param abstol absolute convergence tolerance (only for \code{method = "rev"})
#' @param maxit maximum number of iterations (only for \code{method = "rev"})
#'
#' @details The estimates are implemented mainly for comparison with the Bayesian sampling approach implemented in \code{\link{stationary}}, which quantify estimation uncertainty (i.e., posterior SD) of the posterior model probability estimates.
#'
#' @return a vector with posterior model probability estimates
#' @seealso \code{\link{stationary}}
#'
#' @examples
#' P <- matrix(c(.1,.5,.4,
#'               0,.5,.5,
#'               .9,.1,0), ncol = 3, byrow=TRUE)
#' z <- rmarkov(1000, P)
#' stationary_mle(z)
#'
#' # input: transition frequency
#' tab <- transitions(z)
#' stationary_mle(N = tab)
#' @references
#' Trendelkamp-Schroer, B., Wu, H., Paul, F., & Noé, F. (2015). Estimation and uncertainty of reversible Markov models. The Journal of Chemical Physics, 143(17), 174101. \url{https://doi.org/10.1063/1.4934536}
#' @export
stationary_mle <- function (z, N, labels, method = "rev",
                            abstol = 1e-5, maxit = 1e5){

  if (missing(labels))
    labels <- NULL
  if (!missing(N) && !is.null(N)){
    if (ncol(N) != nrow(N) || any(N<0) )
      stop ("The transition matrix 'N' has negative values.")
    N <- as.matrix(N)
  } else {
    N <- transitions(z, labels = labels)
  }

  method <- match.arg(method, c("iid", "rev", "eigen"))
  if (method == "iid"){
    pi <- colSums(N)/sum(N)

  } else if (method == "rev"){
    M <- ncol(N)
    start <- stationary_mle(N = N, method = "iid")
    pi <- stationary_reversible(start, N, abstol = abstol, maxit = maxit)
    # pi.old <- rep(.5,M)
    # N.row <- rowSums(N)
    # cnt <- 0
    # while(max(abs(pi - pi.old)) > abstol && cnt<maxit){
    #   pi.old <- pi
    #   for(i in 1:M){
    #     pi[i] <-  sum( (N[i,] + N[,i])/(N.row[i]/pi.old[i] + N.row/pi.old))
    #   }
    #   cnt <- cnt + 1
    # }

  } else if (method == "eigen"){
    ev <- Re(eigen(t(N/rowSums(N)))$vectors[,1])
    pi <- ev/sum(ev)
  } else {
    stop ("Method not supported.")
  }

  names(pi) <- colnames(N)
  pi
}
