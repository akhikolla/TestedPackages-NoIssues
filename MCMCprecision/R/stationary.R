#' Precision of stationary distribution for discrete MCMC variables
#'
#' Transdimensional MCMC methods include a discrete model-indicator variable \eqn{z}
#' with a fixed but unknown stationary distribution with probabilities \eqn{\pi}
#' (i.e., the model posterior probabiltiies). The function \code{stationary}
#' draws posterior samples to assess the estimation uncertainty of \eqn{\pi}.
#'
#' @param z MCMC output for the discrete indicator variable with numerical,
#'     character, or factor labels (can also be a \code{\link[coda]{mcmc.list}}
#'     or a matrix with one MCMC chain per column).
#' @param N the observed \code{\link{transitions}} matrix (if supplied, \code{z} is ignored).
#'     A quadratic matrix with sampled transition frequencies
#'     (\code{N[i,j]} = number of switches from \code{z[t]=i} to \code{z[t+1]=j}).
#' @param labels optional: vector of labels for complete set of models
#'     (e.g., models not sampled in the chain \code{z}). If \code{epsilon=0},
#'     this does not affect inferences due to the improper Dirichlet(0,..,0) prior.
#' @param sample number of posterior samples to be drawn for the stationary distribution \eqn{\pi}.
#' @param epsilon prior parameter for the rows of the estimated transition matrix \eqn{P}:
#'     \eqn{P[i,]} ~ Dirichlet\eqn{(\epsilon, ..., \epsilon)}.
#'     The default \code{epsilon="1/M"} (with M = number of sampled models) provides
#'     estimates close to the i.i.d. estimates and is numerically stable.
#'     The alternative \code{epsilon=0} minimizes the impact of the prior and
#'     renders non-sampled models irrelevant.
#'     If \code{method="iid"} (ignores dependencies), a Dirichlet prior is assumed on the stationary
#'     distribution \eqn{\pi} instead of the rows of the transition matrix \eqn{P}.
#' @param summary whether the output should be summarized.
#'     If \code{FALSE}, posterior samples of the stationary probabilities are returned.
#' @param cpu number of CPUs used for parallel sampling.
#'     Will only speed up computations for large numbers of models
#'     (i.e., for large transition matrices).
#' @param progress whether to show a progress bar (not functional for \code{cpu>1})
#' @param method how to compute eigenvectors:
#' \itemize{
#'   \item \code{"arma"} (default): Uses \code{RcppArmadillo::eig_gen}.
#'   \item \code{"base"}: Uses \code{base::\link[base]{eigen}}, which might be more stable,
#'                        but also much slower than \code{"arma"} for small transition matrices.
#'   \item \code{"eigen"}: Uses package \code{RcppEigen::EigenSolver}
#'   \item \code{"armas"}: Uses sparse matrices with \code{RcppArmadillo::eigs_gen},
#'                         which can be faster for very large number of models
#'                         if \code{epsilon=0} (might be numerically unstable).
#'   \item \code{"iid"}: Assumes i.i.d. sampling of the model indicator variable \code{z}.
#'                       This is only implemented as a benchmark, because results cannot
#'                       be trusted if the samples \code{z} are correlated (which is usually
#'                       the case for transdimensional MCMC output)
#'  }
#' @param digits number of digits that are used for checking whether the first
#'    eigenvalue is equal to 1 (any difference must be due to low numerical precision)
#'
#' @details
#' The method draws independent posterior samples of the transition matrix \eqn{P}
#' for the discrete-valued indicator variable \code{z} (usually, a sequence of sampled models).
#' For each row of the transition matrix, a Dirichlet\eqn{(\epsilon,...,\epsilon)}
#' prior is assumed, resulting in a conjugate Dirichlet posterior.
#' For each sample, the eigenvector with eigenvalue 1 is computed and normalized.
#' These (independent) posterior samples can be used to assess the estimation
#' uncertainty in the stationary distribution \eqn{pi} of interest
#' (e.g., the model posterior probabilities) and to estimate the effective sample size
#' (see \code{\link{summary.stationary}}).
#'
#' @return default: a summary for the posterior distribution of the model
#' posterior probabilities (i.e., the fixed but unknown stationary distribution of \code{z}).
#' If \code{summary=FALSE}, posterior samples for \eqn{pi} are returned.
#'
#' @seealso \code{\link{best_models}}, \code{\link{summary.stationary}}
#'
#' @examples
#' # data-generating transition matrix
#' P <- matrix(c(.1,.5,.4,
#'               0, .5,.5,
#'               .9,.1,0), ncol = 3, byrow=TRUE)
#'
#' # input: sequence of sampled models
#' z <- rmarkov(500, P)
#' stationary(z)
#'
#' # input: transition frequencies
#' N <- transitions(z)
#' samples <- stationary(N = N, summary = FALSE)
#'
#' # summaries:
#' best_models(samples, k = 3)
#' summary(samples)
#' @export
stationary <- function (z, N, labels, sample = 1000, epsilon = "1/M",
                        cpu = 1, method = "arma", digits = 6,
                        progress = TRUE, summary = TRUE){

  method <- match.arg(method, c("arma", "armas", "eigen", "base", "iid"))
  if (missing(labels))
    labels <- NULL

  if (!missing(N) && !is.null(N)){
    if (ncol(N) != nrow(N) || any(N<0))
      stop ("The transition matrix 'N' has negative values.")
    tab <- as.matrix(N)
    # postpred: tab2 <- c(array(NA, rep(ncol(tab), 3)))
  } else {
    tab <- transitions(z, labels=labels)
    # postpred: tab2 <- c(transitions(z, labels = labels, order = 2))
  }
  if (method == "armas")
    tab <- Matrix(tab, sparse = TRUE)

  M <- nrow(tab)
  if (epsilon == "1/M"){
    epsilon <- 1/M
  } else if (epsilon <0){
    stop ("'epsilon' must be zero or positive.")
  }

  labels <- rownames(tab)
  if (cpu == 1){
    mcmc <- stationary_samples(tab, sample, epsilon, method, digits, progress)
  } else {
    cl <- makeCluster(cpu)
    sample.cpu <- ceiling(sample/cpu)
    mcmc <- do.call("rbind",
                    clusterCall(cl, stationary_samples, method = method,
                                tab = tab, sample = sample.cpu,
                                epsilon = epsilon, digits = digits, progress = FALSE))
    stopCluster(cl)
  }

  # # posterior predictive checks (only with stationaryArma)
  # # (requires the raw sequence z to compute a second-order transition matrix)
  # if (ncol(mcmc) > M){
  #   x2 <- mcmc[,M + 1:2]
  #   mcmc <- mcmc[,1:M]
  #   postpred <- c("X2.obs" = median(x2[,1], na.rm = TRUE),
  #                 "X2.pred" = median(x2[,2], na.rm = TRUE),
  #                 "p-value" = mean(x2[,1] < x2[,2], na.rm = TRUE))
  # }
  # attr(mcmc, "postpred") <- postpred

  mcmc[mcmc < 0] <- 0
  if (anyNA(mcmc))
    warning (
      "Sampled posterior model probabilities contain NAs/NaNs.\n",
      "  This might be due to low precision of the eigenvalue decomposition.\n",
      "  As a remedy, use 'digits=1' or change the Dirichlet prior (e.g., 'epsilon = .1')")

  if (!all(is.na(mcmc)) && (max(mcmc, na.rm = TRUE) == 1 ))
    warning (
      "Some models are assigned posterior probability of one.\n",
      "  This indicates numerical issues with the eigenvalue decomposition.\n",
      "  As a remedy, the Dirichlet prior can be changed, e.g., 'epsilon = .01'.")

  colnames(mcmc) <- colnames(tab)
  class(mcmc) <- c("stationary", "matrix")
  attr(mcmc, "epsilon") <- epsilon
  attr(mcmc, "method") <- method

  if (summary){
    mcmc <- summary.stationary(mcmc)
    if (method == "iid")
      mcmc$n.eff <- sum(tab) + epsilon*M
  }
  mcmc
}

## get samples
stationary_samples <- function(tab, sample, epsilon, method, digits, progress){
  if (method == "arma"){
    mcmc <- stationaryArma(tab, sample = sample, epsilon = epsilon,
                           digits = digits, progress = progress)

  } else if (method == "armas") {
    if (epsilon != 0)
      stop ("'epsilon=0' required for method='armas' (sparse matrices)")
    mcmc <- stationaryArmaSparse(tab, sample = sample,
                                 digits = digits, progress = progress)

  } else if (method == "iid"){
    mcmc <- rdirichlet(sample, rowSums(tab) + epsilon)

  } else if (method == "base"){
    mcmc <- matrix(NA, sample, ncol(tab))
    if (progress)
      prog <- txtProgressBar(0, sample, style=3, width=50)
    for (i in 1:sample){
      if (progress) setTxtProgressBar(prog, i)
      mcmc[i,] <- posterior.sample(1, tab=tab, epsilon=epsilon, digits = digits)
    }
    if (progress) close(prog)

  } else if (method == "eigen") {
    mcmc <- stationaryEigen(tab, sample = sample, epsilon = epsilon,
                            digits = digits, progress = progress)
    mcmc[mcmc == -99] <- NA

  } else {
    stop ("method not supported.")
  }
  mcmc
}
