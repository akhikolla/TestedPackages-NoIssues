#' Exploratory reduced Reparameterized Unified Model (ErRUM)
#'
#' Obtains samples from posterior distribution for the Exploratory
#' reduced Reparameterized Unified Model (ErRUM).
#'
#' @param y            Binary responses to assessments in `matrix`
#'                     form with dimensions \eqn{N \times J}{N x J}.
#' @param k            Number of Attribute Levels as a positive `integer`.
#' @param burnin       Number of Observations to discard on the chain.
#' @param chain_length Length of the MCMC chain
#' @param verbose      Display estimation progress updates.
#' @param X,v0,v1,cv0,cv1,bnu Additional tuning parameters
#'
#' @return
#'
#' An `errum` object that has:
#' - `PISTAR`
#' - `RSTAR`
#' - `PIs`
#' - `QS`
#' - `m_Delta`
#' - `Delta_biject`
#' - `M2`
#' - `M1`
#' - `NUS`
#'
#' @seealso
#' [simcdm::attribute_bijection()],
#' [simcdm::sim_rrum_items()]
#'
#' @export
#' @examples
#' # Setup Simulation Parameters
#' N = 5
#' K = 3
#' J = 30
#' # Note:
#' # Sample size has been reduced to create a minimally
#' # viable example that can be run during CRAN's automatic check.
#' # Please make sure to have a larger sample size of around 3,000.
#'
#' # Sample true attribute profiles
#' Z         = matrix(rnorm(N * K), N, K)
#' Sig       = matrix(.5, K, K)
#' diag(Sig) = 1
#' theta     = Z %*% chol(Sig)
#'
#' thvals    = matrix(qnorm((1:K) / (K + 1)),
#'                    N, K, byrow = TRUE)
#'
#' Alphas    = 1 * (theta > thvals)
#'
#' # Defining matrix of possible attribute profiles
#' As = as.matrix(expand.grid(c(0, 1), c(0, 1), c(0, 1)))
#' Q = rbind(As[rep(c(2, 3, 5), 4),],
#'           As[rep(c(4, 6, 7), 4),],
#'           As[rep(8, 6),])
#'
#' # Use simulation functions available in simcdm
#' if (requireNamespace("simcdm", quietly = TRUE)) {
#'
#' a = As %*% simcdm::attribute_bijection(K)
#' As = As[a + 1,]
#'
#' # Setting item parameters
#' pistar = rep(.9, J)
#' rstar = matrix(.6, J, K) * Q
#'
#' # Simulate data under rRUM model
#' Y = simcdm::sim_rrum_items(Q, rstar, pistar, Alphas)
#'
#' # Estimation Settings
#' chainLength = 10000  # Run with 20000
#' burnin = chainLength / 2
#'
#' # Gibbs Estimation
#' model = errum(Y, K, burnin, chainLength)
#' }
errum = function(y, k = 3,
                 burnin = 1000, chain_length = 10000,
                 verbose = FALSE,
                 X = matrix(1, nrow = ncol(y)),
                 v0 = 4, v1 = 2,
                 cv0 = .1, cv1 = 10,
                 bnu = 16) {

    message("Beginning the estimation procedure... ")
    o = rRUM_mvnQ_Gibbs(y, k,
                        X, v0, v1, cv0, cv1, bnu,
                        burnin, chain_length, verbose)

    class(o) = c("errum")
    o
}
