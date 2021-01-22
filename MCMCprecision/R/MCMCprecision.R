#' MCMCprecision: Precision of discrete parameters in transdimensional MCMC
#'
#' MCMCprecision estimates the precision of the posterior model probabilities in
#' transdimensional Markov chain Monte Carlo methods (e.g., reversible jump MCMC
#' or product-space MCMC). This is useful for applications of transdimensional
#' MCMC such as model selection, mixtures with varying numbers of components,
#' change-point detection, capture-recapture models, phylogenetic trees,
#' variable selection, and for discrete parameters in MCMC output in general.
#'
#' The main function to assess the estimation uncertainty of discrete MCMC output is
#' is \code{\link{stationary}}.
#'
#' The method is explained in detail in Heck et al. (2019, Statistics & Computing),
#' available in the package by calling: \code{vignette("MCMCprecision")}
#'
#' @author Daniel W. Heck
#' @docType package
#'
#' @importFrom parallel clusterExport makeCluster parSapply clusterEvalQ clusterExport stopCluster clusterCall
#' @importFrom Matrix Matrix
#' @importFrom stats sd quantile rgamma na.omit runif median
#' @importFrom combinat combn
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom Rcpp evalCpp
#' @useDynLib "MCMCprecision", .registration=TRUE
#'
#' @references
#' Heck, D. W., Overstall, A. M., Gronau, Q. F., & Wagenmakers, E.-J. (2019).
#' Quantifying uncertainty in transdimensional Markov chain Monte Carlo
#' using discrete Markov models. Statistics & Computing, 29, 631–643.
#' \url{https://dx.doi.org/10.1007/s11222-018-9828-0}
#'
"_PACKAGE"
