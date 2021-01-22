#' @title slfm
#' @description Bayesian Sparse Latent Factor Model (SLFM) designed for the analysis of coherent patterns in gene expression data matrices. Details about the methodology
#' being applied here can be found in \emph{Duarte and Mayrink (2015)} and \emph{Duarte and Mayrink (2019)}.
#' 
#' @references 
#' 
#' DOI:10.18637/jss.v090.i09 (\href{https://www.jstatsoft.org/article/view/v090i09}{Duarte and Mayrink; 2019})
#' 
#' DOI:10.1007/978-3-319-12454-4_15 (\href{https://link.springer.com/chapter/10.1007/978-3-319-12454-4_15}{Duarte and Mayrink; 2015})
#'  
#' @param x matrix with the pre-processed data.
#' @param a positive shape parameter of the Inverse Gamma prior distribution (default = 2.1). 
#' @param b positive scale parameter of the Inverse Gamma prior distribution (default = 1.1).
#' @param gamma_a positive 1st shape parameter of the Beta prior distribution (default = 1).
#' @param gamma_b positive 2nd shape parameter of the Beta prior distribution (default = 1).
#' @param omega_0 prior variance of the spike mixture component (default = 0.01).
#' @param omega_1 prior variance of the slab mixture component (default = 10).
#' @param sample sample size to be considered for inference after the burn in period (default = 1000).
#' @param burnin size of the burn in period in the MCMC algorithm (default = sample/4).
#' @param lag lag to build the chains based on spaced draws from the Gibbs sampler (defaul = 1).
#' @param degenerate logical argument (default = FALSE) indicating whether to use the degenerate version of the mixture prior for the factor loadings.
#' @return x: data matrix.
#' @return q_star: matrix of MCMC chains for q_star parameter.
#' @return alpha: summary table of MCMC chains for alpha parameter.
#' @return lambda: summary table of MCMC chains for lambda parameter.
#' @return sigma: summary table of MCMC chains for sigma parameter.
#' @return classification: classification of each alpha (`present`, `marginal`, `absent`)
#' @export
#' @importFrom coda as.mcmc
#' @importFrom stats window
#' @importFrom Rcpp evalCpp
#' @useDynLib slfm
#' 
#' @seealso
#' \code{\link{process_matrix}}, \code{\link{plot_matrix}}
#'
#' @examples
#' mat <- matrix(rnorm(2000), nrow = 20)
#' slfm(mat, sample = 1000)
slfm <- function(
  x, a = 2.1, b = 1.1, gamma_a = 1, gamma_b = 1,
  omega_0 = 0.01, omega_1 = 10, sample = 1000, burnin = round(0.25*sample), lag = 1, degenerate = FALSE) {
  
  ite <- sample + burnin

  # Convert the x input to numeric matrix
  x <- data.matrix(x)

  if(degenerate) {
    res <- slfm_MDN(x, a, b, gamma_a, gamma_b, omega_1, burnin, lag, sample)
  } else {
    res <- slfm_MNN(x, a, b, gamma_a, gamma_b, omega_0, omega_1, burnin, lag, sample)
  }

  after_burnin <- (burnin + 1):ite

  q_star_matrix <- coda::as.mcmc(res[["qstar"]][after_burnin,])
  q_star_matrix <- window(q_star_matrix, thin = lag)
  hpds_q_star <- coda::HPDinterval(q_star_matrix)
  stats_q_star <- summary(q_star_matrix)$statistics
  alpha_clas <- format_classification(class_interval(hpds_q_star))
  alpha_clas_mean <- class_mean(stats_q_star)

  z_matrix <- coda::as.mcmc(res[["z"]][after_burnin,])
  z_matrix <- window(z_matrix, thin = lag)

  alpha_matrix <- coda::as.mcmc(res[["alpha"]][after_burnin,])
  alpha_matrix <- window(alpha_matrix, thin = lag)

  table_alpha <- alpha_estimation(res[["alpha"]][after_burnin,], alpha_clas_mean, res[["z"]][after_burnin,])

  lambda_matrix <- coda::as.mcmc(res[["lambda"]][after_burnin,])
  lambda_matrix <- window(lambda_matrix, thin = lag)
  stats_lambda <- summary(lambda_matrix)$statistics
  hpds_lambda <- coda::HPDinterval(lambda_matrix)
  table_lambda <- cbind(stats_lambda, hpds_lambda)[,-4]
  colnames(table_lambda)[4:5] = c("Upper HPD", "Lower HPD")

  sigma_matrix <- coda::as.mcmc(res[["sigma2"]][after_burnin,])
  sigma_matrix <- window(sigma_matrix, thin = lag)
  stats_sigma <- summary(sigma_matrix)$statistics
  hpds_sigma <- coda::HPDinterval(sigma_matrix)
  table_sigma <- cbind(stats_sigma, hpds_sigma)[,-4]
  colnames(table_sigma)[4:5] = c("Upper HPD", "Lower HPD")

  obj <- list(
    x = x,
    alpha = table_alpha,
    lambda = table_lambda,
    sigma2 = table_sigma,
    alpha_matrix = alpha_matrix,
    lambda_matrix = lambda_matrix,
    sigma2_matrix = sigma_matrix,
    q_star = q_star_matrix,
    z_matrix = z_matrix,
    classification = alpha_clas)
  class(obj) <- "slfm"
  obj
}

print.slfm <- function(x) {
  cat("SLFM object", "\n")
  cat("\n")
  cat("Dimensions","\n")
  cat("- alpha:", nrow(x$alpha),"\n")
  cat("- lambda:", nrow(x$lambda),"\n")
  cat("\n")
  cat("Classification:","\n")
  print(x$classification)
}

alpha_estimation <- function(x, alpha_clas, z_matrix) {
  table_list <- lapply(1:length(alpha_clas), function(i) {
    chain_indicator <- z_matrix[, i] == 1
    if(alpha_clas[i]) {
      chain <- x[chain_indicator, i]
    } else {
      chain <- x[!chain_indicator, i]
    }
    chain.mcmc <- coda::as.mcmc(chain)
    stats <- summary(chain.mcmc)$statistics
    hpds <- coda::HPDinterval(chain.mcmc)
    table <- c(stats, hpds)[-4]
  })
  table <- do.call(rbind, table_list)
  colnames(table)[4:5] = c("Upper HPD", "Lower HPD")
  table
}