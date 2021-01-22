#' @title slfm_list
#' @description Function to fit the Bayesian Sparse Latent Factor Model to a group of data matrices within a directory. All matrices are supposed to have values
#' representing the gene expression observed for different genes (rows) and different samples (columns). 
#' @param path path to the directory where the target data matrices are located.
#' @param recursive logical argument (default = TRUE) indicating whether the function should look recursively 
#' inside folders.
#' @param a positive shape parameter of the Inverse Gamma prior distribution (default = 2.1).
#' @param b positive scale parameter of the Inverse Gamma prior distribution (default = 1.1).
#' @param gamma_a positive 1st shape parameter of the Beta prior distribution (default = 1).
#' @param gamma_b positive 2nd shape parameter of the Beta prior distribution (default = 1).
#' @param omega_0 prior variance of the spike mixture component (default = 0.01).
#' @param omega_1 prior variance of the slab mixture component (default = 10).
#' @param sample sample size to be considered for inference after the burn in period (default = 1000).
#' @param burnin size of the burn in period in the MCMC algorithm (default = sample/4).
#' @param lag lag to build the chains based on spaced draws from the Gibbs sampler (default = 1).
#' @param degenerate logical argument (default = FALSE) indicating whether to use the degenerate version of 
#' the mixture prior for the factor loadings.
#' @importFrom coda HPDinterval
#' @importFrom tools file_path_sans_ext
#' @importFrom utils read.table
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @seealso
#' \code{\link{slfm}}, \code{\link{process_matrix}}, \code{\link{plot_matrix}}
#' @export
slfm_list <- function(
  path = ".", recursive = TRUE,
  a = 2.1, b = 1.1, gamma_a = 1, gamma_b = 1,
  omega_0 = 0.01, omega_1 = 10, sample = 1000, burnin = round(0.25*sample), lag = 1, degenerate = FALSE) {

  files_list <- list.files(path, recursive = recursive, full.names = T)

  message("* |Press Ctrl + C to cancel...")
  pb <- txtProgressBar(0, length(files_list), style = 3, title = "")

  MATRIX_CLASSIFICATION <- c("Present", "Marginal", "Absent")
  names(MATRIX_CLASSIFICATION) <- c("S", "I", "N")

  results_list <- list()
  for(i in 1:length(files_list)) {
    file_name <- files_list[i]
    mat <- read.table(file_name)
    res <- slfm(mat, a, b, gamma_a, gamma_b, omega_0, omega_1, sample, burnin, lag, degenerate)
    clas_table <- table(res$classification)
    final_clas <- MATRIX_CLASSIFICATION[names(which.max(clas_table))]
    countS = clas_table["S"]
    if(is.na(countS)==TRUE){countS = 0}
    freq <- format(round(countS/sum(clas_table), 4), nsmall = 4)
    results_list[[i]] <- c(name = basename(tools::file_path_sans_ext(file_name)), clas = final_clas, frequency = freq)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  ret <- as.data.frame(do.call(rbind, results_list))
  names(ret) <- c("File", "Classification", "Significant %")
  ret
}