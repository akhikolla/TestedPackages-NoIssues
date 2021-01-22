# Purpose: Basic score test
# Updated: 2020/10/03

#' Partition Data
#' 
#' Partition y and X according to the missingness pattern of g. 
#' 
#' @param e Numeric residual vector.
#' @param g Genotype vector.
#' @param X Model matrix of covariates.
#' @return List containing:
#' \itemize{
#'   \item "g_obs", observed genotype vector.
#'   \item "X_obs", covariates for subjects with observed genotypes.
#'   \item "X_mis", covariates for subjects with missing genotypes. 
#'   \item "e_obs", residuals for subjects with observed genotypes.
#' }

PartitionData <- function(e, g, X) {
  
  # Ensure matrix formatting. 
  g <- matrix(g, ncol = 1)
  e <- matrix(e, ncol = 1)
  
  # Adjust for missingness.
  is_obs <- !is.na(g)
  any_miss <- (sum(!is_obs) > 0)
  if (any_miss) {
    g_obs <- g[is_obs, , drop = FALSE]
    X_obs <- X[is_obs, , drop = FALSE]
    X_mis <- X[!is_obs, , drop = FALSE]
    e_obs <- e[is_obs, , drop = FALSE]
  } else {
    g_obs <- g
    X_obs <- X
    X_mis <- NULL
    e_obs <- e
  }
  
  # Output.
  out <- list(
    "g_obs" = g_obs,
    "X_obs" = X_obs,
    "X_mis" = X_mis,
    "e_obs" = e_obs
  )
  return(out)
}


# -----------------------------------------------------------------------------

#' Score Statistics
#' 
#' @param e Numeric residual vector.
#' @param g Genotype vector.
#' @param X Model matrix of covariates.
#' @param v Residual variance.
#' @return Numeric vector containing the "score" statistic, standard error "se",
#'   "z", and "p" value.
#'   
#' @importFrom stats pchisq

ScoreStat <- function(e, g, X, v) {
  
  # Split data.
  split_data <- PartitionData(e = e, g = g, X = X)
  
  # Information components.
  info_gg <- matIP(split_data$g_obs, split_data$g_obs)
  info_gx <- matIP(split_data$g_obs, split_data$X_obs)
  info_xx <- matIP(split_data$X_obs, split_data$X_obs)
  
  # Efficient info.
  eff_info <- as.numeric(SchurC(info_gg, info_xx, info_gx))
  
  # Score.
  score <- as.numeric(matIP(split_data$g_obs, split_data$e_obs)) / v
  
  # SE.
  se <- sqrt(eff_info / v)
  
  # Z statistic.
  z_stat <- score / se
  
  # Chi statistic.
  chi_stat <- z_stat^2
  
  # p-value.
  p <- pchisq(q = chi_stat, df = 1, lower.tail = FALSE)
  
  # Output.
  out <- c(
    "score" = score, 
    "se" = se, 
    "z" = z_stat, 
    "p" = p
  )
  return(out)
}


#' Basic Association Score Test
#' 
#' @param y Numeric phenotype vector.
#' @param G Genotype matrix with observations as rows, SNPs as columns.
#' @param X Model matrix of covariates.
#' @return Numeric matrix, with 1 row per SNP, containing these columns:
#' \itemize{
#'   \item "score", the score statistic.
#'   \item "se", its standard error.
#'   \item "z", the Z statistic.
#'   \item "p", the p-value. 
#' }
#' 
#' @importFrom plyr aaply

BAT.ScoreTest <- function(y, G, X) {
  
  # Fit null model.
  fit0 <- fitOLS(y = y, X = X)
  
  # Extract model components.
  e <- matrix(fit0$Resid, ncol = 1)
  v <- fit0$V
  
  # Calculate Score Statistic.
  out <- aaply(.data = G, .margins = 2, .fun = function(g) {
    ScoreStat(e = e, g = g, X = X, v = v)
  })
  
  return(out)
}


# -----------------------------------------------------------------------------

#' Basic Association Score Test
#' 
#' @param y Numeric phenotype vector.
#' @param g Genotype vector.
#' @param X Model matrix of covariates.
#' @return Numeric matrix, with 1 row per SNP, containing these columns:
#' \itemize{
#'   \item "score", the score statistic.
#'   \item "se", its standard error.
#'   \item "z", the Z statistic.
#'   \item "p", the p-value. 
#' }
#' 
#' @importFrom stats pchisq

WaldStat <- function(y, g, X) {
  
  # Split data. 
  split_data <- PartitionData(e = y, g = g, X = X)
  
  # Fit cull model.
  fit1 <- fitOLS(y = split_data$e_obs, X = cbind(split_data$g_obs, split_data$X_obs))
  
  # Coefficient.
  bg <- fit1$Beta[1]
  
  # Variance.
  eff_info_inv <- as.numeric(matInv(fit1$Ibb)[1, 1])
    
  # Standard error.
  se <- sqrt(eff_info_inv)
  
  # Z statistic.
  z_stat <- bg / se
  
  # Chi statistic.
  chi_stat <- z_stat^2
  
  # p-value.
  p <- pchisq(q = chi_stat, df = 1, lower.tail = FALSE)
  
  # Output.
  out <- c(
    "wald" = bg, 
    "se" = se, 
    "z" = z_stat, 
    "p" = p
  )
  return(out)
}


#' Basic Association Wald Test
#' 
#' @param y Numeric phenotype vector.
#' @param G Genotype matrix with observations as rows, SNPs as columns.
#' @param X Model matrix of covariates.
#' @return Numeric matrix, with 1 row per SNP, containing these columns:
#' \itemize{
#'   \item "score", the score statistic.
#'   \item "se", its standard error.
#'   \item "z", the Z statistic.
#'   \item "p", the p-value. 
#' }
#' 
#' @importFrom plyr aaply

BAT.WaldTest <- function(y, G, X) {
  out <- aaply(.data = G, .margins = 2, .fun = function(g) {
    WaldStat(y = y, g = g, X = X)
  })
  return(out)
}


# -----------------------------------------------------------------------------

#' Basic Input Checks
#' 
#' @param y Numeric phenotype vector.
#' @param G Genotype matrix with observations as rows, SNPs as columns.
#' @param X Covariate matrix.

BasicInputChecks <- function(y, G, X) {
  
  # Ensure y is a numeric vector.
  if (!is.vector(y)) {
    stop("A numeric vector is expected for y.")
  }
  
  # Ensure G is a numeric matrix.
  if (!is.matrix(G)) {
    stop("A numeric matrix is expected for G.")
  }
  
  # Ensure X is a numeric matrix.
  if (!is.matrix(X)) {
    stop("A numeric matrix is expected for X.")
  }
  
  # Ensure y and X are complete.
  y_or_x_miss <- sum(is.na(y)) + sum(is.na(X))
  if (y_or_x_miss > 0) {
    stop("Please exclude observations missing phenotype or covariate information.")
  }
}


# -----------------------------------------------------------------------------

#' Basic Association Test
#'
#' Conducts tests of association between the loci in \code{G} and the
#' untransformed phenotype \code{y}, adjusting for the model matrix \code{X}.
#'
#' @param y Numeric phenotype vector.
#' @param G Genotype matrix with observations as rows, SNPs as columns.
#' @param X Model matrix of covariates and structure adjustments. Should include
#'   an intercept. Omit to perform marginal tests of association.
#' @param test Either Score or Wald.
#' @param simple Return the p-values only?
#' @return If \code{simple = TRUE}, returns a vector of p-values, one for each column
#'   of \code{G}. If \code{simple = FALSE}, returns a numeric matrix, including the
#'   Wald or Score statistic, its standard error, the Z-score, and the p-value.
#'   
#' @export
#' @seealso
#' \itemize{
#'   \item Direct INT \code{\link{DINT}}
#'   \item Indirect INT \code{\link{IINT}}
#'   \item Omnibus INT \code{\link{OINT}}
#' }
#'
#' @examples
#' set.seed(100)
#' # Design matrix
#' X <- cbind(1, rnorm(1e3))
#' # Genotypes
#' G <- replicate(1e3, rbinom(n = 1e3, size = 2, prob = 0.25))
#' storage.mode(G) <- "numeric"
#' # Phenotype
#' y <- as.numeric(X %*% c(1, 1)) + rnorm(1e3)
#' # Association test
#' p <- BAT(y = y, G = G, X = X)

BAT <- function(y, G, X = NULL, test = "Score", simple = FALSE) {
  
  # Generate X is omitted.
  if (is.null(X)) {
    X <- array(1, dim = c(length(y), 1))
  }
  
  # Input check.
  BasicInputChecks(y, G, X)

  # Association testing.
  if (test == "Score") {
    out <- BAT.ScoreTest(y = y, G = G, X = X)
  } else if (test == "Wald") {
    out <- BAT.WaldTest(y = y, G = G, X = X)
  } else {
    stop("Select test from among: Score, Wald.")
  }
  if (!is.matrix(out)) {out <- matrix(out, nrow = 1)}
  
  # Check for genotype names.
  gnames <- colnames(G)
  if (is.null(gnames)) {
    gnames <- seq_len(ncol(G))
  }
  
  # Format output.
  if (simple) { 
    out <- out[, 4]
    names(out) <- gnames
  } else {
    colnames(out) <- c(test, "SE", "Z", "P")
    rownames(out) <- gnames
  }
  return(out)
}
