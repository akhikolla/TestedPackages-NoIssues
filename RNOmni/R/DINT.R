# Purpose: Direct INT-based test
# Updated: 2020/10/03

#' Direct-INT
#'
#' Applies the rank-based inverse normal transformation (\code{\link{RankNorm}})
#' to the phenotype \code{y}. Conducts tests of association between the loci in
#' \code{G} and transformed phenotype, adjusting for the model matrix \code{X}.
#'
#' @param y Numeric phenotype vector.
#' @param G Genotype matrix with observations as rows, SNPs as columns.
#' @param X Model matrix of covariates and structure adjustments. Should include
#'   an intercept. Omit to perform marginal tests of association.
#' @param k Offset applied during rank-normalization. See
#'   \code{\link{RankNorm}}.
#' @param test Either Score or Wald.
#' @param simple Return the p-values only?
#' @return If \code{simple = TRUE}, returns a vector of p-values, one for each column
#'   of \code{G}. If \code{simple = FALSE}, returns a numeric matrix, including the
#'   Wald or Score statistic, its standard error, the Z-score, and the p-value.
#'
#' @export
#' @seealso
#' \itemize{
#'   \item Basic association test \code{\link{BAT}}.
#'   \item Indirect INT test \code{\link{IINT}}.
#'   \item Omnibus INT test \code{\link{OINT}}.
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
#' y <- exp(as.numeric(X %*% c(1, 1)) + rnorm(1e3))
#' # Association test
#' p <- DINT(y = y, G = G, X = X)

DINT <- function(y, G, X = NULL, k = 0.375, test = "Score", simple = FALSE) {
  
  # Generate X is omitted.
  if (is.null(X)) {
    X <- array(1, dim = c(length(y), 1))
  }
  
  # Input check.
  BasicInputChecks(y, G, X)

  # Transform phenotype.
  z <- RankNorm(u = y, k = k)
  
  # Apply basic association test to transformed phenotype.
  out <- BAT(y = z, G = G, X = X, test = test, simple = simple)
  return(out)
}
