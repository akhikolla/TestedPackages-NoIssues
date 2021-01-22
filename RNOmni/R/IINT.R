# Purpose: Indirect INT-based method
# Updated: 2020/10/04

#' Basic Association Score Test
#' 
#' @param y Numeric phenotype vector.
#' @param G Genotype matrix with observations as rows, SNPs as columns.
#' @param X Model matrix of covariates.
#' @param k Offset applied during rank-normalization.
#' @return Numeric matrix, with 1 row per SNP, containing these columns:
#' \itemize{
#'   \item "score", the score statistic.
#'   \item "se", its standard error.
#'   \item "z", the Z statistic.
#'   \item "p", the p-value. 
#' }
#' @importFrom plyr aaply

IINT.ScoreTest <- function(y, G, X, k) {
  
  # Fit null model.
  fit0 <- fitOLS(y = y, X = X)
  
  # Extract model components.
  e <- matrix(
    RankNorm(u = as.numeric(fit0$Resid), k = k), 
    ncol = 1
  )
  v <- fit0$V
  
  # Calculate Score Statistic.
  out <- aaply(.data = G, .margins = 2, .fun = function(g) {
    ScoreStat(e = e, g = g, X = X, v = 1)
  })
  
  return(out)
}


# -----------------------------------------------------------------------------

#' Indirect-INT
#' 
#' Two-stage association testing procedure. In the first stage, phenotype 
#' \code{y} and genotype \code{G} are each regressed on the model matrix
#' \code{X} to obtain residuals. The phenotypic residuals are transformed
#' using \code{\link{RankNorm}}. In the next stage, the INT-transformed
#' residuals are regressed on the genotypic residuals. 
#' 
#' @param y Numeric phenotype vector.
#' @param G Genotype matrix with observations as rows, SNPs as columns.
#' @param X Model matrix of covariates and structure adjustments. Should include
#'   an intercept. Omit to perform marginal tests of association. 
#' @param k Offset applied during rank-normalization. See
#'   \code{\link{RankNorm}}.
#' @param simple Return the p-values only? 
#' @return If \code{simple = TRUE}, returns a vector of p-values, one for each column
#'   of \code{G}. If \code{simple = FALSE}, returns a numeric matrix, including the
#'   Wald or Score statistic, its standard error, the Z-score, and the p-value.
#' 
#' @export
#' @seealso
#' \itemize{
#'   \item Basic association test \code{\link{BAT}}.
#'   \item Direct INT test \code{\link{DINT}}.
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
#' y <- exp(as.numeric(X %*% c(1,1)) + rnorm(1e3))
#' # Association test
#' p <- IINT(y = y, G = G, X = X)

IINT <- function(y, G, X = NULL, k = 0.375, simple = FALSE) {

  # Generate X is omitted.
  if (is.null(X)) {
    X <- array(1, dim = c(length(y), 1))
  }
  
  # Input check.
  BasicInputChecks(y, G, X)

  # Score statistics
  out <- IINT.ScoreTest(y = y, G = G, X = X, k = k)
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
    if (!is.matrix(out)) {out <- matrix(out, nrow = 1)}
    colnames(out) <- c("Score", "SE", "Z", "P")
    rownames(out) <- gnames
  }
  return(out)
}
