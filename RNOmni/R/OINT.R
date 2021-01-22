# Purpose: Rank Normal Omnibus test.
# Updated: 2020/10/03

# -----------------------------------------------------------------------------

#' Convert P-value to Cauchy Random
#' 
#' @param p Numeric p-value.
#' @return Numeric Cauchy random variable. 

PtoCauchy <- function(p) {
  if (p < 1e-10) {
    out <- 1 / (p * pi)
  } else {
    out <- tanpi(0.5 - p)
  }
  return(out)
}


#' Convert Cauchy Random Variable to P
#' 
#' @param z Numeric Cauchy random variable. 
#' @return Numeric p-value.
#' 
#' @importFrom stats pcauchy

CauchyToP <- function(z) {
  if (z > 1e10) {
    out <- (1 / z) / pi
  } else {
    out <- pcauchy(q = z, lower.tail = FALSE)
  }
  return(out)
}


# -----------------------------------------------------------------------------

#' Omnibus P-value.
#' 
#' Obtains an omnibus p-value from a vector of potentially dependent p-values 
#' using the method of Cauchy combination. The p-values are converted to Cauchy
#' random deviates then averaged. The distribution of the average of these 
#' deviates is well-approximated by a Cauchy distribution in the tails. See
#' <https://doi.org/10.1080/01621459.2018.1554485>.
#'
#' @param p Numeric vector of p-values.
#' @return OINT p-value.
#'
#' @export

OmniP <- function(p) {
  # Check input
  m <- min(p)
  M <- max(p)
  
  # Cases
  if ((m < 0) | (M > 1)) {
    stop("Cannot have p-values < 0 or > 1.")
  }
  if ((m == 0) & (M == 1)) {
    stop("Cannot have p-values of 0 and 1 simultaneously.")
  }
  if ((m == 0) & (M < 1)) {
    return(0)
  }
  if ((m > 0) & (M == 1)) {
    return(1)
  }
  
  # Convert to Cauchy.
  q <- sapply(p, PtoCauchy)
  
  # Mean of Cauchy random variables, which is again Cauchy.
  z <- mean(q)

  # Invert to obtain p-value.
  p <- CauchyToP(z)
  return(p)
}


# -----------------------------------------------------------------------------

#' Omnibus-INT
#'
#' Association test that synthesizes the \code{\link{DINT}} and
#' \code{\link{IINT}} tests. The first approach is most powerful for traits that
#' could have arisen from a rank-preserving transformation of a latent normal
#' trait. The second approach is most powerful for traits that are linear in
#' covariates, yet have skewed or kurtotic residual distributions. During the
#' omnibus test, the direct and indirect tests are separately applied, then the
#' p-values are combined via the Cauchy combination method.
#'
#' @param y Numeric phenotype vector.
#' @param G Genotype matrix with observations as rows, SNPs as columns.
#' @param X Model matrix of covariates and structure adjustments. Should include
#'   an intercept. Omit to perform marginal tests of association.
#' @param k Offset applied during rank-normalization. See
#'   \code{\link{RankNorm}}.
#' @param simple Return the OINT p-values only?
#' @return A numeric matrix of p-values, three for each column of \code{G}.
#'
#' @importFrom plyr aaply
#' @export
#' @seealso
#' \itemize{
#'   \item Basic association test \code{\link{BAT}}.
#'   \item Direct INT test \code{\link{DINT}}.
#'   \item Indirect INT test \code{\link{IINT}}.
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
#' # Omnibus
#' p <- OINT(y = y, G = G, X = X, simple = TRUE)

OINT <- function(y, G, X = NULL, k = 0.375, simple = FALSE) {
  
  # Generate X is omitted.
  if (is.null(X)) {
    X <- array(1, dim = c(length(y), 1))
  }
  
  # Input check.
  BasicInputChecks(y, G, X)

  # Association testing.
  # Calculate D-INT p-values.
  p_dint <- DINT(y = y, G = G, X = X, k = k, simple = TRUE)
  
  # Calculate I-INT p-values.
  p_iint <- IINT(y = y, G = G, X = X, k = k, simple = TRUE)

  # P Matrix
  p_mat <- cbind(p_dint, p_iint)

  # Omnibus p-values
  p_omni <- aaply(.data = p_mat, .margins = 1, .fun = OmniP)

  # Output

  # Check for genotype names.
  gnames <- colnames(G)
  if (is.null(gnames)) {
    gnames <- seq_len(ncol(G))
  }

  # Format
  if (simple) {
    out <- p_omni
    names(out) <- gnames
  } else {
    out <- cbind(
      "DINT-p" = p_dint,
      "IINT-p" = p_iint, 
      "OINT-p" = p_omni
    )
    rownames(out) <- gnames
  }
  return(out)
}
