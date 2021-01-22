#' Matrix constructor functions
#'
#' These functions construct different parts of matrix components. They are used
#' internally. If you are interested in the weights of a model fitted using
#' \link{rfh} please try to use \link[saeRobust]{weights.fitrfh} on that object.
#'
#' @param .V (Matrix) variance matrix
#' @param .nDomains (integer) number of domains
#' @param .nTime (integer) number of time periods
#'
#' @details
#'
#' \code{matU} computes U. U is the matrix containing only the diagonal
#'   elements of V. This function returns a list of functions which can be
#'   called to compute specific transformations of U.
#'
#' @references
#'
#' Warnholz, S. (2016): "Small Area Estimaiton Using Robust Extension to Area
#' Level Models". Not published (yet).
#'
#' @rdname varianceMatrices
#' @export
#'
#' @examples
#' data("grapes", package = "sae")
#' data("grapesprox", package = "sae")
#'
#' fitRFH <- rfh(
#'   grapehect ~ area + workdays - 1,
#'   data = grapes,
#'   samplingVar = "var"
#' )
#'
#' matV <- variance(fitRFH)
#'
#' # matU:
#' matU(matV$V())$U()
#' matU(matV$V())$sqrt()
#' matU(matV$V())$sqrtInv()
#'
#' # matB (and matA + matW accordingly):
#' matB(
#'   fitRFH$y,
#'   fitRFH$x,
#'   fitRFH$coefficients,
#'   fitRFH$re,
#'   matV,
#'   function(x) psiOne(x, k = fitRFH$k)
#' )
#'
#' matBConst(
#'   fitRFH$y,
#'   fitRFH$x,
#'   fitRFH$coefficients,
#'   matV,
#'   function(x) psiOne(x, k = fitRFH$k)
#' )(fitRFH$re)
#'
#' # construcors for 'Z' in linear mixed models
#' matTZ(2, 3)
#' matTZ1(2, 3)
matU <- function(.V) {

    .diag <- function(x) Diagonal(x = x)

    U <- getter(diag(.V), .diag)
    sqrt <- getter(base::sqrt(diag(.V)), .diag)
    sqrtInv <- getter(1 / base::sqrt(diag(.V)), .diag)
    retList()

}

#' @param x ([m|M]atrix) a matrix
#'
#' @details \code{matTrace} computes the trace of a matrix.
#'
#' @rdname varianceMatrices
#' @export
matTrace <- function(x) {
    sum(diag(x))
}

#' @param y (numeric) response
#' @param X (Matrix) design matrix
#' @param beta (numeric) vector of regression coefficients
#' @param re (numeric) vector of random effects
#' @param reblup (numeric) vector with robust best linear unbiased predictions
#' @param matV (list of functions) see \code{matVFH} for an example
#' @param psi (function) the influence function
#' @param W (Matrix) the weighting matrix
#' @param samplingVar (numeric) the vector of sampling variances
#' @param c (numeric) scalar
#'
#' @details \code{matB} computes the matrix B which is used to compute the
#'   weights in the pseudo linearised representation of the REBLUP.
#'
#' @rdname varianceMatrices
#' @export
matB <- function(y, X, beta, re, matV, psi) {
    matBConst(y, X, beta, matV, psi)(re)
}

#' @details \code{matBConst} returns a function with one argument, u, to compute
#'   the matrix B. This function is used internally to compute B in the fixed
#'   point algorithm.
#'
#' @rdname varianceMatrices
#' @export
matBConst <- function(y, X, beta, matV, psi) {

  Ue <- matU(matV$Ve())
  Uu <- matU(matV$Vu())

  # Helper functions
  resid <- function(u) as.numeric(memResid - as.matrix(matV$Z()) %*% u)

  w2 <- function(u) {
    resids <- resid(u) * diag(Ue$sqrtInv())
    psi(resids) / resids
  }

  w3 <- function(u) {
    resids <- u * diag(Uu$sqrtInv())
    psi(resids) / resids
  }

  # Precalculations - they only have to be done once
  memXB <- X %*% beta
  memResid <- y - memXB
  memZVeU <- crossprod(matV$Z(), matV$VeInv())

  function(re) {
    W2 <- Diagonal(x = w2(re))
    W3 <- Diagonal(x = w3(re))
    memPart1 <- memZVeU %*% W2
    memPart2 <- matV$VuInv() %*% W3
    solve(memPart1 %*% matV$Z() + memPart2) %*% memPart1
  }
}

#' @details \code{matA} computes the matrix A which is used to compute the
#'   weights in the pseudo linearized representation of the REBLUP.
#'
#' @rdname varianceMatrices
#' @export
matA <- function(y, X, beta, matV, psi) {
    matAConst(y, X, matV, psi)(beta)
}

#' @details \code{matAConst} returns a function with one argument, beta, to
#'   compute the matrix A. This function is used internally to compute A in the
#'   fixed point algorithm for beta.
#'
#' @rdname varianceMatrices
#' @export
matAConst <- function(y, X, matV, psi) {
    force(y); force(psi)

    # Helper functions
    resid <- function(beta) as.numeric(U$sqrtInv() %*% (y - X %*% beta))

    w <- function(beta) {
        resids <- resid(beta)
        psi(resids) / resids
    }

    # Precalculations - they only have to be done once
    U <- matU(matV$V())
    memP0 <- crossprod(X, matV$VInv())

    function(beta) {
        W1 <- Diagonal(x = w(beta))
        solve(memP0 %*% W1 %*% X) %*% memP0 %*% W1
    }
}

#' @details \code{matW} returns a matrix containing the weights as they are
#'   defined for the pseudo linear form, such that \code{matW \%*\% y} is the
#'   REBLUP.
#'
#' @rdname varianceMatrices
#' @export
matW <- function(y, X, beta, re, matV, psi) {
  A <- matA(y, X, beta, matV, psi)
  B <- matB(y, X, beta, re, matV, psi)
  XA <- X %*% A
  XA + matV$Z() %*% B %*% (Diagonal(length(y)) - XA)
}

#' @details \code{matWbc} returns a matrix containing the weights as they are
#'   defined for the pseudo linear form, such that \code{matWbc \%*\% y} is the
#'   bias-corrected REBLUP. \code{c} is a multiplyer for the standard deviation.
#'
#' @rdname varianceMatrices
#' @export
matWbc <- function(y, reblup, W, samplingVar, c = 1) {
  samplingVar <- sqrt(samplingVar) # we continue with sd not var
  Wbc <- W
  w1 <- 1 - c * samplingVar / y
  w2 <- 1 + c * samplingVar / y
  cond1 <- reblup < (y - c * samplingVar)
  cond2 <- reblup > (y + c * samplingVar)

  # How can we vectorize this: So long we need operate on each line of W.
  for (i in seq_along(samplingVar)) {
    I <- ifelse(seq_along(samplingVar) == i, 1, 0) # 1(j = i)
    w1_ <- w1[i] * I
    w2_ <- w2[i] * I
    w <- W[i, ]
    (if (cond1[i]) w1_
    else if (cond2[i]) w2_
    else w) -> Wbc[i, ]
  }
  Wbc
}

#' @details \code{matTZ} constructs the Z matrix in a linear mixed model with
#'   autocorrelated random effects.
#'
#' @export
#' @rdname varianceMatrices
matTZ <- function(.nDomains, .nTime) {
  #Construct Z-Matrix for temporal mixed model
  Z1 <- matTZ1(.nDomains, .nTime)
  Z2 <- Diagonal(.nDomains * .nTime)
  cbind(Z1, Z2)
}

#' @details \code{matTZ1} constructs the Z1 matrix in a linear mixed model with
#'   autocorrelated random effects.
#'
#' @export
#' @rdname varianceMatrices
matTZ1 <- function(.nDomains = 10, .nTime = 10) {
  I1 <- Diagonal(.nDomains)
  I2 <- rep(1, .nTime)
  I1 %x% I2
}
