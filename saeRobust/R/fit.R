#' Fitting Precedures
#'
#' Several fitting procedures. The arguments can be passed to these functions
#' using the interface in \link{rfh}. The functions here listed are the low
#' level implementations and are not intended for interactive use.
#'
#' @inheritParams rfh
#' @param y (numeric) response vector
#' @param x ([m|M]atrix) the design matrix
#' @param samplingVar (numeric) vector with sampling variances
#' @param ... arguments passed to \code{fitGenericModel}
#' @param matVFun (function) a function with one argument - the variance
#'   parameters - constructing something like \link{variance}
#' @param fixedPointParam (function) a function with one argument. The vector of
#'   model parameters. Returns a list of results of the next iteration in the
#'   overall algorithm.
#' @param k (numeric) tuning constant
#' @param K (numeric) scaling constant
#' @param x0Coef (numeric) starting values for regression coefficients
#' @param x0Var (numeric) starting values for variance parameters
#' @param x0Re (numeric) starting values for random effects
#' @param tol (numeric) numerical toloerance to be used during optimisation
#' @param maxIter (integer) the maximum number of iterations for model parameters.
#' @param maxIterParam (integer) the maximum number of iterations for each
#'   parameter in each overall iteration
#' @param maxIterRe (integer) the maximum number of iterations for fitting the
#'   random effects
#' @param convCrit (function) a function defining the stopping rule
#' @param W (matrix) proximity matrix
#' @param nTime (integer) number of time periods
#'
#' @details
#' \code{fitrfh} implements the robust Fay-Herriot model; \code{fitrsfh} the
#' spatial, \code{fitrtfh} the temporal, and \code{fitrstfh} the spatio-temporal
#' extension to this model type. See \link{rfh} how to fit such models.
#' \code{fitGenericModel} is used by all these implementations and can be used
#' for possible extensions of the framework.
#'
#' @rdname fit
#' @export
#'
#' @examples
#' data(milk, package = "sae")
#' x <- matrix(1, nrow = NROW(milk))
#' y <- milk$yi
#' samplingVar <- milk$SD^2
#' fitrfh(y, x, samplingVar)
fitrfh <- function(y, x, samplingVar, ...) {
  # Non interactive fitting function for robust FH

  fixedPointParam <- function(parent = parent.frame()) {
    # To connect the returned function to the calling environment:
    enclosingEnv <- environment()
    parent.env(enclosingEnv) <- parent
    function(param) {
      # Defines one Iteration in algorithm:
      param <- as.numeric(param)
      beta <- param[-length(param)]
      sigma2 <- param[length(param)]

      fpBeta <- addHistory(fixedPointRobustBeta(
        y, x, matVFun(sigma2), psi = psi
      ))

      beta <- fixedPoint(fpBeta, beta, addMaxIter(convCrit, maxIterParam))

      fpSigma2 <- fixedPointRobustDelta(
        y, x, beta, matVFun, psi, K
      ) %>% addConstraintMin(0.00001) %>% addHistory

      sigma2 <- fixedPoint(fpSigma2, sigma2, addMaxIter(convCrit, maxIterParam))

      list(beta, sigma2)
    }
  }

  matVFun <- . %>% matVFH(.samplingVar = samplingVar)
  out <- fitGenericModel(y, x, matVFun, fixedPointParam, ...)
  names(out$iterations) <- c("coefficients", "variance", "re")
  names(out$variance) <- "variance"

  stripSelf(retList(public = c("samplingVar"), super = out))

}

#' @rdname fit
#' @export
fitrsfh <- function(y, x, samplingVar, W, x0Var = c(0.01, 1), ...) {
  # Non interactive fitting function for robust spatial FH

  fixedPointParam <- function(parent = parent.frame()) {
    # To connect the returned function to the calling environment:
    enclosingEnv <- environment()
    parent.env(enclosingEnv) <- parent
    function(param) {
      # Defines one Iteration in algorithm:
      param <- as.numeric(param)
      beta <- param[-c(length(param)-1, length(param))]
      sigma2 <- param[length(param)]
      rho <- param[length(param) - 1]

      fpBeta <- addHistory(fixedPointRobustBeta(
        y, x, matVFun(c(rho, sigma2)), psi = psi
      ))

      beta <- fixedPoint(fpBeta, beta, addMaxIter(convCrit, maxIterParam))

      fpSigma2 <- fixedPointRobustDelta(
        y, x, beta, function(x) matVFun(c(rho, x)), psi, K, "sigma2"
      ) %>%
        addConstraintMin(0.00001) %>%
        addHistory

      sigma2 <- fixedPoint(fpSigma2, sigma2, addMaxIter(convCrit, maxIterParam))

      fpRho <- fixedPointNumericDelta(
        y, x, beta, function(x) matVFun(c(x, sigma2)), psi, K, "rho", 0.01, -0.99999
      ) %>%
        addConstraintMin(-0.99999) %>% addConstraintMax(0.99999) %>%
        addHistory

      rho <- fixedPoint(fpRho, rho, addMaxIter(convCrit, maxIterParam))

      list(beta, rho, sigma2)
    }
  }

  matVFun <- function(x, ...) {
    matVSFH(.rho = x[1], .sigma2 = x[2], .W = W, .samplingVar = samplingVar, ...)
  }
  out <- fitGenericModel(y, x, matVFun, fixedPointParam, x0Var = x0Var, ...)
  names(out$iterations) <- c("coefficients", "correlation", "variance", "re")
  names(out$variance) <- c("correlation", "variance")

  stripSelf(retList(
    c("fitrsfh"),
    public = c("samplingVar", "W"), super = out)
  )

}

#' @rdname fit
#' @export
fitrtfh <- function(y, x, samplingVar, nTime, x0Var = c(0.01, 1, 1), ...) {
  # Non interactive fitting function for robust temporal FH

  fixedPointParam <- function(parent = parent.frame()) {
    # To connect the returned function to the calling environment:
    enclosingEnv <- environment()
    parent.env(enclosingEnv) <- parent
    function(param) {
      # Defines one Iteration in algorithm:
      param <- as.numeric(param)
      beta <- param[-c(length(param) - (0:2))]
      sigmas <- param[length(param) - (1:0)]
      rho <- param[length(param) - 2]

      matV <- matVFun(c(rho, sigmas))

      # Regression Coefficients:
      fpBeta <- addHistory(fixedPointRobustBeta(
        y, x, matV, psi = psi
      ))

      beta <- fixedPoint(fpBeta, beta, addMaxIter(convCrit, maxIterParam))

      # Variance Parameters:
      .fpSigma2 <- function(sigmas) {
        # For performance this one is implemented in c++:
        fixedPointSigma(
          y, as.matrix(x), beta, sigmas, rho, as.matrix(matV$Z1()),
          as.matrix(matV$Omega1()), as.matrix(matV$Z()),
          get(".samplingVar", envir = attr(matV, ".self")), k, K
        )
      }

      fpSigma2 <- .fpSigma2 %>% addConstraintMin(0.00001) %>% addHistory

      sigmas <- fixedPoint(
        fpSigma2,
        sigmas,
        addMaxIter(convCrit, maxIterParam)
      )

      fpRho <- fixedPointNumericDelta(
        y, x, beta, function(x) matVFun(c(x, sigmas)), psi, K, "rho", 0.01, -0.99999
      ) %>%
        addConstraintMin(-0.99999) %>% addConstraintMax(0.99999) %>%
        addHistory

      rho <- fixedPoint(fpRho, rho, addMaxIter(convCrit, maxIterParam))

      list(beta, rho, sigmas)
    }
  }

  matVFun <- function(x, ...) {
    matVTFH(.rho = x[1], .sigma2 = x[2:3], .nTime = nTime,
            .samplingVar = samplingVar, ...)
  }

  out <- fitGenericModel(y, x, matVFun, fixedPointParam, x0Var = x0Var, ...)
  names(out$iterations) <- c("coefficients", "correlation", "variance", "re")
  names(out$variance) <- c("correlation", "variance1", "varianceAR")

  stripSelf(retList(
    c("fitrtfh"),
    public = c("samplingVar", "nTime"), super = out)
  )

}

#' @rdname fit
#' @export
fitrstfh <- function(y, x, samplingVar, W, nTime, x0Var = c(0.01, 0.01, 1, 1), ...) {
  # Non interactive fitting function for robust temporal FH

  fixedPointParam <- function(parent = parent.frame()) {
    # To connect the returned function to the calling environment:
    enclosingEnv <- environment()
    parent.env(enclosingEnv) <- parent
    function(param) {
      # Defines one Iteration in algorithm:
      param <- as.numeric(param)
      beta <- param[-c(length(param) - (0:3))]
      sigmas <- param[length(param) - (1:0)]
      rho1 <- param[length(param) - (3)]
      rho2 <- param[length(param) - (2)]

      matV <- matVFun(c(rho1, rho2, sigmas))

      # Regression Coefficients:
      fpBeta <- addHistory(fixedPointRobustBeta(
        y, x, matV, psi = psi
      ))

      beta <- fixedPoint(fpBeta, beta, addMaxIter(convCrit, maxIterParam))

      # Variance Parameters:
      .fpSigma2 <- function(sigmas) {
        # For performance this one is implemented in c++:
        fixedPointSigma(
          y, as.matrix(x), beta, sigmas, rho2, as.matrix(matV$Z1()),
          as.matrix(matV$Omega1()), as.matrix(matV$Z()),
          get(".samplingVar", envir = attr(matV, ".self")), k, K
        )
      }

      fpSigma2 <- .fpSigma2 %>% addConstraintMin(0.00001) %>% addHistory

      sigmas <- fixedPoint(
        fpSigma2,
        sigmas,
        addMaxIter(convCrit, maxIterParam)
      )

      # AR:
      fpRho <- fixedPointNumericDelta(
        y, x, beta, function(x) matVFun(c(rho1, x, sigmas)), psi, K, "rho2", 0.01, -0.99999
      ) %>%
        addConstraintMin(-0.99999) %>% addConstraintMax(0.99999) %>%
        addHistory

      rho2 <- fixedPoint(fpRho, rho2, addMaxIter(convCrit, maxIterParam))

      # SAR
      fpRho <- fixedPointNumericDelta(
        y, x, beta, function(x) matVFun(c(x, rho2, sigmas)), psi, K, "rho1", 0.01, -0.99999
      ) %>%
        addConstraintMin(-0.99999) %>% addConstraintMax(0.99999) %>%
        addHistory

      rho1 <- fixedPoint(fpRho, rho1, addMaxIter(convCrit, maxIterParam))

      list(beta, rho1, rho2, sigmas)
    }
  }

  matVFun <- function(x, ...) {
    matVSTFH(.rho = x[1:2], .sigma2 = x[3:4], .nTime = nTime, .W = W,
            .samplingVar = samplingVar, ...)
  }

  out <- fitGenericModel(y, x, matVFun, fixedPointParam, x0Var = x0Var, ...)
  names(out$iterations) <- c("coefficients", "SAR", "AR", "variance", "re")
  names(out$variance) <- c("SAR", "AR", "varianceSAR", "varianceAR")

  stripSelf(retList(
    c("fitrstfh"),
    public = c("samplingVar", "W", "nTime"), super = out)
  )

}

#' @export
#' @rdname fit
fitGenericModel <- function(
  y, x, matVFun, fixedPointParam,
  k = 1.345, K = getK(k),
  x0Coef = NULL, x0Var = 1, x0Re = NULL,
  tol = 1e-6, maxIter = 100, maxIterParam = 10, maxIterRe = 100, convCrit = convCritRelative(tol), ...) {
  # Non interactive fitting function for robust FH Models
  # ... are ignored

  psi <- . %>% psiOne(k)

  # Fitting Model Parameter:
  if (is.null(x0Coef)) {
    x0Coef <- as.numeric(fitCoefStartingValue(y, x, matVFun(x0Var)))
  }
  out <- fixedPoint(
    addStorage(fixedPointParam()),
    c(x0Coef, x0Var),
    addMaxIter(convCrit, maxIter)
  )
  coefficients <- out[1:length(x0Coef)]
  names(coefficients) <- colnames(x)
  variance <- out[(length(x0Coef) + 1):length(out)]

  # Fitting Random Effects
  matV <- matVFun(variance)
  if (is.null(x0Re)) {
    x0Re <- fitReStartingValues(y, x, coefficients, matV, psi)
  }
  re <- fitRe(
    y, x, coefficients, matV,
    x0Re, psi, addMaxIter(convCrit, maxIterRe)
  )

  # Iterations
  iterations <- storage$reformat(attr(out, "storage"))
  iterations <- c(iterations, re["iterations"])

  # Reformats
  re <- re$re
  fitted <- as.numeric(x %*% coefficients)
  reblup <- as.numeric(fitted + matV$Z() %*% re)
  residuals <- y - reblup

  stripSelf(retList("fitrfh", public = c(
    "coefficients", "variance", "iterations",
    "tol", "maxIter", "maxIterParam", "maxIterRe",
    "k", "K",
    "y", "x", "re", "reblup", "residuals", "fitted"
  )))

}

fitCoefStartingValue <- function(y, x, matV) {
  xPrimeV <- crossprod(x, matV$VInv())
  solve(xPrimeV %*% x) %*% xPrimeV %*% y
}

fitReStartingValues <- function(y, x, beta, matV, psi) {
  # y: (numeric) response
  # x: (Matrix) Design-Matrix
  # beta: (numeric) estimated fixed effects
  # matV: (list) list of functions see e.g. matVFH
  U <- matU(matV$V())
  resids <- psi(U$sqrtInv() %*% (y - x %*% beta))
  as.numeric(
    tcrossprod(matV$Vu(), matV$Z()) %*% matV$VInv() %*% U$sqrt() %*% psi(resids)
  )
}

fitRe <- function(y, x, beta, matV, x0, psi, convCrit) {
  # y: (numeric) response
  # x: (Matrix) Design-Matrix
  # beta: (numeric) estimated fixed effects
  # matV: (list) list of functions see e.g. matVFH
  # psi: (function) influence function
  # convCrit: (function) convergence criterion

  fpFun <- fixedPointRobustRandomEffect(
    y, x, beta, matV, psi
  )

  re <- fixedPoint(addHistory(fpFun), x0, convCrit)
  iterRe <- attr(re, "history")

  list(re = as.numeric(re), iterations = cbind(iterRe, i = 1:NROW(iterRe)))

}
