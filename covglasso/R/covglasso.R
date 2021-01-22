covglasso <- function(data = NULL,
                      S = NULL,
                      n = NULL,
                      lambda = NULL,
                      rho = NULL,
                      duplicated = TRUE,
                      L = 10,
                      crit = c("bic", "ebic"),
                      gamma = 0.8,
                      penalize.diag = FALSE,
                      start = NULL,
                      ctrl = control(),
                      path = FALSE)
{
  if ( all(is.null(data), is.null(S)) ) stop("We need some data! Please input 'data' or 'S' and 'n'")
  if ( is.null(S) & !is.null(data) ) {
    n <- nrow(data)
    S <- cov(data)*(n-1)/n
  } else if ( is.null(n) & is.null(data) ) stop("You need to provide the sample size 'n' in input if 'data' not supplied")
  V <- ncol(S)
  varnames <- colnames(S)
  crit <- match.arg( crit, choices = eval(formals(covglasso)$crit) )
  # if ( is.null(start) ) start <- S
  if ( is.null(start) ) start <- diag(diag(S))


  if ( is.null(lambda) ) {
    lambdaVec <- NULL
    R <- abs( cov2cor(S) )
    SS <- abs(S) + sqrt(.Machine$double.eps)
    if ( is.null(rho) ) {
      rho <- quantile( R[upper.tri(R)], seq(0,1, 1/L) )[-1]
      # rho <- seq(0,1, 1/L)
    } else {
      if ( max(rho) > 1 | min(rho) < 0 ) stop("Vector 'rho' must contain values between 0 and 1")
      L <- length(rho)
    }
    lambda <- lapply( 1:L, function(l) {
      temp <- ifelse(R < rho[l], (1/SS), 0)
      if ( !penalize.diag ) diag(temp) <- 0
      return(temp) } )
    lambda <- array( unlist(lambda), dim = c(V,V,L) )
    if ( duplicated ) {
      rho <- rho[ !duplicated(lambda, MARGIN = 3) ]
      lambda <- unique(lambda, MARGIN = 3)
    }
    LL <- dim(lambda)[3]

  } else {
    rho <- NULL
    if ( is.vector(lambda) ) {
      LL <- length(lambda)
      lambdaVec <- lambda
      lambda <- lapply(lambdaVec, function(l) {
        temp <- matrix(l, V, V)
        if ( !penalize.diag ) diag(temp) <- 0
        return(temp) } )
      lambda <- array( unlist(lambda), dim = c(V,V,LL) )
    } else {
      lambdaVec <- NULL
      LL <- dim(lambda)[3]
    }
  }

  fit <- covglassopath_bic(S, lambda, start = start, n = n, L = LL,
                           tolout = ctrl$tol.out, tolin = ctrl$tol.in,
                           iterout = ctrl$iter.out, iterin = ctrl$iter.in)

  bic <- 2*fit$loglik - fit$npar*log(n) - (crit == "ebic")*4*fit$npar*gamma*log(V)
  best <- which.max(bic)

  sel <- fit$out[[best]]
  dimnames(sel$sigma) <- dimnames(sel$omega) <- list(varnames, varnames)
  out <- list(sigma = sel$sigma, omega = sel$omega,
              loglik = sel$loglik, npar = fit$npar[best], penalty = sel$pen,
              bic = bic[best], BIC = bic, path = if (path) fit$out else NULL,
              rho = rho, lambda = lambdaVec)
  return(out)
}
