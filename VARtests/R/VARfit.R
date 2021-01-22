VARfit <- function(y, p = 1, const = TRUE, trend = FALSE, exogen = NULL, univariate = FALSE){
  
  currentTime <- Sys.time()
  
  # checks that the arguments were correcly entered and possibly casts 'y' and 'exogen' to matrices:
  .checkArgs.VARfit()
  
  K <- ncol(y); N <- nrow(y); NlessP <- N - p
  NnonLagVar <- const + trend; if(!is.null(exogen)) NnonLagVar <- NnonLagVar + ncol(exogen)
  
  # creates a matrix of the lagged regressors: 
  Z <- embed(y, dimension = p + 1)[ , -(1:K)] # ([N - p] x Kp)
  
  # binds the const, trend, and dummies, if applicable, to make Z a [N - p] x [Kp + numberOf(const, trend, exogen)] matrix:
  if(!is.null(exogen)) Z <- cbind(exogen[-(1:p), ], Z)
  if(trend) Z <- cbind((p + 1):N, Z)
  if(const) Z <- cbind(rep(1, NlessP), Z) 
  
  # LS estimates:
  qrTemp <- qr(x = Z)
  coef <- as.matrix(qr.coef(qr = qrTemp, y = y[-(1:p), ])) # LS coefs (K x K)
  resid <- as.matrix(qr.resid(qr = qrTemp, y = y[-(1:p), ])) # residauls ([N - p] x K)
  
  # univariate LS estimates:
  uniCoef <- matrix(nrow = nrow(coef), ncol = ncol(coef))
  uniResid <- matrix(nrow = nrow(resid), ncol = ncol(resid))
  for(i in 1:K){
    Zuni <- NULL
    if(NnonLagVar > 0) Zuni <- Z[, 1:NnonLagVar]
    Zuni <- cbind(Zuni,
                  Z[, NnonLagVar + i + (0:(p - 1)) * K])
    qrTemp <- qr(x = Zuni)
    
    uniCoefTemp <- qr.coef(qr = qrTemp, y = y[-(1:p), i])
    uniResid[, i] <- qr.resid(qr = qrTemp, y = y[-(1:p), i]) 
    
    if(NnonLagVar > 0) uniCoef[1:NnonLagVar, i] <- uniCoefTemp[1:NnonLagVar]
    uniCoef[NnonLagVar + i + (0:(p - 1)) * K, i] <- uniCoefTemp[(NnonLagVar + 1):length(uniCoefTemp)]
  }
  
  # creates the names of the lagged variables. If the columns in 'y' had no names, the names are set to y1,...yK:
  if(is.null(colnames(y))) colnames(y) <- paste0("y", 1:K)
  if(!is.null(exogen)){
    if(is.null(colnames(exogen))){
      exogenNames <- paste0("exogen", 1:ncol(exogen))
      colnames(exogen) <- exogenNames
    } else exogenNames <- colnames(exogen)
  }
  varLagNames <- NULL
  for (j in 1:p) varLagNames <- c(varLagNames, paste0(colnames(y), "[-", j, "]"))
  rownames(coef) <- c(if(const) "const",
                      if(trend) "trend",
                      if(!is.null(exogen)) exogenNames,
                      varLagNames)
  colnames(coef) <- colnames(y)

  rownames(uniCoef) <- rownames(coef)
  colnames(uniCoef) <- colnames(coef) 
  
  returnValue <- list(y = y,
                      p = p,
                      N = N,
                      K = K,
                      const = const,
                      trend = trend,
                      exogen = exogen,
                      Z = Z,
                      call = match.call(),
                      coef = coef,
                      resid = resid,
                      uniCoef = uniCoef,
                      uniResid = uniResid,
                      univariate = univariate,
                      NnonLagVar = NnonLagVar,
                      description = paste("Estimated at", currentTime, "by user", Sys.info()[["user"]]),
                      time = difftime(Sys.time(), currentTime, units = "secs"))
  
  class(returnValue) <-  c("VARfit", class(returnValue))
  
  print(returnValue)
  invisible(returnValue)
}

.checkArgs.VARfit <- function(){
  # this function checks the arguments in the parent function 'VARfit'
  # and possibly casts 'y' and 'exogen' to matrices
  
  VARfitEnv <- parent.frame()
  
  # list of error messages:
  errs <- list()
  
  # checks 'y':
  if(is.null(VARfitEnv$y)) stop("'y' must be provided!")
  if(!is.matrix(VARfitEnv$y)) VARfitEnv$y <- as.matrix(VARfitEnv$y)
  if(!is.numeric(VARfitEnv$y)){
    errs <- append(errs, "'y' must be a numeric matix, vector, or data.frame")
    if(anyNA(VARfitEnv$y)) errs <- append(errs, "'y' cant contain any NA elements")
  }
  
  # checks 'p':
  if(!is.numeric(VARfitEnv$p)){
    errs <- append(errs, "'p' must be a positive integer")
  } else if(VARfitEnv$p %% 1 != 0 || VARfitEnv$p < 1){
    errs <- append(errs, "'p' must be a positive integer")
  } 
  
  # checks 'const', 'trend', and 'univariate':
  if(!is.logical(VARfitEnv$const)) errs <- append(errs, "'const' must be logical")
  if(!is.logical(VARfitEnv$trend)) errs <- append(errs, "'trend' must be logical")
  if(!is.logical(VARfitEnv$univariate)) errs <- append(errs, "'univariate' must be logical")
  
  # checks 'exogen':
  if(!is.null(VARfitEnv$exogen)){
    
    if(!is.matrix(VARfitEnv$exogen)) VARfitEnv$exogen <- as.matrix(VARfitEnv$exogen)
    if(!is.numeric(VARfitEnv$exogen)){
      errs <- append(errs, "'exogen' must be a matix, vector, or data.frame")
      if(anyNA(VARfitEnv$exogen)) errs <- append(errs, "'exogen' cant contain any NA elements")
    }
    if(nrow(VARfitEnv$exogen) != nrow(VARfitEnv$exogen)) errs <- append(errs, "the number of rows of 'exogen' and 'y' must be the same")
  } 
  
  K <- ncol(VARfitEnv$y); N <- nrow(VARfitEnv$y); p <- VARfitEnv$p
  NnonLagVar <- VARfitEnv$const + VARfitEnv$trend; if(!is.null(VARfitEnv$exogen)) NnonLagVar <- NnonLagVar + ncol(VARfitEnv$exogen)
  
  if(length(errs) > 0){
    
    msg <- "\nInput errors in 'VARfit()':"
    
    for(i in 1:length(errs)) msg <- paste0(msg, "\n", i, ". ", errs[[i]])
    
    stop(msg)
    
  } else{
    
    # checks the number of rows and columns of 'y':
    if(N < (K * p + NnonLagVar + 1)) stop("too few rows in 'y'. Perhaps 'y' should have been transposed?")
    
  }
}

print.VARfit <- function(x, ...){
  
  cat("Estimated VAR(", x$p, "):\n\n", sep = "")
  
  print(x$coef)
  
  if(x$univariate == "TRUE"){
    
    cat("\nUnivariate AR(", x$p, ") estimates:\n\n", sep = "")
    
    print(x$uniCoef, na.print = "")
  }
}

residuals.VARfit <- function(object, ...){
  
  return(object$resid)
}

coef.VARfit <- function(object, ...){
  
  return(object$coef)
}