ACtest <- function(fit, h = 4, HCtype = c("LM", "HC0", "HC1", "HC2", "HC3"), univariate = FALSE){
  
  currentTime <- Sys.time()
  
  .checkArgs.ACtest()
  
  inputType <- ""
  if("VARfit" %in% class(fit)){
    inputType <- "VARfit"
  } else if("varest" %in% class(fit)){ 
    inputType <- "varest"
  }
  
  if(inputType == "VARfit"){
    
    Z <- fit$Z
    K <- fit$K
    p <- fit$p
    N <- fit$N
    NnonLagVar <- fit$NnonLagVar
    const <- fit$const
    trend <- fit$trend
    resi <- fit$resid
    
  } else if(inputType == "varest"){
    
    K <- fit$K
    p <- fit$p
    N <- fit$totobs
    
    NnonLagVar <- ncol(fit$datamat) - K * (p + 1)
    
    if(fit$type == "const") {const <- TRUE; trend <- FALSE}
    if(fit$type == "trend") {const <- FALSE; trend <- TRUE}
    if(fit$type == "both") {const <- TRUE; trend <- TRUE}
    if(fit$type == "none") {const <- FALSE; trend <- FALSE}
    
    exogen <- NULL
    if(fit$type != "none") exogen <- as.matrix(fit$datamat[ , -(1:((p + 1) * K + const + trend))])
    
    Z <- NULL
    if(fit$type != "none") Z <- cbind(Z, 
                                      as.matrix(fit$datamat[ , K * (p + 1) + 1:(const + trend)]))
    if(!is.null(exogen)) Z <- cbind(Z,
                                    exogen)
    Z <- cbind(Z,
               as.matrix(fit$datamat[ , K + 1:(p * K)]))
    
    resi <- matrix(nrow = N - p, ncol = K)
    for(i in 1:K){
      
      resi[ , i] <- fit$varresult[[i]]$residuals
    }
    
    
  } else stop("wrong 'fit' argument class")
  

  #LMresults <- .LMtest(z = Z, epshat = resi, h = h, K = K, N = N, p = p, HCtype = HCtype, univariate = univariate)
  
  epshat <- rbind(matrix(0, nrow = h, ncol = K), resi)
  z_e <- cbind(Z, embed(epshat, (h + 1))[, -(1:K)])
  LMresults <- .ACtestCpp(z_e, Z, resi, h, univariate, 
                          "LM" %in% HCtype, "HC0" %in% HCtype, 
                          "HC1" %in% HCtype, "HC2" %in% HCtype, "HC3" %in% HCtype)
  
  Q <- pValues <- NULL
  if(!(univariate == "only")){
    
    Q <- LMresults[[1]]
    pValues <- 1 - pchisq(q = Q, df = h * K^2)
    rownames(pValues) <- "PV"
  }
  
  uniQ <- unipValues <- NULL
  if(univariate == TRUE || univariate == "only"){
    
    uniQ <- LMresults[[2]]
    unipValues <- 1 - pchisq(q = uniQ, df = h)
    
    colnames(uniQ) <- colnames(unipValues) <- colnames(Q)
    rownames(uniQ) <- rownames(unipValues) <- colnames(fit$y)
  }

  returnValue <- list(fit = fit,
                      inputType = inputType,
                      HCtype = HCtype,
                      h = h,
                      pValues = pValues,
                      Q = Q,
                      unipValues = unipValues,
                      uniQ = uniQ,
                      univariate = univariate,
                      description = paste("Test run at", currentTime, "by user", Sys.info()[["user"]]),
                      time = difftime(Sys.time(), currentTime, units = "secs"),
                      call = match.call()
                      )
  
  
  class(returnValue) <-  c("ACtest", class(returnValue))
  
  print(returnValue)
  invisible(returnValue)
}

.checkArgs.ACtest <- function(){
  # this function checks the arguments in the parent function 'ACtest'
  
  ACtestEnv <- parent.frame()
  
  # list of error messages:
  errs <- list()
  
  if(!any(c("VARfit", "varest") %in% class(ACtestEnv$fit))){
    stop("'fit' must be an object of class 'VARfit', as returned by VARfit(),
         or an object of class 'varest', as returned by VAR() in package 'vars'")
  } 
  
  
  # checks 'h':
  if(!is.numeric(ACtestEnv$h)){
    errs <- append(errs, "'h' must be a positive integer")
  } else if(ACtestEnv$h %% 1 != 0 || ACtestEnv$h < 1){
    errs <- append(errs, "'h' must be a positive integer")
  } 
  
  # provides the possibility of entering truncated and/or case mismatched arguments:
  ACtestEnv$HCtype <- match.arg(toupper(ACtestEnv$HCtype), c("LM", "HC0", "HC1", "HC2", "HC3"), several.ok = TRUE)
  # sorts the 'HCtype' arguments:
  ACtestEnv$HCtype <- c("LM", "HC0", "HC1", "HC2", "HC3")[c("LM", "HC0", "HC1", "HC2", "HC3") %in% ACtestEnv$HCtype]
  
  # checks 'univariate':
  if(!is.logical(ACtestEnv$univariate) && !(ACtestEnv$univariate) == "only") 
    errs <- append(errs, "'univariate' must be either FALSE, TRUE, or \"only\"")
  if("varest" %in% class(ACtestEnv$fit)){
    if(ACtestEnv$univariate == TRUE || ACtestEnv$univariate == "only"){
      warning("the 'univariate' version of the test only works for 'fit' of class 'VARfit'")
    } 
  }
  
  # possibly prints the errors and stops the function:
  if(length(errs) > 0){
    
    msg <- "\nInput errors in 'ACtest()':"
    
    for(i in 1:length(errs)) msg <- paste0(msg, "\n", i, ". ", errs[[i]])
    
    stop(msg)
  }
}

print.ACtest <- function(x, ...){
  
  if(x$univariate != "only"){
    cat("------ Multivariate test for error autocorrelations (AC) ------\n")
    cat("h:", x$h, "\n\n")
    
    mat <- matrix(NA, length(x$HCtype), 3)
    mat[ , 1] <- x$Q[c("LM", "HC0", "HC1", "HC2", "HC3") %in% x$HCtype]
    mat[ , 2] <- x$h * x[[1]]$K^2
    mat[ , 3] <- x$pValues[c("LM", "HC0", "HC1", "HC2", "HC3") %in% x$HCtype]
    rownames(mat) <- x$HCtype
    colnames(mat) <- c("Q", "df", "P.values")
    
    printCoefmat(mat, P.values = TRUE, has.Pvalue = TRUE)
    cat("--------------------------------------------------------------\n\n")
  }
  
  if(isTRUE(x$univariate) || x$univariate == "only"){
    
    cat("---------------- Univariate error AC tests -------------------\n")
    cat("h:", x$h, "\n")
    
    for(i in 1:x[[1]]$K){
      cat("         ", colnames(x[[1]]$y)[i], ":\n", sep = "")
      
      mat <- matrix(NA, length(x$HCtype), 3)
      mat[ , 1] <- x$uniQ[i, c("LM", "HC0", "HC1", "HC2", "HC3") %in% x$HCtype]
      mat[ , 2] <- x$h
      mat[ , 3] <- x$unipValues[i, c("LM", "HC0", "HC1", "HC2", "HC3") %in% x$HCtype]
      rownames(mat) <- x$HCtype
      colnames(mat) <- c("Q", "df", "P.values")
      
      printCoefmat(mat, P.values = TRUE, has.Pvalue = TRUE)
      
    }
    
    cat("--------------------------------------------------------------\n")
    
  }
  
  cat("Test time (secs):", x$time)
}