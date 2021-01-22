wildBoot <- function(test, WBtype = c("recursive", "fixed"), B = 199,  WBdist = c("rademacher", "normal", "mammen"), 
                     HCtype = c("LM", "HC0", "HC1", "HC2", "HC3"), univariate = FALSE){
  
  .checkArgs.wildBoot()
  
  if(test$inputType == "VARfit"){
    
    p <- test$fit$p
    N <- test$fit$N
    const <- test$fit$const
    trend <- test$fit$trend
    exogen <- test$fit$exogen
    h <- test$h
    Z <- test$fit$Z
    y <- test$fit$y
    resi <- test$fit$resid 
    coef <- test$fit$coef
    K <- test$fit$K
    
  } else if(test$inputType == "varest"){
    
    p = test$fit$p
    K = test$fit$K
    N = test$fit$totobs
    h <- test$h
    
    NnonLagVar <- ncol(test$fit$datamat) - K * (p + 1)
    
    if(test$fit$type == "const") {const <- TRUE; trend <- FALSE}
    if(test$fit$type == "trend") {const <- FALSE; trend <- TRUE}
    if(test$fit$type == "both") {const <- TRUE; trend <- TRUE}
    if(test$fit$type == "none") {const <- FALSE; trend <- FALSE}
    
    exogen <- NULL
    if(NnonLagVar - const - trend > 0) exogen <- as.matrix(test$fit$datamat[ , -(1:((p + 1) * K + const + trend))])
    
    Z <- NULL
    if(test$fit$type != "none") Z <- cbind(Z, as.matrix(test$fit$datamat[ , K * (p + 1) + 1:(const + trend)]))
    if(!is.null(exogen)) Z <- cbind(Z, exogen)
    Z <- cbind(Z, as.matrix(test$fit$datamat[ , K + 1:(p * K)]))
    
    # adds 0's to the first p rows of 'exogen':
    if(NnonLagVar - const - trend > 0) exogen <- rbind(matrix(0, nrow = p, ncol = ncol(exogen)) , exogen)
    
    resi <- matrix(nrow = N - p, ncol = K)
    for(i in 1:K) resi[ , i] <- test$fit$varresult[[i]]$residuals
    
    coef <- matrix(nrow = ncol(test$fit$datamat) - K, ncol = K)
    for(i in 1:K){
      
      if(NnonLagVar > 0) coef[1:NnonLagVar , i] <- test$fit$varresult[[i]]$coefficients[-(1:(K * p))]
      coef[-(1:NnonLagVar) , i] <- test$fit$varresult[[i]]$coefficients[1:(K * p)]
    }
    
    y <- test$fit$y
    
  } else stop("wrong inputType in the 'test' object")
  
  WBr.Q <- WBf.Q <- matrix(nrow = B, ncol = 5)
  colnames(WBr.Q) <- colnames(WBf.Q) <- c("LM", "HC0", "HC1", "HC2", "HC3")
  
  if(univariate != FALSE){
    
    uniList <- list()
    for(i in 1:K){
      
      if("recursive" %in% WBtype){
        uniList <- append(x = uniList,
                          values = list(matrix(nrow = B, ncol = 5)))
        names(uniList)[length(uniList)] <- paste0("uni", i, "WBr.Q")
      }
      if("fixed" %in% WBtype){
        
        uniList <- append(x = uniList,
                          values = list(matrix(nrow = B, ncol = 5)))
        names(uniList)[length(uniList)] <- paste0("uni", i, "WBf.Q")
      }
    }
  }
  
  startTime <- Sys.time()
  
  WBr.pv <- WBf.pv <- NULL
  numberOfErrors <- 0
  numberOfNA <- 0
  
  
  cat("\nWild Bootstrap simulations started at", format(startTime), "\n")
  
  percDone <- 10
  
  # runs the first two bootstrap simulation (repeats again if it failes):
  b <- 1
  while(b <= 2){
    wbtestTemp <- try(.WBtest(), 
                      silent = FALSE)
    
    # reruns if there were any errors:
    if(inherits(wbtestTemp, "try-error")){
      
      numberOfErrors <- numberOfErrors + 1
      next
      
    } else{
      if(univariate != "only"){
        if(anyNA(wbtestTemp[[1]][c("recursive" %in% WBtype, "fixed" %in% WBtype), 
                                 c("LM", "HC0", "HC1", "HC2", "HC3") %in% HCtype])){
          numberOfNA <- numberOfNA + 1
          next
        }
      }
      if(univariate != FALSE){
        if(anyNA(wbtestTemp[[2]][rep("recursive" %in% WBtype, K), 
                                 c("LM", "HC0", "HC1", "HC2", "HC3") %in% HCtype]) ||
           anyNA(wbtestTemp[[3]][rep("fixed" %in% WBtype, K), 
                                 c("LM", "HC0", "HC1", "HC2", "HC3") %in% HCtype])){
          numberOfNA <- numberOfNA + 1
          next
        }
      }   
    }
  
    if(univariate != "only"){
      WBr.Q[b, ] <- wbtestTemp[[1]][1, ]
      WBf.Q[b, ] <- wbtestTemp[[1]][2, ]
    } 
    if(univariate != FALSE){
      
      for(i in 1:K){
        if("recursive" %in% WBtype) uniList[[paste0("uni", i, "WBr.Q")]][b, ] <- wbtestTemp[[2]][i, ]
        if("fixed" %in% WBtype) uniList[[paste0("uni", i, "WBf.Q")]][b, ] <- wbtestTemp[[3]][i, ]
      }
    }
    
    if(b == 1) startTime2ndB <- Sys.time()
    b <- b + 1
  }
  
  # estimates the time needed to perform the B bootstrap tests:
  timeEst <- difftime(Sys.time(), startTime2ndB)
  timeEst <- round(difftime(startTime2ndB + timeEst * B, startTime2ndB), 1)
  cat("\nEstimated time to complete the", B, "bootstrap simulations:", format(timeEst), "\n")
  cat("Running Bootstrap: ")
  
  # runs the rest of the B bootstrap simulations:
  b <- 3
  while(b <= B){
    
    wbtestTemp <- try(.WBtest(), 
                      silent = FALSE)
    # reruns if there were any errors:
    if(inherits(wbtestTemp, "try-error")){
      
      numberOfErrors <- numberOfErrors + 1
      next
      
    } else{
      if(univariate != "only"){
        if(anyNA(wbtestTemp[[1]][c("recursive" %in% WBtype, "fixed" %in% WBtype), 
                                 c("LM", "HC0", "HC1", "HC2", "HC3") %in% HCtype])){
          numberOfNA <- numberOfNA + 1
          next
        }
      }
      if(univariate != FALSE){
        if(anyNA(wbtestTemp[[2]][rep("recursive" %in% WBtype, K), 
                                 c("LM", "HC0", "HC1", "HC2", "HC3") %in% HCtype]) ||
           anyNA(wbtestTemp[[3]][rep("fixed" %in% WBtype, K), 
                                 c("LM", "HC0", "HC1", "HC2", "HC3") %in% HCtype])){
          numberOfNA <- numberOfNA + 1
          next
        }
      }   
    }
    
    if(univariate != "only"){
      WBr.Q[b, ] <- wbtestTemp[[1]][1, ]
      WBf.Q[b, ] <- wbtestTemp[[1]][2, ]
    } 
    if(univariate != FALSE){
      
      for(i in 1:K){
        if("recursive" %in% WBtype) uniList[[paste0("uni", i, "WBr.Q")]][b, ] <- wbtestTemp[[2]][i, ]
        if("fixed" %in% WBtype) uniList[[paste0("uni", i, "WBf.Q")]][b, ] <- wbtestTemp[[3]][i, ]
      }
    }
    
    if((b / B) >= (percDone / 100)) {cat(percDone, "% ", sep = ""); percDone = percDone + 10}
    b <- b + 1
    
  }
  
  cat("\n\n")
  
  if("recursive" %in% WBtype){
    WBr.pv <- (colSums(WBr.Q >= matrix(test$Q, nrow = nrow(WBr.Q), ncol = 5, byrow = TRUE)) + 1) / (B + 1)
    names(WBr.pv) <- c("LM", "HC0", "HC1", "HC2", "HC3")
  }
  if("fixed" %in% WBtype){
    WBf.pv <- (colSums(WBf.Q >= matrix(test$Q, nrow = nrow(WBf.Q), ncol = 5, byrow = TRUE)) + 1) / (B + 1)
    names(WBf.pv) <- c("LM", "HC0", "HC1", "HC2", "HC3")
  }
  
  if(univariate != FALSE){
    
    for(i in 1:K){
      
      if("recursive" %in% WBtype){
        uniWBr.pv <- (colSums(uniList[[paste0("uni", i, "WBr.Q")]] >= matrix(test$uniQ[i, ], nrow = B, ncol = 5, byrow = TRUE)) + 1) / (B + 1)
        names(uniWBr.pv) <- c("LM", "HC0", "HC1", "HC2", "HC3")
        uniList <- append(x = uniList,
                          values = list(uniWBr.pv))
        names(uniList)[length(uniList)] <- paste0("uni", i, "WBr.pv")
      }
      if("fixed" %in% WBtype){
        uniWBf.pv <- (colSums(uniList[[paste0("uni", i, "WBf.Q")]] >= matrix(test$uniQ[i, ], nrow = B, ncol = 5, byrow = TRUE)) + 1) / (B + 1)
        names(uniWBf.pv) <- c("LM", "HC0", "HC1", "HC2", "HC3")
        uniList <- append(x = uniList,
                          values = list(uniWBf.pv))
        names(uniList)[length(uniList)] <- paste0("uni", i, "WBf.pv")
      }
    }
  }
  
  time = difftime(Sys.time(), startTime)

  matchCall <- match.call()
  returnValue <- .make.wildBoot()
  
  print(returnValue)
  invisible(returnValue)
  
}

.make.wildBoot <- function(){
  
  pEnv <- parent.frame()
  
  returnValue <- list(test = pEnv$test,
                      WBtype = pEnv$WBtype,
                      B = pEnv$B,
                      WBdist = pEnv$WBdist,
                      HCtype = pEnv$HCtype,
                      univariate = pEnv$univariate,
                      description = paste("Estimated at", pEnv$startTime, "by user", Sys.info()[["user"]]),
                      time = pEnv$time,
                      call = pEnv$matchCall,
                      numberOfErrors = pEnv$numberOfErrors,
                      numberOfNA = pEnv$numberOfNA)
  
  if("recursive" %in% pEnv$WBtype){
    returnValue[[length(returnValue) + 1]] <- pEnv$WBr.Q
    names(returnValue)[[length(returnValue)]] <- "WBr.Q"
    returnValue[[length(returnValue) + 1]] <- pEnv$WBr.pv
    names(returnValue)[[length(returnValue)]] <- "WBr.pv"
  }
  if("fixed" %in% pEnv$WBtype){
    returnValue[[length(returnValue) + 1]] <- pEnv$WBf.Q
    names(returnValue)[[length(returnValue)]] <- "WBf.Q"
    returnValue[[length(returnValue) + 1]] <- pEnv$WBf.pv
    names(returnValue)[[length(returnValue)]] <- "WBf.pv"
  }
  if(pEnv$univariate != FALSE){
    returnValue[[length(returnValue) + 1]] <- pEnv$uniList
    names(returnValue)[[length(returnValue)]] <- "uniList"
  }
  
  class(returnValue) <-  c("wildBoot", class(returnValue))
  return(returnValue)
}

.checkArgs.wildBoot <- function(){
  # this function checks the arguments in the parent function 'wildBoot'
  
  wildBootEnv <- parent.frame()
  
  # list of error messages:
  errs <- list()
  
  # checks 'test':
  if(!("ACtest" %in% class(wildBootEnv$test))) stop("'test' object must be of class \"ACtest\", as returned by the ACtest() function")
  
  
  # checks 'B':
  if(!is.numeric(wildBootEnv$B)){
    errs <- append(errs, "'B' must be a positive integer")
  } else if(wildBootEnv$B %% 1 != 0 || wildBootEnv$B < 1){
    errs <- append(errs, "'B' must be a positive integer")
  } 
  
  # provides the possibility of entering truncated and/or case mismatched arguments:
  wildBootEnv$WBtype <- match.arg(tolower(wildBootEnv$WBtype), c("recursive", "fixed"), several.ok = TRUE)
  wildBootEnv$WBdist <- match.arg(tolower(wildBootEnv$WBdist), c("rademacher", "normal", "mammen"), several.ok = FALSE)
  wildBootEnv$HCtype <- match.arg(toupper(wildBootEnv$HCtype), c("LM", "HC0", "HC1", "HC2", "HC3"), several.ok = TRUE)
  # sorts the 'HCtype' arguments:
  wildBootEnv$HCtype <- c("LM", "HC0", "HC1", "HC2", "HC3")[c("LM", "HC0", "HC1", "HC2", "HC3") %in% wildBootEnv$HCtype]
  
  if(is.null(wildBootEnv$HCtype)) wildBootEnv$HCtype <- wildBootEnv$test$HCtype
  HCtypeNotInTest <- wildBootEnv$HCtype[!(wildBootEnv$HCtype %in% wildBootEnv$test$HCtype)]
  
  if(length(HCtypeNotInTest) != 0) errs <- append(errs,
                                                  paste0("The HCtypes ", 
                                                         toString(HCtypeNotInTest), 
                                                         " was not used in the 'test' object"))
  
  # checks 'univariate':
  if(!is.logical(wildBootEnv$univariate) && !(wildBootEnv$univariate) == "only"){ 
    errs <- append(errs, "'univariate' must be either FALSE, TRUE, or \"only\"")
  } else if(wildBootEnv$univariate %in% c(TRUE, "only") && wildBootEnv$test$univariate == FALSE){
    errs <- append(errs, "'univariate' was not used in the 'test' object")
  }
  
  
  # possibly prints the errors and stops the function:
  if(length(errs) > 0){
    
    msg <- "\nInput errors in 'wildBoot()':"
    
    for(i in 1:length(errs)) msg <- paste0(msg, "\n", i, ". ", errs[[i]])
    
    stop(msg)
  }
}

print.wildBoot <- function(x, ...){
  
  if(x$univariate != "only"){
    
    cat("------ Multivariate test for error autocorrelations (AC) ------\n")
    
    cat("\nh:", x$test$h, " B:", x$B, "\n\n")
    
    if("recursive" %in% x$WBtype){
      
      cat("Recursive Wild Bootstrap:\n\n")
      
      mat <- matrix(NA, length(x$HCtype), 4)
      mat[ , 1] <- x$test$Q[c("LM", "HC0", "HC1", "HC2", "HC3") %in% x$HCtype]
      mat[ , 2] <- x$test$h * x$test$fit$K^2
      mat[ , 3] <- x$test$pValues[c("LM", "HC0", "HC1", "HC2", "HC3") %in% x$HCtype]
      mat[ , 4] <- x$WBr.pv[c("LM", "HC0", "HC1", "HC2", "HC3") %in% x$HCtype]
      rownames(mat) <- x$HCtype
      colnames(mat) <- c("Q", "df", "Asy.PV", "WB.PV")
      
      printCoefmat(mat, P.values = TRUE, has.Pvalue = TRUE)
    }
    if("fixed" %in% x$WBtype){
      
      cat("\nFixed Wild Bootstrap:\n\n")
      
      mat <- matrix(NA, length(x$HCtype), 4)
      mat[ , 1] <- x$test$Q[c("LM", "HC0", "HC1", "HC2", "HC3") %in% x$HCtype]
      mat[ , 2] <- x$test$h * x$test$fit$K^2
      mat[ , 3] <- x$test$pValues[c("LM", "HC0", "HC1", "HC2", "HC3") %in% x$HCtype]
      mat[ , 4] <- x$WBf.pv[c("LM", "HC0", "HC1", "HC2", "HC3") %in% x$HCtype]
      rownames(mat) <- x$HCtype
      colnames(mat) <- c("Q", "df", "Asy.PV", "WB.PV")
      
      printCoefmat(mat, P.values = TRUE, has.Pvalue = TRUE)
    }
    cat("--------------------------------------------------------------\n\n")
    
  }
  
  if(x$univariate != FALSE){
    
    cat("---------------- Univariate error AC tests -------------------\n")
    cat("h:", x$test$h, " B:", x$B, "\n")
    
    for(i in 1:x$test$fit$K){
      
      if("recursive" %in% x$WBtype){
        
        cat("\n[", colnames(x$test$fit$y)[i], "] Recursive Wild Bootstrap\n")
        
        mat <- matrix(NA, length(x$HCtype), 4)
        mat[ , 1] <- x$test$uniQ[i, ][c("LM", "HC0", "HC1", "HC2", "HC3") %in% x$HCtype]
        mat[ , 2] <- x$test$h
        mat[ , 3] <- x$test$unipValues[i, ][c("LM", "HC0", "HC1", "HC2", "HC3") %in% x$HCtype]
        mat[ , 4] <- x$uniList[[paste0("uni", i, "WBr.pv")]][c("LM", "HC0", "HC1", "HC2", "HC3") %in% x$HCtype]
        rownames(mat) <- x$HCtype
        colnames(mat) <- c("Q", "df", "Asy.PV", "WB.PV")
        
        printCoefmat(mat, P.values = TRUE, has.Pvalue = TRUE)
      }
      if("fixed" %in% x$WBtype){
        
        cat("\n[", colnames(x$test$fit$y)[i], "] Fixed Wild Bootstrap\n")
        
        mat <- matrix(NA, length(x$HCtype), 4)
        mat[ , 1] <- x$test$uniQ[i, ][c("LM", "HC0", "HC1", "HC2", "HC3") %in% x$HCtype]
        mat[ , 2] <- x$test$h
        mat[ , 3] <- x$test$unipValues[i, ][c("LM", "HC0", "HC1", "HC2", "HC3") %in% x$HCtype]
        mat[ , 4] <- x$uniList[[paste0("uni", i, "WBf.pv")]][c("LM", "HC0", "HC1", "HC2", "HC3") %in% x$HCtype]
        rownames(mat) <- x$HCtype
        colnames(mat) <- c("Q", "df", "Asy.PV", "WB.PV")
        
        printCoefmat(mat, P.values = TRUE, has.Pvalue = TRUE)
      }
    }
    
    cat("--------------------------------------------------------------\n")
    
  }
  
  cat("\n", x$numberOfErrors + x$numberOfNA, " of the bootstrap simulations failed and had to be resimulated", sep = "")
    
  cat("\nTest time:", round(x$time, 1), units(x$time))
}



.WBtest <- function(){
  
  pEnv <- parent.frame()
  
  e <- .rwb(N = pEnv$N - pEnv$p, WBdist = pEnv$WBdist) * pEnv$resi
  
  returnValue <- list(matrix(nrow = 2, ncol = 5),
                      matrix(nrow = pEnv$K, ncol = 5),
                      matrix(nrow = pEnv$K, ncol = 5))
  
  names(returnValue) <- c("multivariate" ,"univariateRecursive", "univariateFixed")
  
  rownames(returnValue[[1]]) <- c("recursive", "fixed")
  rownames(returnValue[[2]]) <- rownames(returnValue[[3]]) <- colnames(pEnv$y)
  
  colnames(returnValue[[1]]) <- colnames(returnValue[[2]]) <- colnames(returnValue[[3 ]]) <- 
    c("LM", "HC0", "HC1", "HC2", "HC3")
  
  if("recursive" %in% pEnv$WBtype){
    # ySim.r <- .varmodel(x0 = matrix(pEnv$y[1:pEnv$p, ], nrow = pEnv$p), e = e, pihat = pEnv$coef,
    #                     p = pEnv$p, const = pEnv$const, trend = pEnv$trend, 
    #                     exogen = if(is.null(pEnv$exogen)) NULL else as.matrix(pEnv$exogen[-(1:pEnv$p), ]))
    
    ySim.r <- .makeVar(matrix(pEnv$y[1:pEnv$p, ], nrow = pEnv$p), e, pEnv$coef,
                       pEnv$p, pEnv$const, pEnv$trend, !is.null(pEnv$exogen), 
                      if(is.null(pEnv$exogen)) matrix() else as.matrix(pEnv$exogen[-(1:pEnv$p), ]))
    
    # appends the lagged ySim.r with the non-endogenous variables:
    if(ncol(pEnv$Z) > pEnv$K * pEnv$p){
      Zsim.r <- cbind(pEnv$Z[, -(ncol(pEnv$Z) + 1 - 1:(pEnv$K * pEnv$p))],
                      embed(ySim.r, pEnv$p + 1)[, -(1:pEnv$K)])
    } else{
      Zsim.r <- embed(ySim.r, pEnv$p + 1)[, -(1:pEnv$K)]
    }
    
    eWb.r <- qr.resid(qr(Zsim.r), ySim.r[-(1:pEnv$p), ])
    
    epshat <- rbind(matrix(0, nrow = pEnv$h, ncol = pEnv$K), eWb.r)
    z_e <- cbind(Zsim.r, embed(epshat, (pEnv$h + 1))[, -(1:pEnv$K)])
    LMtemp <- .ACtestCpp(z_e, Zsim.r, eWb.r, pEnv$h, pEnv$univariate, 
                            "LM" %in% pEnv$HCtype, "HC0" %in% pEnv$HCtype, 
                            "HC1" %in% pEnv$HCtype, "HC2" %in% pEnv$HCtype, "HC3" %in% pEnv$HCtype)
    
    
    # LMtemp <- .LMtest(Zsim.r, eWb.r, h = pEnv$h, K = pEnv$K, N = pEnv$N, p = pEnv$p, 
    #                   HCtype = pEnv$HCtype, univariate = pEnv$univariate)
    
    if(pEnv$univariate != "only") returnValue[[1]][1, ] <- LMtemp[[1]]
    if(pEnv$univariate != FALSE) returnValue[[2]] <- LMtemp[[2]]
      
  }
  if("fixed" %in% pEnv$WBtype){
    
    ySim.f <- pEnv$Z %*% pEnv$coef + e
    eWb.f <- qr.resid(qr(pEnv$Z), ySim.f)
    
    epshat <- rbind(matrix(0, nrow = pEnv$h, ncol = pEnv$K), eWb.f)
    z_e <- cbind(pEnv$Z, embed(epshat, (pEnv$h + 1))[, -(1:pEnv$K)])
    LMtemp <- .ACtestCpp(z_e, pEnv$Z, eWb.f, pEnv$h, pEnv$univariate, 
                         "LM" %in% pEnv$HCtype, "HC0" %in% pEnv$HCtype, 
                         "HC1" %in% pEnv$HCtype, "HC2" %in% pEnv$HCtype, "HC3" %in% pEnv$HCtype)
    
    # LMtemp <- .LMtest(pEnv$Z, eWb.f, h = pEnv$h, K = pEnv$K, N = pEnv$N, p = pEnv$p,
    #                             HCtype = pEnv$HCtype, univariate = pEnv$univariate)
    
    
    if(pEnv$univariate != "only") returnValue[[1]][2, ] <- LMtemp[[1]]
    if(pEnv$univariate != FALSE) returnValue[[3]] <- LMtemp[[2]]
    
  }
  
  return(returnValue)
}

.rwb <- function(N, WBdist){
  # this function returns 'N' random samples from either the "rademacher", "normal" or "mammen" distribution
  
  if(WBdist == "rademacher") return(sample(x = c(-1, 1), 
                                           size = N, 
                                           replace = TRUE))
  
  if(WBdist == "normal") return(rnorm(N))
  
  if(WBdist == "mammen") return(sample(x = c(-(sqrt(5) - 1) / 2,
                                             (sqrt(5) + 1) / 2), 
                                       prob = c((sqrt(5) + 1) / (2 * sqrt(5)),
                                                (sqrt(5) - 1) / (2 * sqrt(5))),
                                       size = N,
                                       replace = TRUE))
  
}