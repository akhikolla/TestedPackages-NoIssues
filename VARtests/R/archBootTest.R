archBootTest <- function(fit, h = 2, B = 499, CA = TRUE, ET = TRUE, MARCH = TRUE, 
                         dist = "norm", skT.param = c(0, 1, 0, 5)){
  
  .checkArgs.archBootTest
  
  currentTime <- Sys.time()
  
  
  if("VARfit" %in% class(fit)){
    inputType <- "VARfit"
  } else if("varest" %in% class(fit)){
    inputType <- "varest"
  } 
  
  if(inputType == "VARfit"){
    
    K <- fit$K
    p <- fit$p
    N <- fit$N
    Nresi <- N - p
    const <- fit$const
    trend <- fit$trend
    exogen <- fit$exogen
    
    resi <- fit$resid
    Z <- fit$Z 
    coef <- fit$coef
    
  } else if(inputType == "varest"){
    
    K <- fit$K
    p <- fit$p
    N <- fit$totobs
    Nresi <- N - p
    
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
    
    coef <- matrix(nrow = ncol(fit$datamat) - K, ncol = K)
    
    
    for(i in 1:K){
      
      if(NnonLagVar > 0) coef[1:NnonLagVar , i] <- fit$varresult[[i]]$coefficients[-(1:(K * p))]
      coef[-(1:NnonLagVar) , i] <- fit$varresult[[i]]$coefficients[1:(K * p)]
      
    }
  } else stop("wrong 'fit' argument class")
  
  # Standardizes the residuals:
  Su <- chol(crossprod(resi)/ Nresi)
  W <- resi %*% solve(Su)
  
  # computes r2 and LMi from auxiliary regression:
  if(CA){
    CA_LMi <- R2i <- numeric(K)
    for (i2 in 1:K) {
      
      ei <- W[, i2]^2
      ei <- embed(ei, (h + 1))
      epsl <- ei[, -1]
      ey <- ei[, 1]
      err <- qr.resid(qr(cbind(1, epsl)), ey)
      Zv <- crossprod(err)
      Z0 <- crossprod(ey - mean(ey))
      
      R2i[i2] <- 1 - Zv/Z0
      CA_LMi[i2] <- R2i[i2] * Nresi
    }
  }
  
  CA_LM <- pchisq(q = max(CA_LMi), df = h, lower.tail = TRUE)
  
  # runs the MARCH test:
  if(MARCH){
    MARCH_LM <- computeMARCH(W, h)
    MARCH_PV <- pchisq(q = MARCH_LM, df = K^2 * (K + 1)^2 * h / 4, lower.tail = FALSE)
  } 
  
  # runs the ET test:
  if(ET){
    ET_LM <- as.numeric(computeET_LM(W, h))
    ET_PV <- pchisq(q = ET_LM, df = K * h, lower.tail = FALSE)
  } 
  
  # sets up the function 'drawWj' to draw random errors:
  if(is.function(dist)){
    drawWj <- dist
  } else if(dist == "norm"){
    drawWj <- .rnormMatrix
  } else if(dist == "skT") {
    drawWj <- .rmstMatrix
  } else(stop("'dist' must be \"norm\", \"skT\", or a function"))
  
  CA_LMijStar <- PijStar <- matrix(nrow = B, ncol = K)
  MARCH_LMStar <- ET_LMStar <- numeric(B)
  Yhat <- Z %*% coef
  
  percDone <- 10
  
  ### start bootstrap:  
  for(b in 1:B) {
    
    # checks that the supplied error function works as intended:
    if(b == 1 && is.function(dist)){
      Wj0 <- try(drawWj()) # simulates the matrix of (standardized) error terms
      
      if(inherits(Wj0, "try-error")) stop("the supplied function in 'dist' failed")
      
      if(is.matrix(Wj0) || is.data.frame(Wj0)){
        
        if(length(Wj0) != Nresi * K) stop("the supplied function in 'dist' doesn't return the correct number of elements")
        if(ncol(Wj0) != K){
          drawWj <- function() matrix(dist(), ncol = K)
          Wj0 <- matrix(Wj0, ncol = K)
          warning("the supplied function in 'dist' returned the wrong dimensions, but was cast into the correct dimensions automatically")
        }
        
      } else if(is.vector(Wj0)){
        if(length(Wj0) != Nresi * K) stop("the supplied function in 'dist' doesn't return the correct number of elements")
        drawWj <- function() matrix(dist(), ncol = K)
        Wj0 <- matrix(Wj0, ncol = K)
        
      } else stop("the supplied function in 'dist' must return a matrix, data.frame or vector of correct size.")
      
      
    } else Wj0 <- drawWj() # simulates the matrix of (standardized) error terms
    
    # Creates the bootstrap VAR series and computes its residuals from its 
    # regression on the designmatrix Z 
    Yj <- Yhat + Wj0 %*% Su
    Uj <- qr.resid(qr(Z), Yj)
    
    # standerdizes the residuals:
    Suj <- chol(crossprod(Uj)/ Nresi)
    Wj <- Uj %*% solve(Suj)
    
    # computes LMi from auxiliary regression:
    if(CA){
      for (k in 1:K) {
        
        ei <- Wj[, k]^2
        ei <- embed(ei, (h + 1))
        epsl <- ei[, -1]
        ey <- ei[, 1]
        err <- qr.resid(qr(cbind(1, epsl)), ey)
        Zv <- crossprod(err)
        Z0 <- crossprod(ey - mean(ey))
        Rv2 <- 1 - Zv/Z0
        
        CA_LMijStar[b, k] <- Nresi * Rv2
        PijStar[b, k] <- pchisq(q = CA_LMijStar[b, k], df = h, lower.tail = FALSE)
      }
    }
    
    # runs the MARCH test:
    if(MARCH) MARCH_LMStar[b] <- computeMARCH(Wj, h)
    
    # runs the ET test:
    if(ET) ET_LMStar[b] <- computeET_LM(Wj, h)
    
    # prints infromation about the progress:
    if(b == 1) startTime2ndB <- Sys.time()
    if(b == 2){
      timeEst <- difftime(Sys.time(), startTime2ndB)
      timeEst <- round(difftime(startTime2ndB + timeEst * B, startTime2ndB), 1)
      cat("\nEstimated time to complete the", B, "bootstrap simulations:", format(timeEst), "\n")
      cat("Running Bootstrap: ")
    }
    if(b > 2){
      if ((b / B) >= (percDone / 100)) {cat(percDone, "% ", sep = ""); percDone = percDone + 10}
    }
    
  } 
  ### end bootstrap
  
  CA_uniBootPV <- numeric(K)
  for(k in 1:K){
    CA_uniBootPV[k] <- (sum(CA_LMijStar[, k] >= CA_LMi[k]) + 1) / (B + 1)
  }
  
  LMjStar <- 1 - apply(PijStar, MARGIN = 1, FUN = min)
  CA_bootPV <- (sum(LMjStar >= CA_LM) + 1) / (B + 1)
  
  if(MARCH) MARCH_bootPV <- (sum(MARCH_LMStar >= MARCH_LM) + 1) / (B + 1)
  if(ET) ET_bootPV <- (sum(ET_LMStar >= ET_LM) + 1) / (B + 1)
  
  returnValue <- list(fit = fit,
                      inputType = inputType,
                      h = h,
                      B = B,
                      K = K,
                      CA = CA, 
                      ET = ET, 
                      MARCH = MARCH,
                      dist = dist,
                      standardizedResi = W)
               
  if(CA){
    returnValue[[length(returnValue) + 1]] <- CA_LM
    names(returnValue)[[length(returnValue)]] <- "CA_LM"
    returnValue[[length(returnValue) + 1]] <- CA_bootPV
    names(returnValue)[[length(returnValue)]] <- "CA_bootPV"
    returnValue[[length(returnValue) + 1]] <- CA_LMi
    names(returnValue)[[length(returnValue)]] <- "CA_LMi"
    returnValue[[length(returnValue) + 1]] <- CA_LMijStar
    names(returnValue)[[length(returnValue)]] <- "CA_LMijStar"
    returnValue[[length(returnValue) + 1]] <- CA_uniBootPV
    names(returnValue)[[length(returnValue)]] <- "CA_uniBootPV"
  }
  
  if(MARCH){
    returnValue[[length(returnValue) + 1]] <- MARCH_LM
    names(returnValue)[[length(returnValue)]] <- "MARCH_LM"
    returnValue[[length(returnValue) + 1]] <- MARCH_PV
    names(returnValue)[[length(returnValue)]] <- "MARCH_PV"
    returnValue[[length(returnValue) + 1]] <- MARCH_bootPV
    names(returnValue)[[length(returnValue)]] <- "MARCH_bootPV"
    returnValue[[length(returnValue) + 1]] <- MARCH_LMStar
    names(returnValue)[[length(returnValue)]] <- "MARCH_LMStar"
  }
  
  if(ET){
    returnValue[[length(returnValue) + 1]] <- ET_LM
    names(returnValue)[[length(returnValue)]] <- "ET_LM"
    returnValue[[length(returnValue) + 1]] <- ET_PV
    names(returnValue)[[length(returnValue)]] <- "ET_PV"
    returnValue[[length(returnValue) + 1]] <- ET_bootPV
    names(returnValue)[[length(returnValue)]] <- "ET_bootPV"
    returnValue[[length(returnValue) + 1]] <- ET_LMStar
    names(returnValue)[[length(returnValue)]] <- "ET_LMStar"
  }
  
  returnValue[[length(returnValue) + 1]] <- paste("Test run at", currentTime, "by user", Sys.info()[["user"]])
  names(returnValue)[[length(returnValue)]] <- "description"
  returnValue[[length(returnValue) + 1]] <- difftime(Sys.time(), currentTime, units = "secs")
  names(returnValue)[[length(returnValue)]] <- "time"
  returnValue[[length(returnValue) + 1]] <- match.call()
  names(returnValue)[[length(returnValue)]] <- "call"
  
  class(returnValue) <-  c("archBootTest", class(returnValue))
  
  print(returnValue)
  invisible(returnValue)
  
}

.checkArgs.archBootTest <- function(){
  # this function checks the arguments in the parent function 'archBootTest'
  
  archBootTestEnv <- parent.frame()
  
  # list of error messages:
  errs <- list()
  
  if(!any(c("VARfit", "varest") %in% class(archBootTestEnv$fit))){
    stop("'fit' must be an object of class 'VARfit', as returned by VARfit(),
         or an object of class 'varest', as returned by VAR() in package 'vars'")
  } 
  
  
  # checks 'h':
  if(!is.numeric(archBootTestEnv$h)){
    errs <- append(errs, "'h' must be a positive integer")
  } else if(archBootTestEnv$h %% 1 != 0 || archBootTestEnv$h < 1){
    errs <- append(errs, "'h' must be a positive integer")
  } 
  
  # checks 'B':
  if(!is.numeric(archBootTestEnv$B)){
    errs <- append(errs, "'B' must be a positive integer")
  } else if(archBootTestEnv$B %% 1 != 0 || archBootTestEnv$B < 1){
    errs <- append(errs, "'B' must be a positive integer")
  } 
  
  # checks 'dist':
  
  
  # checks 'dist':
  if(!is.function(archBootTestEnv$dist) && archBootTestEnv$dist %in% c("norm", "skT")) 
    errs <- append(errs, "'dist' must be \"norm\", \"skT\", or a function")
  
  # checks 'skT.param':
  if(length(archBootTestEnv$skT.param) != 4){
    errs <- append(errs, "'skT.param' must be a numeric vector of length four")
  } else if(!is.numeric(archBootTestEnv$skT.param)){
    errs <- append(errs, "'skT.param' must be numeric")
  }
  
  # possibly prints the errors and stops the function:
  if(length(errs) > 0){
    
    msg <- "\nInput errors in 'archBootTest()':"
    
    for(i in 1:length(errs)) msg <- paste0(msg, "\n", i, ". ", errs[[i]])
    
    stop(msg)
  }
  
  }

print.archBootTest <- function(x, ...){
  
  cat("\nTest for ARCH errors in VAR models:\n\n")
  
  cat("h:", x$h)
  cat("\nB:", x$B)
  if(!is.function(x$dist)){
    cat("\nError distribution:", x$dist,"\n\n")
  } else{
    cat("\nUser supplied error distribution function.\n\n")
  }
  
  mat <- matrix(nrow = ifelse(x$CA, 1, 0) + ifelse(x$ET, 1, 0) + ifelse(x$MARCH, 1, 0), ncol = 4)
  
  row <- 1
  rownames_text <- NULL
  
  if(x$CA){
    
    mat[row , 1] <- x$CA_LM
    mat[row , 4] <- x$CA_bootPV
    
    row <- row + 1
    rownames_text <- c(rownames_text, "CA")
  }
  if(x$ET){
    
    mat[row , 1] <- x$ET_LM
    mat[row , 2] <- x$K * x$h
    mat[row , 3] <- x$ET_PV
    mat[row , 4] <- x$ET_bootPV
    
    row <- row + 1
    rownames_text <- c(rownames_text, "ET")
  }
  if(x$MARCH){
    
    mat[row , 1] <- x$MARCH_LM
    mat[row , 2] <- x$K^2 * (x$K + 1)^2 * x$h / 4
    mat[row , 3] <- x$MARCH_PV
    mat[row , 4] <- x$MARCH_bootPV
    
    rownames_text <- c(rownames_text, "MARCH")
  }
  
  rownames(mat) <- rownames_text
  colnames(mat) <- c("LM", "df", "Asy.PV", "Boot.PV")
  
  cat("-------------------- Multivariate tests: ---------------------\n")
  printCoefmat(mat, P.values = TRUE, has.Pvalue = TRUE)
  cat("--------------------------------------------------------------\n\n\n")
  
  if(x$CA){  
    mat <- matrix(nrow = x$fit$K, ncol = 4)
    mat[ , 1] <- x$CA_LMi
    mat[ , 2] <- x$h
    mat[ , 3] <- pchisq(x$CA_LMi, df = x$h, lower.tail = FALSE)
    mat[ , 4] <- x$CA_uniBootPV
    
    rownames(mat) <- colnames(x$fit$y)
    colnames(mat) <- c("LMi", "df", "Asy.PV", "Boot.PV")
    
    cat("---------- Catani and Ahlgren equation by equation: ----------\n")
    printCoefmat(mat, P.values = TRUE, has.Pvalue = TRUE)
    cat("--------------------------------------------------------------\n")
  }
  
  cat("\nTest time (secs):", x$time, "\n")
}

.rnormMatrix <- function(){
  archBootTestEnv <- parent.frame()
  matrix(rnorm(archBootTestEnv$Nresi * archBootTestEnv$K), ncol = archBootTestEnv$K)
} 

.rmstMatrix <- function(){
  archBootTestEnv <- parent.frame()
  
  matrix(sn::rmst(n = archBootTestEnv$Nresi * archBootTestEnv$K,
                  xi = archBootTestEnv$skT.param[1],
                  Omega = archBootTestEnv$skT.param[2],
                  alpha = archBootTestEnv$skT.param[3],
                  nu = archBootTestEnv$skT.param[4]), 
         ncol = archBootTestEnv$K)
} 
