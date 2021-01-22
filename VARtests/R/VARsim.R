VARsim <- function(N = 200, K = 2, p = 1, const = TRUE, trend = FALSE, exogen = NULL,
                   coef = NULL, dist = "normal", Ystart = NULL, errors = NULL, fittedModel = NULL){
  
  
  if(!is.null(fittedModel)){
    
    if("VARfit" %in% class(fittedModel)){
      
      if(!hasArg(N)) N <- fittedModel$N
      K <- fittedModel$K
      p <- fittedModel$p
      
      const <- fittedModel$const
      trend <- fittedModel$trend
      if(is.null(exogen)) exogen <- fittedModel$exogen
      
      coef <- fittedModel$coef
      
      if(!hasArg(Ystart)) Ystart <- matrix(fittedModel$y[1:p, ], nrow = p, ncol = K)
      
    } else if("varest" %in% class(fittedModel)){
      
      if(!hasArg(N)) N <- fittedModel$totobs
      K <- fittedModel$K
      p <- fittedModel$p
      
      if(fittedModel$type == "const") {const <- TRUE; trend <- FALSE}
      if(fittedModel$type == "trend") {const <- FALSE; trend <- TRUE}
      if(fittedModel$type == "both") {const <- TRUE; trend <- TRUE}
      if(fittedModel$type == "none") {const <- FALSE; trend <- FALSE}
      
      NnonEndo <- ncol(fittedModel$datamat) - K * (p + 1)
      coef <- matrix(nrow = ncol(fittedModel$datamat) - K, ncol = K)
      
      for(i in 1:K){
        if(NnonEndo > 0){ 
          coef[1:NnonEndo , i] <- fittedModel$varresult[[i]]$coefficients[-(1:(K * p))]
          coef[-(1:NnonEndo) , i] <- fittedModel$varresult[[i]]$coefficients[1:(K * p)]
        } else{
          coef[, i] <- fittedModel$varresult[[i]]$coefficients
        }
      }
      
      if(!hasArg(Ystart)) Ystart <- matrix(fittedModel$y[1:p, ], nrow = p, ncol = K)
      
      if(is.null(exogen) && NnonEndo > const + trend) 
        exogen <- as.matrix(fittedModel$datamat[ , -(1:((p + 1) * K + const + trend))])
      
    } else stop("wrong class of 'fittedModel'")
    
  }
  
  #check args:
  .checkArgs.VARsim()
  
  if(is.null(Ystart)) Ystart <- matrix(data = 0, nrow = p, ncol = K)
  if(is.null(errors) && dist == "normal") errors <- matrix(data = rnorm((N - p) * K), nrow = N - p, ncol = K)
  
  y <- .makeVar(Ystart, errors, coef,
               p, const, trend, !is.null(exogen), 
               if(is.null(exogen)) matrix() else as.matrix(exogen))
  
  return(y)
}

.checkArgs.VARsim <- function(){
  # this function checks the arguments in the parent function 'VARsim'
  
  env <- parent.frame()
  
  # list of error messages:
  errs <- list()
  
  # checks 'N':
  if(!is.numeric(env$N)){
    errs <- append(errs, "'N' must be a positive integer")
  } else if(env$N %% 1 != 0 || env$N < 1){
    errs <- append(errs, "'N' must be a positive integer")
  } 
  
  # checks 'K':
  if(!is.numeric(env$K)){
    errs <- append(errs, "'K' must be a positive integer")
  } else if(env$K %% 1 != 0 || env$K < 1){
    errs <- append(errs, "'K' must be a positive integer")
  } 
  
  # checks 'p:
  if(!is.numeric(env$p)){
    errs <- append(errs, "'p' must be a positive integer")
  } else if(env$p %% 1 != 0 || env$p < 1){
    errs <- append(errs, "'p' must be a positive integer")
  } else if(env$p > env$N - 1){
    errs <- append(errs, "'p' must be smaller than 'N'")
  }
  
  # checks logical flags:
  if(!(is.logical(env$const) || is.numeric(env$const))){
    errs <- append(errs, "'const' must be TRUE of FALSE")
  }
  if(!(is.logical(env$trend) || is.numeric(env$trend))){
    errs <- append(errs, "'trend' must be TRUE of FALSE")
  }
  
  .stopAndPrint(errs)
  
  
  NexoCol <- 0
  # checks 'exogen':
  if(!is.null(env$exogen)){
    if(!is.matrix(env$exogen) || !is.numeric(env$exogen)){ 
      errs <- append(errs, "'exogen' must be a numerical matrix")
      .stopAndPrint(errs)
    }
    else{ 
      
      # check number of rows:
      if(nrow(env$exogen) == env$N){
        env$exogen <- as.matrix(env$exogen[(env$p + 1):env$N, ])
      } else if(nrow(env$exogen) != env$N - env$p){
        errs <- append(errs, "the number of rows of 'exogen' must be either N or N - p (see documentation)")
      }
      
      NexoCol <- ncol(env$exogen)
    } 
  }
  
  # checks 'coef':
  if(is.null(env$coef) && is.null(env$fittedModel)){
    errs <- append(errs, "'coef' must be provided when 'fittedModel' is missing")
    .stopAndPrint(errs)
  }
  if(!is.null(env$coef)){
    if(!is.matrix(env$coef) || !is.numeric(env$coef)){ 
      errs <- append(errs, "'coef' must be a numerical matrix")
      .stopAndPrint(errs)
    }
    else{ 
      
      # check number of rows:
      if(nrow(env$coef) != env$K * env$p + isTRUE(env$const) + isTRUE(env$trend) + NexoCol){
        errs <- append(errs, "the number of rows of 'coef' are wrong (see documentation)")
      }
      
      # check number of columns:
      if(ncol(env$coef) != env$K) 
        errs <- append(errs, "the number of columns of 'coef' must be K")
    } 
  }
  
  # checks 'Ystart':
  if(!is.null(env$Ystart)){
    if(!is.matrix(env$Ystart) || !is.numeric(env$Ystart)){ 
      errs <- append(errs, "'Ystart' must be a numerical matrix")
    }
    else{ 
      
      # check number of rows:
      if(nrow(env$Ystart) != env$p){
        errs <- append(errs, "the number of rows of 'Ystart' must be 'p' (see documentation)")
      }
      
      # check number of columns:
      if(ncol(env$Ystart) != env$K) 
        errs <- append(errs, "the number of columns of 'Ystart' must be K")
    } 
  }
  
  
  # checks 'errors':
  if(!is.null(env$errors)){
    if(!is.matrix(env$errors) || !is.numeric(env$errors)){ 
      errs <- append(errs, "'errors' must be a numerical matrix")
    }
    else{ 
      
      # check number of rows:
      if(nrow(env$errors) != env$p){
        errs <- append(errs, "the number of rows of 'errors' must be N - p (see documentation)")
      }
      
      # check number of columns:
      if(ncol(env$errors) != env$K) 
        errs <- append(errs, "the number of columns of 'errors' must be K")
    } 
  }
  
  .stopAndPrint(errs)
}

.stopAndPrint <- function(errs){
  
  if(length(errs) > 0){
    msg <- "\nInput errors in 'VARsim()':"
    
    for(i in 1:length(errs)) msg <- paste0(msg, "\n", i, ". ", errs[[i]])
    
    stop(msg)
  }
}