cointBootTest <- function(y, r = "sequence", p, model = 1, signif = .05, dummies = NULL, B = 999, 
                          boot_type = c("B", "WB"), WB_dist = c("rademacher", "normal", "mammen")){
  
  # checks arguments
  .checkArgs.cointBootTest()
  
  k <- ncol(y)
  
  if(r[[1]] == "sequence"){ 
    r_seq <- 0:(ncol(y) - 1)
  } else{
    r_seq <- r
  }
  
  # 'B.Q' and 'WB.Q' are the matrices with the bootstrap Q statistic 
  #   on the rows and null hypothesis r on the columns:
  if("B" %in% boot_type){ 
    B.Q <- matrix(NA, nrow = B, ncol = length(r_seq))
    B.pv <- rep(NA, length(r_seq))
  }
  if("WB" %in% boot_type){ 
    WB.Q <- matrix(NA, nrow = B, ncol = length(r_seq))
    WB.pv <- rep(NA, length(r_seq))
  }
  
  Q <- rep(NA, length(r_seq))
  B.r <- WB.r <- NA
  B.errors <- WB.errors <- 0
  repeate <- FALSE
  
  
  percDone <- 10
  n_simulations_done <- 0
  n_total_simulations <- B * length(r_seq) * length(boot_type)
  first_b_completed <- FALSE
  
  alpha <- beta <- gamma <- rho <- phi <- dummy_coefs <- residuals <- eigen_values <- list()
  
  for(r0 in r_seq){
    
    vecm_est <- try(.estimate_VECM(y = y, r = r0, p = p, model = model, dummy = dummies))
    
    if(inherits(vecm_est, "try-error")){
      stop("failed to estimate VECM")
    }
    
    # gets the eigen values of the VAR(1) companion matrix :
    eigen_values[[length(eigen_values) + 1]] <- .vecm_check_eigen(vecm_est$alpha, vecm_est$beta, vecm_est$gamma, r0, k, p)$eigval

    pos <- which(r0 == r_seq)
    Q[pos] <- vecm_est$Q
    
    if(n_simulations_done == 0){
      cat("\nBootstrap simulations started at", format(startTime <- Sys.time()), "\n")
    } 
    
    if("B" %in% boot_type){
      
      for(b in 1:B){
        
        repeate <- TRUE
        while(repeate){
          # draws the bootstrap errors:
          e <- .rb(vecm_est$residuals, dummies)
          
          # simulates the bootstrap series:
          y_temp <- try(.make_VECM(Ystart = y[1:p, ], e = e, alpha = vecm_est$alpha,
                                   beta = vecm_est$beta, gamma = vecm_est$gamma, rho = vecm_est$rho, phi = vecm_est$phi, model = model))
          
          repeate <- FALSE
          if(inherits(y_temp, "try-error")){
            B.errors <- B.errors + 1
            repeate <- TRUE
          }
          
          # computes the bootstrap Q statistic:
          vecm_est_temp <- .estimate_VECM(y = y_temp, r = r0, p = p, model = model, dummy = matrix(nrow = 0, ncol = 0))
          
          repeate <- FALSE
          if(inherits(vecm_est_temp, "try-error")){
            WB.errors <- WB.errors + 1
            repeate <- TRUE
          } else(
            B.Q[b, pos] <- vecm_est_temp$Q
          )
        }
        
        # estimates the time needed to perform the tests:
        n_simulations_done <- n_simulations_done + 1
        if(n_simulations_done == 50){  
          timeEst <- difftime(Sys.time(), startTime) / 50
          timeEst <- round(difftime(startTime + timeEst * n_total_simulations, startTime), 1)
          cat("\nEstimated time to complete the simulations:", format(timeEst), "\n")
          cat("Running simulations: ")
        }
        # prints the progress:
        if(((n_simulations_done - 1) / n_total_simulations) >= (percDone / 100)) {cat(percDone, "% ", sep = ""); percDone = percDone + 10}
        
      }
      
      B.pv[pos] <- (sum(B.Q[, pos] >= Q[pos]) + 1) / (B + 1)
      if(r[[1]] == "sequence" && B.pv[pos] > signif && is.na(B.r)) B.r <- r0
      
    }
    if("WB" %in% boot_type){
      
      for(b in 1:B){
        
        repeate <- TRUE
        while(repeate){
          # draws the wild bootstrap errors:
          e <- .rwb(N = nrow(vecm_est$residuals), WBdist = WB_dist) * vecm_est$residuals 
          
          # simulates the bootstrap series:
          y_temp <- try(.make_VECM(Ystart = y[1:p, ], e = e, alpha = vecm_est$alpha,
                                   beta = vecm_est$beta, gamma = vecm_est$gamma, rho = vecm_est$rho, phi = vecm_est$phi, model = model))
          repeate <- FALSE
          if(inherits(y_temp, "try-error")){
            numberOfErrors <- numberOfErrors + 1
            repeate <- TRUE
          }
          
          # computes the bootstrap Q statistic:
          vecm_est_temp <- .estimate_VECM(y = y_temp, r = r0, p = p, model = model, dummy = matrix(nrow = 0, ncol = 0))
          
          repeate <- FALSE
          if(inherits(vecm_est_temp, "try-error")){
            numberOfErrors <- numberOfErrors + 1
            repeate <- TRUE
          } else(
            WB.Q[b, pos] <- vecm_est_temp$Q
          )
        }
        
        # estimates the time needed to perform the tests:
        n_simulations_done <- n_simulations_done + 1
        if(n_simulations_done == 50){ 
          timeEst <- difftime(Sys.time(), startTime) / 50
          timeEst <- round(difftime(startTime + timeEst * n_total_simulations, startTime), 1)
          cat("\nEstimated time to complete the simulations:", format(timeEst), "\n")
          cat("Running simulations: ")
        }
        # prints the progress:
        if(((n_simulations_done - 1) / n_total_simulations) >= (percDone / 100)) {cat(percDone, "% ", sep = ""); percDone = percDone + 10}
        
      }
      
      WB.pv[pos] <- (sum(WB.Q[, pos] >= Q[pos]) + 1) / (B + 1)
      if(r[[1]] == "sequence" && WB.pv[pos] > signif && is.na(WB.r)) WB.r <- r0
    }
    
    
    alpha[[length(alpha) + 1]] <- vecm_est$alpha
    beta[[length(beta) + 1]] <- vecm_est$beta
    gamma[[length(gamma) + 1]] <- vecm_est$gamma
    rho[[length(rho) + 1]] <- vecm_est$rho
    phi[[length(phi) + 1]] <- vecm_est$phi
    dummy_coefs[[length(dummy_coefs) + 1]] <- vecm_est$dummy_coefs
    residuals[[length(residuals) + 1]] <- vecm_est$residuals
  }
  
  cat("100%\n")
  
  if("B" %in% boot_type){
    if(r[[1]] == "sequence" && is.na(B.r)) B.r <- ncol(y)
  } 
  if("WB" %in% boot_type){
    if(r[[1]] == "sequence" && is.na(WB.r)) WB.r <- ncol(y)
  }
  
  # estimates with r = k to get the full alpha, beta and rho parameter matrices:
  vecm_est <- try(.estimate_VECM(y = y, r = ncol(y), p = p, model = model, dummy = dummies))
  if(inherits(vecm_est, "try-error")){
    stop("failed to estimate VECM")
  }
  
  # fills the list of return values:
  return_list <- list()
  
  return_list[[length(return_list) + 1]] <- vecm_est$eigval_dbl
  names(return_list)[[length(return_list)]] <- "eigen_val"
  #return_list[[length(return_list) + 1]] <- vecm_est$eigvec_dbl %*% diag(1/diag(vecm_est$eigvec_dbl))
  return_list[[length(return_list) + 1]] <- vecm_est$eigvec_dbl
  names(return_list)[[length(return_list)]] <- "eigen_vec"

  return_list[[length(return_list) + 1]] <- vecm_est$alpha
  names(return_list)[[length(return_list)]] <- "alpha"
  return_list[[length(return_list) + 1]] <- vecm_est$beta
  names(return_list)[[length(return_list)]] <- "beta"
  return_list[[length(return_list) + 1]] <- gamma
  names(return_list)[[length(return_list)]] <- "gamma"
  return_list[[length(return_list) + 1]] <- vecm_est$rho
  names(return_list)[[length(return_list)]] <- "rho"
  return_list[[length(return_list) + 1]] <- phi
  names(return_list)[[length(return_list)]] <- "phi"
  return_list[[length(return_list) + 1]] <- dummy_coefs
  names(return_list)[[length(return_list)]] <- "dummy_coefs"
  return_list[[length(return_list) + 1]] <- residuals
  names(return_list)[[length(return_list)]] <- "residuals"
  
  return_list[[length(return_list) + 1]] <- y
  names(return_list)[[length(return_list)]] <- "y"
  return_list[[length(return_list) + 1]] <- r
  names(return_list)[[length(return_list)]] <- "r"
  return_list[[length(return_list) + 1]] <- p
  names(return_list)[[length(return_list)]] <- "p"
  return_list[[length(return_list) + 1]] <- model
  names(return_list)[[length(return_list)]] <- "model"
  return_list[[length(return_list) + 1]] <- signif
  names(return_list)[[length(return_list)]] <- "signif"
  return_list[[length(return_list) + 1]] <- dummies
  names(return_list)[[length(return_list)]] <- "dummies"
  return_list[[length(return_list) + 1]] <- B
  names(return_list)[[length(return_list)]] <- "B"
  return_list[[length(return_list) + 1]] <- boot_type
  names(return_list)[[length(return_list)]] <- "boot_type"
  return_list[[length(return_list) + 1]] <- WB_dist
  names(return_list)[[length(return_list)]] <- "WB_dist"
  
  return_list[[length(return_list) + 1]] <- Q
  names(return_list)[[length(return_list)]] <- "Q"
  
  
  if("B" %in% boot_type){
    return_list[[length(return_list) + 1]] <- B.Q
    names(return_list)[[length(return_list)]] <- "B.Q"
    return_list[[length(return_list) + 1]] <- B.pv
    names(return_list)[[length(return_list)]] <- "B.pv"
    return_list[[length(return_list) + 1]] <- B.r
    names(return_list)[[length(return_list)]] <- "B.r"
  } 
  if("WB" %in% boot_type){
    return_list[[length(return_list) + 1]] <- WB.Q
    names(return_list)[[length(return_list)]] <- "WB.Q"
    return_list[[length(return_list) + 1]] <- WB.pv
    names(return_list)[[length(return_list)]] <- "WB.pv"
    return_list[[length(return_list) + 1]] <- WB.r
    names(return_list)[[length(return_list)]] <- "WB.r"
  }
  
  return_list[[length(return_list) + 1]] <- r_seq
  names(return_list)[[length(return_list)]] <- "r_seq"
  return_list[[length(return_list) + 1]] <- B.errors
  names(return_list)[[length(return_list)]] <- "B.errors"
  return_list[[length(return_list) + 1]] <- WB.errors
  names(return_list)[[length(return_list)]] <- "WB.errors"
  
  return_list[[length(return_list) + 1]] <- eigen_values
  names(return_list)[[length(return_list)]] <- "companion_eigen"

  class(return_list) <-  c("cointBootTest", class(return_list))
  return(return_list)
}


.checkArgs.cointBootTest <- function(){
  # this function checks the arguments in the parent function 'wildBoot'
  
  parent_env <- parent.frame()
  
  # list of error messages:
  errs <- list()
  
  # checks 'B':
  if(!is.numeric(parent_env$B)){
    errs <- append(errs, "'B' must be a positive integer")
  } else if(parent_env$B %% 1 != 0 || parent_env$B < 1){
    errs <- append(errs, "'B' must be a positive integer")
  }
  
  # checks 'p':
  if(!is.numeric(parent_env$p)){
    stop("'p' must be a positive integer")
  } else if(parent_env$p %% 1 != 0 || parent_env$p < 1){
    stop("'p' must be a positive integer")
  }
  
  # checks 'signif':
  if(!is.numeric(parent_env$signif)){
    errs <- append(errs, "'signif' must be numeric")
  } else if(parent_env$signif < 0 || parent_env$signif > 1){
    errs <- append(errs, "'signif' must be in [0,1]")
  }
  
  # checks 'model':
  if(!parent_env$model %in% c(1, 2, 3)){
    errs <- append(errs, "'model' must be either 1, 2, or 3")
  }
  
  # matches argument for 'boot_type' and 'WB_dist':
  parent_env$boot_type <- match.arg(toupper(parent_env$boot_type), c("B", "WB"), several.ok = TRUE)
  parent_env$WB_dist <- match.arg(tolower(parent_env$WB_dist), c("rademacher", "normal", "mammen"), several.ok = FALSE)
  
  # checks 'y':
  if(is.null(parent_env$y)) stop("'y' must be provided!")
  if(!is.matrix(parent_env$y)) parent_env$y <- try(as.matrix(parent_env$y))
  if(inherits(parent_env$y, "try-error")) stop("'y' is not a matrix and cant be transformed to one")
  if(!is.numeric(parent_env$y)) stop("'y' must be a numeric matix, vector, or data.frame")
  if(anyNA(parent_env$y)) errs <- append(errs, "'y' cant contain any NA elements")
  
  # checks 'dummies':
  if(!is.null(parent_env$dummies)){
    if(!is.matrix(parent_env$dummies)) parent_env$dummies <- try(as.matrix(parent_env$dummies))
    if(inherits(parent_env$dummies, "try-error")) stop("'dummies' is not a matrix and cant be transformed to one")
    if(!is.numeric(parent_env$dummies)) stop("'dummies' must be a numeric matix, vector, or data.frame")
    if(anyNA(parent_env$dummies)) errs <- append(errs, "'dummies' cant contain any NA elements")
    if(nrow(parent_env$dummies) != nrow(parent_env$y)) errs <- append(errs, "'dummies' must have the same number of rows as 'y")
  } else{
    parent_env$dummies <- matrix(0, nrow = 0, ncol = 0)
  }
  
  # checks 'r':
  if(is.numeric(parent_env$r)){
    if(any(parent_env$r < 0) || any(parent_env$r > ncol(parent_env$y) - 1)){
      errs <- append(errs, "'r' can't be below 0 or above k - 1")
    }
  } else if(!parent_env$r[[1]] == "sequence"){
    errs <- append(errs, "'r' must either be \"sequence\" or a numeric vector")
  }
  
  
  # possibly prints the errors and stops the function:
  if(length(errs) > 0){
    
    msg <- "\nInput errors in 'cointBootTest()':"
    
    for(i in 1:length(errs)) msg <- paste0(msg, "\n", i, ". ", errs[[i]])
    
    stop(msg)
  }
}


print.cointBootTest <- function(x, ...){
  
  cat("\n-------------------Test of cointegration rank-----------------\n\n")
  
  cat("B = ", x$B, ", N = ", nrow(x$y), ", k = ", ncol(x$y), ", p = ", x$p, "\n", sep = "")
  model_choice <- switch(x$model,
                         `1` = "no deterministic component (1)",
                         `2` = "restricted constant (2)",
                         `3` = "restricted linear trend (3)")
  cat("model: ", model_choice, "\n", sep = "")
  
  
  cat("\nTesting rank r against k:\n", sep = "")
  
  if("B" %in% x$boot_type){
    
    cat("-------------------------iid bootstrap------------------------\n")
    mat <- matrix(NA, length(x$r_seq), 3)
    mat[ , 1] <- x$r_seq
    mat[ , 2] <- x$Q
    mat[ , 3] <- x$B.pv
    rownames(mat) <- rep("", length(x$r_seq))
    colnames(mat) <- c("r", "Q", "boot.pv")
    
    printCoefmat(mat, P.values = TRUE, has.Pvalue = TRUE)
    cat("--------------------------------------------------------------\n")
    
  } 
  if("WB" %in% x$boot_type){
    
    cat("------------------------wild bootstrap------------------------\n")
    mat <- matrix(NA, length(x$r_seq), 3)
    mat[ , 1] <- x$r_seq
    mat[ , 2] <- x$Q
    mat[ , 3] <- x$WB.pv
    rownames(mat) <- rep("", length(x$r_seq))
    colnames(mat) <- c("r", "Q", "wb.pv")
    
    printCoefmat(mat, P.values = TRUE, has.Pvalue = TRUE)
    cat("--------------------------------------------------------------\n")
  } 
  
  if(x$r[[1]] == "sequence"){
    cat("\nSequential test conclusion (signif = ", x$signif, "):\n", sep = "")
    if("B" %in% x$boot_type){
      cat("Bootstrap selected rank:", x$B.r, "\n")
    } 
    if("WB" %in% x$boot_type){
      cat("Wild bootstrap selected rank:", x$WB.r, "\n")
    }
  }
  
  if(x$B.errors + x$WB.errors == 0){
    cat("\nAll of the simulations were successful\n")
  } else if(x$B.errors > 0){
    cat("\n", x$B.errors , " of the bootstrap simulations failed and had to be resimulated\n", sep = "")
  } else if(x$WB.errors > 0) {
    cat("\n", x$WB.errors , " of the wild bootstrap simulations failed and had to be resimulated\n", sep = "")
  }
  
  rs_violated <- NULL
  if(x$r[[1]] == "sequence"){
    for(i in seq_along(x$companion_eigen)){
      agien_abs <- sapply(x$companion_eigen[[i]], Mod)
      if(sum(abs(agien_abs - 1) < 1e-06) != ncol(x$y) - i + 1 || any(agien_abs > 1 + 1e-06)){
        rs_violated <- c(rs_violated, i)
      }
    }
  } else{
    for(i in seq_along(x$r)){
      agien_abs <- sapply(x$companion_eigen[[i]], Mod)
      if(sum(abs(agien_abs - 1) < 1e-06) != ncol(x$y) - x$r[i] || any(agien_abs > 1 + 1e-06)){
        rs_violated <- c(rs_violated, x$r[i])
      }
    }
  }
  
  
  if(!is.null(rs_violated)){
    if(length(rs_violated) > 1) rs_violated <- paste(paste(rs_violated[-length(rs_violated)], ", ", sep = ""), paste(rs_violated[length(rs_violated)], ".", sep = ""), sep = "")
    else rs_violated <- paste(rs_violated, ".", sep = "")
    
    cat("\nWarning: the conditions for the roots of the characteristic equation (see step 2 of the bootstrap algorithm in the .pdf help page) were violated for r = ",
        rs_violated, " The eigenvalues can be seen in $companion_eigen of the returned object.\n", sep = "")
  }

}

.rb <- function(resi, dummies){
  # This function draws bootstrap errors from the rows of the residuals 'resi', with replacement.
  #  Rows in 'resi' that corresponds to rows where 'dummies' are 1, are not included in the draws
  
  n <- nrow(resi)
  
  # index vector for rows without dummies:
  if(nrow(dummies) > 0){
    if(ncol(dummies) > 1)  available_rows <- (1:n)[rowSums(dummies[-(1:(nrow(dummies) - n)), ]) == 0]
    else  available_rows <- (1:n)[dummies[-(1:(nrow(dummies) - n)), ] == 0]
  } else{
    available_rows <- 1:n
  } 
  
  # samples from the available rows:
  sampled_rows <- sample(x = available_rows, size = n, replace = TRUE)
  
  # returns the selected rows:
  return(resi[sampled_rows, ])
}
