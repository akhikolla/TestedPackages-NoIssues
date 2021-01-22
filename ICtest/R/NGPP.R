NGPP <-
function(X, k, nl = c("skew", "pow3"), alpha = 0.8, method="symm", eps = 1e-6, verbose = FALSE, maxiter = 100){
  
  # Parse non-linearities to integers
  nl_int <- match(nl, c("skew", "pow3", "tanh", "gauss"))
  
  if(length(alpha) != length(nl_int)){
    alpha <- c(alpha, 1 - sum(alpha))
  }
  
  # Whiten data
  MU <- colMeans(X)
  p <- ncol(X)
  cov_x <- cov(X)
  X <- whiten(X)
  
  # Estimate the initial k rows of U with FOBI
  FOBI_U <- FOBI(X)$W
  FOBI_S <- X%*%t(FOBI_U)
  FOBI_obj <- computeObj_C(FOBI_S, nl_int, alpha)
  U <- FOBI_U[(order(FOBI_obj, decreasing=TRUE)[1:k]), , drop=FALSE]
  
  # Deflation-based estimation
  if(method == "defl"){
    
    orth <- diag(p)
    
    for(i in 1:k){
      if(i > 1){
        orth <- orth - U[(i-1), ]%*%t(U[(i-1), ])
        
      }
      u <- U[i, ]
      
      crit <- 1
      iter <- 0
      while(crit > eps){
        
        # Orthogonalization
        Tvec <- computeTVec_C(u, X, nl_int, alpha)
        u2 <- orth%*%Tvec
        u2 <- u2/sqrt(sum(u2^2))
        
        # Criterion evaluation
        crit <- sqrt(min(sum((u2 - u)^2), sum((u2 + u)^2)))
        u <- u2
        
        iter <- iter + 1
        if(iter > maxiter){
          stop("too many iterations")
        } 
      }
      
      if(verbose == TRUE) print(iter)
      U[i, ] <- t(u)
    }
  }
  
  
  # Symmetric estimation
  if(method == "symm"){
    crit <- 1
    iter <- 0
    
    while(crit > eps){
      
      # Symmetric orthogonalization
      Tmat <- matrix(0, k, p)
      for(i in 1:k){
        Tmat[i, ] <- computeTVec_C(U[i, ], X, nl_int, alpha)
      }
      U2 <- symmetricPower_C(tcrossprod(Tmat), -0.5)%*%Tmat
      
      # Criterion evaluation
      if(k == 1){
        J <- as.matrix(sign(U2%*%t(U)))
      }
      else{
        J <- diag(sign(diag(U2%*%t(U))))
      }
      crit <- sqrt(sum((U2 - J%*%U)^2))
      U <- U2
      
      iter <- iter + 1
      if(iter > maxiter){
        stop("too many iterations")
      } 
      
    }
    
    if(verbose == TRUE){
      print(iter)
    }
  }
  
  # Ordering the found components in decreasing order w.r.t. objective function value
  S <- X%*%t(U)
  obj <- computeObj_C(S, nl_int, alpha)
  obj_order <- order(obj, decreasing=TRUE)
  
  S_order <- S[, obj_order]
  colnames(S_order) <- sapply(1:k, function(i) paste("IC.", i, sep=""))

  res <- list(W = U[obj_order, ]%*%symmetricPower_C(cov_x, -0.5),
              S = S_order,
              D = c(obj[obj_order]),
              MU = MU)
  class(res) <- "bss"
  
  return(res)
}
