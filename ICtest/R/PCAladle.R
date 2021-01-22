# Function to return M matrix as in paper
# input data matrix X

MPCA <- function(X)
{
  n <- nrow(X)
  MEAN <- colMeans(X)
  X.C <- sweep(X,2,MEAN,"-")
  COV <- crossprod(X.C)/(n-1)
}



# function to compute the qunatities for one bootstrap sample
# input
# X - bootstrapsample
# EVdata - eigenvectors for the M matrix for the observed data
# r - the value for the number of components to consider
PCAladleboot <- function(X, EVdata, r)
{
  Mboot <- MPCA(X)
  EVboot <- eigen(Mboot, symmetric = TRUE)$vectors
  fi(EVboot, EVdata, r)
}


# main function for ladle for FOBI
# input
# X - data matrix
# n.boot - number of bootstrapping samples to be taken
# cutDenom - rule how to compute "r"
PCAladle <- function(X, n.boot=200, ncomp=ifelse(ncol(X) > 10, floor(ncol(X)/log(ncol(X))), ncol(X)-1))
{
  data.name <-  deparse(substitute(X))
  method <- "PCA"                        
  
  p <- ncol(X)
  n <- nrow(X)
  
  MEAN <- colMeans(X)
  X.C <- sweep(X, 2, MEAN, "-")
  Mdata <- crossprod(X.C)/(n-1)
  
  EV.Mdata <- eigen(Mdata, symmetric = TRUE)
  EVdata <- EV.Mdata$vectors
  
  fis <- replicate(n.boot, PCAladleboot(X[sample(1:n, n, replace=TRUE),], EVdata, ncomp))
  
  fn0 <- c(0,rowMeans(fis))
  fn <- fn0 / (1+sum(fn0))
  phin <- EV.Mdata$values[1:(ncomp+1)] / (1+sum(EV.Mdata$values[1:(ncomp+1)]))
  gn <- fn + phin
  est.k <- which.min(gn)-1
  
  W <- EVdata
  #S <- tcrossprod(X.C, W)
  S <- X.C %*% W
  colnames(S) <- paste0("PC.", 1:p)
  
  
  RES <- list(method = method, k = est.k, fn = fn, phin = phin, gn = gn, 
              lambda = EV.Mdata$values[1:(ncomp+1)], W = W, S = S, MU = MEAN,
              data.name = data.name)
  class(RES) <- "ladle"
  RES
}
