# Function to return M matrix as in paper
# input data matrix X

MSIR <- function(X,y, h, ...)
  {
    n <- nrow(X)
    MEAN <- colMeans(X)
    X.C <- sweep(X,2,MEAN,"-")
    COV <- crossprod(X.C)/(n-1)

    EVD.COV <- eigen(COV, symmetric = TRUE)

    COV.inv.sqrt <- EVD.COV$vectors %*% tcrossprod(diag((1/EVD.COV$values)^0.5),
                                                   EVD.COV$vectors)
    Z <- tcrossprod(X.C, COV.inv.sqrt)
    M <-  covSIR(Z, y=y, h = h, ...)
    M
  }



# function to compute the qunatities for one bootstrap sample
# input
# X - bootstrapsample
# EVdata - eigenvectors for the M matrix for the observed data
# r - the value for the number of components to consider
SIRladleboot <- function(ind, X, y,h, EVdata, r,...)
  {
  Mboot <- MSIR(X[ind, ], y[ind], h,...)
  EVboot <- eigen(Mboot, symmetric = TRUE)$vectors
  fi(EVboot, EVdata, r)
  }


# main function for ladle for SIR
# input
# X - data matrix
# y - response
# h - number of slices
# n.boot - number of bootstrapping samples to be taken
# cutDenom - rule how to compute "r"
# ... - passed on to covSIR
SIRladle <- function(X,y, h=10, n.boot=200, ncomp=ifelse(ncol(X) > 10, floor(ncol(X)/log(ncol(X))), ncol(X)-1),...)
  {
  data.name <-  deparse(substitute(X))
  method <- "SIR"

  p <- ncol(X)
  n <- nrow(X)

  MEAN <- colMeans(X)
  X.C <- sweep(X, 2, MEAN, "-")
  COV <- crossprod(X.C)/(n-1)

  EVD.COV <- eigen(COV, symmetric = TRUE)

  COV.inv.sqrt <- EVD.COV$vectors %*% tcrossprod(diag((1/EVD.COV$values)^0.5),
                                                 EVD.COV$vectors)
  Z <- tcrossprod(X.C, COV.inv.sqrt)
  
  Mdata <- covSIR(Z,y=y,h=h, ...)

  EV.Mdata <- eigen(Mdata, symmetric = TRUE)
  EVdata <- EV.Mdata$vectors

  fis <- replicate(n.boot, SIRladleboot(ind = sample(1:n, n, replace=TRUE), X=X, y=y, h=h, EVdata, ncomp, ...))

  fn0 <- c(0,rowMeans(fis))
  fn <- fn0 / (1+sum(fn0))
  phin <- EV.Mdata$values[1:(ncomp+1)] / (1+sum(EV.Mdata$values[1:(ncomp+1)]))
  gn <- fn + phin
  est.k <- which.min(gn)-1

  W <- crossprod(EVdata,  COV.inv.sqrt)
  S <- tcrossprod(X.C, W)
  colnames(S) <- paste0("SIC.", 1:p)


  RES <- list(method = method, k = est.k, fn = fn, phin = phin, gn = gn,
              lambda = EV.Mdata$values[1:(ncomp+1)], W = W, S = S, MU = MEAN,
              data.name = data.name)
  class(RES) <- "ladle"
  RES
}
