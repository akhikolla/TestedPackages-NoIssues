SIL.F <- function (Xca, U, alpha, distance = FALSE)
  {
    if (missing(Xca))
      stop("The data set must be given")
    if (is.null(Xca))
      stop("The data set Xca is empty")
    n=nrow(Xca)
    Xca=as.matrix(Xca)
    if (any(is.na(Xca)))
      stop("The data set Xca must not contain NA values")
    if (!is.numeric(Xca))
      stop("The data set Xca is not a numeric data.frame or matrix")
    if (is.null(rownames(Xca)))
      rn=paste("Obj",1:n,sep=" ")
    else
      rn=rownames(Xca)
    if (missing(U))
      stop("The membership degree matrix U must be given")
    if (is.null(U))
      stop("The membership degree matrix U is empty")
    U=as.matrix(U)
    if (any(is.na(U)))
      stop("The membership degree matrix U must not contain NA values")
    if (!is.numeric(U))
      stop("The membership degree matrix U is not numeric")
    if (ncol(U)==1)
      stop("There is only k=1 cluster: the fuzzy silhouette index is not computed")
    if (nrow(U)!=nrow(Xca))
      stop("The numbers of rows of U and Xca must be the same")
    if (missing(alpha))
    {
      alpha=1
      cat("The default value alpha=1 has been set ",fill=TRUE)
    }
    if (!is.numeric(alpha))
    {
      alpha=1
      cat("The weighting coefficient alpha is not numeric: the default value alpha=1 will be used ",fill=TRUE)
    }
    if (alpha<0)
    {
      alpha=1
      cat("The number of clusters k must be non negative: the value alpha=1 will be used ",fill=TRUE)
    }


    sil.f = silhouetteFuzzy(X = Xca,U = U,alpha = alpha,p = ncol(Xca),k = ncol(U),n = n, distance = distance)
    return(sil.f)
  }
