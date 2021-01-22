SIL <- function (Xca, U, distance = FALSE)
  {
    if (missing(Xca))
      stop("The data set must be given")
    if (is.null(Xca))
      stop("The data set in Xca is empty")
    n=nrow(Xca)
    Xca=as.matrix(Xca)
    if (any(is.na(Xca)))
      stop("The data set in Xca must not contain NA values")
    if (!is.numeric(Xca))
      stop("The data set in Xca is not a numeric data.frame or matrix")
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
      stop("There is only k=1 cluster: the silhouette index is not computed")
    if (nrow(U)!=nrow(Xca))
      stop("The numbers of rows of U and Xca must be the same")

    out = silhouette(X = Xca,U = U,p = ncol(Xca),k = ncol(U),n = n, distance = distance)
    sil.obj = as.vector(out$sil.obj)
    names(sil.obj)=rn
    sil=out$sil
    out=list(sil.obj = sil.obj,
             sil = sil)
    return(out)
  }
