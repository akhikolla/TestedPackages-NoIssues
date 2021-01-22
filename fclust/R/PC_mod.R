PC <- function (U)
  {
    if (missing(U))
      stop("The membership degree matrix U must be given")
    if (is.null(U))
      stop("The membership degree matrix U is empty")
    U=as.matrix(U)
    if (any(is.na(U)))
      stop("The membership degree matrix U must not contain NA values")
    if (!is.numeric(U))
      stop("The membership degree matrix U is not numeric")
    part.coeff=partCoef(U = U,n = nrow(U))
    return(part.coeff)
  }
