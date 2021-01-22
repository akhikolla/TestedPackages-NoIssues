#' Estimate Mean Covariance Matrix
#'
#' For a given 3-dimensional array where symmetric positive definite (SPD) matrices are stacked slice
#' by slice, it estimates Frechet mean on an open cone of SPD matrices under corresponding metric/distance
#' measure.
#'
#' @param A a \eqn{(p\times p\times N)} 3d array of \eqn{N} SPD matrices.
#' @param method the type of distance measures to be used;
#' \describe{
#' \item{\code{"A"}}{(AIRM) Affine Invariant Riemannian Metric}
#' \item{\code{"C"}}{(Cholesky) Cholesky difference in Frobenius norm}
#' \item{\code{"E"}}{(Euclidean) naive Frobenius norm as distance}
#' \item{\code{"L"}}{(LERM) Log Euclidean Riemannian Metric}
#' \item{\code{"PS"}}{(Procrustes.SS) Procrustes Size and Shape measure}
#' \item{\code{"PF"}}{(Procrustes.Full) Procrustes analysis with scale}
#' \item{\code{"PE"}}{(PowerEuclidean) weighted eigenvalues by some exponent}
#' \item{\code{"RE"}}{(RootEuclidean) matrix square root}
#' }
#' @param power a non-zero number for PowerEuclidean distance.
#' @return a \eqn{(p\times p)} mean covariance matrix estimated.
#'
#' @examples
#' \dontrun{
#' ## generate 50 sample covariances of size (10-by-10).
#' pdim    = 10
#' samples = samplecovs(50,pdim)
#'
#' ## compute means of first 50 sample covariances from data under Normal(0,Identity).
#' mA = CovMean(samples, method="A")
#' mC = CovMean(samples, method="C")
#' mE = CovMean(samples, method="E")
#' mL = CovMean(samples, method="L")
#' mPS = CovMean(samples, method="PS")
#' mPF = CovMean(samples, method="PF")
#' mPE = CovMean(samples, method="PE")
#' mRE = CovMean(samples, method="RE")
#'
#' #' ## visualize
#' opar <- par(mfrow=c(3,3), pty="s")
#' image(diag(pdim)[,pdim:1], main="true covariance")
#' image(mA[,pdim:1], main="AIRM")
#' image(mC[,pdim:1], main="Cholesky")
#' image(mE[,pdim:1], main="Euclidean")
#' image(mL[,pdim:1], main="LERM")
#' image(mPS[,pdim:1], main="Procrustes.SS")
#' image(mPF[,pdim:1], main="Procrustes.Full")
#' image(mPE[,pdim:1], main="PowerEuclidean")
#' image(mRE[,pdim:1], main="RootEuclidean")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{dryden_non-euclidean_2009}{CovTools}
#'
#' @export
CovMean <- function(A,method=c("A","C","E","L","PS","PF","PE","RE"),power=1.0){
  ## PREPROCESSING
  ## 1) 3d array, 2) square, 3) symmetric, 4) sequentially check PDness
  if (length(dim(A))!=3){
    stop("* CovMean : 3d array should be used as input.")
  }
  if (dim(A)[1]!=dim(A)[2]){
    stop("* CovMean : A should be stack of square matrices.")
  }
  SymApply = apply(A, 3, isSymmetric, tol=sqrt(.Machine$double.eps))
  if (any(SymApply)==FALSE){
    SymIdx = which(!SymApply)
    if (length(SymIdx)==1){
      stop(paste("* CovMean : slice number",SymApply,"is not Symmetric."))
    } else {
      stop(paste("* CovMean : multiple of slices are not Symmetric."))
    }
  }
  if (any(apply(A, 3, isSymmetric, tol=sqrt(.Machine$double.eps))==FALSE)){
    stop("* CovMean : each slice of A should be symmetric matrix.")
  }
  for (i in 1:dim(A)[3]){
    if (!check_pd(A[,,i])){
      stop(paste(" CovMean : slice number",i,"is not Positive Definite."))
    }
  }

  ## Main Computation with Switch Argument
  if (missing(method)){
    mymethod = "a"
  } else {
    mymethod = tolower(match.arg(method))
  }
  if (all(mymethod=="pe")){
    if (power==0){
      stop("* CovMean : 'power' should be a nonzero element. Suggests > 0.")
    }
    power = as.double(power)
  }
  outmean= switch(mymethod,
                  a = meancov.Riemannian(A),
                  c = meancov.Cholesky(A),
                  e = meancov.Euclidean(A),
                  l = meancov.LogEuclidean(A),
                  ps = meancov.Procrustes.SS(A),
                  pf = meancov.Procrustes.Full(A),
                  pe = meancov.PowerEuclidean(A,power),
                  re = meancov.RootEuclidean(A))
  return(outmean)
}
