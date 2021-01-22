#' Compute Pairwise Distance for Symmetric Positive-Definite Matrices
#'
#' For a given 3-dimensional array where symmetric positive definite (SPD) matrices are stacked slice
#' by slice, it computes pairwise distance using various popular measures. Some of measures
#' are \emph{metric} as they suffice 3 conditions in mathematical context; nonnegative definiteness,
#' symmetry, and triangle inequalities. Other non-metric measures represent \emph{dissimilarities} between
#' two SPD objects.
#'
#' @param A a \eqn{(p\times p\times N)} 3d array of \eqn{N} SPD matrices.
#' @param method the type of distance measures to be used;
#' \describe{
#' \item{\code{"A"}}{(AIRM) Affine Invariant Riemannian Metric}
#' \item{\code{"B"}}{(Bhattacharyya) Bhattacharyya distance based on normal model}
#' \item{\code{"C"}}{(Cholesky) Cholesky difference in Frobenius norm}
#' \item{\code{"E"}}{(Euclidean) naive Frobenius norm as distance}
#' \item{\code{"H"}}{(Hellinger) Hellinger distance based on normal model}
#' \item{\code{"J"}}{(JBLD) Jensen-Bregman Log Determinant Distance}
#' \item{\code{"K"}}{(KLDM) symmetrized Kullback-Leibler Distance Measure}
#' \item{\code{"L"}}{(LERM) Log Euclidean Riemannian Metric}
#' \item{\code{"PS"}}{(Procrustes.SS) Procrustes Size and Shape measure}
#' \item{\code{"PF"}}{(Procrustes.Full) Procrustes analysis with scale}
#' \item{\code{"PE"}}{(PowerEuclidean) weighted eigenvalues by some exponent}
#' \item{\code{"RE"}}{(RootEuclidean) matrix square root}
#' }
#' @param power a non-zero number for PowerEuclidean distance.
#' @param as.dist a logical; \code{TRUE} to return a \code{dist} object, \code{FALSE} otherwise.
#'
#' @return an \eqn{(N\times N)} symmetric matrix of pairwise distances or a \code{dist} object via \code{as.dist} option.
#'
#' @examples
#' ## generate 50 SPD matrices of size (5-by-5)
#' samples = samplecovs(50,5)
#'
#' ## get pairwise distance for several methods
#' distA = CovDist(samples, method="A")
#' distB = CovDist(samples, method="B")
#' distC = CovDist(samples, method="C")
#'
#' ## dimension reduction using MDS
#' ssA = stats::cmdscale(distA, k=2)
#' ssB = stats::cmdscale(distB, k=2)
#' ssC = stats::cmdscale(distC, k=2)
#'
#' ## visualize
#' opar <- par(mfrow=c(1,3), pty="s")
#' plot(ssA ,main="project with AIRM", pch=19)
#' plot(ssB ,main="project with Bhattacharyya", pch=19)
#' plot(ssC ,main="project with Cholesky", pch=19)
#' par(opar)
#'
#' @references
#' \insertRef{arsigny_log-euclidean_2006}{CovTools}
#'
#' \insertRef{dryden_non-euclidean_2009}{CovTools}
#'
#' @export
CovDist <- function(A,method=c("A","B","C","E","H","J","K","L","PS","PF","PE","RE"),
                    power=1.0, as.dist=FALSE){
  ## PREPROCESSING
  ## 1) 3d array, 2) square, 3) symmetric, 4) sequentially check PDness
  if (length(dim(A))!=3){
    stop("* CovDist : 3d array should be used as input.")
  }
  if (dim(A)[1]!=dim(A)[2]){
    stop("* CovDist : A should be stack of square matrices.")
  }
  SymApply = apply(A, 3, isSymmetric, tol=sqrt(.Machine$double.eps))
  if (any(SymApply)==FALSE){
    SymIdx = which(!SymApply)
    if (length(SymIdx)==1){
      stopmessage = paste("* CovDist : slice number",SymApply,"is not Symmetric.")
    } else {
      stopmessage = paste("* CovDist : multiple of slices are not Symmetric.")
    }
    stop(stopmessage)
  }
  if (any(apply(A, 3, isSymmetric, tol=sqrt(.Machine$double.eps))==FALSE)){
    stop("* CovDist : each slice of A should be symmetric matrix.")
  }
  for (i in 1:dim(A)[3]){
    if (!check_pd(A[,,i])){
      stopmessage = paste(" CovDist : slice number",i,"is not Positive Definite.")
      stop(stopmessage)
    }
  }

  ## Main Iteration with Switch Argument
  if (missing(method)){
    mymethod = tolower("A")
  } else {
    mymethod = tolower(match.arg(method))
  }
  if (all(mymethod=="pe")){
    if (power==0){
      stop("* CovDist : 'power' should be a nonzero element. Suggests > 0.")
    }
    power = as.double(power)
  }

  # \item{\code{"AIRM"}}{Affine Invariant Riemannian Metric}
  # \item{\code{"Bhattacharyya"}}{Bhattacharyya distance based on normal model}
  # \item{\code{"Cholesky"}}{Cholesky difference in Frobenius norm}
  # \item{\code{"Euclidean"}}{naive Frobenius norm as distance}
  # \item{\code{"Hellinger"}}{Hellinger distance based on normal model}
  # \item{\code{"JBLD"}}{Jensen-Bregman Log Determinant Distance}
  # \item{\code{"KLDM"}}{symmetrized Kullback-Leibler Distance Measure}
  # \item{\code{"LERM"}}{Log Euclidean Riemannian Metric}
  # \item{\code{"Procrustes.SS"}}{Procrustes Size and Shape measure}
  # \item{\code{"Procrustes.Full"}}{Procrustes analysis with scale}
  # \item{\code{"PowerEuclidean"}}{weighted eigenvalues by some exponent}
  # \item{\code{"RootEuclidean"}}{matrix square root}

  outdist= switch(mymethod,
                  a = measure.AIRM.3d(A),
                  b = measure.Bhattacharyya.3d(A),
                  c = measure.Choleksy.3d(A),
                  e = measure.Euclidean.3d(A),
                  h = measure.Hellinger.3d(A),
                  j = measure.JBLD.3d(A),
                  k = measure.KLDM.3d(A),
                  l = measure.LERM.3d(A),
                  ps = measure.Procrustes.SS.3d(A),
                  pf = measure.Procrustes.Full.3d(A),
                  pe = measure.PowerEuclidean.3d(A,power),
                  re = measure.RootEuclidean.3d(A)
                  )
  if (as.dist){
    return(stats::as.dist(outdist))
  } else {
    return(outdist)
  }
}
