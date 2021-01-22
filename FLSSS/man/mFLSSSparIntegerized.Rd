\name{mFLSSSparIntegerized}
\alias{mFLSSSparIntegerized}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An advanced version of \code{mFLSSSpar()}
}
\description{
This function maps a real-value multidimensional Subset Sum problem to the integer domain with minimal precision loss. Those integers are further compressed in 64-bit buffers for dimension reduction and SWAR (SIMD within a register) that could lead to substantial acceleration.
}
\usage{
mFLSSSparIntegerized(
  maxCore = 7L,
  len,
  mV,
  mTarget,
  mME,
  solutionNeed = 1L,
  precisionLevel = integer(ncol(mV)),
  returnBeforeMining = FALSE,
  tlimit = 60,
  dl = ncol(mV),
  du = ncol(mV),
  useBiSrchInFB = FALSE,
  avgThreadLoad = 8L,
  verbose = TRUE
  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{maxCore}{
See \code{maxCore} in \code{mFLSSSpar()}.
}
  \item{len}{
See \code{len} in \code{mFLSSSpar()}.
}
  \item{mV}{
See \code{mV} in \code{mFLSSSpar()}.
}
  \item{mTarget}{
See \code{mTarget} in \code{mFLSSSpar()}.
}
  \item{mME}{
See \code{mME} in \code{mFLSSSpar()}.
}
  \item{solutionNeed}{
See \code{solutionNeed} in \code{mFLSSSpar()}.
}
\item{precisionLevel}{
An integer vector of size equal to the dimensionality of \code{mV}. This argument controls the precision of real-to-integer conversion.

If \code{precisionLevel[i] = 0}, \code{mV[,i]} is shifted, scaled and rounded to the nearest integers such that the maximum becomes no less than \code{nrow(mV) * 8}.

If \code{precisionLevel[i] > 0}, e.g. \code{precisionLevel[i] = 1000}, \code{mV[,i]} is shifted, scaled and rounded to the nearest integers such that the maximum becomes no less than 1000.

If \code{precisionLevel[i] = -1}, \code{mV[,i]} is shifted, scaled and rounded to the nearest integers such that ranks of elements stay the same.

The shift operator contributes no precision loss. It only lowers the number of bits used for storing integers.
}
\item{returnBeforeMining}{
A boolean value. If \code{TRUE}, function returns the integerized \code{mV}, \code{mTarget} and \code{mME}.
}
  \item{tlimit}{
See \code{tlimit} in \code{mFLSSSpar()}.
}
  \item{dl}{
See \code{dl} in \code{mFLSSSpar()}.
}
  \item{du}{
See \code{du} in \code{mFLSSSpar()}.
}
  \item{useBiSrchInFB}{
See \code{useBiSrchInFB} in \code{mFLSSSpar()}.
}
  \item{avgThreadLoad}{
See \code{avgThreadLoad} in \code{mFLSSSpar()}.
}
\item{verbose}{
If TRUE, prints mining progress.
}
}
\value{
A list of two.
\item{Value$solution}{is a list of solution index vectors.}
\item{Value$INT}{is a list of three.}
\item{Value$INT$mV}{is the integerized superset.}
\item{Value$INT$mTarget}{is the integerized subset sum.}
\item{Value$INT$mME}{is the integerized subset sum error threshold.}
\item{Value$INT$compressedDim}{is the dimensionality after integerization.}
}
\note{32-bit architecture unsupported.}
\examples{
if(.Machine$sizeof.pointer == 8L){
# =====================================================================================
# 64-bit architecture required.
# =====================================================================================
rm(list = ls()); gc()
subsetSize = 7L
supersetSize = 60L
dimension = 5L # dimensionality


# Create a supertset at random:
N = supersetSize * dimension
superset = matrix(1000 * (rnorm(N) ^ 3 + 2 * runif(N) ^ 2 + 3 * rgamma(N, 5, 1) + 4),
                  ncol = dimension)
rm(N)


# Plant a subset sum:
solution = sample(1L : supersetSize, subsetSize)
subsetSum = colSums(superset[solution, ])
subsetSumError = abs(subsetSum) * 0.01 # relative error within 1%
rm(solution)


# Mine subsets, dimensions fully bounded
system.time({rst = FLSSS::mFLSSSparIntegerized(
  maxCore = 2, len = subsetSize, mV = superset, mTarget = subsetSum,
  mME = subsetSumError, solutionNeed = 2, dl = ncol(superset),
  du = ncol(superset), tlimit = 2, useBiSrchInFB = FALSE, avgThreadLoad = 8L)})


# Compare the time cost of 'mFLSSSparIntegerized()' and 'mFLSSSpar()'. The
# speed advantage of 'mFLSSSparIntegerized()' may not be pronounced for toy
# examples.
system.time(FLSSS::mFLSSSpar(
  maxCore = 2, len = subsetSize, mV = superset, mTarget = subsetSum,
  mME = subsetSumError, solutionNeed = 2, dl = ncol(superset),
  du = ncol(superset), tlimit = 2, useBiSrchInFB = FALSE, avgThreadLoad = 8L))


# Verify:
cat("Number of solutions = ", length(rst$solution), "\n")
if(length(rst$solution) > 0)
{
  cat("Solutions unique: ")
  cat(length(unique(lapply(rst$solution, function(x)
    sort(x)))) == length(rst$solution), "\n")


  cat("Solutions correct regarding integerized data: ")
  cat(all(unlist(lapply(rst$solution, function(x)
    abs(colSums(rst$INT$mV[x, ]) - rst$INT$mTarget) <= rst$INT$mME))), "\n")


  cat("Solutions correct regarding original data: ")
  boolean = all(unlist(lapply(rst$solution, function(x)
    abs(colSums(superset[x, ]) - subsetSum) <= subsetSumError)))
  cat(boolean, "\n")
  if(!boolean)
  {
    cat("The given error threshold relative to subset sum:\n")
    givenRelaErr = round(abs(subsetSumError / subsetSum), 5)
    cat(givenRelaErr, "\n")


    cat("Solution subset sum relative error:\n")
    tmp = lapply(rst$solution, function(x)
    {
      err = round(abs(colSums(superset[x, ]) / subsetSum -1), 5)
      for(i in 1L : length(err))
      {
        if(givenRelaErr[i] < err[i]) message(paste0(err[i], " "), appendLF = FALSE)
        else cat(err[i], "")
      }
      cat("\n")
    })
    cat("Integerization caused the errors. Future versions of")
    cat("'mFLSSSparIntegerized()' would have a parameter of precision level.\n")
  }
} else
{
  cat("No solutions exist or time ended too soon.\n")
}


# Mine subsets, the first 3 dimensions lower bounded,
# the last 4 dimension upper bounded
rst = FLSSS::mFLSSSparIntegerized(
  maxCore = 2, len = subsetSize, mV = superset, mTarget = subsetSum,
  mME = subsetSumError, solutionNeed = 2, dl = 3L, du = 4L, tlimit = 2,
  useBiSrchInFB = FALSE, avgThreadLoad = 8L)


# Verify:
cat("Number of solutions = ", length(rst$solution), "\n")
if(length(rst$solution) > 0)
{
  cat("Solutions unique: ")
  cat(length(unique(lapply(rst$solution, function(x)
    sort(x)))) == length(rst$solution), "\n")


  cat("Solutions correct regarding integerized data: ")
  cat(all(unlist(lapply(rst$solution, function(x)
  {
    lowerBoundedDim = 1L : 3L
    lowerBounded = all(colSums(rst$INT$mV[x, lowerBoundedDim]) >=
      rst$INT$mTarget[lowerBoundedDim] - rst$INT$mME[lowerBoundedDim])


    upperBoundedDim = (ncol(rst$INT$mV) - 3L) : ncol(rst$INT$mV)
    upperBounded = all(colSums(rst$INT$mV[x, upperBoundedDim]) <=
      rst$INT$mTarget[upperBoundedDim] + rst$INT$mME[upperBoundedDim])


    lowerBounded & upperBounded
  }))), "\n")
} else
{
  cat("No solutions exist or timer ended too soon.\n")
}
# =====================================================================================
# =====================================================================================
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
