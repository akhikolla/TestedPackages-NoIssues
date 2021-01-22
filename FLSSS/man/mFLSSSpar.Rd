\name{mFLSSSpar}
\alias{mFLSSSpar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multithreaded multidimensional Subset Sum given error thresholds
}
\description{
The multidimensional version of \code{FLSSS()}. See \code{decomposeMflsss()} for the multi-process version.
}
\usage{
mFLSSSpar(
  maxCore = 7L,
  len,
  mV,
  mTarget,
  mME,
  solutionNeed = 1L,
  tlimit = 60,
  dl = ncol(mV),
  du = ncol(mV),
  useBiSrchInFB = FALSE,
  avgThreadLoad = 8L
  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{maxCore}{
Maximal threads to invoke. Ideally \code{maxCore} should not surpass the total logical processors on machine.
}
  \item{len}{
An integer as the subset size. See \code{len} in \code{FLSSS()}.
}
  \item{mV}{
A data frame or a matrix as the multidimensional set, columns as dimensions.
}
  \item{mTarget}{
A numeric vector of size \code{ncol(mV)} as the subset sum.
}
  \item{mME}{
A numeric vector of size \code{ncol(mV)} as the subset sum error thresholds.
}
  \item{solutionNeed}{
See \code{solutionNeed} in \code{FLSSS()}.
}
  \item{tlimit}{
See \code{tlimit} in \code{FLSSS()}.
}
  \item{dl}{
An integer no greater than \code{ncol(mV)}. Let \code{sol} be the index vector of a solution. Let \code{dls <- 1L : dl}. The following is true:

\code{colSums(mV[sol, dls]) >= mTarget[dls] - mME[dls]}.
}
  \item{du}{
An integer no greater than \code{ncol(mV)}. Let \code{sol} be the index vector of a solution. Let \code{dus <- (ncol(mV) - du + 1) : ncol(mV)}. The following is true:

\code{colSums(mV[sol, dus]) <= mTarget[dus] + mME[dus]}.
}
  \item{useBiSrchInFB}{
See \code{useBiSrchInFB} in \code{FLSSS()}.
}
\item{avgThreadLoad}{
If \code{mV} is comonotonic, \code{mFLSSSpar()} warms up with a breadth-first search and then spawns at least \code{B} branches for parallelization. \code{B} equals the first power-of-two integer no less than \code{avgThreadLoad} \code{* maxCore}.
}
}
%\details{
%}
\value{
A list of index vectors.
}

\examples{
rm(list = ls()); gc()
subsetSize = 7L
supersetSize = 60L
dimension = 5L # dimensionality


# Create a supertset at random:
N = supersetSize * dimension
superset = matrix(1000 * (rnorm(N) ^ 3 + 2 * runif(N) ^ 2 +
  3 * rgamma(N, 5, 1) + 4), ncol = dimension)
rm(N)


# Plant a subset sum:
solution = sample(1L : supersetSize, subsetSize)
subsetSum = colSums(superset[solution, ])
subsetSumError = abs(subsetSum) * 0.01 # relative error within 1%
rm(solution)


# Mine subsets, dimensions fully bounded
rst = FLSSS::mFLSSSpar(maxCore = 2, len = subsetSize, mV = superset,
                       mTarget = subsetSum, mME = subsetSumError,
                       solutionNeed = 2, dl = ncol(superset), du = ncol(superset),
                       tlimit = 2, useBiSrchInFB = FALSE, avgThreadLoad = 8L)


# Verify:
cat("Number of solutions = ", length(rst), "\n")
if(length(rst) > 0)
{
  cat("Solutions unique: ")
  cat(length(unique(lapply(rst, function(x) sort(x)))) == length(rst), "\n")
  cat("Solutions correct: ")
  cat(all(unlist(lapply(rst, function(x)
    abs(colSums(superset[x, ]) - subsetSum) <= subsetSumError))), "\n")
} else
{
  cat("No solutions exist or timer ended too soon.\n")
}




# Mine subsets, the first 3 dimensions lower bounded,
# the last 4 dimension upper bounded
rst = FLSSS::mFLSSSpar(maxCore = 2, len = subsetSize, mV = superset,
                       mTarget = subsetSum, mME = subsetSumError,
                       solutionNeed = 2, dl = 3L, du = 4L,
                       tlimit = 2, useBiSrchInFB = FALSE, avgThreadLoad = 8L)


# Verify:
cat("Number of solutions = ", length(rst), "\n")
if(length(rst) > 0)
{
  cat("Solutions unique: ")
  cat(length(unique(lapply(rst, function(x) sort(x)))) == length(rst), "\n")
  cat("Solutions correct: ")
  cat(all(unlist(lapply(rst, function(x)
  {
    lowerBoundedDim = 1L : 3L
    lowerBounded = all(colSums(superset[x, lowerBoundedDim]) >=
      subsetSum[lowerBoundedDim] - subsetSumError[lowerBoundedDim])


    upperBoundedDim = (ncol(superset) - 3L) : ncol(superset)
    upperBounded = all(colSums(superset[x, upperBoundedDim]) <=
      subsetSum[upperBoundedDim] + subsetSumError[upperBoundedDim])


    lowerBounded & upperBounded
  }))), "\n")
} else
{
  cat("No solutions exist or timer ended too soon.\n")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
