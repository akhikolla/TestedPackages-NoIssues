\name{mFLSSSparImposeBounds}
\alias{mFLSSSparImposeBounds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multithreaded multidimensional Subset Sum in bounded solution space given error thresholds}
\description{
For comparison, function \code{mFLSSSpar()} puts no bounds on the solution space so it sorts \code{mV} internally in a special order for mining accerlation.
}
\usage{
mFLSSSparImposeBounds(
  maxCore = 7L,
  len,
  mV,
  mTarget,
  mME,
  LB = 1L : len,
  UB = (nrow(mV) - len + 1L) : nrow(mV),
  solutionNeed = 1L,
  tlimit = 60,
  dl = ncol(mV),
  du = ncol(mV),
  targetsOrder = NULL,
  useBiSrchInFB = FALSE,
  avgThreadLoad = 8L
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
  \item{LB}{
See \code{LB} in \code{FLSSS()}.
}
  \item{UB}{
See \code{UB} in \code{FLSSS()}.
}
  \item{solutionNeed}{
See \code{solutionNeed} in \code{mFLSSSpar()}.
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
  \item{targetsOrder}{
This argument is mainly for research and unrecommended for use. Depending on the structure of \code{mV}, \code{mFLSSSpar()} or \cr\code{mFLSSSparImposeBounds()} would break the mining task into a collection of no more than \code{len * (nrow(mV) - len) + 1} independent subtasks. Threads work on these subtasks, sequentially coordinated by an atomic counter [1]. Different subtasks have different probabilities of yielding a qualified subset, thus the order of subtasks matters to the mining speed. \code{targetsOrder} is an index vector of size \code{len} \code{*} \code{(nrow(mV)} \code{-} \code{len)} \code{+} \code{1} for ordering the subtasks. \code{targetsOrder <- NULL} makes a special order, and is implicitly the choice in \code{mFLSSSpar()}. This order is empirically optimal based on simulations. See the \href{https://arxiv.org/abs/1612.04484}{package documentation} for details.
}
  \item{useBiSrchInFB}{
See \code{useBiSrchInFB} in \code{mFLSSSpar()}.
}
  \item{avgThreadLoad}{
See \code{avgThreadLoad} in \code{mFLSSSpar()}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list of index vectors.
}
\references{
[1] Atomic template class in \href{https://software.intel.com/content/www/us/en/develop/tools/threading-building-blocks.html}{Intel TBB}. An atomic counter is used to coordinate heterogeneous subtasks to avoid idle threads. The atomic operation overhead is negligible compared to the time cost of the lightest subtask.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
rm(list = ls()); gc()
subsetSize = 7L
supersetSize = 60L
dimension = 5L # dimensionality

\donttest{
# Create a supertset at random:
N = supersetSize * dimension
superset = matrix(1000 * (rnorm(N) ^ 3 + 2 * runif(N) ^ 2 +
                  3 * rgamma(N, 5, 1) + 4), ncol = dimension)
rm(N)


# Make up the lower and upper bounds for the solution space:
tmp = sort(sample(1L : supersetSize, subsetSize))
tmp2 = sort(sample(1L : supersetSize, subsetSize))
lowerBounds = pmin(tmp, tmp2)
upperBounds = pmax(tmp, tmp2)
rm(tmp, tmp2)


# Exclude elements not covered by 'lowerBounds' and 'upperBounds':
remainIndex = unique(unlist(apply(cbind(lowerBounds, upperBounds), 1,
                                  function(x) x[1] : x[2])))
lowerBounds = match(lowerBounds, remainIndex)
upperBounds = match(upperBounds, remainIndex)
superset = superset[remainIndex, ]


# Plant a subset sum:
solution = apply(rbind(lowerBounds, upperBounds), 2, function(x)
  sample(x[1] : x[2], 1))
subsetSum = colSums(superset[solution, ])
subsetSumError = abs(subsetSum) * 0.01 # relative error within 1%
rm(solution)


rst = FLSSS::mFLSSSparImposeBounds(
  maxCore = 2L, len = subsetSize, mV = superset, mTarget = subsetSum,
  mME = subsetSumError, LB = lowerBounds, UB = upperBounds,
  solutionNeed = 1, tlimit = 2, dl = ncol(superset), du = ncol(superset),
  targetsOrder = NULL, useBiSrchInFB = FALSE, avgThreadLoad = 8L)


# Verify:
cat("Number of solutions = ", length(rst), "\n")
if(length(rst) > 0)
{
  cat("Solutions unique: ")
  cat(length(unique(lapply(rst, function(x) sort(x)))) == length(rst), "\n")
  cat("Solution in bounded space: ")
  cat(all(unlist(lapply(rst, function(x)
    sort(x) <= upperBounds & sort(x) >= lowerBounds))), "\n")
  cat("Solutions correct: ")
  cat(all(unlist(lapply(rst, function(x)
    abs(colSums(superset[x, ]) - subsetSum) <= subsetSumError))), "\n")
} else
{
  cat("No solutions exist or timer ended too soon.\n")
}
}
}
