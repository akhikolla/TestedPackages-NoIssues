\name{FLSSSmultiset}
\alias{FLSSSmultiset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multi-Subset Sum given error threshold
}
\description{
Find a subet of a given size for each of multiple supersets such that all the subsets sum in a given range.
}
\usage{
FLSSSmultiset(
  len,
  buckets,
  target,
  ME,
  solutionNeed = 1L,
  tlimit = 60,
  useBiSrchInFB = FALSE,
  useFloat = FALSE
  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{len}{
A positive integer vector as the subset sizes for the supersets.
}
  \item{buckets}{
A list of the supersets. \code{buckets[[i]]} is an unsorted numeric vector of size \code{len[i]}.
}
  \item{target}{
See \code{target} in \code{FLSSS()}.
}
  \item{ME}{
See \code{ME} in \code{FLSSS()}.
}
  \item{solutionNeed}{
See \code{solutionNeed} in \code{FLSSS()}.
}
  \item{tlimit}{
See \code{tlimit} in \code{FLSSS()}.
}
  \item{useBiSrchInFB}{
See \code{useBiSrchInFB} in \code{FLSSS()}.
}
  \item{useFloat}{
See \code{useFloat} in \code{FLSSS()}.
}
}
\value{
A list of solutions. Each solution is a list of index vectors. Assume \code{X} is a solution. \code{X[[i]]} indexes the subset of superset \code{buckets[[i]]}.
}
\examples{
# rm(list = ls()); gc()
Nsupersets = 30L
supersetSizes = sample(5L : 20L, Nsupersets, replace = TRUE)
subsetSizes = sapply(supersetSizes, function(x) sample(1L : x, 1))


# Create supersets at random:
supersets = lapply(supersetSizes, function(n)
{
  1000 * (rnorm(n) ^ 3 + 2 * runif(n) ^ 2 + 3 * rgamma(n, 5, 1) + 4)
})
str(supersets) # see the structure


# Give a subset sum
solution = mapply(function(n, l) sample(1L : n, l), supersetSizes, subsetSizes)
str(solution) # See structure
subsetsSum = sum(mapply(function(x, s) sum(x[s]), supersets, solution, SIMPLIFY = TRUE))
subsetsSumError = abs(subsetsSum) * 1e-7 # relative error within 0.00001%
rm(solution)


# Mine subsets:
rst = FLSSS::FLSSSmultiset(len = subsetSizes, buckets = supersets, target = subsetsSum,
                           ME = subsetsSumError, solutionNeed = 3, tlimit = 4)
cat("Number of solutions =", length(rst), "\n")


# Verify:
ver = all(unlist(lapply(rst, function(sol)
{
  S = sum(unlist(mapply(function(x, y) sum(x[y]), supersets, sol)))
  abs(S - subsetsSum) <= subsetsSumError
})))
cat("All subsets are qualified:", ver)
}














