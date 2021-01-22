\name{decomposeMflsss}
\alias{decomposeMflsss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
mFLSSS decomposition
}
\description{
Decompose an \code{mFLSSS} instance into sub-problems for distributed computing.
}
\usage{
decomposeMflsss(
  len,
  mV,
  mTarget,
  mME,
  solutionNeed = 1L,
  dl = ncol(mV),
  du = ncol(mV),
  useBiSrchInFB = FALSE,
  approxNinstance = 50000L
  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{len}{
See the same argument in \code{mFLSSSpar()}.
}
  \item{mV}{
See the same argument in \code{mFLSSSpar()}.
}
  \item{mTarget}{
See the same argument in \code{mFLSSSpar()}.
}
  \item{mME}{
See the same argument in \code{mFLSSSpar()}.
}
  \item{solutionNeed}{
See the same argument in \code{mFLSSSpar()}.
}
  \item{dl}{
See the same argument in \code{mFLSSSpar()}.
}
  \item{du}{
See the same argument in \code{mFLSSSpar()}.
}
  \item{useBiSrchInFB}{
See the same argument in \code{mFLSSSpar()}.
}
  \item{approxNinstance}{
Approximately how many instances should the problem be decomposed into.
}
}
\details{
This function and \code{mFLSSSobjRun()} constitute a multi-process counterpart of \code{mFLSSSpar()}. It decomposes a multidimensional subset sum problem into numerous independent instances that can be submitted to any computing resource of CPU threads, where \code{mFLSSSobjRun()} receives and solves the instance.

For example, if 1000 threads are available, either on a computing cluster integrated by LAN or on 1000 PCs connected by the Internet, one could (i) decompose the problem of interest into 100000 instances using \code{decomposeMflsss()}, (ii) transmit each instance and an R script that calls \code{mFLSSSobjRun()} on the instance to any available thread, (iii) collect the results from all threads. A good practice is to decompose the problem of interest into much more instances than the threads, provided that an automatic queueing system exists, so there would be less chance of having idling threads during computation --- if the number of instances equals the number of threads, some threads may finish much earlier than others due to the heterogeneous nature of the instances, thus the waste.


The pair \code{decomposeMflsss()} and \code{mFLSSSobjRun()} is designed for exploiting distrbuted resource to solve large and hard multidimensional subset sum instances.
}
\value{
A list of two:

\code{$mflsssObjects}: a list. Each element is an \code{mFLSSS} object that would be supplied to \code{mFLSSSobjRun()}.

\code{$solutionsFound}: solutions found during decomposition.
}

\examples{
N = 30L # Superset size.
len = 6L # Subset size.
dimen = 5L # Dimension.
set.seed(8120)
v = matrix(runif(N * dimen) * 1000, nrow = N) # Superset.
sol = sample(N, len)
target = colSums(v[sol, ]) # Target sum.
ME = target * 0.03 # Error threshold.
approxNinstance = 1000


validate = function(len, v, target, ME, result)
{
  all(unlist(lapply(result, function(x)
    all(abs(colSums(v[x, ]) - target) <= ME))))
}


decompedFlsss = FLSSS::decomposeMflsss(
  len = len, mV = v, mTarget = target, mME = ME, solutionNeed = 1e6,
  approxNinstance = approxNinstance)


str(decompedFlsss$solutionsFound) # See if the agent already found
# some solutions and validate them.
if(length(decompedFlsss$solutionsFound) > 0)
  print(validate(len, v, target, ME, decompedFlsss$solutionsFound))


length(decompedFlsss$mflsssObjects) # Number of independent small jobs.
someOtherSolutions = FLSSS::mFLSSSobjRun(
  decompedFlsss$mflsssObjects[[620]], tlimit = 3, solutionNeed = 1e6)


if(length(someOtherSolutions) > 0) # Validate solutions.
{
  print(someOtherSolutions)
  print(validate(len, v, target, ME, someOtherSolutions))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
