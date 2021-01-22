\name{mFLSSSobjRun}
\alias{mFLSSSobjRun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run an \code{mFLSSS} instance
}
\description{
Run a multidimensional subset sum instance decomposed from \code{decomposeMflsss()}.
}
\usage{
mFLSSSobjRun(
  mflsssObj,
  solutionNeed = 1,
  tlimit = 60
  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mflsssObj}{
An \code{mFLSSS} object.
}
  \item{solutionNeed}{
See the same argument in \code{mFLSSSpar()}.
}
  \item{tlimit}{
See the same argument in \code{mFLSSSpar()}.
}
}
\details{
See the details about \code{decomposeMflsss()}.
}
\value{
See the value of \code{mFLSSSpar()}.
}

\examples{
# See the example for decomposeMflsss().
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
