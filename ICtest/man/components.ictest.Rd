\name{components}
\alias{components}
\alias{components.ictest}
\alias{components.ladle}

\title{
Generic Components Extraction Function
}
\description{
Function to extract components from an object. If the object is of class ictest or ladle the user can choose if all components
are extracted or only those which were interesting under the null hypothesis.
}
\usage{
components(x, ...)
\method{components}{ictest}(x, which = "all", ...)
\method{components}{ladle}(x, which = "all", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object which has a components method, like for example an ictest object.}  
  \item{which}{for an object of class ictest. If \code{"all"}, then all components S in the ictest object are extracted. If \code{"k"}, then only the first k components are extracted,
   where the value of \code{k} is taken from the ictest object. This is only meaningful if \code{k} was at least 1.}
  \item{...}{arguments passed on to other methods.}
}
\value{
a matrix with the components.
}

\author{
Klaus Nordhausen
}



\examples{
n <- 200
X <- cbind(rnorm(n, sd = 2), rnorm(n, sd = 1.5), rnorm(n), rnorm(n), rnorm(n))

TestCov <- PCAasymp(X, k = 2)
head(components(TestCov))
head(components(TestCov, which = "k"))
}
   

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods	}
