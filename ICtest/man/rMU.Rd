\name{rMU}
\alias{rMU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Greek Letter mu Shaped Bivariate Data Generation
}
\description{
A function to generate bivariate data with the scatterplot resembling the greek letter mu.
}
\usage{
rMU(n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{the sample size.}
}


\value{
A \code{n} times \code{2} matrix
}

\author{
Klaus Nordhausen, Joni Virta
}

\examples{
x <- rMU(1000)

plot(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}% use one of  RShowDoc("KEYWORDS")

