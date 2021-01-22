\name{rorth}
\alias{rorth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Random Orthogonal Matrix Creation Uniform WRT the Haar Measure.
}
\description{
A function to create a random orthogonal matrix uniformly distributed with respect to the Haar measure.
}
\usage{
rorth(k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{the desired numer of columns (=rows) of the orthogonal matrix.}
}

\details{
The function fills a \code{k}x\code{k} matrix with N(0,1) random variables and perfroms then a QR decompoistion using \code{\link{qr}}. If the diagonal elements of R are all positive the resulting orthogonal matrix Q is uniform distributed wrt to the Haar measure. Note that the function currently does not check if
all diagonal measurements are indeed positive (however this will happen with probability 1 in theory). 
}

\value{
An orthogonal \code{k} times \code{k} matrix
}
\references{
\cite{Stewart, G.W. (1980), The efficient generation of random orthogonal matrices with an application to condition estimators, \emph{SIAM Journal on Numerical Analysis}, \bold{17},  403--409. <doi:10.1137/0717034>.}

}
\author{
Klaus Nordhausen
}

\examples{
Orth <- rorth(4)

crossprod(Orth)
tcrossprod(Orth)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}% use one of  RShowDoc("KEYWORDS")

