\name{Tran}
\alias{Tran}
\title{
Preprocessing - Tran
}
\description{
This function computes the transition matrix corresponding to a specific matrix of latent states.
}
\usage{
Tran(xi)
}
\arguments{
  \item{xi}{
Matrix of estimated latent states
}
}
\details{
The matrix of latent states must consider only 4 states:
\enumerate{
\item Loss
\item Neutral
\item Gain
\item Amplification
}
}
\value{
A 4 by 4 matrix with transition probabilities derived from the input matrix of latent states.
}
\references{
Guha S, Li Y, Neuberg D. Bayesian hidden Markov modelling of array cgh data. JASA. 2008;103(482):485-497.
}
\author{
Alberto Cassese
}

\examples{
data(NCI_60)

X <- NCI_60$aCGH
xi <- InitXi(X)
tran <- Tran(xi)
}
\keyword{Preprocessing}
