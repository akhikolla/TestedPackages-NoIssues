\name{InferenceR}
\alias{InferenceR}
\title{
Postprocessing - Inference on the association matrix
}
\description{
This function performs posterior inference on the association matrix, returning the matrix of posterior probabilities of inclusion for each association gene expression - CGH. It also returns the binary matrix of significant links that exceed a threshold given as argument.
}
\usage{
InferenceR(listR, G, M, niter, burnin, threshold = 0.5)
}
\arguments{
  \item{listR}{
First niter objects of the output of the main function
}
  \item{G}{
Number of gene expression probes
}
  \item{M}{
Number of aCGH probes
}
  \item{niter}{
Number of Monte Carlo Markov Chains iterations
}
  \item{burnin}{
Burn-in
}
  \item{threshold}{
Threshold on the posterior probability of inclusion
}
}
\value{
A list made of the following items
\item{FreqMat }{Matrix of posterior probabilities of inclusion}
\item{final }{Binary association matrix}
}
\references{
Cassese A, Guindani M, Tadesse M, Falciani F, Vannucci M. A hierarchical Bayesian model for inference of copy number variants and their association to gene expression. Annals of Applied Statistics, 8(1), 148-175.\cr
Cassese A, Guindani M, Vannucci M. A Bayesian integrative model for genetical genomics with spatially informed variable selection. Cancer Informatics.
}
\author{
Alberto Cassese
}

\seealso{
See Also as \code{\link{Inference}}
}
\examples{
##See Inference
}
\keyword{Postprocessing}
