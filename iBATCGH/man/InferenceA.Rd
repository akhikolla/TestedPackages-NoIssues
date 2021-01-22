\name{InferenceA}
\alias{InferenceA}
\title{
Postprocessing - Inference on the transition matrix
}
\description{
This function returns a matrix obtained as the elementwise mean of the association matrices, output of the Monte Carlo Markov Chain iterations, after removing burn-in.
}
\usage{
InferenceA(listA, niter, burnin)
}
\arguments{
  \item{listA}{
Second niter objects of the output of the main function
}
  \item{niter}{
Number of Monte Carlo Markov Chain iterations
}
  \item{burnin}{
Burn-in
}
}
\value{
Estimated transition matrix.
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
