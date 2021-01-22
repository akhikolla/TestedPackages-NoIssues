\name{InferenceXi}
\alias{InferenceXi}
\title{
Postprocessing - Inference on the latent states
}
\description{
This function returns the modal latent states.
}
\usage{
InferenceXi(listXi, niter, burnin)
}
\arguments{
  \item{listXi}{
Last three objects of the output of the main function
}
  \item{niter}{
Number of Monte Carlo Markov Chain iterations
}
  \item{burnin}{
Burn-in
}
}
\details{
Must use the same burn-in as in the main function. 
}
\value{
Matrix of modal latent states, i.e. estimated Copy Number Variants. A four class classification is considered:
\enumerate{
\item Loss
\item Neutral
\item Gain
\item Amplification
}
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
