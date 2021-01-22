\name{InitXi}
\alias{InitXi}
\title{
Preprocessing - Initialize matrix of latent states
}
\description{
This function takes a matrix of CGH data as the only argument and returns a crude estimate of the corresponding latent copy number states.
}
\usage{
InitXi(X, bounds = c(-0.5, 0.29, 0.79))
}
\arguments{
  \item{X}{
Matrix of aCGH data
}
  \item{bounds}{
Vector of threshold used to estimate the latent states
}
}
\details{
Given as argument a vector of threshold bounds the function simply applies the thresholding to the data and groups them into four subsets. Each subset is associated to a specific latent state.
}
\value{
Return a matrix of estimated latent states, that could be used as input of the main function.
}
\references{
Cassese A, Guindani M, Tadesse M, Falciani F, Vannucci M. A hierarchical Bayesian model for inference of copy number variants and their association to gene expression. Annals of Applied Statistics, 8(1), 148-175.\cr
Cassese A, Guindani M, Vannucci M. A Bayesian integrative model for genetical genomics with spatially informed variable selection. Cancer Informatics.
}
\author{
Alberto Cassese
}

\seealso{
See Also \code{\link{Tran}}
}
\examples{
data(NCI_60)
X <- NCI_60$aCGH
xi <- InitXi(X)
}
