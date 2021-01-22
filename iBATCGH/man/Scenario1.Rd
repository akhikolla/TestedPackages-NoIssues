\name{Scenario1}
\alias{Scenario1}
\title{
Simulated data - Scenario 1
}
\description{
Simulates the data as described in the reference provided below (Scenario 1).
}
\usage{
Scenario1(sigmak = 0.1)
}

\arguments{
  \item{sigmak}{
Standard deviation of the error term
}
}
\value{
Return a list made of the following items
\item{Y }{Matrix of simulated gene expression}
\item{X }{Matrix of simulated CGH}
\item{Xi }{True matrix of hidden states}
\item{A }{Empirical transition matrix}
\item{mu }{True vector of state specific mean}
\item{Sd }{True vector of state specific sd}
\item{coeff }{True matrix of association coefficients between gene expression and CGH probes}
\item{distance }{Vector of distance between CGH probes}
\item{disfix }{Length of the chromosome}
}
\references{
Cassese A, Guindani M, Tadesse M, Falciani F, Vannucci M. A hierarchical Bayesian model for inference of copy number variants and their association to gene expression. Annals of Applied Statistics, 8(1), 148-175.\cr
Cassese A, Guindani M, Vannucci M. A Bayesian integrative model for genetical genomics with spatially informed variable selection. Cancer Informatics.
}
\author{
Alberto Cassese
}
\examples{
data <- Scenario1(sigmak = 0.1) 
}
\keyword{Simulated data}
