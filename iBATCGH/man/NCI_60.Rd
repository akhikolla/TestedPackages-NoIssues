\name{NCI_60}
\alias{NCI_60}
\docType{data}
\title{
NCI-60 cancer cell lines data
}
\description{
Processed and filtered NCI-60 cancer cell lines data, as described in the reference.
}
\usage{data(NCI_60)}
\format{
The format is a list of 3 objects\cr
 $ aCGH\cr
 $ Affy\cr
 $ distance
}

\source{
Full matrix of data downloaded from \url{discover.nci.nih.gov/cellminer}.
}
\references{
Cassese A, Guindani M, Tadesse M, Falciani F, Vannucci M. A hierarchical Bayesian model for inference of copy number variants and their association to gene expression. Annals of Applied Statistics, 8(1), 148-175.
}
\examples{
data <- data(NCI_60)
}
\keyword{datasets}
