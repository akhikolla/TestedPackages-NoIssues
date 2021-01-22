\name{Center}
\alias{Center}
\title{
Preprocessing - Center gene expression data
}
\description{
This function takes as argument a matrix of gene expression measurements Y, and returns the matrix obtained after centering each column with respect to its mean.
}
\usage{
Center(Y)
}
\arguments{
  \item{Y}{
The matrix of data to be centered.
}
}
\details{
This function center each column of the gene expression matrix.
}
\value{
The matrix of data obtained after centering each column with respect to its mean.
}

\author{
Alberto Cassese
}

\examples{
data(NCI_60)
Y <- NCI_60$Affy
Y<-Center(Y)
}
\keyword{Preprocessing}
