\name{synt.data}
\alias{synt.data}
\docType{data}
\title{Synthetic data}
\description{
  Synthetic dataset with 2 non-spherical clusters.
}
\usage{data(synt.data)}
\format{
 A matrix with 302 rows and 2 columns.
}
\details{
Although two clusters are clearly visible, fuzzy k-means fails to discover them.
The Gustafson and Kessel-like fuzzy k-means should be used for finding the known-in-advance clusters.}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{Fclust}}, \code{\link{FKM}}, \code{\link{FKM.gk}}, \code{\link{plot.fclust}}}
\examples{
\dontrun{
## synthetic data
data(synt.data)
plot(synt.data)
## fuzzy k-means
syntFKM=FKM(synt.data)
## Gustafson and Kessel-like fuzzy k-means
syntFKM.gk=FKM.gk(synt.data)
## plot of cluster structures from fuzzy k-means and Gustafson and Kessel-like fuzzy k-means
par(mfcol = c(2,1))
plot(syntFKM)
plot(syntFKM.gk)
}
}
\keyword{data}
\keyword{multivariate}
