\name{synt.data2}
\alias{synt.data2}
\docType{data}
\title{Synthetic data}
\description{Synthetic dataset with 2 non-spherical clusters.}
\usage{data(synt.data2)}
\format{
  A matrix with 240 rows and 2 columns.
}
\details{
Although three clusters are clearly visible, Gustafson and Kessel - like fuzzy \emph{k}-means clustering algorithm \code{FKM.gk} fails due to singularity of some covariance matrix.
The Gustafson, Kessel and Babuska - like fuzzy \emph{k}-means clustering algorithm \code{FKM.gkb} should be used to avoid singularity problem.
}
\references{
 Gustafson E.E., Kessel W.C., 1978. Fuzzy clustering with a fuzzy covariance matrix. Proceedings of the IEEE Conference on Decision and Control, pp. 761-766.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{Fclust}}, \code{\link{FKM.gk}}, \code{\link{FKM.gkb}}, \code{\link{plot.fclust}}}


\examples{
data(synt.data2)
plot(synt.data2)

## Gustafson and Kessel-like fuzzy k-means
syntFKM.gk=FKM.gk(synt.data2, k = 3, RS = 1, seed = 123)
## Gustafson, Kessel and Babuska-like fuzzy k-means
syntFKM.gkb=FKM.gkb(synt.data2, k = 3, RS = 1, seed = 123)

}
\keyword{datasets}
