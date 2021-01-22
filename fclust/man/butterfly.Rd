\name{butterfly}
\alias{butterfly}
\docType{data}
\title{Butterfly data}
\description{
  Synthetic dataset with 2 clusters and some outliers.
}
\usage{data(butterfly)}
\format{
 A matrix with 17 rows and 2 columns.
}
\details{
The butterfly data motivate the need for the fuzzy approach to clustering.
\cr The presence of outliers can be handled using fuzzy \emph{k}-means with noise cluster. In fact, differently from fuzzy \emph{k}-means, the membership degrees of the outliers are low for all the clusters.}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{Fclust}}, \code{\link{FKM}}, \code{\link{FKM.noise}}}
\examples{
## butterfly data
data(butterfly)
plot(butterfly,type='n')
text(butterfly[,1],butterfly[,2],labels=rownames(butterfly),cex=0.7,lwd=2)
## membership degree matrix using fuzzy k-means (rounded)
round(FKM(butterfly)$U,2)
## membership degree matrix using fuzzy k-means with noise cluster (rounded)
round(FKM.noise(butterfly,delta=3)$U,2)
}
\keyword{data}
\keyword{multivariate}
