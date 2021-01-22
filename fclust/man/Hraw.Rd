\name{Hraw}
\alias{Hraw}
\title{Raw prototypes}
\description{
 Produces prototypes using the original units of measurement of \emph{X} (useful if the clustering algorithm is run using standardized data).
}
\usage{
 Hraw (X, H)
}
\arguments{
 \item{X}{Matrix or data.frame}
 \item{H}{Prototype matrix}
}
\value{
 \item{Hraw}{Prototypes matrix using the original units of measurement of \code{X}}
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{Fclust}}, \code{\link{unemployment}}}
\examples{
## example n.1 (k-means case)
## unemployment data
data(unemployment)
## fuzzy k-means
unempFKM=FKM(unemployment,k=3,stand=1)
## standardized prototypes
unempFKM$H
## prototypes using the original units of measurement
unempFKM$Hraw=Hraw(unempFKM$X,unempFKM$H)
## example n.2  (k-medoids case)
## unemployment data
data(unemployment)
## fuzzy k-medoids
\dontrun{
## It may take more than a few seconds
unempFKM.med=FKM.med(unemployment,k=3,RS=10,stand=1)
## prototypes using the original units of measurement:
## in fuzzy k-medoids one can equivalently use
unempFKM.med$Hraw1=Hraw(unempFKM.med$X,unempFKM.med$H)
unempFKM.med$Hraw2=unempFKM.med$X[unempFKM.med$medoid,]}
}
\keyword{multivariate}
