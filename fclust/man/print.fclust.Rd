\name{print.fclust}
\alias{print.fclust}
\title{Printing fuzzy clustering output}
\description{
Print method for class \code{fclust}.
}
\usage{
 ## S3 method for class fclust
 \method{print}{fclust} (x, ...)
}
\arguments{
 \item{x}{Object of class \code{fclust}}
  \item{\dots}{Additional arguments for \code{\link[base:print]{print}}}
}
\details{
The function displays the number of objects, the number of clusters, the closest hard clustering partition (objects assigned to the clusters with the highest membership degree) and the membership degree matrix (rounded).
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{Fclust}}, \code{\link{summary.fclust}}, \code{\link{plot.fclust}}, \code{\link{unemployment}}}
\examples{
## unemployment data
data(unemployment)
## fuzzy k-means
unempFKM=FKM(unemployment,k=3,stand=1)
unempFKM
}
\keyword{multivariate}
