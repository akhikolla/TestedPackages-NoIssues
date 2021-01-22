\name{summary.fclust}
\alias{summary.fclust}
\title{Summarizing fuzzy clustering output}
\description{
Summary method for class \code{fclust}.
}
\usage{
 ## S3 method for class fclust
 \method{summary}{fclust} (object, ...)
}
\arguments{
 \item{object}{Object of class \code{fclust}}
  \item{\dots}{Additional arguments for \code{\link[base:summary]{summary}}}
}
\details{
The function displays the number of objects, the number of clusters, the cluster sizes, the closest hard clustering partition (objects assigned to the clusters with the highest membership degree), the cluster memberships (using the closest hard clustering partition), the number of objects with unclear assignment (when the maximal membership degree is lower than 0.5), the objects with unclear assignment and the cluster sizes without unclear assignments (only if objects with unclear assignment are present), the cluster summary (for every cluster: size, minimal membership degree, maximal membership degree, average membership degree, number of objects with unclear assignment) and the Euclidean distance matrix for the cluster prototypes.}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{Fclust}}, \code{\link{print.fclust}}, \code{\link{plot.fclust}}, \code{\link{unemployment}}}
\examples{
## unemployment data
data(unemployment)
## fuzzy k-means
unempFKM=FKM(unemployment,k=3,stand=1)
summary(unempFKM)
}
\keyword{multivariate}
