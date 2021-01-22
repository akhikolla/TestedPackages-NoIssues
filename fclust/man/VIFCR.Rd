\name{VIFCR}
\alias{VIFCR}
\title{Visual inspection of fuzzy clustering results}
\description{
Plots for validation of fuzzy clustering results. Three plots (selected by \code{which}) are available.
}
\usage{
 VIFCR (fclust.obj, which)
}
\arguments{
 \item{fclust.obj}{Object of class \code{fclust}}
 \item{which}{If a subset of the plots is required, specify a subset of the numbers \code{1:3} (default: \code{1:3})}.
}
\details{
Plot 1 (\code{which=1}). Histogram of the membership degrees setting \code{breaks=seq(from=0,to=1,by=0.1)}. The frequencies are scaled so that the heights of the first and the latter rectangles are the same in the ideal case of crisp (non-fuzzy) memberships. The fuzzy clustering solution should be such that the heights of the first and the latter rectangles are high and those of the rectangles in the middle are low. High heights of rectangles in the middle denote the presence of ambiguous membership degrees. This is an indicator for a non-optimal clustering result.
\cr Plot 2 (\code{which=2}). Scatter plot of the objects at the co-ordinates (u1,u2). For each object, u1 and u2 denote, respectively, the highest and the second highest membership degrees. All points lie within the triangle with vertices (0,0), (0.5,0.5) and (1,0). In the ideal case of (almost) crisp membership degrees all points are near the vertex (1,0). Points near the vertex (0.5,0.5) highlight ambiguous objects shared by two clusters. Points near the vertex (0,0) are usually outliers characterized by low membership degrees to all clusters (provided that the noise approach is considered).
\cr Plot 3 (\code{which=3}). For each cluster, scatter plot of the of the objects at the co-ordinates (dc,uc). For each object, dc is the squared Euclidean distance between the object and the cluster prototype and uc is the membership degree of the object to the cluster. The ideal case is such that points are in the upper left area or in the lower right area. In fact, this highlights high membership degrees for small distances and low membership degrees for large distances.
}
\references{
Klawonn F., Chekhtman V., Janz E., 2003. Visual inspection of fuzzy clustering results. In Benitez J.M., Cordon O., Hoffmann, F., Roy R. (Eds.):  Advances in Soft Computing - Engineering Design and Manufacturing, pp. 65-76. Springer, London.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{plot.fclust}}, \code{\link{VAT}}, \code{\link{VCV}}, \code{\link{VCV2}}, \code{\link{unemployment}}}
\examples{
## unemployment data
data(unemployment)
## fuzzy k-means
unempFKM=FKM(unemployment,k=3,stand=1)
## all plots
VIFCR(unempFKM)
## plots 1 and 3
VIFCR(unempFKM,c(1,3))
}
\keyword{multivariate}
