\name{Fclust.index}
\alias{Fclust.index}
\title{Cluster validity indexes}
\description{
 Performs some cluster validity indexes for choosing the optimal number of clusters \emph{k}.
}
\usage{
 Fclust.index (fclust.obj, index, alpha)
}
\arguments{
 \item{fclust.obj}{Object of class \code{fclust}}
 \item{index}{Cluster validity indexes to select the number of clusters: \code{PC} (partition coefficient), \code{PE} (partition entropy), \code{MPC} (modified partition coefficient), \code{SIL} (silhouette), \code{SIL.F} (fuzzy silhouette), \code{XB} (Xie and Beni),\code{ALL} for all the indexes (default: "ALL")}
 \item{alpha}{Weighting coefficient for the fuzzy silhouette index \code{SIL.F} (default: 1)}
}
\value{
 \item{out.index}{Vector containing the index values}
}
\details{
 \code{index} is not case-sensitive.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{PC}}, \code{\link{PE}}, \code{\link{MPC}}, \code{\link{SIL}}, \code{\link{SIL.F}}, \code{\link{XB}}, \code{\link{Fclust}}, \code{\link{Mc}}}
\examples{
## McDonald's data
data(Mc)
names(Mc)
## data normalization by dividing the nutrition facts by the Serving Size (column 1)
for (j in 2:(ncol(Mc)-1))
Mc[,j]=Mc[,j]/Mc[,1]
## removing the column Serving Size
Mc=Mc[,-1]
## fuzzy k-means
## (excluded the factor column Type (last column))
clust=FKM(Mc[,1:(ncol(Mc)-1)],k=6,m=1.5,stand=1)
## cluster validity indexes
all.indexes=Fclust.index(clust)
## Xie and Beni cluster validity index
XB.index=Fclust.index(clust,'XB')
}
\keyword{multivariate}
