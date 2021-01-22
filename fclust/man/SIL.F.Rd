\name{SIL.F}
\alias{SIL.F}
\title{Fuzzy silhouette index}
\description{
 Produces the fuzzy silhouette index. The optimal number of clusters \emph{k} is  is such that the index takes the maximum value.
}
\usage{
 SIL.F (Xca, U, alpha, distance)
}
\arguments{
 \item{Xca}{Matrix or data.frame}
 \item{U}{Membership degree matrix}
 \item{alpha}{Weighting coefficient (default: 1)}
 \item{distance}{If \code{distance=TRUE}, \code{Xca} is assumed to contain distances/dissimilarities (default: FALSE)}
}
\value{
 \item{sil.f}{Value of the fuzzy silhouette index}
}
\details{
 \code{Xca} should contain the same dataset used in the clustering algorithm, i.e., if the clustering algorithm is run using standardized data, then \code{SIL.F} should be computed using the same standardized data.
 \cr Set \code{distance=TRUE} if \code{Xca} is a distance/dissimilarity matrix.
}
\references{
 Campello R.J.G.B., Hruschka E.R., 2006. A fuzzy extension of the silhouette width criterion for cluster analysis. Fuzzy Sets and Systems, 157, 2858-2875.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{PC}}, \code{\link{PE}}, \code{\link{MPC}}, \code{\link{SIL}}, \code{\link{XB}}, \code{\link{Fclust}}, \code{\link{Mc}}}
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
## fuzzy silhouette index
sil.f=SIL.F(clust$Xca,clust$U)
}
\keyword{multivariate}
