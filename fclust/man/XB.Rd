\name{XB}
\alias{XB}
\title{Xie and Beni index}
\description{
 Produces the Xie and Beni index. The optimal number of clusters \emph{k} is  is such that the index takes the minimum value.
}
\usage{
 XB (Xca, U, H, m)
}
\arguments{
 \item{Xca}{Matrix or data.frame}
 \item{U}{Membership degree matrix}
 \item{H}{Prototype matrix}
 \item{m}{Parameter of fuzziness (default: 2)}
}
\value{
 \item{xb}{Value of the Xie and Beni index}
}
\details{
 \code{Xca} should contain the same dataset used in the clustering algorithm, i.e., if the clustering algorithm is run using standardized data, then \code{XB} should be computed using the same standardized data.
 \cr \code{m} should be the same parameter of fuzziness used in the clustering algorithm.
}
\references{
 Xie X.L., Beni G. (1991). A validity measure for fuzzy clustering, IEEE Transactions on Pattern Analysis and Machine Intelligence, 13, 841-847.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{PC}}, \code{\link{PE}}, \code{\link{MPC}}, \code{\link{SIL}}, \code{\link{SIL.F}}, \code{\link{Fclust}}, \code{\link{Mc}}}
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
## Xie and Beni index
xb=XB(clust$Xca,clust$U,clust$H,clust$m)
}
\keyword{multivariate}
