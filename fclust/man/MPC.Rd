\name{MPC}
\alias{MPC}
\title{Modified partition coefficient}
\description{
 Produces the modified partition coefficient index. The optimal number of clusters \emph{k} is such that the index takes the maximum value.
}
\usage{
 MPC (U)
}
\arguments{
 \item{U}{Membership degree matrix}
}
\value{
 \item{mpc}{Value of the modified partition coefficient index}
}
\references{
 Dave' R.N., 1996. Validating fuzzy partitions obtained through \emph{c}-shells clustering. Pattern Recognition Letters, 17, 613-623.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{PC}}, \code{\link{PE}}, \code{\link{SIL}}, \code{\link{SIL.F}}, \code{\link{XB}}, \code{\link{Fclust}}, \code{\link{Mc}}}
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
## modified partition coefficient
mpc=MPC(clust$U)
}
\keyword{multivariate}
