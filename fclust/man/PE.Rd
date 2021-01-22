\name{PE}
\alias{PE}
\title{Partition entropy}
\description{
 Produces the partition entropy index. The optimal number of clusters \emph{k} is  is such that the index takes the minimum value.
}
\usage{
 PE (U, b)
}
\arguments{
 \item{U}{Membership degree matrix}
 \item{b}{Logarithmic base (default: exp(1))}
}
\value{
 \item{pe}{Value of the partition entropy index}
}
\references{
 Bezdek J.C., 1981. Pattern Recognition with Fuzzy Objective Function Algorithms. Plenum Press, New York.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{PC}}, \code{\link{MPC}}, \code{\link{SIL}}, \code{\link{SIL.F}}, \code{\link{XB}}, \code{\link{Fclust}}, \code{\link{Mc}}}
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
## partition entropy index
pe=PE(clust$U)
}
\keyword{multivariate}
