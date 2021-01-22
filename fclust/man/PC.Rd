\name{PC}
\alias{PC}
\title{Partition coefficient}
\description{
 Produces the partition coefficient index. The optimal number of clusters \emph{k} is  is such that the index takes the maximum value.
}
\usage{
 PC (U)
}
\arguments{
 \item{U}{Membership degree matrix}
}
\value{
 \item{pc}{Value of the partition coefficient index}
}
\references{
 Bezdek J.C., 1974. Cluster validity with fuzzy sets. Journal of Cybernetics, 3, 58-73.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{PE}}, \code{\link{MPC}}, \code{\link{SIL}}, \code{\link{SIL.F}}, \code{\link{XB}}, \code{\link{Fclust}}, \code{\link{Mc}}}
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
## partition coefficient
pc=PC(clust$U)
}
\keyword{multivariate}
