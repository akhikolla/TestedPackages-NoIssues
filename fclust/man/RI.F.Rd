\name{RI.F}
\alias{RI.F}
\title{
Fuzzy Rand index}
\description{
Produces the fuzzy version of the Rand index between a hard (reference) partition and a fuzzy partition.
}
\usage{
RI.F(VC, U, t_norm)
}
\arguments{
\item{VC}{Vector of class labels}
\item{U}{Fuzzy membership degree matrix or data.frame}
\item{t_norm}{Type of the triangular norm: "minimum" (minimum triangular norm), "triangular product" (product norm) (default: "minimum")}
}
\value{
\code{ri.f}{Value of the fuzzy adjusted Rand index}
}
\references{
Campello, R.J., 2007. A fuzzy extension of the Rand index and other related indexes for clustering and classification assessment. Pattern Recognition Letters, 28, 833-841.
\cr Rand, W.M., 1971. Objective criteria for the evaluation of clustering methods. Journal of the American Statistical Association, 66, 846-850.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{
\code{\link{ARI.F}}, \code{\link{JACCARD.F}}, \code{\link{Fclust.compare}}
}
\examples{
\dontrun{
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
## fuzzy Rand index
ri.f=RI.F(VC=Mc$Type,U=clust$U)
}
}
\keyword{multivariate}
