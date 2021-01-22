\name{ARI.F}
\alias{ARI.F}
\title{
Fuzzy adjusted Rand index}
\description{
Produces the fuzzy version of the adjusted Rand index between a hard (reference) partition and a fuzzy partition.
}
\usage{
ARI.F(VC, U, t_norm)
}
\arguments{
\item{VC}{Vector of class labels}
\item{U}{Fuzzy membership degree matrix or data.frame}
\item{t_norm}{Type of the triangular norm: "minimum" (minimum triangular norm), "triangular product" (product norm) (default: "minimum")}
}
\value{
\code{ari.f}{Value of the fuzzy adjusted Rand index}
}
\references{
Campello, R.J., 2007. A fuzzy extension of the Rand index and other related indexes for clustering and classification assessment. Pattern Recognition Letters, 28, 833-841.
\cr Hubert, L., Arabie, P., 1985. Comparing partitions. Journal of Classification, 2, 193-218.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{
\code{\link{RI.F}}, \code{\link{JACCARD.F}}, \code{\link{Fclust.compare}}
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
## fuzzy adjusted Rand index
ari.f=ARI.F(VC=Mc$Type,U=clust$U)
}
}
\keyword{multivariate}
