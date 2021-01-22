\name{Fclust.compare}
\alias{Fclust.compare}
\title{
Similarity between partitions}
\description{
Performs some measures of similarity between a hard (reference) partition and a fuzzy partition.}
\usage{
Fclust.compare(VC, U, index, tnorm)
}
\arguments{
\item{VC}{Vector of class labels}
\item{U}{Fuzzy membership degree matrix or data.frame}
\item{index}{Measures of similarity: "ARI.F" (fuzzy version of the adjuster Rand index), "RI.F" (fuzzy version of the Rand index), "JACCARD.F" (fuzzy version of the Jaccard index), "ALL" for all the indexes (default: "ALL")}
\item{tnorm}{Type of the triangular norm: "minimum" (minimum triangular norm), "triangular product" (product norm) (default: "minimum")}
}
\value{
 \code{out.index}{Vector containing the similarity measures}
}
\details{
\code{index} is not case-sensitive. All the measures of similarity share the same properties of their non-fuzzy counterpart.
}
\references{
Campello, R.J., 2007. A fuzzy extension of the Rand index and other related indexes for clustering and classification assessment. Pattern Recognition Letters, 28, 833-841.
\cr Hubert, L., Arabie, P., 1985. Comparing partitions. Journal of Classification, 2, 193-218.
\cr Jaccard, P., 1901. Étude comparative de la distribution florale dans une portion des Alpes et des Jura. Bulletin de la Société Vaudoise des Sciences Naturelles, 37, 547-579.
\cr Rand, W.M., 1971. Objective criteria for the evaluation of clustering methods. Journal of the American Statistical Association, 66, 846-850.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{
\code{\link{RI.F}}, \code{\link{ARI.F}}, \code{\link{JACCARD.F}}
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
## all measures of similarity
all.indexes=Fclust.compare(VC=Mc$Type,U=clust$U)
## fuzzy adjusted Rand index
Fari.index=Fclust.compare(VC=Mc$Type,U=clust$U,index="ARI.F")
}
}
\keyword{multivariate}
