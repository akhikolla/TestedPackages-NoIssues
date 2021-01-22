\name{VAT}
\alias{VAT}
\title{Visual Assessment of (Cluster) Tendency}
\description{Digital intensity image to inspect the number of clusters}
\usage{
 VAT (Xca)
}
\arguments{
 \item{Xca}{Matrix or data.frame (usually data to be used in the clustering algorithm)}
}
\details{
Each cell refers to a dissimilarity between a pair of objects. Small dissimilarities are represented by dark shades and large dissimilarities are represented by light shades. In the plot the dissimilarities are reorganized in such a way that, roughly speaking, (darkly shaded) diagonal blocks correspond to clusters in the data. Therefore, \emph{k} dark blocks along its main diagonal suggest that the data contain \emph{k} (as yet unfound) clusters and the size of each block represents the approximate size of the cluster.
}
\references{
Bezdek J.C., Hathaway, R.J., 2002. VAT: a tool for visual assessment of (cluster) tendency. Proceedings of the IEEE International Joint Conference on Neural Networks, , pp. 2225?2230.
\cr Hathaway R.J., Bezdek J.C., 2003. Visual cluster validity for prototype generator clustering models. Pattern Recognition Letters, 24, 1563?1569.
\cr Huband J.M., Bezdek J.C., 2008. VCV2 ? Visual Cluster Validity. In Zurada J.M., Yen G.G., Wang J. (Eds.): Lecture Notes in Computer Science, 5050, pp. 293?308. Springer-Verlag, Berlin Heidelberg.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{plot.fclust}}, \code{\link{VIFCR}}, \code{\link{VCV}}, \code{\link{VCV2}}, \code{\link{Mc}}}
\examples{
## McDonald's data
data(Mc)
names(Mc)
## data normalization by dividing the nutrition facts by the Serving Size (column 1)
for (j in 2:(ncol(Mc)-1))
Mc[,j]=Mc[,j]/Mc[,1]
## data standardization (after removing the column Serving Size)
Mc=scale(Mc[,1:(ncol(Mc)-1)],center=TRUE,scale=TRUE)[,]
## plot of VAT
VAT(Mc)
}
\keyword{multivariate}
