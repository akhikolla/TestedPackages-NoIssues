\name{VCV}
\alias{VCV}
\title{Visual Cluster Validity}
\description{Digital intensity image generated using the prototype matrix (and the membership degree matrix) to do cluster validation. The function also plots the VAT image.}
\usage{
 VCV (Xca, U, H, which)
}
\arguments{
 \item{Xca}{Matrix or data.frame (usually data used in the clustering algorithm)}
 \item{U}{Membership degree matrix}
 \item{H}{Prototype matrix}
 \item{which}{If a subset of the plots is required, specify a subset of the numbers \code{1:2} (default: \code{1:2})}.

}
\details{
Plot 1 (\code{which=1}): VAT. Each cell refers to a dissimilarity between a pair of objects. Small dissimilarities are represented by dark shades and large dissimilarities are represented by light shades. In the plot the dissimilarities are reorganized in such a way that, roughly speaking, (darkly shaded) diagonal blocks correspond to clusters in the data. Therefore, \emph{k} dark blocks along its main diagonal suggest that the data contain \emph{k} (as yet unfound) clusters and the size of each block represents the approximate size of the cluster.
\cr Plot 2 (\code{which=2}): VCV. Each cell refers to a dissimilarity between a pair of objects computed with respect to the cluster prototypes. Small dissimilarities are represented by dark shades and large dissimilarities are represented by light shades. In the plot the dissimilarities are organized by reordering the clusters (the original first cluster is the first reordered cluster and the remaining clusters are reordered so that (new) cluster \emph{c}+1 is the nearest of the remaining clusters to (newly indexed) cluster \emph{c}) and the objects (in accordance with decreasing membership degrees). If \emph{k} dark blocks along its main diagonal are visible, then a \emph{k}-cluster structure is revealed. Note that the actual number of clusters can be revealed even when a larger number of clusters is used. This suggests that the correct value of \emph{k} can sometimes be found by running the algorithm with a large value of \emph{k}, and then ascertaining its correct value from the visual evidence in the VCV image.
}
\references{
Bezdek J.C., Hathaway, R.J., 2002. VAT: a tool for visual assessment of (cluster) tendency. Proceedings of the IEEE International Joint Conference on Neural Networks, , pp. 2225?2230.
\cr Hathaway R.J., Bezdek J.C., 2003. Visual cluster validity for prototype generator clustering models. Pattern Recognition Letters, 24, 1563?1569.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{plot.fclust}}, \code{\link{VIFCR}}, \code{\link{VAT}}, \code{\link{VCV2}}, \code{\link{Mc}}}
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
## plots of VAT and VCV
VCV(clust$Xca,clust$U,clust$H)
## plot of VCV
VCV(clust$Xca,clust$U,clust$H, 2)
}
\keyword{multivariate}
