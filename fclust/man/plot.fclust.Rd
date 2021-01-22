\name{plot.fclust}
\alias{plot.fclust}
\title{Plotting fuzzy clustering output}
\description{
Plot method for class \code{fclust}. The function creates a scatter plot visualizing the cluster structure. The objects are represented by points in the plot using observed variables or principal components.
}
\usage{
 ## S3 method for class fclust
 \method{plot}{fclust} (x, v1v2, colclus, umin, ucex, pca, ...)
}
\arguments{
 \item{x}{Object of class \code{fclust}}
 \item{v1v2}{Vector with two elements specifying the numbers of the variables (or of the principal components) to be plotted (default: \code{1:2}); in case of relational data, the argument is ignored}
 \item{colclus}{Vector specifying the color palette for the clusters (default: \code{palette(rainbow(k))})}
 \item{umin}{Lowest maximal membership degree such that an object is assigned to a cluster (default: 0)}
 \item{ucex}{Logical value specifying if the points are magnified according to the maximal membership degree (if \code{ucex=TRUE}) (default: \code{ucex=FALSE})}
 \item{pca}{Logical value specifying if the objects are represented using principal components (if \code{pca=TRUE}) (default: \code{pca=FALSE}); in case of relational data, the argument is ignored}
 \item{\dots}{Additional arguments arguments for \code{\link[graphics:plot]{plot}}}
}
\details{
In the scatter plot the objects are represented by circles (\code{pch=16}) and the prototypes by stars (\code{pch=8}) using observed variables (if \code{pca=FALSE}) or principal components (if \code{pca=TRUE}), the numbers of which are specified in \code{v1v2}. Their colors differ for every cluster according to \code{colclus}. Objects such that their maximal membership degrees are lower than \code{umin} are in black. The sizes of the circles depends on the maximal membership degrees of the corresponding objects if \code{ucex=TRUE}. Also note that principal components are extracted using standardized data.
\cr In case of relational data, the first two components resulting from Non-metric Multidimensional Scaling performed using the package \pkg{MASS} are used.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{VIFCR}}, \code{\link{VAT}}, \code{\link{VCV}}, \code{\link{VCV2}}, \code{\link{Fclust}}, \code{\link{print.fclust}}, \code{\link{summary.fclust}}, \code{\link{Mc}}}
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
## Scatter plot of Calories vs Cholesterol (mg)
names(Mc)
plot(clust,v1v2=c(1,5))
## Scatter plot of Calories vs Cholesterol (mg) using gray levels for the clusters
plot(clust,v1v2=c(1,5),colclus=gray.colors(6))
## Scatter plot of Calories vs Cholesterol (mg)
## coloring in black objects with maximal membership degree lower than 0.5
plot(clust,v1v2=c(1,5),umin=0.5)
## Scatter plot of Calories vs Cholesterol (mg)
## coloring in black objects with maximal membership degree lower than 0.5
## and magnifying the points according to the maximal membership degree
plot(clust,v1v2=c(1,5),umin=0.5,ucex=TRUE)
## Scatter plot using the first two principal components and
## coloring in black objects with maximal membership degree lower than 0.3
plot(clust,v1v2=1:2,umin=0.3,pca=TRUE)
}
\keyword{multivariate}
