\name{Fclust}
\alias{Fclust}
\title{Fuzzy clustering}
\description{
 Performs fuzzy clustering by using the algorithms available in the package.
}
\usage{
 Fclust (X, k, type, ent, noise, stand, distance)
}
\arguments{
 \item{X}{Matrix or data.frame}
 \item{k}{An integer value specifying the number of clusters (default: 2)}
 \item{type}{Fuzzy clustering algorithm: \code{"standard"} (standard algorithms: FKM - type if \code{distance=TRUE}, NEFRC - type if if \code{distance=FALSE}), \code{"polynomial"} (algorithms with the polynomial fuzzifier), \code{"gk"} (Gustafson and Kessel - like algorithms), \code{"gkb"} (Gustafson, Kessel and Babuska - like algorithms), \code{"medoids"} (Medoid - based algorithms) (default: "standard")}
 \item{ent}{If \code{ent=TRUE}, the entropy regularization variant of the algorithm is run (default: FALSE)}
 \item{noise}{If \code{noise=TRUE}, the noise cluster variant of the algorithm is run (default: FALSE)}
 \item{stand}{Standardization: if \code{stand=1}, the clustering algorithm is run using standardized data (default: no standardization)}
 \item{distance}{If \code{distance=TRUE}, \code{X} is assumed to be a distance/dissimilarity matrix (default: FALSE)}
}
\value{
 \item{clust}{Object of class \code{fclust}}
}
\details{
 The clustering algorithms are run by using default options.
 \cr To specify different options, use the corresponding function.
} 
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{print.fclust}}, \code{\link{summary.fclust}}, \code{\link{plot.fclust}}, \code{\link{FKM}}, \code{\link{FKM.ent}}, \code{\link{FKM.gk}}, \code{\link{FKM.gk.ent}}, \code{\link{FKM.gkb}}, \code{\link{FKM.gkb.ent}}, \code{\link{FKM.med}}, \code{\link{FKM.pf}}, \code{\link{FKM.noise}}, \code{\link{FKM.ent.noise}}, \code{\link{FKM.gk.noise}}, \code{\link{FKM.gkb.ent.noise}},  \code{\link{FKM.gkb.noise}}, \code{\link{FKM.gk.ent.noise}},\code{\link{FKM.med.noise}}, \code{\link{FKM.pf.noise}}, \code{\link{NEFRC}}, \code{\link{NEFRC.noise}}, \code{\link{Fclust.index}}, \code{\link{Fclust.compare}}}
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
clust=Fclust(Mc[,1:(ncol(Mc)-1)],k=6,type="standard",ent=FALSE,noise=FALSE,stand=1,distance=FALSE)
## fuzzy k-means with polynomial fuzzifier 
## (excluded the factor column Type (last column))
clust=Fclust(Mc[,1:(ncol(Mc)-1)],k=6,type="polynomial",ent=FALSE,noise=FALSE,stand=1,distance=FALSE)
## fuzzy k-means with entropy regularization
## (excluded the factor column Type (last column))
clust=Fclust(Mc[,1:(ncol(Mc)-1)],k=6,type="standard",ent=TRUE,noise=FALSE,stand=1,distance=FALSE)
## fuzzy k-means with noise cluster
## (excluded the factor column Type (last column))
clust=Fclust(Mc[,1:(ncol(Mc)-1)],k=6,type="standard",ent=FALSE,noise=TRUE,stand=1,distance=FALSE)
}
}
\keyword{multivariate}