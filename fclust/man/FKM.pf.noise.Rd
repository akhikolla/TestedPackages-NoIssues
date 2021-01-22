\name{FKM.pf.noise} 
\alias{FKM.pf.noise}
\title{Fuzzy k-means with polynomial fuzzifier and noise cluster}
\description{
 Performs the fuzzy \emph{k}-means clustering algorithm with polynomial fuzzifier function and noise cluster.
 \cr The polynomial fuzzifier creates areas of crisp membership degrees around the prototypes while, outside of these areas of crisp membership degrees, fuzzy membership degrees are given. Therefore, the polynomial fuzzifier produces membership degrees equal to one for objects clearly assigned to clusters, that is, very close to the cluster prototypes.
 \cr The noise cluster is an additional cluster (with respect to the \emph{k} standard clusters) such that objects recognized to be outliers are assigned to it with high membership degrees.
} 
\usage{ 
 FKM.pf.noise (X, k, b, delta, RS, stand, startU, index, alpha, conv, maxit, seed)
}
\arguments{
 \item{X}{Matrix or data.frame}
 \item{k}{An integer value or vector specifying the number of clusters for which the \code{index} is to be calculated (default: 2:6)}
 \item{b}{Parameter of the polynomial fuzzifier (default: 0.5)}
 \item{delta}{Noise distance (default: average Euclidean distance between objects and prototypes from \code{FKM.pf} using the same values of \code{k} and \code{m})}
 \item{RS}{Number of (random) starts (default: 1)}
 \item{stand}{Standardization: if \code{stand=1}, the clustering algorithm is run using standardized data (default: no standardization)}
 \item{startU}{Rational start for the membership degree matrix \code{U} (default: no rational start)}
\item{index}{Cluster validity index to select the number of clusters: \code{PC} (partition coefficient), \code{PE} (partition entropy), \code{MPC} (modified partition coefficient), \code{SIL} (silhouette), \code{SIL.F} (fuzzy   silhouette), \code{XB} (Xie and Beni) (default: "SIL.F")}
 \item{alpha}{Weighting coefficient for the fuzzy silhouette index \code{SIL.F} (default: 1)}
 \item{conv}{Convergence criterion (default: 1e-9)}
 \item{maxit}{Maximum number of iterations (default: 1e+6)}
 \item{seed}{Seed value for random number generation (default: NULL)}
}
\value{
 Object of class \code{fclust}, which is a list with the following components:
 \item{U}{Membership degree matrix}
 \item{H}{Prototype matrix}
 \item{F}{Array containing the covariance matrices of all the clusters (\code{NULL} for \code{FKM.pf.noise})}
 \item{clus}{Matrix containing the indexes of the clusters where the objects are assigned (column 1) and the associated membership degrees (column 2)}
 \item{medoid}{Vector containing the indexes of the medoid objects (\code{NULL} for \code{FKM.pf.noise})}
 \item{value}{Vector containing the loss function values for the \code{RS} starts}
 \item{criterion}{Vector containing the values of the cluster validity index}
 \item{iter}{Vector containing the numbers of iterations for the \code{RS} starts}
 \item{k}{Number of clusters}
 \item{m}{Parameter of fuzziness (\code{NULL} for \code{FKM.pf.noise})}
 \item{ent}{Degree of fuzzy entropy (\code{NULL} for \code{FKM.pf.noise})}
 \item{b}{Parameter of the polynomial fuzzifier}
 \item{vp}{Volume parameter (\code{NULL} for \code{FKM.pf.noise})} 
 \item{delta}{Noise distance}
 \item{gam}{Weighting parameter for the fuzzy covariance matrices (\code{NULL} for \code{FKM.pf.noise})}
 \item{mcn}{Maximum condition number for the fuzzy covariance matrices (\code{NULL} for \code{FKM.pf.noise})}
 \item{stand}{Standardization (Yes if \code{stand=1}, No if \code{stand=0})}
 \item{Xca}{Data used in the clustering algorithm (standardized data if \code{stand=1})}
 \item{X}{Raw data}
  \item{D}{Dissimilarity matrix (\code{NULL} for \code{FKM.pf.noise})}
 \item{call}{Matched call}
} 
\details{
 If \code{startU} is given, the argument \code{k} is ignored (the number of clusters is \code{ncol(startU)}).
 \cr If \code{startU} is given, the first element of \code{value}, \code{cput} and \code{iter} refer to the rational start.
} 
\references{
Dave' R.N., 1991. Characterization and detection of noise in clustering. Pattern Recognition Letters, 12, 657-664.
\cr Winkler R., Klawonn F., Hoeppner F., Kruse R., 2010. Fuzzy cluster analysis of larger data sets. In: Scalable Fuzzy Algorithms for Data Management and Analysis: Methods and Design IGI Global, pp. 302-331. IGI Global, Hershey.
\cr Winkler R., Klawonn F., Kruse R., 2011. Fuzzy clustering with polynomial fuzzifier function in connection with M-estimators. Applied and Computational Mathematics, 10, 146-163.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{FKM.pf}}, \code{\link{Fclust}}, \code{\link{Fclust.index}}, \code{\link{print.fclust}}, \code{\link{summary.fclust}}, \code{\link{plot.fclust}}, \code{\link{Mc}}}
\examples{
## McDonald's data
data(Mc)
names(Mc)
## data normalization by dividing the nutrition facts by the Serving Size (column 1)
for (j in 2:(ncol(Mc)-1))
Mc[,j]=Mc[,j]/Mc[,1]
## removing the column Serving Size
Mc=Mc[,-1]
## fuzzy k-means with polynomial fuzzifier and noise cluster, fixing the number of clusters 
## (excluded the factor column Type (last column))
clust=FKM.pf.noise(Mc[,1:(ncol(Mc)-1)],k=6,stand=1)
## fuzzy k-means with polynomial fuzzifier and noise cluster, selecting the number of clusters 
## (excluded the factor column Type (last column))
clust=FKM.pf.noise(Mc[,1:(ncol(Mc)-1)],k=2:6,stand=1)
}
\keyword{multivariate}
