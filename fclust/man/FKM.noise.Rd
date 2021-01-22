\name{FKM.noise} 
\alias{FKM.noise}
\title{Fuzzy k-means with noise cluster}
\description{
 Performs the fuzzy \emph{k}-means clustering algorithm with noise cluster. 
 \cr The noise cluster is an additional cluster (with respect to the \emph{k} standard clusters) such that objects recognized to be outliers are assigned to it with high membership degrees.
} 
\usage{ 
 FKM.noise (X, k, m, delta, RS, stand, startU, index, alpha, conv, maxit, seed)
}
\arguments{
 \item{X}{Matrix or data.frame}
 \item{k}{An integer value or vector specifying the number of clusters for which the \code{index} is to be calculated (default: 2:6)}
 \item{m}{Parameter of fuzziness (default: 2)}
 \item{delta}{Noise distance (default: average Euclidean distance between objects and prototypes from \code{FKM} using the same values of \code{k} and \code{m})}
 \item{RS}{Number of (random) starts (default: 1)}
 \item{stand}{Standardization: if \code{stand=1}, the clustering algorithm is run using standardized data (default: no standardization)}
 \item{startU}{Rational start for the membership degree matrix \code{U} (default: no rational start)}
\item{index}{Cluster validity index to select the number of clusters: \code{PC} (partition coefficient), \code{PE} (partition entropy), \code{MPC} (modified partition coefficient), \code{SIL} (silhouette), \code{SIL.F} (fuzzy silhouette), \code{XB} (Xie and Beni) (default: "SIL.F")}
 \item{alpha}{Weighting coefficient for the fuzzy silhouette index \code{SIL.F} (default: 1)}
 \item{conv}{Convergence criterion (default: 1e-9)}
 \item{maxit}{Maximum number of iterations (default: 1e+6)}
 \item{seed}{Seed value for random number generation (default: NULL)}
}
\value{
 Object of class \code{fclust}, which is a list with the following components:
 \item{U}{Membership degree matrix}
 \item{H}{Prototype matrix}
 \item{F}{Array containing the covariance matrices of all the clusters (\code{NULL} for \code{FKM.noise})}
 \item{clus}{Matrix containing the indexes of the clusters where the objects are assigned (column 1) and the associated membership degrees (column 2)}
 \item{medoid}{Vector containing the indexes of the medoid objects (\code{NULL} for \code{FKM.noise})}
 \item{value}{Vector containing the loss function values for the \code{RS} starts}
 \item{criterion}{Vector containing the values of the cluster validity index}
 \item{iter}{Vector containing the numbers of iterations for the \code{RS} starts}
 \item{k}{Number of clusters}
 \item{m}{Parameter of fuzziness}
 \item{ent}{Degree of fuzzy entropy (\code{NULL} for \code{FKM.noise})}
 \item{b}{Parameter of the polynomial fuzzifier (\code{NULL} for \code{FKM.noise})}
 \item{vp}{Volume parameter (\code{NULL} for \code{FKM.noise})} 
 \item{delta}{Noise distance}
 \item{gam}{Weighting parameter for the fuzzy covariance matrices (\code{NULL} for \code{FKM.noise})}
 \item{mcn}{Maximum condition number for the fuzzy covariance matrices (\code{NULL} for \code{FKM.noise})}
 \item{stand}{Standardization (Yes if \code{stand=1}, No if \code{stand=0})}
 \item{Xca}{Data used in the clustering algorithm (standardized data if \code{stand=1})}
 \item{X}{Raw data}
  \item{D}{Dissimilarity matrix (\code{NULL} for \code{FKM.noise})}
 \item{call}{Matched call}
} 
\details{
 If \code{startU} is given, the argument \code{k} is ignored (the number of clusters is \code{ncol(startU)}).
 \cr If \code{startU} is given, the first element of \code{value}, \code{cput} and \code{iter} refer to the rational start.
} 
\references{
 Dave' R.N., 1991. Characterization and detection of noise in clustering. Pattern Recognition Letters, 12, 657-664.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{FKM}}, \code{\link{Fclust}}, \code{\link{Fclust.index}}, \code{\link{print.fclust}}, \code{\link{summary.fclust}}, \code{\link{plot.fclust}}, \code{\link{butterfly}}}
\examples{
## butterfly data
data(butterfly)
## fuzzy k-means with noise cluster, fixing the number of clusters
clust=FKM.noise(butterfly, k = 2, RS=5,delta=3)
## fuzzy k-means with noise cluster, selecting the number of clusters
clust=FKM.noise(butterfly,RS=5,delta=3)
}
\keyword{multivariate}