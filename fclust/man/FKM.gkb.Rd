\name{FKM.gkb}
\alias{FKM.gkb}
\title{Gustafson, Kessel and Babuska - like fuzzy k-means}
\description{
 Performs the Gustafson, Kessel and Babuska - like fuzzy \emph{k}-means clustering algorithm.
 \cr Differently from fuzzy \emph{k}-means, it is able to discover non-spherical clusters.
 \cr The Babuska et al. variant improves the computation of the fuzzy covariance matrices in the standard Gustafson and Kessel clustering algorithm.
}
\usage{
 FKM.gkb (X, k, m, vp, gam, mcn, RS, stand, startU, index, alpha, conv, maxit, seed)
}
\arguments{
 \item{X}{Matrix or data.frame}
 \item{k}{An integer value or vector specifying the number of clusters for which the \code{index} is to be calculated (default: 2:6)}
 \item{m}{Parameter of fuzziness (default: 2)}
 \item{vp}{Volume parameter (default: rep(1,k))}
 \item{gam}{Weighting parameter for the fuzzy covariance matrices (default: 0)}
 \item{mcn}{Maximum condition number for the fuzzy covariance matrices (default: 1e+15)}
 \item{RS}{Number of (random) starts (default: 1)}
 \item{stand}{Standardization: if \code{stand=1}, the clustering algorithm is run using standardized data (default: no standardization)}
 \item{startU}{Rational start for the membership degree matrix \code{U} (default: no rational start)}
\item{index}{Cluster validity index to select the number of clusters: \code{PC} (partition coefficient), \code{PE} (partition entropy), \code{MPC} (modified partition coefficient), \code{SIL} (silhouette), \code{SIL.F} (fuzzy silhouette), \code{XB} (Xie and Beni) (default: "SIL.F")}
 \item{alpha}{Weighting coefficient for the fuzzy silhouette index \code{SIL.F} (default: 1)}
 \item{conv}{Convergence criterion (default: 1e-9)}
 \item{maxit}{Maximum number of iterations (default: 1e+2)}
 \item{seed}{Seed value for random number generation (default: NULL)}
}
\value{
 Object of class \code{fclust}, which is a list with the following components:
 \item{U}{Membership degree matrix}
 \item{H}{Prototype matrix}
 \item{F}{Array containing the covariance matrices of all the clusters}
 \item{clus}{Matrix containing the indexes of the clusters where the objects are assigned (column 1) and the associated membership degrees (column 2)}
 \item{medoid}{Vector containing the indexes of the medoid objects (\code{NULL} for \code{FKM.gkb})}
 \item{value}{Vector containing the loss function values for the \code{RS} starts}
 \item{criterion}{Vector containing the values of clustering index}
 \item{iter}{Vector containing the numbers of iterations for the \code{RS} starts}
 \item{k}{Number of clusters}
 \item{m}{Parameter of fuzziness}
 \item{ent}{Degree of fuzzy entropy (\code{NULL} for \code{FKM.gkb})}
 \item{b}{Parameter of the polynomial fuzzifier (\code{NULL} for \code{FKM.gkb})}
 \item{vp}{Volume parameter}
 \item{delta}{Noise distance (\code{NULL} for \code{FKM.gkb})}
 \item{gam}{Weighting parameter for the fuzzy covariance matrices}
 \item{mcn}{Maximum condition number for the fuzzy covariance matrices}
 \item{stand}{Standardization (Yes if \code{stand=1}, No if \code{stand=0})}
 \item{Xca}{Data used in the clustering algorithm (standardized data if \code{stand=1})}
 \item{X}{Raw data}
 \item{D}{Dissimilarity matrix (\code{NULL} for \code{FKM.gkb})}
 \item{call}{Matched call}
}
\details{
 If \code{startU} is given, the argument \code{k} is ignored (the number of clusters is \code{ncol(startU)}).
 \cr If \code{startU} is given, the first element of \code{value}, \code{cput} and \code{iter} refer to the rational start.
 \cr If a cluster covariance matrix becomes singular, then the algorithm stops and the element of \code{value} is {NaN}.
}
\references{
 Babuska R., van der Veen P.J., Kaymak U., 2002.  Improved covariance estimation for Gustafson-Kessel clustering. Proceedings of the IEEE International Conference on Fuzzy Systems (FUZZ-IEEE), 1081-1085.
 \cr Gustafson E.E., Kessel W.C., 1978. Fuzzy clustering with a fuzzy covariance matrix. Proceedings of the IEEE Conference on Decision and Control, pp. 761-766.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{FKM.gk}}, \code{\link{Fclust}}, \code{\link{Fclust.index}}, \code{\link{print.fclust}}, \code{\link{summary.fclust}}, \code{\link{plot.fclust}}, \code{\link{unemployment}}}
\examples{
\dontrun{
## unemployment data
data(unemployment)
## Gustafson, Kessel and Babuska-like fuzzy k-means, fixing the number of clusters
clust=FKM.gkb(unemployment,k=3,RS=10)
## Gustafson, Kessel and Babuska-like fuzzy k-means, selecting the number of clusters
clust=FKM.gkb(unemployment,k=2:6,RS=10)}
}
\keyword{multivariate}
