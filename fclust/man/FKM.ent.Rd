\name{FKM.ent} 
\alias{FKM.ent}
\title{Fuzzy k-means with entropy regularization}
\description{
 Performs the fuzzy \emph{k}-means clustering algorithm with entropy regularization.
 \cr The entropy regularization allows us to avoid using the artificial fuzziness parameter \emph{m}. This is replaced by the degree of fuzzy entropy \emph{ent}, related to the concept of temperature in statistical physics.
 An interesting property of the fuzzy \emph{k}-means with entropy regularization is that the prototypes are obtained as weighted means with weights equal to the membership degrees (rather than to the membership degrees at
 the power of \emph{m} as is for the fuzzy k-means).
} 
\usage{ 
 FKM.ent (X, k, ent, RS, stand, startU, index, alpha, conv, maxit, seed)
}
\arguments{
 \item{X}{Matrix or data.frame}
 \item{k}{An integer value or vector specifying the number of clusters for which the \code{index} is to be calculated (default: 2:6)}
 \item{ent}{Degree of fuzzy entropy (default: 1)}
 \item{RS}{Number of (random) starts (default: 1)}
 \item{stand}{Standardization: if \code{stand=1}, the clustering algorithm is run using standardized data (default: no standardization)}
 \item{startU}{Rational start for the membership degree matrix \code{U} (default: no rational start)}
 \item{index}{Cluster validity index to select the number of clusters: \code{"PC"} (partition coefficient), \code{"PE"} (partition entropy), \code{"MPC"} (modified partition coefficient), \code{"SIL"} (silhouette), \code{"SIL.F"} (fuzzy silhouette), \code{"XB"} (Xie and Beni) (default: "SIL.F")}
 \item{alpha}{Weighting coefficient for the fuzzy silhouette index \code{SIL.F} (default: 1)}
 \item{conv}{Convergence criterion (default: 1e-9)}
 \item{maxit}{Maximum number of iterations (default: 1e+6)}
 \item{seed}{Seed value for random number generation (default: NULL)}
}
\value{
 Object of class \code{fclust}, which is a list with the following components:
 \item{U}{Membership degree matrix}
 \item{H}{Prototype matrix}
 \item{F}{Array containing the covariance matrices of all the clusters (\code{NULL} for \code{FKM.ent})}
 \item{clus}{Matrix containing the indexes of the clusters where the objects are assigned (column 1) and the associated membership degrees (column 2)}
 \item{medoid}{Vector containing the indexes of the medoid objects (\code{NULL} for \code{FKM.ent})}
 \item{value}{Vector containing the loss function values for the \code{RS} starts}
 \item{criterion}{Vector containing the values of the cluster validity index}
 \item{iter}{Vector containing the numbers of iterations for the \code{RS} starts}
 \item{k}{Number of clusters}
 \item{m}{Parameter of fuzziness (\code{NULL} for \code{FKM.ent})}
 \item{ent}{Degree of fuzzy entropy}
 \item{b}{Parameter of the polynomial fuzzifier (\code{NULL} for \code{FKM.ent})}
 \item{vp}{Volume parameter (\code{NULL} for \code{FKM.ent})} 
 \item{delta}{Noise distance (\code{NULL} for \code{FKM.ent})}
 \item{gam}{Weighting parameter for the fuzzy covariance matrices (\code{NULL} for \code{FKM.ent})}
 \item{mcn}{Maximum condition number for the fuzzy covariance matrices (\code{NULL} for \code{FKM.ent})}
 \item{stand}{Standardization (Yes if \code{stand=1}, No if \code{stand=0})}
 \item{Xca}{Data used in the clustering algorithm (standardized data if \code{stand=1})}
 \item{X}{Raw data}
  \item{D}{Dissimilarity matrix (\code{NULL} for \code{FKM.ent})}
 \item{call}{Matched call}
} 
\details{
 If \code{startU} is given, the argument \code{k} is ignored (the number of clusters is \code{ncol(startU)}).
 \cr If \code{startU} is given, the first element of \code{value}, \code{cput} and \code{iter} refer to the rational start.
 \cr The default value for \code{ent} is in general not reasonable if \code{FKM.ent} is run using raw data. 
 \cr The update of the membership degrees requires the computation of exponential functions. In some cases, this may produce \code{NaN} values and the algorithm stops. Such a problem is usually solved by running \code{FKM.ent} using standardized data (\code{stand=1}).
} 
\references{
 Li R., Mukaidono M., 1995. A maximum entropy approach to fuzzy clustering. Proceedings of the Fourth IEEE Conference on Fuzzy Systems (FUZZ-IEEE/IFES '95), pp. 2227-2232.
 \cr Li R., Mukaidono M., 1999. Gaussian clustering method based on maximum-fuzzy-entropy interpretation. Fuzzy Sets and Systems, 102, 253-258.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{FKM.ent.noise}}, \code{\link{Fclust}}, \code{\link{Fclust.index}}, \code{\link{print.fclust}}, \code{\link{summary.fclust}}, \code{\link{plot.fclust}}, \code{\link{Mc}}}
\examples{
## McDonald's data
data(Mc)
names(Mc)
## data normalization by dividing the nutrition facts by the Serving Size (column 1)
for (j in 2:(ncol(Mc)-1))
Mc[,j]=Mc[,j]/Mc[,1]
## removing the column Serving Size
Mc=Mc[,-1]
## fuzzy k-means with entropy regularization, fixing the number of clusters
## (excluded the factor column Type (last column))
clust=FKM.ent(Mc[,1:(ncol(Mc)-1)],k=6,ent=3,RS=10,stand=1)
## fuzzy k-means with entropy regularization, selecting the number of clusters
## (excluded the factor column Type (last column))
clust=FKM.ent(Mc[,1:(ncol(Mc)-1)],k=2:6,ent=3,RS=10,stand=1)
}

\keyword{multivariate}