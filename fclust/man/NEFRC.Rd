\name{NEFRC}
\alias{NEFRC}
\title{
Non-Euclidean Fuzzy Relational Clustering}
\description{
Performs the Non-Euclidean Fuzzy Relational data Clustering algorithm.}
\usage{
NEFRC(D, k, m, RS, startU, index, alpha, conv, maxit, seed)
}
\arguments{
 \item{D}{Matrix or data.frame containing distances/dissimilarities}
 \item{k}{An integer value or vector specifying the number of clusters for which the \code{index} is to be calculated (default: 2:6)}
  \item{m}{Parameter of fuzziness (default: 2)}
 \item{RS}{Number of (random) starts (default: 1)}
 \item{startU}{Rational start for the membership degree matrix \code{U} (default: no rational start)}
 \item{conv}{Convergence criterion (default: 1e-9)}
 \item{index}{Cluster validity index to select the number of clusters: \code{"PC"} (partition coefficient), \code{"PE"} (partition entropy), \code{"MPC"} (modified partition coefficient), \code{"SIL"} (silhouette), \code{"SIL.F"} (fuzzy silhouette) (default: "SIL.F")}
 \item{alpha}{Weighting coefficient for the fuzzy silhouette index \code{SIL.F} (default: 1)}
 \item{maxit}{Maximum number of iterations (default: 1e+6)}
 \item{seed}{Seed value for random number generation (default: NULL)}
}

\value{
 Object of class \code{fclust}, which is a list with the following components:
 \item{U}{Membership degree matrix}
 \item{H}{Prototype matrix (\code{NULL} for \code{NEFRC})}
 \item{F}{Array containing the covariance matrices of all the clusters (\code{NULL} for \code{NEFRC}).}
 \item{clus}{Matrix containing the indexes of the clusters where the objects are assigned (column 1) and the associated membership degrees (column 2)}
 \item{medoid}{Vector containing the indexes of the medoid objects (\code{NULL} for \code{NEFRC}).}
 \item{value}{Vector containing the loss function values for the \code{RS} starts}
 \item{criterion}{Vector containing the values of the cluster validity index}
 \item{iter}{Vector containing the numbers of iterations for the \code{RS} starts}
 \item{k}{Number of clusters}
 \item{m}{Parameter of fuzziness}
 \item{ent}{Degree of fuzzy entropy (\code{NULL} for \code{NEFRC})}
 \item{b}{Parameter of the polynomial fuzzifier (\code{NULL} for \code{NEFRC})}
 \item{vp}{Volume parameter (\code{NULL} for \code{NEFRC})}
 \item{delta}{Noise distance (\code{NULL} for \code{NEFRC})}
 \item{stand}{Standardization (Yes if \code{stand=1}, No if \code{stand=0}) (\code{NULL} for \code{NEFRC})}
 \item{Xca}{Data used in the clustering algorithm (\code{NULL} for \code{NEFRC}, \code{D} is used)}
 \item{X}{Raw data (\code{NULL} for \code{NEFRC})}
 \item{D}{Dissimilarity matrix}
 \item{call}{Matched call}
}
\details{
 If \code{startU} is given, the argument \code{k} is ignored (the number of clusters is \code{ncol(startU)}).
 \cr If \code{startU} is given, the first element of \code{value}, \code{cput} and \code{iter} refer to the rational start.
}
\references{
Davé, R. N., & Sen, S. 2002. Robust fuzzy clustering of relational data. IEEE Transactions on Fuzzy Systems, 10(6), 713-727.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}

\seealso{
\code{\link{NEFRC.noise}}, \code{\link{print.fclust}}, \code{\link{summary.fclust}}, \code{\link{plot.fclust}}}
\examples{
\dontrun{
require(cluster)
data("houseVotes")
X <- houseVotes[,-1]
D <- daisy(x = X, metric = "gower")
clust.NEFRC <- NEFRC(D = D, k = 2:6, m = 2, index = "SIL.F")
summary(clust.NEFRC)
plot(clust.NEFRC)
}
}
\keyword{multivariate}

