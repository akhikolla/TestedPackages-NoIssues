\encoding{utf-8}
\name{ClusterStability-package}
\alias{ClusterStability-package}
\docType{package}

\title{Assessment of the stability of individual objects, clusters and a whole clustering solution based on repeated runs of a clustering algorithm.}

\description{The ClusterStability package uses a probabilistic framework and some well-known clustering
 criteria (e.g. Calinski-Harabasz, Silhouette, Dunn and Davies-Bouldin) to compute the stability scores (\emph{ST})
  of each individual object (i.e., \emph{element}) in the clustering solution provided by the K-means and K-medoids
  partitioning algorithms.
}

\details{
	\tabular{ll}{
		Package: \tab ClusterStability\cr
		Type: \tab Package\cr
		Version: \tab 1.0.2\cr
		Date: \tab 2015-10-14\cr
		License: \tab GPL-2\cr
		Maintainer: \tab Etienne Lord <m.etienne.lord@gmail.com>, \cr \tab Vladimir Makarenkov <makarenkov.vladimir@uqam.ca>
	}

	Function \code{\link{ClusterStability}} computes the individual and global stability scores (\emph{ST}) for a partitioning solution using either K-means or K-medoids (the approximate solution is provided).
	
	Function \code{\link{ClusterStability_exact}} is similar to the \code{\link{ClusterStability}} function but uses the Stirling numbers of the second kind to compute the exact stability scores (but is limited to a small number of objects).

	Function \code{\link{Kcombination}} computes the \emph{k}-combination of a set of numbers for a given \emph{k}.
	
	Function \code{\link{Reorder}} returns the re-ordered partitioning of a series of clusters.

	Function \code{\link{Stirling2nd}} computes the Stirling numbers of the second kind.

}

\author{Etienne Lord, François-Joseph Lapointe and Vladimir Makarenkov}



%\references{}

\keyword{package}
\keyword{stability score}
\keyword{partitioning criteria}

\seealso{
\code{\link{ClusterStability}},
\code{\link{ClusterStability_exact}},
\code{\link{Kcombination}},
\code{\link{Reorder}},
\code{\link{Stirling2nd}}
}

%\examples{}
