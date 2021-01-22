\encoding{ISO-8859-2}
\name{ClusterStability_exact}
\alias{ClusterStability_exact}
\title{Calculates the exact stability score (\emph{ST}) for individual objects in a clustering solution.}
\description{This function will return the exact individual stability score \emph{ST} and the exact global score \emph{STglobal} using either the K-means or K-medoids algorithm and four different clustering indices: Calinski-Harabasz, Silhouette, Dunn or Davies-Bouldin. \strong{Variable overflow errors are possible for large numbers of objects}.}
\usage{ClusterStability_exact(dat, k, replicate, type) }
\arguments{
	\item{dat}{the input dataset: either a matrix or a dataframe.}
    \item{k}{the number of classes for the K-means or K-medoids algorithm (default=3).}
    \item{replicate}{the number of replicates to perform (default=1000).}
	\item{type}{the algorithm used in the partitioning: either 'kmeans' or 'kmedoids' algorithm (default=kmeans).}
}

\value{Returns the exact individual (\emph{ST}) and global (\emph{ST_global}) stability scores for the four clustering indices: Calinski-Harabasz (\emph{ch}), Silhouette (\emph{sil}), Dunn (\emph{dunn}) or Davies-Bouldin (\emph{db}).}

\examples{
   ## Calculate the stability scores of individual objects of the Iris dataset
   ## using K-means, 100 replicates (random starts) and k=3
   ClusterStability_exact(dat=iris[1:4],k=3,replicate=100,type='kmeans');
}

\keyword{Stability score,ST,individual,global,exact}

