\name{cl.memb.t}
\alias{cl.memb.t}
\title{Cluster membership}
\description{
 Produces a summary of the membership degree information according to a threshold.
}
\usage{
 cl.memb.t (U, t)
}
\arguments{
 \item{U}{Membership degree matrix}
 \item{t}{Threshold in [0,1] (default: 0)}
}
\value{
 \item{info.U}{Matrix containing the indexes of the clusters where the objects are assigned (row 1) and the associated membership degrees (row 2)}
}
\details{
 An object is assigned to a cluster according to the maximal membership degree provided that such a maximal membership degree is >= \code{t}, otherwise it is assumed that an object is not assigned to any cluster (denoted by cluster index = 0 in row 1).
 The function can be useful to select the subset of objects clearly assigned to clusters (objects with maximal membership degrees >= \code{t}).
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{cl.memb}}, \code{\link{cl.memb.H}}}
\examples{
n=20
k=3
## randomly generated membership degree matrix
U=matrix(runif(n*k,0,1), nrow=n, ncol=k)
U=U/apply(U,1,sum)
## threshold t=0.6
info.U=cl.memb.t(U,0.6)
## objects clearly assigned to clusters
rownames(info.U[info.U[,1]!=0,])
}
\keyword{multivariate}
