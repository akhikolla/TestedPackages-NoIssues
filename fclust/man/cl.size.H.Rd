\name{cl.size.H}
\alias{cl.size.H}
\title{Cluster size}
\description{
 Produces the sizes of the clusters in the hard clustering sense (objects are considered to be assigned to clusters only if the corresponding membership degree are >=0.5).
}
\usage{
 cl.size.H (U)
}
\arguments{
 \item{U}{Membership degree matrix}
}
\value{
 \item{clus.size}{Vector containing the sizes of the clusters}
}
\details{
 An object is assigned to a cluster according to the maximal membership degree provided that such a maximal membership degree is >=0.5, otherwise it is assumed that an object is not assigned to any cluster.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{cl.size}}}
\examples{
n=20
k=3
## randomly generated membership degree matrix
U=matrix(runif(n*k,0,1), nrow=n, ncol=k)
U=U/apply(U,1,sum)
## cluster size in the hard clustering sense
clus.size=cl.size.H(U)
}
\keyword{multivariate}
