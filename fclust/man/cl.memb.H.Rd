\name{cl.memb.H}
\alias{cl.memb.H}
\title{Cluster membership}
\description{
 Produces a summary of the membership degree information in the hard clustering sense (objects are considered to be assigned to clusters only if the corresponding membership degree are >=0.5).
}
\usage{
 cl.memb.H (U)
}
\arguments{
 \item{U}{Membership degree matrix}
  }
\value{
 \item{info.U}{Matrix containing the indexes of the clusters where the objects are assigned (row 1) and the associated membership degrees (row 2)}
}
\details{
 An object is assigned to a cluster according to the maximal membership degree provided that such a maximal membership degree is >=0.5, otherwise it is assumed that an object is not assigned to any cluster (denoted by cluster index = 0 in row 1).
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{cl.memb}}, \code{\link{cl.memb.t}}}
\examples{
n=20
k=3
## randomly generated membership degree matrix
U=matrix(runif(n*k,0,1), nrow=n, ncol=k)
U=U/apply(U,1,sum)
info.U=cl.memb.H(U)
## objects assigned to clusters in the hard clustering sense
rownames(info.U[info.U[,1]!=0,])
}
\keyword{multivariate}
