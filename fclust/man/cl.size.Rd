\name{cl.size}
\alias{cl.size}
\title{Cluster size}
\description{
 Produces the sizes of the clusters.
}
\usage{
 cl.size (U)
}
\arguments{
 \item{U}{Membership degree matrix}
}
\value{
 \item{clus.size}{Vector containing the sizes of the clusters}
}
\details{
 An object is assigned to a cluster according to the maximal membership degree.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{cl.size.H}}}
\examples{
n=20
k=3
## randomly generated membership degree matrix
U=matrix(runif(n*k,0,1), nrow=n, ncol=k)
U=U/apply(U,1,sum)
clus.size=cl.size(U)
}
\keyword{multivariate}
