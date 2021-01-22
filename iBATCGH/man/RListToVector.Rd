\name{RListToVector}
\alias{RListToVector}
\title{
Internal function
}
\description{
Internal function
}
\usage{
RListToVector(xList, G, T)
}
\arguments{
  \item{xList}{
List of included associations (C format)
}
  \item{G}{
Number of gene expression probes
}
  \item{T}{
Number of CGH probes
}
}
\details{
This is an internal function.
}
\value{
Returns a vector whose elements are the number of times, across the MCMC iterations, the corresponding position of the transition matrix has been set to one.
}
\author{
Alberto Cassese
}
