\name{simdata1}
\alias{simdata1}
\alias{simdata2}
\alias{simdata3}
\alias{simdata4}
\alias{simdata5}
\docType{data}
\title{
Simulated data.
}
\description{
Five simulated data sets are included. 
}
\usage{
data("simdata1")
data("simdata2")
data("simdata3")
data("simdata4")
data("simdata5")
}
\format{
Each data set contains a list that comprises a tree (called "tree") and data (called "data") as its components. The tree is in the ape package phylo format. The data component consists of a matrix of discrete character patterns with the different patterns as the rows and the taxa as the columns. 
}
\details{
These data were simulated with either the geiger or the phangorn packages.
}
\references{
Schliep K.P. 2011. phangorn: phylogenetic analysis in R. Bioinformatics, 27(4)  592--593. R package version 1.99-12.

Harmon Luke J, Jason T Weir, Chad D Brock, Richard E Glor, and Wendell Challenger. 2008. GEIGER: investigating evolutionary radiations. Bioinformatics 24:129-131. R package version 2.0.3.
}
\examples{
data(simdata1)
data(simdata2)
data(simdata3)
data(simdata4)
data(simdata5)
}
\keyword{datasets}