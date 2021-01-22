\name{genotype2}
\alias{genotype2}
\docType{data}

\title{A genotype data frame compatible with \code{\link{ClustMMDD}} main functions.
}

\description{
  This data frame can be obtained using \code{\link{cutEachCol}}\code{(genotype1[, -11], ploidy = 2)} (see \code{\link{genotype1}}).
}

\usage{data(genotype2)}

\format{
  The format is:
 chr [1:1000, 1:20] "109" "107" "105" "101" "106" "106" "107" ..., representing observed alleles for the considered \eqn{10} loci, \eqn{2} 
 column per locus.
}
\source{
  Simulated data.
}

\details{
  \code{ploidy = 2} for diploid individual.
}

\author{
  Wilson Toussile
}

\seealso{
  \code{\link{genotype1}}.
}

\examples{
data(genotype2)
head(genotype2)
data(genotype1)
genotype3 = cutEachCol(genotype1[,-11], ploidy = 2)
head(genotype3)
}
\keyword{Genotype}
\keyword{Genotype data}
\keyword{Alleles}
