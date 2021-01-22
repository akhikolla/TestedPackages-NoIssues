\name{DSM_SingularValues}
\alias{DSM_SingularValues}
\docType{data}
\title{
  Typical Singular Values of a Term-Context Matrix (wordspace)
}
\description{

  Typical singular values of a term-document matrix based on encyclopedia articles.
  
}
\usage{
  
DSM_SingularValues

}
\format{

  A numeric vector of length 2623.
  
}
\details{

  The data were obtained by singular value decomposition of a term-document matrix representing 100,000 Wikipedia articles, with 2623 target terms from a Basic English vocabulary.  Articles were truncated to the first ca. 500 words.  Occurrence frequencies of the target terms were log-scaled and rows of the matrix were L2-normalized before applying the SVD.
  
}
\examples{

\dontrun{
plot(DSM_SingularValues, type="h", xaxs="i", yaxs="i")
}

}
\keyword{datasets}
