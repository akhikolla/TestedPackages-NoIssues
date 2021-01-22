\name{as.dsm}
\alias{as.dsm}
\title{
  Create DSM Object From Various R Data Structures (wordspace)
}
\description{

  Convert co-occurrence data from various in-memory formats to DSM object.

}
\usage{

as.dsm(obj, \dots)

}
\arguments{

  \item{obj}{an object of a suitable class (for which a method implementation is available)}
  
  \item{...}{additional arguments passed on to the method implementation
    (see respective manpages for details)}
}
\value{

  An object of class \code{\link{dsm}}.

}
\seealso{

  Currently available implementations: \code{\link{as.dsm.TermDocumentMatrix}}, \code{\link{as.dsm.DocumentTermMatrix}}
  
}
%% \examples{
%% }
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
\keyword{ methods }
