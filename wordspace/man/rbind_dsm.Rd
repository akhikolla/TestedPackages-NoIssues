\name{rbind.dsm}
\alias{rbind.dsm}
\alias{cbind.dsm}
\title{
  Combine DSM Objects by Rows and Columns (wordspace)
}
\description{

  Combine conformable DSM matrices by rows or columns.
  Additional information in the DSM objects (such as marginal frequencies)
  is checked for consistency and updated automatically.

  \bold{Warning: these functions are experimental and may be removed or modified in a future release of \code{wordspace}}
  
}
\usage{
  \method{rbind}{dsm}(\dots, term.suffix=NULL, deparse.level = 1)
  
  \method{cbind}{dsm}(\dots, term.suffix=NULL, deparse.level = 1)
  
}
\arguments{

  \item{\dots}{
    one or more objects of class \code{dsm}, which must have the same feature dimensions (\code{rbind}) or target terms (\code{cbind})
  }

  \item{term.suffix}{
    optional character vector specifying one suffix string for each DSM, which will be appended to row (\code{rbind}) or column (\code{cbind}) labels in order to make them unique
  }
  
  \item{deparse.level}{
    ignored
  }

}
%% \details{
%% }
\value{
  if \code{term.suffix} is specified, row information of returned DSM object will be extended with variables \code{orig.term} specifying the original terms and \code{orig.part} specifying the original component model (identified by the corresponding entry from \code{term.suffix}) 
}
%% \seealso{
%% }
%% \examples{
%% }
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
