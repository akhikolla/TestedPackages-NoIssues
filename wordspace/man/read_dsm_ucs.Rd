\name{read.dsm.ucs}
\alias{read.dsm.ucs}
\title{
  Load Raw DSM Data from Disk Files in UCS Export Format (wordspace)
}
\description{

  This function loads raw DSM data -- a cooccurrence frequency matrix and tables of marginal frequencies -- in \bold{UCS} export format. The data are read from a directory containing several text files with predefined names, which can optionally be compressed (see \sQuote{File Format} below for details).

}
\usage{

read.dsm.ucs(filename, encoding = getOption("encoding"), verbose = FALSE)

}
\arguments{

  \item{filename}{
    the name of a directory containing files with the raw DSM data.
  }

  \item{encoding}{
    character encoding of the input files, which will automatically be converted to R's internal representation if possible.  See \sQuote{Encoding} in \code{\link{file}} for details.
  }

  \item{verbose}{
    if \code{TRUE}, a few progress and information messages are shown
  }
}
%% \details{
%%   \bold{TODO}
%% }
\value{

  An object of class \code{\link{dsm}} containing a dense or sparse DSM.

  Note that the information tables for target terms (field \code{rows}) and feature terms (field \code{cols}) include the correct marginal frequencies from the UCS export files.  Nonzero counts for rows are and columns are added automatically unless they are already present in the disk files.  Additional fields from the information tables as well as all global variables are preserved with their original names.

}
\section{File Format}{

  The UCS export format is a directory containing the following files with the specified names:

  \itemize{
    \item \file{M} \strong{or} \file{M.mtx}

    cooccurrence matrix (dense, plain text) or sparse matrix (MatrixMarket format)
  
    \item \file{rows.tbl}
  
    row information (labels \code{term}, marginal frequencies \code{f})
  
    \item \file{cols.tbl}
  
    column information (labels \code{term}, marginal frequencies \code{f})

    \item \file{globals.tbl}

    table with single row containing global variables; must include variable \code{N} specifying sample size
  }

  Each individual file may be compressed with an additional filename extension \code{.gz}, \code{.bz2} or \code{.xz}; \code{read.dsm.ucs} automatically decompresses such files when loading them.

}
\seealso{

  \code{\link{dsm}}, \code{\link{read.dsm.triplet}}
  
}
\references{

  The UCS toolkit is a software package for collecting and manipulating co-occurrence data available from \url{http://www.collocations.de/software.html}.

  UCS relies on compressed text files as its main storage format.  They can be exported as a DSM with \code{ucs-tool export-dsm-matrix}.
  
}
%% \examples{
%% }
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
\keyword{file}
