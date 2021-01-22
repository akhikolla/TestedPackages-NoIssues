\encoding{UTF-8}
\name{read.dsm.triplet}
\alias{read.dsm.triplet}
\title{
  Load DSM Data from Triplet Representation (wordspace)
}
\description{

  This function loads a sparse distributional semantic model in triplet representation -- (target label, feature label, score) -- from a disk file or a pipe.  Such a triplet file usually represents a pre-scored DSM, but it can also be used to read raw co-occurrence frequencies.  In this case, marginals and sample size can either be derived from the co-occurrence matrix (for syntactic and term-context models) or provided in separate TAB-delimited tables (for surface and textual co-occurrence, or if frequency thresholds have been applied).

}
\usage{

read.dsm.triplet(filename, freq = FALSE, value.first = FALSE, tokens = FALSE,
                 rowinfo = NULL, rowinfo.header = NULL,
                 colinfo = NULL, colinfo.header = NULL,
                 N = NA, span.size = 1,
                 sep = "\t", quote = "", nmax = -1, sort = FALSE, 
                 encoding = getOption("encoding"), verbose = FALSE) 

}
\arguments{

  \item{filename}{
    the name of a file containing the triplet data (see \sQuote{File Format} below for details), which may be compressed (\file{.gz}, \file{.bz2}, \file{.xz}).  If \code{filename} ends in \code{|}, it is opened as a Unix pipe for reading.
  }

  \item{freq}{
    whether values are raw co-occurrence frequencies (\code{TRUE}) or pre-computed scores (\code{FALSE})
  }

  \item{value.first}{
    if \code{TRUE}, triplets are given as (score, row label, column label) instead of the default (row label, column label, score)
  }

  \item{tokens}{
    if \code{TRUE}, the input file contains pair \emph{tokens}, i.e. row and column labels without score/frequency values.
    Co-occurrence frequencies will automatically be calculated, but this input format should only be used for small samples up to a few millon tokens.
  }

  \item{rowinfo}{
    the name of an optional TAB-delimited table file with additional information about the target terms (see \sQuote{File Format} below for details), which may be compressed (\file{.gz}, \file{.bz2}, \file{.xz}).
  }
  
  \item{rowinfo.header}{
    if the \code{rowinfo} file does not start with a header row, specify its column names as a character vector here
  }

  \item{colinfo}{
    the name of an optional TAB-delimited table file with additional information about the feature terms or contexts (see \sQuote{File Format} below for details), which may be compressed (\file{.gz}, \file{.bz2}, \file{.xz}).
  }

  \item{colinfo.header}{
    if the \code{colinfo} file does not start with a header row, specify its column names as a character vector here
  }

  \item{N}{
    sample size to assume for the distributional model (see \sQuote{Details} below)
  }

  \item{span.size}{
    if marginal frequencies are provided externally for surface co-occurrence, they need to be adjusted for span size.
    If this hasn't been taken into account in data extraction, it can be approximated by specifying the total number of tokens in a span here (see \sQuote{Details} below).
  }

  \item{sep, quote}{
    specify field separator and the types of quotes used by the disk file (see the \code{\link{scan}} documentation for details).  By default, a TAB-delimited file without quotes is assumed.
  }

  \item{nmax}{
    if the number of entries (= text lines) in the triplet file is known, it can be specified here in order to make loading faster and more memory-efficient.  Caution: If \code{nmax} is smaller than the number of lines in the disk file, the extra lines will silently be discarded.
  }

  \item{sort}{
    if \code{TRUE}, the rows and columns of the co-occurrence matrix will be sorted alphabetically according to their labels (i.e. the target and feature terms); otherwise they are listed as encountered in the triplet representation
  }

  \item{encoding}{
    character encoding of the input files, which will automatically be converted to \R's internal representation if possible.  See \sQuote{Encoding} in \code{\link{file}} for details.
  }

  \item{verbose}{
    if \code{TRUE}, a few progress and information messages are shown
  }

}
\details{

  The function \code{read.dsm.triplet} can be used to read triplet representations of three different types of DSM.
  
  \subsection{1. A pre-scored DSM matrix}{
    If \code{freq=FALSE} and \code{tokens=FALSE}, the triplet file is assumed to contain pre-scored entries of the DSM matrix.
    Marginal frequencies are not required for such a model, but additional information about targets and features can be provided in separate \code{rowinfo=} and \code{colinfo=} files.
  }

  \subsection{2. Raw co-occurrence frequencies (syntactic or term-context)}{
    If the triplet file contains syntactic co-occurrence frequencies or term-document frequency counts, specify \code{freq=TRUE}. For small data sets, frequencies can also be aggregated directly in \R from co-occurrence tokens; specify \code{tokens=TRUE}.
    
    Unless high frequency thresholds or other selective filters have been applied to the input data, the marginal frequencies of targets and features as well as the sample size can automatically be derived from the co-occurrence matrix.  \emph{Do not specify \code{rowinfo=} or \code{colinfo=} in this case!}

    Evert (2008) explains the differences between syntactic, textual and surface co-occurrence.
  }

  \subsection{3. Raw co-occurrence frequencies with explicit marginals}{
    For surface and textual co-occurrence data, the correct marginal frequencies cannot be derived automatically and have to be provided in auxiliary table files specified with \code{rowinfo=} and \code{colinfo}.  These files must contain a column \code{f} with the marginal frequency data.  In addition, the total sample size (which cannot be derived from the marginals) has to be passed in the argument \code{N=}.  Of course, it is still necessary to specify \code{freq=TRUE} (or \code{token=TRUE}) in order to indicate that the input data aren't pre-computed scores.
    
    The computation of consistent marginal frequencies is particulary tricky for surface co-occurrence (Evert 2008, p. 1233f) and specialized software should be used for this purpose.  As an approximation, simple corpus frequencies of target and feature terms can be corrected by a factor corresponding to the total size of the collocational span (e.g. \code{span.size=8} for a symmetric L4/R4 span, cf. Evert 2008, p. 1225).  The \code{read.dsm.triplet} function applies this correction to the row marginals.
    
    Explicit marginals should also be provided if syntactic co-occurrence data or text-context frequencies have been filtered, either individually with a frequency threshold or by selecting a subset of the targets and features.  See the examples below for an illustration.
  }

}

\value{

  An object of class \code{\link{dsm}} containing a sparse DSM.
  
  For a model of type 1 (pre-scored) it will include the score matrix \code{$S} but no co-occurrence frequency data. Such a DSM object cannot be passed to \code{\link{dsm.score}}, except with \code{score="reweight"}.  For models of type 2 and 3 it will include the matrix of raw co-occurrence frequencies \code{$M}, but no score matrix.

}
\section{File Format}{

  \subsection{Triplet files}{
    The triplet file must be a plain-text table with two or three TAB-delimited columns and no header.  It may be compressed in \code{.gz}, \code{.bz2} or \code{.xz} format.
  
    For \code{tokens=TRUE}, each line represents a single pair token with columns
    \enumerate{
      \item target term
      \item feature term / context
    }
  
    For \code{tokens=FALSE}, each line represents a pair type (i.e. a unique cell of the co-occurrence matrix) with columns:
    \enumerate{
      \item target term
      \item feature term / context
      \item score (\code{freq=FALSE}) \emph{or} co-occurrence frequency (\code{freq=TRUE})
    }
    If \code{value.first=TRUE}, the score entry is expected in the first column:
    \enumerate{
      \item score \emph{or} co-occurrence frequency
      \item target term
      \item feature term / context
    }
    Note that the triplet file may contain multiple entries for the same cell, whose values will automatically be added up.
    This might not be very sensible for pre-computed scores.
  }
  
  \subsection{Row and column information}{
    Additional information about target terms (matrix rows) and feature terms / contexts (matrix columns) can be provided in additional TAB-delimited text tables, optionally compressed in \code{.gz}, \code{.bz2} or \code{.xz} format.
    
    Such tables can have an arbitrary number of columns whose data types are inferred from the first few rows of the table.
    Tables should start with a header row specifying the column labels; otherwise they must be passed in the \code{rowinfo.header} and \code{colinfo.header} arguments.
    
    Every table must contain a column \code{term} listing the target terms or feature terms / contexts.  Their ordering need not be the same as in the main co-occurrence matrix, and redundant entries will silently be dropped.
    
    If \code{freq=TRUE} or \code{tokens=TRUE}, the tables must also contain marginal frequencies in a column \code{f}.  Nonzero counts for rows and columns of the matrix are automatically added unless a column \code{nnzero} is already present.
  }
  
}
\seealso{

  \code{\link{dsm}}, \code{\link{read.dsm.ucs}}
  
}
\references{
  
  Evert, Stefan (2008). Corpora and collocations.
  In A. Lüdeling and M. Kytö (eds.), \emph{Corpus Linguistics. An International Handbook}, chapter 58, pages 1212--1248. Mouton de Gruyter, Berlin, New York.

}
\examples{
## this helper function displays the cooccurrence matrix together with marginals
with.marginals <- function (x) {
  y <- x$M
  rownames(y) <- with(x$rows, sprintf("\%-8s | \%6d", term, f))
  colnames(y) <- with(x$cols, sprintf("  \%s | \%d", term, f))
  y
}

## we will read this term-context DSM from a triplet file included in the package
with.marginals(DSM_TermContext)

## the triplet file with term-document frequencies
triplet.file <- system.file("extdata", "term_context_triplets.gz", package="wordspace")
cat(readLines(triplet.file), sep="\n") # file format

## marginals incorrect because matrix covers only subset of targets & features
TC1 <- read.dsm.triplet(triplet.file, freq=TRUE)
with.marginals(TC1) # marginal frequencies far too small

## TAB-delimited file with marginal frequencies and other information
marg.file <- system.file("extdata", "term_context_marginals.txt.gz", package="wordspace")
cat(readLines(marg.file), sep="\n") # notice the header row with "term" and "f"

## single table with marginals for rows and columns, but has to be specified twice
TC2 <- read.dsm.triplet(triplet.file, freq=TRUE, 
                        rowinfo=marg.file, colinfo=marg.file, N=108771103)
with.marginals(TC2) # correct marginal frequencies

## marginals table without header: specify column lables separately
no.hdr <- system.file("extdata", "term_context_marginals_noheader.txt", 
                      package="wordspace")
hdr.names <- c("term", "f", "df", "type")
TC3 <- read.dsm.triplet(triplet.file, freq=TRUE, 
                        rowinfo=no.hdr, rowinfo.header=hdr.names,
                        colinfo=no.hdr, colinfo.header=hdr.names, N=108771103)
all.equal(TC2, TC3, check.attributes=FALSE) # same result
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
\keyword{file}
