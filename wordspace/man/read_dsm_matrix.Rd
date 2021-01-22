\name{read.dsm.matrix}
\alias{read.dsm.matrix}
\title{
  Load DSM Matrix from File (wordspace)
}
\description{

  This function loads a DSM matrix from a disk file in the specified format (see section sQuote(Formats) for details).

}
\usage{

read.dsm.matrix(file, format = c("word2vec"),
                encoding = "UTF-8", batchsize = 1e6, verbose=FALSE)

}
\arguments{

  \item{file}{
    either a character string naming a file or a \code{\link{connection}} open for writing (in text mode)
  }

  \item{format}{
    input file format (see section sQuote(Formats)). The input file format cannot be guessed automatically.
  }
  
  \item{encoding}{
    character encoding of the input file (ignored if \code{file} is a connection)
  }

  \item{batchsize}{
    for certain input formats, the matrix is read in batches of \code{batchsize} cells each in order to limit memory overhead
  }

  \item{verbose}{
    if \code{TRUE}, show progress bar when reading in batches
  }
}
\details{

  In order to read text formats from a compressed file, pass a \code{\link{gzfile}}, \code{\link{bzfile}} or \code{\link{xzfile}} connection with appropriate \code{encoding} in the argument \code{file}.  Make sure not to open the connection before passing it to \code{read.dsm.matrix}.

}
\section{Formats}{

  Currently, the only supported file format is \code{word2vec}.

  \describe{
    \item{\code{word2vec}}{
      This widely used text format for word embeddings is only suitable for a dense matrix. Row labels must be unique and may not contain whitespace.
      Values are usually rounded to a few decimal digits in order to keep file size manageable.
      
      The first line of the file lists the matrix dimensions (rows, columns) separated by a single blank.
      It is followed by one text line for each matrix row, starting with the row label. The label and are cells are separated by single blanks, so row labels cannot contain whitespace.
    }
  }

}
\seealso{
  \code{\link{write.dsm.matrix}}, \code{\link{read.dsm.triplet}}, \code{\link{read.dsm.ucs}}
}
\examples{
fn <- system.file("extdata", "word2vec_hiero.txt", package="wordspace")
read.dsm.matrix(fn, format="word2vec")
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
\keyword{file}
