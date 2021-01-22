\name{write.dsm.matrix}
\alias{write.dsm.matrix}
\title{
  Export DSM Matrix to File (wordspace)
}
\description{

  This function exports a DSM matrix to a disk file in the specified format (see section sQuote(Formats) for details).

}
\usage{

write.dsm.matrix(x, file, format = c("word2vec"), round=FALSE,
                 encoding = "UTF-8", batchsize = 1e6, verbose=FALSE)

}
\arguments{

  \item{x}{
    a dense or sparse matrix representing a DSM, or an object of class \code{dsm}
  }

  \item{file}{
    either a character string naming a file or a \code{\link{connection}} open for writing (in text mode)
  }

  \item{format}{
    desired output file format. See section sQuote(Formats) for a list of available formats and their limitations.
  }
  
  \item{round}{
    for some output formats, numbers can be rounded to the specified number of decimal digits in order to reduce file size
  }

  \item{encoding}{
    character encoding of the output file (ignored if \code{file} is a connection)
  }

  \item{batchsize}{
    for certain output formats, the matrix is written in batches of \code{batchsize} cells each in order to limit memory overhead
  }
  
  \item{verbose}{
    if \code{TRUE}, show progress bar when writing in batches
  }
}
\details{

  In order to save text formats to a compressed file, pass a \code{\link{gzfile}}, \code{\link{bzfile}} or \code{\link{xzfile}} connection with appropriate \code{encoding} in the argument \code{file}.  Make sure not to open the connection before passing it to \code{write.dsm.matrix}.  See section \sQuote{Examples} below.

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
  \code{\link{read.dsm.matrix}}
}
\examples{
model <- dsm.score(DSM_TermTerm, score="MI", normalize=TRUE) # a typical DSM

# save in word2vec text format (rounded to 3 digits)
fn <- tempfile(fileext=".txt")
write.dsm.matrix(model, fn, format="word2vec", round=3)
cat(readLines(fn), sep="\n")

# save as compressed file in word2vec format
fn <- tempfile(fileext=".txt.gz")
fh <- gzfile(fn, encoding="UTF-8") # need to set file encoding here
write.dsm.matrix(model, fh, format="word2vec", round=3)
# write.dsm.matrix() automatically opens and closes the connection
cat(readLines(gzfile(fn)), sep="\n")
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
\keyword{file}
