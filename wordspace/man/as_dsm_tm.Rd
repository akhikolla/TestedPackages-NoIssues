\name{as.dsm.tm}
\alias{as.dsm.TermDocumentMatrix}
\alias{as.dsm.DocumentTermMatrix}
\title{
  Create DSM Object From \code{tm} Package (wordspace)
}
\description{

  Convert a \bold{tm} term-document or document-term matrix into a \code{wordspace} DSM object.

}
\usage{

\method{as.dsm}{TermDocumentMatrix}(obj, \dots, verbose=FALSE)
\method{as.dsm}{DocumentTermMatrix}(obj, \dots, verbose=FALSE)

}
\arguments{

  \item{obj}{
    an term-document or document-term matrix from the \bold{tm} package, i.e.
    an object of a class \code{\link[tm]{TermDocumentMatrix}} or \code{\link[tm]{DocumentTermMatrix}}.
  }
  
  \item{...}{
    additional arguments are ignored
  }
  
  \item{verbose}{
    if \code{TRUE}, a few progress and information messages are shown
  }
}
\value{

  An object of class \code{\link{dsm}}.

}
\seealso{

  \code{\link{as.dsm}} and the documentation of the \bold{tm} package

}
\examples{

\dontrun{
library(tm) # tm package needs to be installed
data(crude) # news messages on crude oil from Reuters corpus

cat(as.character(crude[[1]]), "\n") # a text example

corpus <- tm_map(crude, stripWhitespace) # some pre-processing
corpus <- tm_map(corpus, content_transformer(tolower))
corpus <- tm_map(corpus, removePunctuation)
corpus <- tm_map(corpus, removeWords, stopwords("english"))

cat(as.character(corpus[[1]]), "\n") # pre-processed text

dtm <- DocumentTermMatrix(corpus) # document-term matrix
inspect(dtm[1:5, 90:99])   # rows = documents

wordspace_dtm <- as.dsm(dtm, verbose=TRUE) # convert to DSM
print(wordspace_dtm$S[1:5, 90:99]) # same part of dtm as above

wordspace_tdm <- t(wordspace_dtm) # convert to term-document matrix
print(wordspace_tdm)
}

}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
