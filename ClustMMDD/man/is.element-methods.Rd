\name{is.element-methods}
\docType{methods}

\alias{is.element-methods}
\alias{is.element}

\title{
Check if a \code{\linkS4class{modelKS}} object is in a set of such objects.
}

\arguments{
  \item{el}{An instance of \code{\linkS4class{modelKS}} class.}

  \item{set}{A set of instances of \code{\linkS4class{modelKS}} class.}
}

\section{Methods}{
\describe{
  \item{\code{signature(el = c("modelKS"), set = c("modelKS"))}}{The two arguments must be vectors (see examples)}
}
}

\value{
  TRUE if the object \code{el} belongs to a given \code{set} of \code{\linkS4class{modelKS}}.
}

\description{
  Return TRUE if an instance of \code{\linkS4class{modelKS}} belongs to a set.
}

\author{
  Wilson Toussile
}

\examples{
data(exModelKS)
is.element(c(exModelKS), c(exModelKS))
is.element(c(exModelKS, 1, c(1:5)), c(exModelKS))
is.element(c(exModelKS), c(exModelKS, 1, list(1:5, 0)))
}

