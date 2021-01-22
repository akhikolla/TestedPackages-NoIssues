\name{is.modelKS-methods}
\docType{methods}

\alias{is.modelKS-methods}
\alias{is.modelKS}

\title{
Is an object from class \code{\linkS4class{modelKS}}?
}

\description{
Function to test inheritance relationships between an object and
     a class \code{\linkS4class{modelKS}}.
}

\arguments{
  \item{object}{Any \code{R} object.}
}

\section{Methods}{
\describe{

\item{\code{signature(object = "modelKS")}}{Is an \code{object} from class \code{\linkS4class{modelKS}}?}
}}


\value{
  TRUE if \code{object} is from class \code{\linkS4class{modelKS}}, and FALSE if not.
}

\author{
  Wilson Toussile
}

\examples{
data(exModelKS)
is.modelKS(exModelKS)
is.modelKS(1:7)
}
