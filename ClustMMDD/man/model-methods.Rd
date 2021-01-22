\name{model-methods}
\docType{methods}

\alias{model-methods}
\alias{model}

\title{
  Retrieve a list of model  \eqn{\left(K,S\right)} from a \code{\linkS4class{modelKS}} object.
}

\section{Methods}{
\describe{
\item{\code{signature(object = "modelKS")}}{Retrieve a list of model  \eqn{\left(K,S\right)} from a \code{\linkS4class{modelKS}} object.}
}}

\description{
Recall that a model is defined by \eqn{\left(K,S\right)} where \eqn{K} is the number of clusters and \eqn{S} that indicates the clustering variable. This method retrieves a list of model  \eqn{\left(K,S\right)} from a \code{\linkS4class{modelKS}} object.
}


\seealso{
  \code{\linkS4class{modelKS}}, \code{\link{slotNames}}, \code{\link{new}}, \code{\link{methods}}, \code{\link{show}}
}

\author{
  Wilson Toussile.
}

\examples{
data(exModelKS)
showClass("modelKS")
slotNames("modelKS")
exModelKS
exModelKS["K"]
exModelKS["S"]
model(exModelKS)
}

\keyword{methods}

