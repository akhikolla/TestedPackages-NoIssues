\name{==-methods}
\docType{methods}

\alias{==-methods}
\alias{==}
\title{ Methods for Function \code{==} }

\description{
  Check if two objects of class \code{\linkS4class{modelKS}} are equal.
}

\section{Methods}{
\describe{

\item{\code{signature(e1 = "modelKS", e2 = "modelKS")}}{
  
}
}}

\author{
Wilson Toussile.
}

\seealso{
  \code{\link{slotNames}}, \code{\link{new}}, \code{\link{methods}}, \code{\link{show}}
}

\examples{
showClass("modelKS")
slotNames("modelKS")
data(exModelKS)
exModelKS
exModelKS == exModelKS
}


\keyword{methods}
