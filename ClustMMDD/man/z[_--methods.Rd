\name{[<--methods}
\docType{methods}

\alias{[<--methods}
\alias{[<-}

\title{ 
  Get or set a slot from \code{\linkS4class{modelKS}}.
}


\description{
 Get or set a slot from \code{\linkS4class{modelKS}}.
}

\section{Methods}{
\describe{
\item{\code{signature(x = "modelKS")}}{See examples}.
}
}

\author{
  Wilson Toussile.
}

\seealso{
  \code{\linkS4class{modelKS}}, \code{\link{slotNames}}, \code{\link{new}}, \code{\link{methods}}, \code{\link{show}}
}


\examples{
data(exModelKS)
slotNames(exModelKS)
exModelKS["K"]
exModelKS["S"]
}

\keyword{methods}
