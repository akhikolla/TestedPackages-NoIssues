\name{SW}
\Rdversion{7.1}
\docType{class}
\alias{SW-class}
\alias{summary,SW-method}

\title{Class to contain the results from function S. Weight.}
\description{
 The function S.Weight returns an object of class SW that contains the optimal weights and the first and second derivatives of the criterion.
}
\section{Objects from the Class}{
new("SW",Opt.W=new("matrix"),First.C=new("vector"),Second.C=new("vector"))
}
\section{Slots}{
  \describe{
    \item{\code{Opt.W}:}{A numeric matrix showing the optimal weights.}
    \item{\code{First.C}:}{A numeric vector showing the first derivative of the criterion.}
    \item{\code{Second.C}:}{A numeric vector showing the second derivative of the criterion.}
  }
}
\section{Methods}{
  \describe{
    \item{summary}{\code{signature(object = "SW")}: ... }
	 }
}
\author{
Seung Won Hyun, Weng Kee Wong, and Yarong Yang
}
\examples{
showClass("SW")
}
\keyword{classes}
