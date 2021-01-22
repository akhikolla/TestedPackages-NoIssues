\name{PAR}
\Rdversion{7.1}
\docType{class}
\alias{PAR-class}

\title{Class to contain the variables for the verification plots from function MOPT, ceff1, ceff2, and Deff.}
\description{
For function MOPT, ceff1, ceff2, and Deff, the plot of the sensitivity function over the given
dose interval can be used to verify the optimality of the obtained design. Class PAR is used to contain the variables for the plot.
}
\section{Objects from the Class}{
new("PAR",fid="character",LB="numeric",UB="numeric",grid="numeric",ds="vector")
}
\section{Slots}{
  \describe{
    \item{\code{fid}:}{A character showing the name of the function from where the four arguments are.}
    \item{\code{LB}:}{A numeric number showing the predetermined lower bound of the dose range for the log dose.}
    \item{\code{UB}:}{A numeric number showing the predetermined upper bound of the dose range for the log dose.}
    \item{\code{grid}:}{A numeric number showing the grid density to discretize the predetermined dose interval.}
    \item{\code{ds}:}{A numeric vector represents the values of criterion (an absolute value of the difference between the sensitivity function and 1) computed at each iteration. }
  }
}
\author{
Seung Won Hyun, Weng Kee Wong, and Yarong Yang
}
\examples{
showClass("PAR")
}
\keyword{classes}


