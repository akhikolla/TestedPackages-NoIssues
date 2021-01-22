\name{OPT}
\Rdversion{7.1}
\docType{class}
\alias{OPT-class}
\alias{plot,OPT-method}
\alias{summary,OPT-method}

\title{Class to contain the results from function MOPT, ceff1, ceff2, and Deff.}
\description{
The function MOPT, ceff1, ceff2, Deff return object of class OPT that contains the variables for the verification plot, the optimal matrix, and the efficiency.
}
\section{Objects from the Class}{
new("OPT",Par="PAR",OPT=new("matrix"),Eff=new("numeric"))
}
\section{Slots}{
  \describe{
    \item{\code{Par}:}{An object of class PAR containing the variables for the verification plot.}
    \item{\code{Opt}:}{A numeric matrix showing the multiple-objective optimal design for estimating model parameters, the ED50, and the MED (for function MOPT), or c-optimal design for ED50 (for function ceff1),
   or c-optimal design for MED (for function ceff2), or D-optimal design (for function Deff).}
    \item{\code{Eff}:}{A numeric number showing the c-efficiency for ED50 (for function ceff1), or c-efficiency for MED (for function ceff2), or D-efficiency (for function Deff).}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "OPT")}: ... }
    \item{summary}{\code{signature(object = "OPT")}: ... }
	 }
}
\author{
Seung Won Hyun, Weng Kee Wong, and Yarong Yang
}
\examples{
showClass("OPT")
}
\keyword{classes}
