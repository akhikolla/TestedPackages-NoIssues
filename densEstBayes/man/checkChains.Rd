\name{checkChains}
\alias{checkChains}
\title{Check Markov chain Monte Carlo samples}
\description{Facilitates a graphical check of the Markov chain Monte Carlo samples ("chains") corresponding to vertical slices of the log-density function estimate at five representative abscissa values.}
\usage{
checkChains(fitObject,locations = "equally-spaced")
}
\arguments{
\item{fitObject}{\code{densEstBayes()} fit object.}
\item{locations}{Character string that specifies the locations of the vertical slices for the chains to be displayed. The
options are:\cr
"equally-spaced" for 5 equally-spaced locations over the range of the data.\cr 
"hexiles" for the locations being the sample hexiles of the data.\cr
The default value is "equally-spaced".
}
}
\details{A graphic is produced that summarises the Markov chain Monte Carlo samples ("chains") corresponding to vertical slices of the log-density function estimate at five representative abscissae. If the location is specified to be "equally-spaced" and the range.x array of the density estimate is the interval [a,b] then the abscissae are a + j(b-a)/6 for j=1,2,3,4,5. If the location is specified to be "hexiles" then the abscissae are the sample hexiles of the data. The columns of the graphic are the following summaries of each chain: (1) trace (time series) plot, (2) lag-1 plot in which each chain value is plotted against its previous value and (3) sample autocorrelation function plot as produced by the R function acf().}
\author{Matt P. Wand \email{matt.wand@uts.edu.au}}
\examples{
library(densEstBayes) ; data(OldFaithful2011)

# Obtain a density estimate for the `OldFaithful2011' data:
\donttest{
dest <- densEstBayes(OldFaithful2011)

# Obtain a graphic for checking the chains:

checkChains(dest)
}
}

