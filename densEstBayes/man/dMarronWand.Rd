\name{dMarronWand}
\alias{dMarronWand}
\title{Marron and Wand density function}
\description{Returns ordinates of a member of the family of Normal Mixture density functions devised in Marron and Wand (1992).}
\usage{
dMarronWand(x,densNum,drv)
}
\arguments{
\item{x}{Vector of abscissae values.}
\item{densNum}{An integer between 1 and 15 that specifies the density function according to Table 1 of Marron and Wand (1992).}
\item{drv}{Either -1,0,1,2 depending on the order of the derivative required, with \code{drv = -1} corresponding to the anti-derivative. The default value is \code{drv = 0}.}
}
\author{Matt P. Wand \email{matt.wand@uts.edu.au}}
\references{Marron, J.S. and Wand, M.P. (1992). Exact mean integrated squared error. \emph{The Annals of Statistics}, \bold{20}, 712-736.}
\examples{
library(densEstBayes)
xg <- seq(-3,3,length = 1001)
densg <- dMarronWand(xg,8)
plot(xg,densg,type = "l")
}

