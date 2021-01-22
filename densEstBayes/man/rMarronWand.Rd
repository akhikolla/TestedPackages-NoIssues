\name{rMarronWand}
\alias{rMarronWand}
\title{Marron and Wand random sample}
\description{Returns a random sample from a member of the family of Normal Mixture density functions devised in Marron and Wand (1992).}
\usage{
rMarronWand(n,densNum)
}
\arguments{
\item{n}{The sample size, which is a positive integer.}
\item{densNum}{An integer between 1 and 15 that specifies the density function according to Table 1 of Marron and Wand (1992).}
}
\author{Matt P. Wand \email{matt.wand@uts.edu.au}}
\references{Marron, J.S. and Wand, M.P. (1992). Exact mean integrated squared error. \emph{The Annals of Statistics}, \bold{20}, 712-736.}
\examples{
library(densEstBayes)
x <- rMarronWand(1000,8)
hist(x,breaks = 50,col = "gold")
dest <- densEstBayes(x,method = "SMFVB")
plot(dest) ; rug(x,col = "dodgerblue",quiet = TRUE)
}
