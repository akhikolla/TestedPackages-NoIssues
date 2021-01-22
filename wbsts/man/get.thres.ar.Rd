\name{get.thres.ar}
\alias{get.thres.ar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Selection of thresholds by fitting an AR(p) model
}
\description{
The function returns data-driven thresholds and it is described in Korkas and Fryzlewicz (2015) where
it is referred as Bsp1. See also the supplementary material for this work.
}
\usage{
get.thres.ar(y, q=.95, r=100, scales=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
The time series.
}
  \item{q}{
The quantile of the r simulations.
}
  \item{r}{
Number of simulations.
}
  \item{scales}{
The wavelet periodogram scales to be used. If NULL (DEFAULT) then this is selected as described in the main text.
}
}

\references{
K. Korkas and P. Fryzlewicz (2017), Multiple change-point detection for non-stationary time series using Wild Binary Segmentation. Statistica Sinica, 27, 287-311. (http://stats.lse.ac.uk/fryzlewicz/WBS_LSW/WBS_LSW.pdf)

K. Korkas and P. Fryzlewicz (2017), Supplementary material: Multiple change-point detection for non-stationary time series using Wild Binary Segmentation.
}
\author{
K. Korkas and P. Fryzlewicz
}
\examples{
#cps=seq(from=100,to=1200,by=350)
#y=sim.pw.arma(N =1200,sd_u = c(1,1.5,1,1.5,1),
#b.slope=rep(0.99,5),b.slope2 = rep(0.,5), mac = rep(0.,5),br.loc = cps)[[2]]
#C_i=get.thres.ar(y=y, q=.95, r=100, scales=NULL)
#wbs.lsw(y,M=1, C_i = C_i)$cp.aft


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{change-point detection }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
