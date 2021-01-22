\name{wbsts-package}
\alias{wbsts-package}
\docType{package}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multiple change-point detection for nonstationary time series
}
\description{Implements the Wild Binary Segmentation method of Fryzlewicz (2014) for nostationary time series as described in Korkas and Fryzlewicz (2017). Its purpose is the estimation of the number and locations of the change-points in a time series utilising the wavelet periodogram.

}
\references{
P. Fryzlewicz (2014), Wild Binary Segmentation for multiple change-point detection. Annals of
Statistics, 42, 2243-2281. (http://stats.lse.ac.uk/fryzlewicz/wbs/wbs.pdf)


K. Korkas and P. Fryzlewicz (2017), Multiple change-point detection for non-stationary time series using Wild Binary Segmentation. Statistica Sinica, 27, 287-311. (http://stats.lse.ac.uk/fryzlewicz/WBS_LSW/WBS_LSW.pdf)
}
\author{
K. Korkas and P. Fryzlewicz
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
#### Generate a highly persistent time series with changing variance and of length 5,000
###Location of the change-points
#cps=seq(from=1000,to=2800,by=200)
#y=sim.pw.arma(N =3000,sd_u = c(1,1.5,1,1.5,1,1.5,1,1.5,1,1.5,1),
#b.slope=rep(0.99,11),b.slope2 = rep(0.,11), mac = rep(0.,11),br.loc = cps)[[2]]
###Estimate the change points via Binary Segmentation
#wbs.lsw(y,M=1)$cp.aft
###Estimate the change points via Wild Binary Segmentation
#wbs.lsw(y,M=0)$cp.aft

}

