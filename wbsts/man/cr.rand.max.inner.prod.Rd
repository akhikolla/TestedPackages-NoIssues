\name{cr.rand.max.inner.prod}
\alias{cr.rand.max.inner.prod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The value that maximises the random CUSUM statistic across all the scales

}
\description{
The function finds the value which yields the maximum inner product  with the input time series (CUSUM) located between \eqn{100(1-p)\%} and \eqn{100p\%} of their support across all the wavelet periodogram scales.
}
\usage{
cr.rand.max.inner.prod(XX,Ts,C_i,epp,M = 0,Plot = FALSE,cstar=0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{XX}{
The wavelet periodogram.
}
  \item{Ts}{
The sample size of the series.
}
  \item{C_i}{
The CUSUM threshold.
}

  \item{epp}{
A minimum adjustment for the bias present in \eqn{E^{(i)}_{t,T}}.
}

  \item{M}{
Number of random CUSUM to be generated.
}

  \item{Plot}{
Plot the threhsold CUSUM statistics across the wavelet scales.
}

  \item{cstar}{
A scalar in (0.67,1]
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{1}{Candidate change point}
\item{2}{The maximum CUSUM value}
\item{3}{The starting point \eqn{s} of the favourable draw}
\item{4}{The ending point \eqn{e} of the favourable draw}

%% ...
}
\references{
K. Korkas and P. Fryzlewicz (2017), Multiple change-point detection for non-stationary time series using Wild Binary Segmentation. Statistica Sinica, 27, 287-311. (http://stats.lse.ac.uk/fryzlewicz/WBS_LSW/WBS_LSW.pdf)
}
\author{
K. Korkas and P. Fryzlewicz
}
\examples{ 
#cps=seq(from=1000,to=2000,by=200)
#y=sim.pw.arma(N =3000,sd_u = c(1,1.5,1,1.5,1,1.5,1),
#b.slope=rep(0.99,7),b.slope2 = rep(0.,7), mac = rep(0.,7),br.loc = cps)[[2]]
#z=ews.trans(y,scales=c(11,9,8,7,6))
#out=cr.rand.max.inner.prod(z, Ts = length(y),C_i = tau.fun(y), 
#epp = rep(32,5), M = 2000, cstar = 0.75, Plot = 1)
#abline(v=cps,col="red")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ piecewise ARMA }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
