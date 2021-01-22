\name{uh.wbs}
\alias{uh.wbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Wild Binary Segmentation algorithm
}
\description{
The function implements the Wild Binary Segmentation method and aggregates the change-points across the wavelet periodogram. Currently only the Method 2 of aggregation is implemented.
}
\usage{
uh.wbs(z,C_i, del=-1, epp, scale,M=0,cstar=0.75)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
The wavelet periodogram matrix.
}
  \item{C_i}{
The CUSUM threshold.
}

  \item{del}{
The minimum allowed size of a segment.
}

  \item{epp}{
A minimum adjustment for the bias present in \eqn{E^{(i)}_{t,T}}.
}
  \item{scale}{
Which wavelet periodogram scales to be used.
}
  \item{M}{
The maximum number of random intervals drawn. If M=0 (DEFAULT) this is selected to be a linear function of the sample size of y. If M=1 then the segmentation is conducted via the Binary segmentation method.
}
  \item{cstar}{
This refers to the unbalanceness parameter \eqn{c_{\star}}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{cp.bef}{Returns the estimated change-points before post-processing}
  \item{cp.aft}{Returns the estimated change-points after post-processing}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
K. Korkas and P. Fryzlewicz (2017), Multiple change-point detection for non-stationary time series using Wild Binary Segmentation. Statistica Sinica, 27, 287-311. (http://stats.lse.ac.uk/fryzlewicz/WBS_LSW/WBS_LSW.pdf)
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{binary segmentation}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
