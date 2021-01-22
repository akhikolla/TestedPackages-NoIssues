\name{wbs.lsw}
\alias{wbs.lsw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Change point detection for a nonstationary process using Wild Binary Segmentation
}
\description{
The function returns the estimated locations of the change-points in a nonstationary time series. Currently only the Method 2 of aggregation is implemented.
}
\usage{
wbs.lsw(y, C_i = tau.fun(y), scales = NULL, M = 0, cstar = 0.75, lambda = 0.75)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
The time series.
}
  \item{C_i}{
A vector of threshold parameters for different scales.
}
  \item{scales}{
The wavelet periodogram scales to be used. If NULL (DEFAULT) then this is selected as described in the main text.
}
  \item{M}{
The maximum number of random intervals drawn. If M=0 (DEFAULT) this is selected to be a linear function of the sample size of y. If M=1 then the segmentation is conducted via the Binary segmentation method.
}
  \item{cstar}{
This refers to the unbalanceness parameter \eqn{c_{\star}}.
}
  \item{lambda}{
This parameter defines the maximum number of the wavelet periodogam scales. This is used if scales = NULL.
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
\author{
K. Korkas and P. Fryzlewicz
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

