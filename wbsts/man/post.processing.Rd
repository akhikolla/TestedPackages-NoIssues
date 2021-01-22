\name{post.processing}
\alias{post.processing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Post-processing of the change-points
}
\description{
A function to control the number of change-points estimated from the WBS algorithm
and to reduce the risk of over-segmentation.
}
\usage{
post.processing(z,br,del=-1,epp=-1,C_i=NULL,scales=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
The wavelet periodogram matrix.
}
  \item{br}{
The change-points to be post-processed.
}
  \item{del}{
The minimum allowed size of a segment.
}
  \item{epp}{
A minimum adjustment for the bias present in \eqn{E^{(i)}_{t,T}}.
}
  \item{C_i}{
The CUSUM threshold.
}
  \item{scales}{
Which wavelet periodogram scales to be used.
}
}

\references{
K. Korkas and P. Fryzlewicz (2017), Multiple change-point detection for non-stationary time series using Wild Binary Segmentation. Statistica Sinica, 27, 287-311. (http://stats.lse.ac.uk/fryzlewicz/WBS_LSW/WBS_LSW.pdf)
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

%\examples{
%#### Generate a highly persistent time series with changing variance and of length 5,000
%###Location of the change-points
%cps=seq(from=1000,to=2800,by=200)
%y=sim.pw.arma(N =3000,sd_u = c(1,1.5,1,1.5,1,1.5,1,1.5,1,1.5,1),
%b.slope=rep(0.99,11),b.slope2 = rep(0.,11), mac = rep(0.,11),br.loc = cps)[[2]]
%###Estimate the change points via Wild Binary Segmentation
%beforeProcessing=wbs.lsw(y,M=0)$cp.bef
%###Post-processing of the change points
%post.processing(z=ews.trans(y,c(11,10,9,8)),br=beforeProcessing,C_i=tau.fun(y),scales=c(11,10,9,8))

%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{change-point detection }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
