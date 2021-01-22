\name{finner_prod_maxp}
\alias{finner_prod_maxp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The function finds the value which yields the maximum inner product
(CUSUM) of a a time series located between \eqn{100(1-p)\%} and \eqn{100p\%} of its support
}
\description{
This function is an internal C++ function wrapped by cr.rand.max.inner.prod.
}
\usage{
finner_prod_maxp(x,p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A time series.
}
  \item{p}{
A scalar in (0.67,1]
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{1}{maximum cusum value}
\item{2}{the location of a point with the maximum cusum value}

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
##balanced
#finner_prod_maxp(z[,1],0.51)
##unbalanced
#finner_prod_maxp(z[,1],0.99)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal }% use one of  RShowDoc("KEYWORDS")
