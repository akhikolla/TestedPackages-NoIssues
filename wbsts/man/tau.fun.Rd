\name{tau.fun}
\alias{tau.fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Universal thresholds
}
\description{
The function returns \eqn{C^{(i)}}. \eqn{C^{(i)}} tends to increase as we move to coarser scales due to
the increasing dependence in the wavelet periodogram sequences. Since the method applies to non-dyadic structures it is reasonable to propose a general rule that will apply in most cases. To accomplish this the \eqn{C^{(i)}} are obtained for \eqn{T=50,100,...,6000}. Then, for each scale \eqn{i} the following regression is fitted 

\eqn{C^{(i)}=c_0^{(i)}+c_1^{(i)} T+ c_2^{(i)} \frac{1}{T} + c_3^{(i)} T^2 +\varepsilon.}

The adjusted \eqn{R^2} was above 90\% for all the scales. Having estimated the values for \eqn{\hat{c}_0^{(i)}, \hat{c}_1^{(i)}, \hat{c}_2^{(i)}, \hat{c}_3^{(i)}} the values can be retrieved for any sample size \eqn{T}.
}
\usage{
tau.fun(y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
A time series
}
 
}

\value{
Thresholds for every wavelet scale
}
\references{
P. Fryzlewicz (2014), Wild Binary Segmentation for multiple change-point detection. Annals of
Statistics, 42, 2243-2281. (http://stats.lse.ac.uk/fryzlewicz/wbs/wbs.pdf)

K. Korkas and P. Fryzlewicz (2017), Multiple change-point detection for non-stationary time series using Wild Binary Segmentation. Statistica Sinica, 27, 287-311. (http://stats.lse.ac.uk/fryzlewicz/WBS_LSW/WBS_LSW.pdf)
}
\author{
K. Korkas and P. Fryzlewicz
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{ 
##not run##
#cps=c(400,470)
#set.seed(101)
#y=sim.pw.ar(N =2000,sd_u = 1,b.slope=c(0.4,-0.6,0.5),br.loc=cps)[[2]]
#tau.fun(y) is the default value for C_i
##Binary segmentation
#wbs.lsw(y,M=1)$cp.aft
##Wild binary segmentation
#wbs.lsw(y,M=3500)$cp.aft

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{threshold}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
