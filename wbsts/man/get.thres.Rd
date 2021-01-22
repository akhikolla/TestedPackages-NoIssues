\name{get.thres}
\alias{get.thres}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Universal thresholds calculation
}
\description{
The function returns universal thresholds and the method is described in Korkas and Fryzlewicz (2017) and Cho and Fryzlewicz (2012). See also the supplementary material for the former work. The function works for any sample size.
}
\usage{
get.thres(n, q=.95, r=100, scales=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
The length of the time series.
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

Cho, H. and Fryzlewicz, P. (2012). Multiscale and multilevel technique for consistent segmentation of nonstationary time series. Statistica Sinica, 22(1), 207-229.
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{change-point detection }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
