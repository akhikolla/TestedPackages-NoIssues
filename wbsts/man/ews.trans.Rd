\name{ews.trans}
\alias{ews.trans}

\title{
Computation of the Evolutionary Wavelet Spectrum (EWS)
}
\description{
The function computes the EWS from a time series of any (non-dyadic) size by utilising the
maximal overlap discrete wavelet transform; see also W. Constantine and D. Percival (2015).
}
\usage{
ews.trans(x,scales=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The time series.
}
  \item{scales}{
The wavelet periodogram scales to compute starting from the finest.
}}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
The evolutionary wavelet spectral estimate of y.
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
	Eric Aldrich (2020), wavelets: Functions for Computing Wavelet Filters, Wavelet Transforms and Multiresolution Analyses.
}
\examples{
ews=ews.trans(rnorm(1000),c(9,8,7))
barplot(ews[,1])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
