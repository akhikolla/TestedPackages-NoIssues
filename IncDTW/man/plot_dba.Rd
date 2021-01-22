\name{plot.dba}
\alias{plot_dba}
\alias{plot.dba}
\alias{plotBary}
\alias{plotM2m}
\alias{plotM2lot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the results from Dynamic Time Warping Barycenter Averaging
}
\description{
Plot function for objects of type \code{dba}, the output of \code{dba()}.}
\usage{
\method{plot}{dba}(x, type = c("barycenter", "m2m", "m2lot"), ...)
# an alias for plot_dba
plot_dba(x, type = c("barycenter", "m2m", "m2lot"), ...)

plotBary(x, ...)

plotM2m(x, ...)

plotM2lot(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{output from \code{dba()}}
 \item{type}{character, one of c('barycenter', 'm2m', 'm2lot')}
\item{...}{Other arguments passed on to methods. Currently not used.}
}
\details{
\itemize{
\item{\code{'barycenter'} plots the iterations of the barycenter per dimension.}
\item{\code{'m2m'} plots the distances (distance method set by \code{iter_dist_method}, see \code{\link{dba}}) of one barycenter-iteration to the former iteration step.}
\item{\code{'m2lot'} plots the distances (if \code{step_pattern == 'symmetric2'} the normalized distances are plotted) of the barycenter to the list of time series per iteration.}
}
 
}


\examples{
# see examples of dba()
}


\seealso{
\code{\link{dba}}
}
