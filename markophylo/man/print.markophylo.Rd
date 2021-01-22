\name{print.markophylo}
\alias{print.markophylo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print summary information for the model fit.
}
\description{
Summary command for use on an object of class "markophylo". The estimated substitution rates for each partition specified along with the estimated \eqn{\alpha}(s) for the gamma distribution (if any) and prior probability of discrete characters at the root (if estimated) are printed. If branch groupings (or clades) were specified, then the rates (and corresponding standard errors) are displayed in a matrix with the columns representing the different branch groupings (ordered by the subsets of x$bg where x is an object of class "idea"). The rows represent the indices of the model matrix.}
\usage{
\method{print}{markophylo}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class "markophylo".}
  \item{\dots}{
Ignore this
}
}
\author{
Utkarsh J. Dang and G. Brian Golding

\email{udang@mcmaster.ca}}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{estimaterates}}.
}