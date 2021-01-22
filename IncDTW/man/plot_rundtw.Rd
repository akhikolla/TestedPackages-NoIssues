\name{plot.rundtw}
\alias{plot_rundtw}
\alias{plot.rundtw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot
}
\description{
Plot the results from \code{\link{rundtw}}.}
\usage{
\method{plot}{rundtw}(x, knn = TRUE, minima = TRUE, 
     scale = c("none", "01", "z"), 
     selDim = 1, lix = 1, Q = NULL, C = NULL, normalize = c("none", "01", "z"),  ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{
 output from \code{\link{rundtw}}
 }
 \item{knn}{
 logical, if TRUE ( = default) and the k nearest neighbors were found by \code{\link{rundtw}}, then they are plotted. See details.
 }
\item{minima}{
logical, if TRUE ( = default) and \code{Q} is either passed or also returned by \code{\link{rundtw}} then the local (with window size of the lengh of \code{Q}) minima of the vector of distances is plotted. See details. 
}
\item{scale}{
character, one of c("none", "01", "z"). If "01" or "z" then the detected minima and knn are normed and plotted.
}
\item{selDim}{
integer vector, default = 1. Set the dimensions to be plotted for multivariate time series \code{Q} and \code{C}.
}
\item{lix}{
list index, integer, default = 1. If \code{C} is a list of time series, set with \code{lix} the list entry of \code{C} to be plotted.
}
\item{Q}{
time series, default = NULL, either passed as list entry of \code{x} (when the parameter \code{return_QC} of \code{\link{rundtw}} is set to TRUE) or passed manually. Necessary for plotting the minima.
}
\item{C}{
time series, default = NULL, either passed as list entry of \code{x} (when the parameter \code{return_QC} or \code{\link{rundtw}} is set to TRUE) or passed manually. Necessary for plotting the minima and knn.
}
\item{normalize}{deprecated, use \code{scale} instead. If \code{normalize} is set, then \code{scale} is overwritten by \code{normalize} for compatibility.
}
\item{...}{
additional arguments passed to ggplot()}
}



\details{
Only for those subsequences for which the calculations were finished by  \code{rundtw}, the distances  are plotted (see the parameters \code{threshold}, \code{k} and \code{early_abandon} of \code{\link{rundtw}}).
}


\examples{

#--- Simulate a query pattern Q and a longer time series C,
# and detect rescaled versions of Q in C
set.seed(123)
Q <- sin(seq(0, 2*pi, length.out = 20))
Q_rescaled <- Q * abs(rnorm(1)) + rnorm(1)
C <- c(rnorm(20), Q_rescaled , rnorm(20))

# Force rundtw to finish all calculations and plot the vector of DTW distances
ret <- rundtw(Q, C, threshold = NULL, lower_bound = FALSE)
ret
plot(ret)

# Allow early abandoning and lower bounding, and also plot C
ret <- rundtw(Q, C, return_QC = TRUE, ws = 5)
ret
plot(ret)

# Get 1 nearest neighbor -> allow early abandon and lower bounding,
# and plot C and also plot the scaled detected nearest neighbors
ret <- rundtw(Q, C, ws = 5, k = 1, return_QC = TRUE)
ret
plot(ret, scale = "01")


#--- See the help page of rundtw() for further examples.

}

\seealso{
\code{\link{rundtw}}
}
