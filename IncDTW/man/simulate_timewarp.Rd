\name{simulate_timewarp}
\alias{simulate_timewarp}
\alias{insert_const}
\alias{insert_linear_interp}
\alias{insert_norm}
\alias{insert_linear_norm}

\title{
Simulate time warp
}
\description{
Simulate a time warp for a given time series. 
}
\usage{
simulate_timewarp(x, stretch = 0, compress = 0, 
                  stretch_method = insert_linear_interp,
                  p_index = "rnorm", p_number = "rlnorm", 
                  p_index_list = NULL, p_number_list = NULL, 
                  preserve_length = FALSE, seed = NULL, ...) 


insert_const(x, ix, N, const = NULL)

insert_linear_interp(x, ix, N)

insert_norm(x, ix, N, mean = 0, sd = 1)

insert_linear_norm(x, ix, N, mean = 0, sd = 1)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{time series, vector or matrix}
\item{stretch}{numeric parameter for stretching the time series. \code{stretch >= 0}, see details}
\item{compress}{numeric parameter for compressing the time series. \code{compress >= 0 & compress < 1}, see details}
\item{stretch_method}{function, either one of (insert_const, insert_linear_interp,
insert_norm, insert_linear_norm), or any user defined function that needs the parameters \code{x} (univariate time series as vector), \code{ix} (index where to insert), \code{N} (number of observations to insert) and any other arguments required for that function. See Details. }
\item{p_index}{string, distribution for simulating the indices where to insert simulated observations, e.g. "rnorm", "runif", etc.}
\item{p_number}{string, distribution for simulating the number of observations to insert per index, e.g. "rnorm", "runif", etc.}
\item{p_index_list}{list of named parameters for the distribution \code{p_index}}
\item{p_number_list}{list of named parameters for the distribution \code{p_number}}
\item{preserve_length}{logical, if TRUE (default = FALSE) then the length of the return time series is the same as before the warping, so the compression and stretching do not change the length of the time series, nevertheless perform local warpings}
\item{seed}{set a seed for reproducible results}
\item{...}{named parameters for the \code{stretch_method}}
\item{ix}{index of x where after which to insert}
\item{N}{number of simulated observations to insert at index ix}
\item{const}{the constant to be inserted, if NULL (default), then \code{const <- x[ix]}}
\item{mean}{mean for \code{rnorm}}
\item{sd}{sd for \code{rnorm}}
}



\details{
The different distributions \code{p_index} and \code{p_number} also determine the behavior of the warp. A uniform distribution for \code{p_number} more likely draws high number than e.g. log-normal distributions with appropriate parameters. So, a uniform distribution more likely simulates fewer, but longer warps, that is points of time where the algorithm inserts simulations.


The algorithm stretches by randomly selecting an index between 1 and the length of the time series. Then a number of observations to be inserted is drawn out of the range 1 to the remaining number of observations to be inserted. These observations are inserted. Then the algorithm starts again with drawing an index, drawing a number of observations to be inserted, and proceeds until the requested time series length is achieved.

The algorithm for compressing works analogous, except it simply omits observations instead of linear interpolation.


The parameter \code{stretch} describes the ratio how much the time series \code{x} is stretched:
e.g. if \code{compress = 0} and ...
\itemize{
\item \code{stretch = 0} then \code{length(x_new) = length(x)}, or
\item \code{stretch = 0.1} then \code{length(x_new) = length(x) * 1.1}, or             
\item \code{stretch = 1} then \code{length(x_new) = length(x) * 2}                        
}

The parameter \code{compress} describes the ratio how much the time series \code{x} is compressed:
e.g. if \code{stretch = 0} and ...
\itemize{
\item \code{compress = 0} then \code{length(x_new) = length(x)}, or
\item \code{compress = 0.2} then \code{length(x_new) = length(x) * 0.8}
}


There are four functions to chose from to insert new simulated observations. You can also define your own function and apply this one. The four functions to chose from are:
\itemize{
\item insert a constant, either a constant defined by the user via the input parameter \code{const}, or if \code{const = NULL}, then the last observation of the time series where the insertion starts is set as \code{const}
\item insert linear interpolated observations (default)
\item insert a constant with gaussian noise 
\item insert linear interpolated observations and add gaussian noise.
}

For the methods with Gaussian noise the parameters \code{mean} and \code{sd} for \code{rnorm} can be set at the function call of \code{simulate_timewarp()}.

}
\value{
\item{}{time warped time series}
}


\examples{
\dontrun{
#--- Simulate a time warped version of a time series x
set.seed(123)
x <- cumsum(rnorm(100))
x_warp <- simulate_timewarp(x, stretch = 0.1, compress = 0.2, seed = 123) 
plot(x, type = "l")
lines(x_warp, col = "red")


#--- Simulate a time warp of a multivariate time series
y <- matrix(cumsum(rnorm(10^3)), ncol = 2)
y_warp <- simulate_timewarp(y, stretch = 0.1, compress = 0.2, seed = 123) 
plot(y[,1], type = "l")
lines(y_warp[,1], col = "red")


#--- Stretchings means to insert at new values at randomly 
# selected points of time. Next the new values are set as constant NA,
# and the points of time simulated uniformly:
y_warp <- simulate_timewarp(y, stretch = 0.2, p_number = "runif", p_index = "runif",
                            stretch_method = insert_const, 
                            const = NA)
matplot(y_warp, type = "l")


# insert NA and simulate the points of time by log normal
y_warp <- simulate_timewarp(y, stretch = 0.2, p_number = "rlnorm", 
                            p_number_list = list(meanlog = 0, sdlog = 1),
                            stretch_method = insert_const, 
                            const = NA)
matplot(y_warp, type = "l")


# insert linear interpolation
y_warp <- simulate_timewarp(y, stretch = 0.2, p_number = "rlnorm", 
                            stretch_method = insert_linear_interp)
matplot(y_warp, type = "l")


# insert random walk with gaussian noise
y_warp <- simulate_timewarp(y, stretch = 0.2, p_number = "rlnorm", 
                            stretch_method = insert_norm,
                            sd = 1, mean = 0)
matplot(y_warp, type = "l")


# insert constant, only 1 observation per random index
y_warp <- simulate_timewarp(y, stretch = 0.2, p_number = "runif", p_index = "runif",
                            p_number_list = list(min = 1, max = 1),   
                            stretch_method = insert_const)
matplot(y_warp, type = "l")


# insert by customized insert function
my_stretch_method <- function(x, ix, N, from, to){
   c(x[1:ix], 
     sin(seq(from = from, to = to, length.out = N)) + x[ix],
     x[(ix + 1):length(x)])
}
y_warp <- simulate_timewarp(y, stretch = 0.5, p_number = "rlnorm", 
                            stretch_method = my_stretch_method,
                            from = 0, to = 4 * pi)
matplot(y_warp, type = "l")

}
}

\keyword{ ts }% __ONLY ONE__ keyword per line

