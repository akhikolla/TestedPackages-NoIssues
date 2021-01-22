\name{scale}
\alias{scale}
\alias{norm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Time Series Scaling
}
\description{
scales a time series per dimension/column. 
}
\usage{
scale(x, type = c('z', '01'), threshold = 1e-5,
     xmean = NULL, xsd = NULL, xmin = NULL, xmax = NULL)

# deprecated
norm(x, type = c('z', '01'), threshold = 1e-5,
     xmean = NULL, xsd = NULL, xmin = NULL, xmax = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{time series as vector or matrix
}
\item{type}{character, describes the method how to scale (or normalize) the time series (per column if x is multivariate). The parameter type is either 'z' for z-scaling or '01' for max-min scaling.
}
\item{threshold}{
double, defines the minimum value of the standard deviation, or difference of minimum and maximum. If this value is smaller than the threshold, then no scaling is performed, only shifting by the mean or minimum, respectively. Default value = \code{1e-5}.
}
\item{xmean}{mean used for z-scaling. See details.
}
\item{xsd}{standard deviation used for z-scaling. See details.
}
\item{xmin}{minimum used for 0-1 scaling. See details.
}
\item{xmax}{maximum used for 0-1 scaling. See details.
}
}


\details{
For a vector \code{x} the z-scaling subtracts the mean and devides by the standard deviation:  of \code{(x-mean(x))/sd(x)}. The min-max scaling performs \code{(x-min(x))/(max(x)-min(x))}.

The parameters \code{xmean}, \code{xsd}, \code{xmin}, can be set \code{xmax} or passed as \code{NULL} (= default value). If these values are NULL, they are calculated based on x.
}

\value{
\item{x}{the scaled vector or matrix}
}

\examples{

# min-max scaling or z-scaling for a vector 
x <- cumsum(rnorm(100, 10, 5))
y <- scale(x, "01")
z <- scale(x, "z")
par(mfrow = c(1, 3))
plot(x, type="l")
plot(y, type="l")
plot(z, type="l")
par(mfrow = c(1, 1))


# columnwise for matrices
x <- matrix(c(1:10, sin(1:10)), ncol = 2)
y <- scale(x, "01")
z <- scale(x, "z")
par(mfrow = c(1, 3))
matplot(x, type="l")
matplot(y, type="l")
matplot(z, type="l")
par(mfrow = c(1, 1))


# IncDTW::scale() and base::scale() perform same z-scaling
x <- cumsum(rnorm(100))
xi <- IncDTW::scale(x, type = "z")
xb <- base::scale(x, TRUE, TRUE)
sum(abs(xi-xb))
}

\keyword{ classif }% use one of  RShowDoc("KEYWORDS")
\keyword{ cluster }% use one of  RShowDoc("KEYWORDS")
\keyword{ ts }% __ONLY ONE__ keyword per line

