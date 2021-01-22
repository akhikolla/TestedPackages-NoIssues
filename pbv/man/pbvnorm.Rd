%% File Name: pbvnorm.Rd
%% File Version: 0.34

\name{pbvnorm}
\alias{pbvnorm}
\alias{dbvnorm}
\alias{pbv_rcpp_dbvnorm}
\alias{pbv_rcpp_dbvnorm0}
\alias{pbv_rcpp_pbvnorm}
\alias{pbv_rcpp_pbvnorm0}


%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Probabilities for Bivariate Normal Distribution
}

\description{
The function \code{pbvnorm} computes probabilities \eqn{\Phi_2(x,y,\rho)} for
the standardized bivariate normal distribution (Drezner & Wesolowsky, 1990;
West, 2004).

The function \code{dbvnorm} computes the corresponding density \eqn{\phi_2(x,y,\rho)}.
}

\usage{
pbvnorm(x, y, rho)

dbvnorm(x, y, rho, log=FALSE)

## exported Rcpp functions
pbv_rcpp_pbvnorm0( h1, hk, r)
pbv_rcpp_pbvnorm( x, y, rho)
pbv_rcpp_dbvnorm0( x, y, rho, use_log)
pbv_rcpp_dbvnorm( x, y, rho, use_log)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Vector of first ordinate
}
  \item{y}{
Vector of second ordinate
}
  \item{rho}{
Vector of correlations
}
\item{log}{Logical indicating whether logarithm of the density
should be calculated}
\item{h1}{Numeric}
\item{hk}{Numeric}
\item{r}{Numeric}
\item{use_log}{Logical}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
A vector
}

\references{
Drezner, Z., & Wesolowsky, G. O. (1990). On the computation of the bivariate normal
integral. \emph{Journal of Statistical Computation and Simulation, 35}(1-2), 101-107.
\doi{10.1080/00949659008811236}

Genz, A. (1992). Numerical computation of multivariate normal probabilities.
\emph{Journal of Computational and Graphical Statistics, 1}(2), 141-149.

West, G. (2005). Better approximations to cumulative normal functions.
\emph{Wilmott Magazine, 9}, 70-76.
}
%% Genz (1992) \doi{10.2307/1390838}

%\author{
%%  ~~who you are~~
%}
\note{
The \pkg{pbv} package can also be used to include \CRANpkg{Rcpp} functions for
computing bivariate probabilities at the C++ level. Numeric and vector versions are
\cr \cr
\code{double pbv::pbv_rcpp_pbvnorm0( double h1, double hk, double r)} \cr \cr
\code{Rcpp::NumericVector pbv::pbv_rcpp_pbvnorm( Rcpp::NumericVector x, }\cr
\code{Rcpp::NumericVector y, Rcpp::NumericVector rho)}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link[pbivnorm:pbivnorm]{pbivnorm::pbivnorm}} in the \pkg{pbivnorm} package
and \code{mnormt::biv.nt.prob} in the \pkg{mnormt} package
for alternative implementations (Genz, 1992).
}



\examples{
#############################################################################
# EXAMPLE 1: Comparison with alternative implementations
#############################################################################

#*** simulate different values of ordinates and correlations
set.seed(9898)
N <- 3000
x <- stats::runif(N,-3,3)
y <- stats::runif(N,-3,3)
rho <- stats::runif(N,-.95,.95)

#*** compute probabilities
res1 <- pbv::pbvnorm(x=x,y=y,rho=rho)

#-- compare results with pbivnorm package
library(pbivnorm)
res2 <- pbivnorm::pbivnorm(x=x, y=y, rho=rho)

summary(abs(res1-res2))

#*** compute density values
log <- TRUE    # logical indicating whether log density should be evaluated
res1 <- pbv::dbvnorm(x=x, y=y, rho=rho, log=log )

\dontrun{
#-- compare results with mvtnorm package
library(mvtnorm)
res2 <- rep(NA, N)
sigma <- diag(2)
for (ii in 1:N){
    sigma[1,2] <- sigma[2,1] <- rho[ii]
    res2[ii] <- mvtnorm::dmvnorm(x=c(x[ii],y[ii]), sigma=sigma, log=log)
}
summary(abs(res1-res2))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{\code{mnormt::biv.nt.prob}}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
