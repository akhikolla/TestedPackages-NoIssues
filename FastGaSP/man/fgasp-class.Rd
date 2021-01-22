\name{fgasp-class}
\docType{class}
\alias{fgasp-class}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fast GaSP class}
\description{S4 class for fast computation of the Gaussian stochastic process (GaSP) model with the Matern kernel function with or without a noise.}

\section{Objects from the Class}{
Objects of this class are created and initialized with the function \code{\link{fgasp}} that computes the calculations needed for setting up the estimation and prediction.}

\section{Slots}{
	 \describe{
    \item{\code{num_obs}:}{object of class \code{integer}. The number of experimental observations.}
    \item{\code{have_noise}:}{object of class \code{logical} to specify whether the the model has a noise or not. "TRUE" means the model contains a noise and "FALSE" means the model does not contain a noise.}
    \item{\code{kernel_type}:}{a \code{character} to specify the type of kernel to use.The current version supports kernel_type to be "matern_5_2" or "exp", meaning that the matern kernel with roughness parameter being 2.5 or 0.5 (exponent kernel), respectively. }
    \item{\code{input}:}{object of class \code{vector} with dimension num_obs x 1 for the sorted input locations.}
        \item{\code{delta_x}:}{object of class \code{vector} with dimension (num_obs-1) x 1 for the differences between the sorted input locations.}
    \item{\code{output}:}{object of class \code{vector} with dimension num_obs x 1 for the observations at the sorted input locations.}
    }
}


\section{Methods}{
  \describe{
    \item{show}{Prints the main slots of the object. }
    \item{predict}{See \code{\link[=predict.fgasp]{predict}}.}
	 }
}


\references{

Hartikainen, J. and Sarkka, S. (2010). \emph{Kalman filtering and smoothing solutions to temporal gaussian process regression models},  \emph{Machine Learning for Signal Processing (MLSP), 2010 IEEE International Workshop},  379-384.

M. Gu, Y. Xu (2017), \emph{Nonseparable Gaussian stochastic process: a unified
view and computational strategy}, arXiv:1711.11501.

M. Gu, X. Wang and J.O. Berger (2018), \emph{Robust Gaussian Stochastic Process Emulation}, \emph{Annals of Statistics}, \bold{46}, 3038-3066.

}

\author{
\packageAuthor{FastGaSP}

Maintainer: \packageMaintainer{FastGaSP}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{fgasp}} for more details about how to create a \code{fgasp} object.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classes}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
