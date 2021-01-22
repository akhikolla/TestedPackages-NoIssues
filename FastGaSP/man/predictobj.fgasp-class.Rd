\name{predictobj.fgasp-class}
\alias{predictobj.fgasp-class}
\alias{predictobj.fgasp}
\docType{class}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Predictive results for the Fast GaSP class }
\description{S4 class for prediction for a Fast GaSP model with or without a noise.}

\section{Objects from the Class}{
Objects of this class are created and initialized with the function \code{\link{predict}} that computes the prediction and the uncertainty quantification.}

\section{Slots}{
	 \describe{
	     \item{\code{num_testing}:}{object of class \code{integer}. Number of testing inputs.}
    \item{\code{testing_input}:}{object of class \code{vector}. The testing input locations.}
      \item{param}{a vector of parameters. The first parameter is the natural logarithm of the inverse range parameter in the kernel function. If the data contain noise, the second parameter is the logarithm of the nugget-variance ratio parameter.
}
    \item{\code{mean}:}{object of class \code{vector}. The predictive mean at testing inputs.}
        \item{\code{var}:}{object of class \code{vector}. The predictive variance at testing inputs. If the \code{var_data} is true, the predictive variance of the data is calculated. Otherwise, the predictive variance of the mean is calculated. }
        \item{\code{var_data}:}{object of class \code{logical}. If the \code{var_data} is true, the predictive variance of the data is calculated for \code{var}. Otherwise, the predictive variance of the mean is calculated for \code{var}. }

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

\seealso{\code{\link{predict.fgasp}} for more details about how to do prediction for a \code{fgasp} object.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classes}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
