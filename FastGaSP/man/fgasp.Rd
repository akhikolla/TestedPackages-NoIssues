\name{fgasp}
\alias{fgasp}
%\alias{show.rgasp}
\alias{fgasp-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Setting up the Fast GaSP model
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Creating an \code{fgasp} class for a GaSP model with matern covariance. 
%The range and noise-variance ratio parameters are given and/or have been estimated.
}
\usage{
  fgasp(input, output, have_noise=TRUE, kernel_type='matern_5_2')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{input}{a vector with dimension num_obs x 1 for the sorted input locations.}
    \item{output}{a  vector with dimension n x 1 for the observations at the sorted input locations.}

  \item{have_noise}{a bool value. If it is true, it means the model contains a noise. 
%%     ~~Describe \code{response} here~~
}

  \item{kernel_type}{a \code{character} to specify the type of kernel to use. The current version supports kernel_type to be "matern_5_2" or "exp", meaning that the matern kernel with roughness parameter being 2.5 or 0.5 (exponent kernel), respectively. 
%%     ~~Describe \code{trend} here~~
}

}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%expand here the details.
%}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\code{fgasp} returns an S4 object of class \code{fgasp} (see \code{fgasp}).

%If there is an emulator (i.e. simul_type is 0), \code{rcalibration} returns a list with 
%\item{rcalibration}{an S4 object of class \code{rcalibration} (see \code{rcalibration-class}.}
%\item{emulator}{an S4 object of class \code{rgasp} produced by RobustGaSP R Package .}


}

\references{
%% ~put references to the literature/web site here ~

Hartikainen, J. and Sarkka, S. (2010). \emph{Kalman filtering and smoothing solutions to temporal gaussian process regression models},  \emph{Machine Learning for Signal Processing (MLSP), 2010 IEEE International Workshop},  379-384.

M. Gu, Y. Xu (2017), \emph{Nonseparable Gaussian stochastic process: a unified
view and computational strategy}, arXiv:1711.11501.

M. Gu, X. Wang and J.O. Berger (2018), \emph{Robust Gaussian Stochastic Process Emulation}, \emph{Annals of Statistics}, \bold{46}, 3038-3066.

}

\author{
\packageAuthor{FastGaSP}

Maintainer: \packageMaintainer{FastGaSP}
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{

library(FastGaSP)

#-------------------------------------
# Example 1: a simple example with noise 
#-------------------------------------

y_R<-function(x){
  cos(2*pi*x)
}

###let's test for 2000 observations
set.seed(1)
num_obs=2000
input=runif(num_obs)

output=y_R(input)+rnorm(num_obs,mean=0,sd=0.1)

##constucting the fgasp.model
fgasp.model=fgasp(input, output)
show(fgasp.model)

#------------------------------------------
# Example 2: a simple example with no noise 
#------------------------------------------

y_R<-function(x){
  sin(2*pi*x)
}


##generate some data without noise
num_obs=50
input=seq(0,1,1/(num_obs-1))

output=y_R(input)


##constucting the fgasp.model
fgasp.model=fgasp(input, output,have_noise=FALSE)

show(fgasp.model)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
