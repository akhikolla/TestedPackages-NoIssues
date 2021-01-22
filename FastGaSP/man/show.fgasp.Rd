\name{show}
\alias{show}
\alias{show.fgasp}
\alias{show.fgasp-class}
\alias{show,fgasp-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Show an \code{fgasp} object.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function to print the code{fgasp} object. 
}
\usage{
\S4method{show}{fgasp}(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of  class \code{fgasp}.}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%Provide here some details.
%}

\references{
Hartikainen, J. and Sarkka, S. (2010). \emph{Kalman filtering and smoothing solutions to temporal gaussian process regression models},  \emph{Machine Learning for Signal Processing (MLSP), 2010 IEEE International Workshop},  379-384.

M. Gu, Y. Xu (2017), \emph{Nonseparable Gaussian stochastic process: a unified
view and computational strategy}, arXiv:1711.11501.

M. Gu, X. Wang and J.O. Berger (2018), \emph{Robust Gaussian Stochastic Process Emulation}, \emph{Annals of Statistics}, \bold{46}, 3038-3066.
}

\author{
%%  ~~who you are~~
\packageAuthor{FastGaSP}

Maintainer: \packageMaintainer{FastGaSP}
}

\examples{

#-------------------------------------
# Example: a simple example with noise 
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


}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
