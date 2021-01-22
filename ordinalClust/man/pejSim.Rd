\name{pejSim}
\alias{pejSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
pejSim
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes the probabilty for a level ej to be sampled from a BOS distribution of parameters (mu,pi), with the number of levels equal to m. It can be used to generate data with a BOS distribution.
}
\usage{
pejSim(ej,m,mu,p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ej}{
%%     ~~Describe \code{x} here~~
Levels to be sampled
}
  \item{m}{
%%     ~~Describe \code{kr} here~~
Number of levels.
}
  \item{mu}{
%%     ~~Describe \code{m} here~~
mu parameter for BOS distribution.
}
  \item{p}{
%%     ~~Describe \code{nbSEM} here~~
pi parameter for BOS distribution.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  Returns the probability of ej to be sampled from a BOS distribution of parameters (mu,pi), with the number of levels equal to m.
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
%%  ~~who you are~~
Margot Selosse, Julien Jacques, Christophe Biernacki.
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
  library(ordinalClust)
  data("dataqol")
  set.seed(5)

  m=7
  nr=10000
  mu=5
  pi=0.5

  probaBOS=rep(0,m)
  for (im in 1:m) probaBOS[im]=pejSim(im,m,mu,pi)
  M <- sample(1:m,nr,prob = probaBOS, replace=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{"ordinal-data"}% use one of  RShowDoc("KEYWORDS")
\keyword{"clustering"}% __ONLY ONE__ keyword per line
