\name{sim.pw.arma}
\alias{sim.pw.arma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of a piecewise constant ARMA(p,q) model for p=2 and q=1
}
\description{
The function simulates a piecewise constant ARMA model with multiple change-points
}
\usage{
sim.pw.arma(N, sd_u, b.slope, b.slope2, mac, br.loc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
Length of the series
}
  \item{sd_u}{
A vector of the innovation standard deviation for every segment
}
  \item{b.slope}{
A vector of the AR(1) coefficients
}
  \item{b.slope2}{
A vector of the AR(2) coefficients
}
  \item{mac}{A vector of the MA(1) coefficients

}
  \item{br.loc}{A vector with the location of the change-points

}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
A simulated series
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\examples{ 
cps=c(125,532,704)
y=sim.pw.arma(N = 1024,sd_u = 1,b.slope=c(0.7,0.3,0.9,0.1),
b.slope2 = c(0,0,0,0), mac = c(0.6,0.3,0,-0.5),br.loc = cps)[[2]]
ts.plot(y)
abline(v=cps,col="red")                    
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ piecewise ARMA }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
