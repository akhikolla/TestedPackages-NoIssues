\name{sim.pw.ar2}
\alias{sim.pw.ar2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of a piecewise constant AR(2) model
}
\description{
The function simulates a piecewise constant AR(2) model with multiple change-points
}
\usage{
sim.pw.ar2(N, sd_u, b.slope, b.slope2, br.loc)
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
  \item{br.loc}{
A vector with the location of the change-points
}
}

\value{
A simulated series
}

\examples{ 
cps=c(512,754)
y=sim.pw.ar2(N =1024,sd_u = 1,b.slope=c(0.9,1.68,1.32),
b.slope2 = c(0.0,-0.81,-0.81),br.loc = cps)[[2]]
ts.plot(y)
abline(v=cps,col="red")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ piecewise AR(2) }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
