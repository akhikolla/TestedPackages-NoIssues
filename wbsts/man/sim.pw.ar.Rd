\name{sim.pw.ar}
\alias{sim.pw.ar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation of a piecewise constant AR(1) model
}
\description{
The function simulates a piecewise constant AR(1) model with multiple change-points
}
\usage{
sim.pw.ar(N, sd_u, b.slope, br.loc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
Length of the series.
}
  \item{sd_u}{
A vector of the innovation standard deviation for every segment.
}
  \item{b.slope}{
A vector of the AR(1) coefficients.
}
  \item{br.loc}{
A vector with the location of the change-points.
}
}
\value{
A simulated series
}
\examples{
cps=c(400,612)
y=sim.pw.ar(N =1024,sd_u = 1,b.slope=c(0.4,-0.6,0.5),br.loc=cps)[[2]]
ts.plot(y)
abline(v=cps,col="red")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ piecewise AR(1) }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
