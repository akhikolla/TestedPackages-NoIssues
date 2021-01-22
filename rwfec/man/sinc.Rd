\name{sinc}
\alias{sinc}
\title{sinc}
\usage{
sinc(x)
}
\description{
sinc(x) = 1 at x = 0, sinc(x) = sin(x)/x for x != 0
}
\arguments{
\item{x}{- input vector}
}
\value{
Response (Double, Numeric vector) of the sinc function applied to the input vector x.
}

\examples{
x <- seq(-6*pi, 6*pi, by = pi/10)
p <- sinc(x)
plot(x/pi,p,ylim=c(-0.3,1.1), pch=19, cex=0.25, ylab="p", xlab="x*pi", type="l")
grid( col = "grey50", lty = "dotted")
}
