\name{cbinom-package}
\alias{cbinom-package}
\alias{cbinom}
\docType{package}
\title{
Continuous Analog of a Binomial Distribution
}
\description{
Implementation of the d/p/q/r family of functions for a continuous analog to the
standard discrete binomial with continuous size parameter and continuous support
with \code{x} in \code{[0, size + 1]}.
}
\details{
Included in the package are functions \code{dcbinom(x, size, prob, log = FALSE)},
\code{pcbinom(q, size, prob, lower.tail = TRUE, log.p = FALSE)},
\code{qcbinom(p, size, prob, lower.tail = TRUE, log.p = FALSE)}, and
\code{rcbinom(n, size, prob}. Usage closely parallels that of the \code{binom}
family of functions in the \code{stats} R package.
}
\author{
Dan Dalthorp <ddalthorp@usgs.gov>
}
\references{
Ilienko, Andreii (2013). Continuous counterparts of Poisson and binomial
distributions and their properties. Annales Univ. Sci. Budapest., Sect. Comp.
39: 137-147. \url{http://ac.inf.elte.hu/Vol_039_2013/137_39.pdf}
}
\keyword{ continuous binomial, R package }
\seealso{
\code{\link[cbinom]{pcbinom}}
}
