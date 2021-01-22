\name{cbinom}
\alias{dcbinom}
\alias{pcbinom}
\alias{qcbinom}
\alias{rcbinom}
\title{The Continuous Binomial Distribution}

\description{Density, distribution function, quantile function and random generation for
a continuous analog to the binomial distribution with parameters \code{size}
and \code{prob}. The usage and help pages are modeled on the d-p-q-r families of
functions for the commonly-used distributions (e.g., \code{\link[stats]{dbinom}})
in the \code{stats} package.

Heuristically speaking, this distribution spreads the standard probability mass
(\code{\link[stats]{dbinom}}) at integer \code{x} to the interval
\code{[x, x + 1]} in a continuous manner. As a result, the distribution looks
like a smoothed version of the standard, discrete binomial but shifted slightly
to the right. The support of the continuous binomial is \code{[0, size + 1]},
and the mean is approximately \code{size * prob + 1/2}.
}

\usage{
  dcbinom(x, size, prob, log = FALSE)
  pcbinom(q, size, prob, lower.tail = TRUE, log.p = FALSE)
  qcbinom(p, size, prob, lower.tail = TRUE, log.p = FALSE)
  rcbinom(n, size, prob)
}

\arguments{
  \item{x, q}{   vector of quantiles.}
  \item{p}{   vector of probabilities.}
  \item{n}{   number of observations. If \code{length(n) > 1}, the length is
   taken to be the number required.}
  \item{size}{   the \code{size} parameter.}
  \item{prob}{   the \code{prob} parameter.}
  \item{log, log.p}{   logical; if TRUE, probabilities p are given as log(p)}
  \item{lower.tail}{   logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]}
}

\details{
The \code{cbinom} package is an implementation of Ilienko's (2013) continuous
binomial distribution.

The continuous binomial distribution with \code{size =} \emph{N} and
\code{prob =} \emph{p} has cumulative distribution function
  \deqn{F(x) = \frac{B(x, N + 1 - x, p)}{B(x, N + 1 - x)}}{%
F(x) = B_p(x, N - x + 1)/B(x, N - x + 1)}
for \code{x} in \code{[0, N + 1]}, where
\deqn{{B(x, N + 1 - x, p) = \int_{p}^{1} t^{x-1}(1-t)^{y-1}dt}}{%
B_p(x, N - x + 1) = integral_0^p (t^(x-1)(1-t)^(y-1))dt} is the incomplete beta
function and \deqn{{B(x, N + 1 - x) = \int_{0}^{1} t^{x-1}(1-t)^{y-1}dt}}{%
B(x, N - x + 1) = integral_0^1 t^(x-1)(1-t)^(y-1)dt} is the beta function (or
\code{beta(x, N - x + 1)} in R). The CDF can be expressed in R as
F(x) = \code{1 - pbeta(prob, x, size - x + 1)} and the mean calculated as
\code{integrate(function(x) pbeta(prob, x, size - x + 1), lower = 0, upper = size + 1)}.

If an element of \code{x} is not in \code{[0, N + 1]}, the result of
\code{dcbinom} is zero. The PDF \code{dcbinom(x, size, prob)} is computed via
numerical differentiation of the CDF = \code{1 - pbeta(prob, x, size - x + 1)}.}

\value{
\code{dcbinom} is the density, \code{pcbinom} is the distribution function,
\code{qcbinom} is the quantile function, and \code{rcbinom} generates random
deviates.

The length of the result is determined by \code{n} for \code{rbinom}, and is the
maximum of the lengths of the numerical arguments for the other functions.

The numerical arguments other than \code{n} are recycled to the length of the
result.
}

\references{
Ilienko, Andreii (2013). Continuous counterparts of Poisson and binomial
distributions and their properties. Annales Univ. Sci. Budapest., Sect. Comp.
39: 137-147. \url{http://ac.inf.elte.hu/Vol_039_2013/137_39.pdf}
}

\examples{
require(graphics)
# Compare continous binomial to a standard binomial
size <- 20
prob <- 0.2
x <- 0:20
xx <- seq(0, 21, length = 200)
plot(x, pbinom(x, size, prob), xlab = "x", ylab = "P(X <= x)")
lines(xx, pcbinom(xx, size, prob))
legend('bottomright', legend = c("standard binomial", "continuous binomial"),
  pch = c(1, NA), lty = c(NA, 1))
mtext(side = 3, line = 1.5, text = "pcbinom resembles pbinom but continuous and shifted")
pbinom(x, size, prob) - pcbinom(x + 1, size, prob)

# Use "log = TRUE" for more accuracy in the tails and an extended range:
n <- 1000
k <- seq(0, n, by = 20)
cbind(exp(dcbinom(k, n, .481, log = TRUE)), dcbinom(k, n, .481))
}