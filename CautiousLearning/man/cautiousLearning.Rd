\name{cautiousLearning}
\alias{cautiousLearning}
\title{
Applications of control charts based on the cautious learning approach
}

\description{
This function applies and, optionally, plots a control chart based on the cautious
learning approach described in Capizzi and Masarotto (2019).
}

\usage{
cautiousLearning(chart, x, mu0, s0, plot = TRUE)
}

\arguments{
  \item{chart}{
    list with the same elements as those returned by
    \code{\link{x.cl}}, \code{\link{ewma.cl}} and \code{\link{cusum.cl}}.
  }
  \item{x}{
      numeric vector containing the Phase II data.
  }
  \item{mu0, s0}{
    estimates of the in-control mean and standard deviation
    obtained by the Phase I reference sample.
  }
  \item{plot}{
    if \code{TRUE} the control statistics and the cautiuos control
    limits are plotted.
  }
}


\value{
The function returns (invisibly when \code{plot==TRUE}) a numeric matrix 
containing 

\item{column 1 for X and EWMA, columns 1-2 for CUSUM}{
      control statistic[s]}
\item{columns 2-4 for X and EWMA, columns 3-5 for CUSUM}{
      central line, lower and upper control limits}
\item{columns 5-7 for X and EWMA, columns 6-8 for CUSUM}{
    "cautious" estimates of the mean, standard deviation and critical
    value, i.e., using the notation in Capizzi and Masarotto (2019),
    \eqn{\hat{\mu}_{i-d_i}}{mu.hat[i-d[i]]}, 
    \eqn{\hat{\sigma}_{i-d_i}}{sigma.hat[i-d[i]]} and
    \eqn{L_{i-d_i}}{L[i-d[i]]}. }
}

\references{
Capizzi, G. and Masarotto, G. (2019) 
"Guaranteed In-Control Control Chart Performance with Cautious Parameter
Learning", accepted for publication in \emph{Journal of Quality Technology},
a copy of the paper can be obtained from the authors.
}

\author{Giovanna Capizzi and Guido Masarotto}

\examples{
## EWMA control chart (nominal ARL=500, 
## initial estimates based on 100 in-control observations)
chart <- list(chart = "EWMA",
              lambda = 0.2,
              limit = c(Linf=3.187, Delta=0.427, A=1.5, B=50, m=100))
## Phase I estimates
set.seed(12345)
xic <- rnorm(100, 12 , 3)
m0 <- mean(xic)
s0 <- sd(xic)
## Phase II observations (one sigma mean shift starting at i=501)
x <- c(rnorm(500, 12, 3), rnorm(50, 15, 3))
## Monitoring
y <- cautiousLearning(chart, x, m0, s0)
head(y)
tail(y)
}
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
