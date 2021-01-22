\name{simulation}
\alias{ruv}
\alias{rcrl}
\title{Estimation errors and conditional run-length simulation}

\description{
    Function \code{ruv} simulates the standardized
    estimation errors at the starting of the monitoring phase
    (see Section 2.3 of Capizzi and Masarotto (2019)). 

    Function \code{rcrl} simulates, under different in-control or
    out-control scenarios, the conditional run-length given
    the standardized estimation errors. When openMP is supported,
    computation can be distribuited on multiple cores.
    See \code{\link{omp}}. 
}
\usage{
ruv(n, m)

rcrl(n, chart, u, v, tau, delta, omega, maxrl = 1000000L)
}

\arguments{
  \item{n}{
      number of simulated values.
    }
  \item{m}{
    number of in-control observations available at the starting of the
    monitoring phase.
  }
  \item{chart}{
     list with the same elements as those returned by
    \code{\link{x.cl}}, \code{\link{ewma.cl}} and \code{\link{cusum.cl}}.
  }

  \item{u, v}{
    values of the estimation errors (scalars).
  }
  \item{tau, delta, omega}{
    when i<tau, observations are distributed as N(mu,sigma^2) random
    variables; when i>=tau, observations are distributed as
    N(mu+delta*sigma, (omega*sigma)^2) random variables.
  }
  \item{maxrl}{
    run-length are truncated at i=maxrl.
  }
}
\value{
  \item{ruv}{numeric matrix of dimension nx2.}

  \item{rcrl}{integer vector of length n.}
}


\references{
Capizzi, G. and Masarotto, G. (2019) 
"Guaranteed In-Control Control Chart Performance with Cautious Parameter
Learning", accepted for publication in \emph{Journal of Quality Technology},
a copy of the paper can be obtained from the authors.
}

\author{Giovanna Capizzi and Guido Masarotto}




\examples{
ruv(5, 100)
## EWMA control chart (nominal ARL=500, 
## initial estimates based on 100 in-control observations)
chart <- list(chart = "EWMA",
              lambda = 0.2,
              limit = c(Linf=3.187, Delta=0.427, A=1.5, B=50, m=100))
rcrl(10, chart, 2, 1, 50, 2, 1) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
