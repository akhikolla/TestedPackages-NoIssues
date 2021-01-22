\name{design}
\alias{x.cl}
\alias{cusum.cl}
\alias{ewma.cl}
\title{Design of control charts based on the cautious learning approach}

\description{
    These functions compute the control limits  
    of X (\code{x.cl}), EWMA (\code{ewma.cl}) and CUSUM (\code{cusum.cl})
    control charts based on the cautious learning approach. 
    The stochastic approximation algorithm, described in 
    the Appendix A of Capizzi and Masarotto (2019), is used. 

    When openMP is supported, computation can be distribuited on multiple cores.
    See \code{\link{omp}}. 
}

\usage{
x.cl(m, arl0, alpha = 0.1, beta = 0.05, H = 200, A = 1.5, B = 50, 
     Ninit = 1000, Nfinal = 30000)

ewma.cl(lambda, m, arl0, alpha = 0.1, beta = 0.05, H = 200, A = 1.5, B = 50, 
        Ninit = 1000, Nfinal = 30000)

cusum.cl(k, m, arl0, alpha = 0.1, beta = 0.05, H = 200, A = 1.5, B = 50, 
         Ninit = 1000, Nfinal = 30000)
}

\arguments{
  \item{lambda}{EWMA smoothing constant.}
  \item{k}{CUSUM reference value.}
  \item{m}{
    number of in-control observations used to estimate the process mean
    and standard deviation at the beginning of the monitoring phase.
  }
  \item{arl0, alpha, beta, H}{
    desired in-control average run-length and constants defining the
    empirical guaranteed in-control performance condition. See
    equations (2) and (6) in Capizzi and Masarotto (2019).
  }
  \item{A, B}{
    constants controlling when the parameters estimate are updated.
    See equation (3) in Capizzi and Masarotto (2019). 
    If \code{A=NA} and \code{B=NA}, the no-learning control limits
    are computed.
  }
  \item{Ninit, Nfinal}{
   number of iterations used in the stochastic approximation algorithm.
   See Capizzi and Masarotto (2019), Appendix A.
  }
}

\value{
A list with the following elements:

    \item{chart}{string describing the control chart ("X", "EWMA" or
    "CUSUM").}
    \item{lambda}{EWMA smoothing constant (only when
        \code{chart=="EWMA"}).}
    \item{k}{CUSUM reference value (only when
        \code{chart=="CUSUM"}).}
    \item{limit}{numeric vector of length equal to five containing the
        constants defining the cautiuos learning control limits, i.e,
        \eqn{L_\infty}{Linf}, \eqn{\Delta}{Delta}, 
        A, B and m (see equation (3) and (4) in Capizzi and
        Masarotto (2019)).}
}

\references{
Capizzi, G. and Masarotto, G. (2019) 
"Guaranteed In-Control Control Chart Performance with Cautious Parameter
Learning", accepted for publication in \emph{Journal of Quality Technology},
a copy of the paper can be obtained from the authors.
}

\author{Giovanna Capizzi and Guido Masarotto}

\examples{
## Only for testing: the number of iterations is reduced 
## to reduce the computing time
Ninit <- 50
Nfinal <- 100
H <- 50
x.cl(100, 500, Ninit=Ninit, Nfinal=Nfinal, H=H)
x.cl(100, 500, A=NA, B=NA, Ninit=Ninit, Nfinal=Nfinal, H=H)
ewma.cl(0.2, 100, 500, Ninit=Ninit, Nfinal=Nfinal, H=H)
cusum.cl(1, 100, 500, Ninit=Ninit, Nfinal=Nfinal, H=H)
\donttest{
## Using the default number of iterations
x.cl(100, 500)
x.cl(100, 500, A=NA, B=NA)
ewma.cl(0.2,100, 500)
cusum.cl(1, 100, 500)
}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
