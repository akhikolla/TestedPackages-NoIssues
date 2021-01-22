\name{NEWS}
\title{durmod news}
\section{Changes in version 1.1-2}{
  \itemize{
    \item Some adjustments to default parameters.
    \item Use median, not mean, for the interval search for new
    points. Use Neumaier compensated addition in computation of gradients.
    \item Some more documentation.
  }
}
\section{Changes in version 1.1-1}{
  \itemize{
    \item Fixed a bug in the interval censored gradient which slowed down
    convergence dramatically.

    \item Fixed a bug which prevents calling \code{mphcrm} from within a
    function with variables defined there. Like in an \code{mclapply} which
    creates the data frame.
    
    \item Fixed a bug which prevented the use of factors when running on
    a cluster.
  }
}
\section{Changes in version 1.1}{
  \itemize{
    \item Remove covariates and factor levels which are constant for a
    transition. They will act as intercepts, we don't want that. If you
    really want such a thing, create it with some neglible variation,
    e.g. as \code{rnorm(N,1,10*.Machine$double.eps)}.
    \item Added generics for \code{vcov()}, \code{coef()}.
  }
}
