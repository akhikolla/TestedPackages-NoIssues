\name{omp}
\alias{omp}
\alias{hasOMP}
\alias{setOMPThreads}
\alias{setSITMOSeeds}

\title{Support for parallel computation}

\description{
    The functions can be used 
        \itemize{
            \item to check if the current system supports the openMP standard;
            \item to control the number of used cores;
            \item to set the seeds of the random number generators.
        }
}

\usage{
hasOMP()

setOMPThreads(nthreads)

setSITMOSeeds(seed)
}

\arguments{
    \item{nthreads}{
        number of OpenMP threads to be used.
    }  
    \item{seed}{
        number between 0 and 1 used to set the seeds of the
        random number generators in each threads.
    }
}

\details{
    Each openMP thread (or the single thread used on systems where openMP is not
    supported) uses a separate \code{sitmo} random number generator.
    See \code{\link[sitmo]{sitmo-package}}.
}

\note{
    When the package is loaded, the following code is automatically executed
    \itemize{
        \item \code{if (hasOMP()) setOMPThreads(parallel::detectCores())}
        \item \code{setSITMOSeeds(runif(1))}
    }
}

\value{
    Function \code{hasOMP} returns TRUE/FALSE if the system
    supports/does not support openMP.

    Functions \code{setOMPThreads} and \code{setSITMOSeeds}
    do not return any value.
}

\author{Giovanna Capizzi and Guido Masarotto}

\keyword{ utilities }% use one of  RShowDoc("KEYWORDS")
