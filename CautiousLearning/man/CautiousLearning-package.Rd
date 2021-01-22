\name{CautiousLearning-package}
\alias{CautiousLearning-package}
\docType{package}
\title{  
Guaranteed In-Control Control Chart Performance with
Cautious Parameter Learning
}
\description{
Functions in this package allow to design, study and apply control
charts based on the cautious parameter learning approach described in 
Capizzi and Masarotto (2019). 

On system where the openMP standard is
supported, these functions can take advantage of the
computing power offered by a multicore workstation.
See \code{\link{omp}} for the default setting.
}

\details{
    The package includes the following functions:
    \itemize{
    \item 
        Computation of the control limits via stochastic
        approximation: \code{\link{x.cl}}, \code{\link{ewma.cl}},
        \code{\link{cusum.cl}};
    \item 
        Estimation errors and conditional run-length simulation: 
        \code{\link{ruv}} and \code{\link{rcrl}};
    \item 
        Application to real data: \code{\link{cautiousLearning}};
    \item
        Controlling the number of used openMP cores and the random number
        generator seeds: \code{\link{hasOMP}}, \code{\link{setOMPThreads}} and
        \code{\link{setSITMOSeeds}}.    
    }
}

\author{
Giovanna Capizzi <giovanna.capizzi@unipd.it> and
Guido Masarotto <guido.masarotto@unipd.it>

Maintainer: Giovanna Capizzi
}



\references{
Capizzi, G. and Masarotto, G. (2019) 
"Guaranteed In-Control Control Chart Performance with Cautious Parameter
Learning", accepted for publication in \emph{Journal of Quality Technology},
a copy of the paper can be obtained from the authors.
}
