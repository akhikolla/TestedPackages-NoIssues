\name{continue}
\alias{continue}
\docType{package}
\title{
Resume MCMC run
}
\description{
Self-defined S3 method that resumes MCMC run from an object of class dpm ,ddp, dpmcomp or ddpcomp.
}
\usage{
continue(previous,iteration=1000,...)
}
\arguments{
\item{previous}{An object of class dpm, ddp, dpmcomp or ddpcomp}
\item{iteration}{The number of iterations to continue sampling. The default is 1000.}
\item{...}{Arguments to be passed to method}
}
\value{
An object of the class of previous.
}
