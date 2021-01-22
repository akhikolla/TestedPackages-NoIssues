\name{summary.predddpcomp}
\alias{summary.predddpcomp}
\docType{package}
\title{
generate summary of the predict output with predddpcomp class.
}
\description{
generate estimated cumulative incidence functions and corresponding credible intervals at 4 time points (1/4, 1/2, 3/4 and 1 of high.pct) for each covaraite configuration.
}
\usage{
\method{summary}{predddpcomp}(object,...)
}
\arguments{
\item{object}{Output from dpweib of predddpcomp class}
\item{...}{Arguments to be passed to method}
}
\value{
an object of class summary.predddpcomp
}
