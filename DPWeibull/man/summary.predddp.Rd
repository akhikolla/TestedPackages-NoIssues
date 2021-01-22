\name{summary.predddp}
\alias{summary.predddp}
\docType{package}
\title{
generate summary of the predict output with predddp class.
}
\description{
generate estimated survival and corresponding credible intervals at 4 time points (1/4, 1/2, 3/4 and 1 of high.pct) for each covaraite configuration.
}
\usage{
\method{summary}{predddp}(object,...)
}
\arguments{
\item{object}{Output from dpweib of predddp class}
\item{...}{Arguments to be passed to method}
}
\value{
an object of class summary.predddp
}
