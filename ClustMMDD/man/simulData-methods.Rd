\name{simulData-methods}
\docType{methods}

\alias{simulData-methods}
\alias{simulData}


\title{
Simulate a dataset from a given set of parameters in an instance of \code{\linkS4class{modelKS}}.
}


\section{Methods}{
\describe{

\item{\code{signature(object = "modelKS", N = "numeric", ploidy = "numeric")}}{Simulate a dataset for a given set of parameters in a \code{\linkS4class{modelKS}} object.}
}}

\description{
Simulate a dataset from a given instance of \code{\linkS4class{modelKS}} containing a set of parameters.
}

\arguments{
  \item{object}{An instance of \code{\linkS4class{modelKS}}.}

  \item{N}{The size of the sample to simulate.}
  
  \item{ploidy}{The number of columns for each variable in the data. For example, \eqn{ploidy = 2} for genotypic data from diploid individual.}
}
\value{
  A list :
  \itemize{
  \item{[["data"]] : }{The simulated dataset.}
  
  \item{[["class"]] : }{The membership class.}
  }
}

\author{
  Wilson Toussile
}

\seealso{
  \code{\linkS4class{modelKS}}, \code{\link{exModelKS}}.
}

\examples{
data(exModelKS)
exModelKS
exData = simulData(exModelKS, 1000, 2)
str(exData)
head(exData$data)
head(exData$class)
}

