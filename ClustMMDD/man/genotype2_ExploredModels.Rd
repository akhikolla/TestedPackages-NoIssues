\name{genotype2_ExploredModels}
\alias{genotype2_ExploredModels}
\docType{data}
\title{
  A data frame of competing models gathered by \code{\link{backward.explorer}}.
}
\description{
  A data frame of competing models gathered by \code{\link{backward.explorer}} for \eqn{Kmax = 10}. Such data file can be used for a final model selection process.
}
\usage{data("genotype2_ExploredModels")}
\format{
  A data frame with 2667 explored models on the following 16 variables.
  \describe{
    \item{\code{N} : }{The size of the data}
    \item{\code{P} : }{The number of variables}
    \item{\code{K} : }{the number of clusters}
    \item{\code{S1} : }{1st variable}
    \item{\code{S2} : }{2nd variable}
    \item{\code{S3} : }{3th variable}
    \item{\code{S4} : }{4th variable}
    \item{\code{S5} : }{5th variable}
    \item{\code{S6} : }{6th variable}
    \item{\code{S7} : }{7th variable}
    \item{\code{S8} : }{8th variable}
    \item{\code{S9} : }{9th variable}
    \item{\code{S10} : }{10th variable}
    \item{\code{logLik} : }{The log-likelihood}
    \item{\code{dim} : }{The dimension = number of free parameters}
    \item{\code{entropy} : }{Entropy.}
  }
}

\details{TODO
}

\source{Wilson Toussile
}

\seealso{
\code{\link{dimJump.R}} and \code{\link{model.selection.R}}.
}

\examples{
data(genotype2_ExploredModels)
head(genotype2_ExploredModels)
plot(genotype2_ExploredModels[, c("dim", "logLik")],
	 col = "blue", xlab = "Dimension", ylab = "Log-likelihood")

# Data-driven calibration of the penalty
dimJump.R(genotype2_ExploredModels, h = 5, N=1000, header=T)
}

\keyword{Backward-Stepwise selection}
\keyword{Model selection in clustering}
\keyword{Population structure}
\keyword{Genotype data}
