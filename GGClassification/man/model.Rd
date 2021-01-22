\name{model}
\alias{model}
\title{Calculates the parameters of a large-margin classifier based on the Gabriel Graph.}
\usage{
model(X, y, normalize=FALSE)
}
\arguments{
  \item{X}{Matrix containing data.}
  
  \item{y}{Vector of labels.}
  
  \item{normalize}{If data should be normalized or not.}
}
\value{

A named list containing the parameters of the classifier.

}
\description{

The function computes the Gabriel Graph and applies a filter on the graph to remove noise in the margin region.

Then a new Graph is calculated considering the filtered data, and the parameters of a
large-margin classifier are obtained and returned as a named list.

}

\references{
Gabriel, K. R., & Sokal, R. R. (1969). A New Statistical Approach to Geographic Variation Analysis. Systematic Zoology, 18(3), 259. doi:10.2307/2412323 

L. C. B. Torres, C. L. Castro, F. Coelho, F. Sill Torres and A. P. Braga, "Distance-based large margin classifier suitable for integrated circuit implementation," in Electronics Letters, vol. 51, no. 24, pp. 1967-1969, 19 11 2015.
}

\examples{
X <- matrix(0.3 * rnorm(100) + 2.5, ncol = 2)
y <- c(rep(0, times = 25), rep(1, times = 25))

mdl <- model(X, y)
}