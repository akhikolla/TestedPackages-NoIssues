\name{GabrielGraph}
\alias{GabrielGraph}
\title{Constructs a Gabriel Graph from data.}
\usage{
GabrielGraph(X)
}
\arguments{
  \item{X}{Matrix containing data.}
}
\value{
The adjacency matrix of the Gabriel Graph, given the data matrix \code{X}.  
}
\description{
Constructs a Gabriel Graph from data and returns its adjacency matrix.
}
\references{

Gabriel, K. R., & Sokal, R. R. (1969). A New Statistical Approach to Geographic Variation Analysis. Systematic Zoology, 18(3), 259. doi:10.2307/2412323 

Matula, David & Sokal, Robert. (1980). Properties of Gabriel Graphs Relevant to Geographic Variation Research and the Clustering of Points in the Plane. Geographical Analysis. 12. 205 - 222. 10.1111/j.1538-4632.1980.tb00031.x. 

}
\examples{
X <- matrix(0.3 * rnorm(100) + 2.5, ncol = 2)
adjacency_matrix = GabrielGraph(X)
}