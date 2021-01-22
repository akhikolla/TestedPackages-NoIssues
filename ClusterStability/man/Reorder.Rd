\encoding{ISO-8859-2}
\name{Reorder}
\alias{Reorder}
\title{This function returns the ordering of a partitioning solution in ascending order.}
\description{This function returns the ordered partition of a set of numbers in ascending order and reorderd to start at one. This is an auxiliary function.}
\usage{Reorder(data)}
\arguments{
	\item{data}{vector of partition numbers to reorder.} 
}
\value{A vector of ordered partition numbers for this data.}

\examples{
  Reorder(c(1,3,4,4,3,1))
  # Expected : 1 2 3 3 2 1
}

%\keyword{ }

