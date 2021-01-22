\encoding{ISO-8859-2}
\name{Kcombination}
\alias{Kcombination}
\title{Kcombination returns the list of all possible combinations of a set of numbers of a given length \emph{k}.}
\description{This function, given a vector of numbers, will return all the possible combinations of a given length \emph{k}.}
\usage{Kcombination(data, k, selector)}
\arguments{
	\item{data}{the vector of numbers (\emph{i.e.} elements) to consider.}
    \item{k}{the length of the returned combination (between 2 and 6 in this version).}
    \item{selector}{if set, returns only the combinations containing this number.}
}

\value{Return a list of all possible combinations for the given vector of numbers.}

\examples{
   ## Returns the k-combination of the list of numbers: 1,2,3 of length=2.
   ## i.e. (1,2), (1,3), (2,3) 
	Kcombination(c(1,2,3),k=2)
   ## Returns only the k-combination containing the number 1.
   ## i.e. (1,2), (1,3)	
	Kcombination(c(1,2,3),k=2,selector=1)
}

\keyword{k-combination}

