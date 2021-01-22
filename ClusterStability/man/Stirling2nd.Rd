\encoding{ISO-8859-2}
\name{Stirling2nd}
\alias{Stirling2nd}
\title{Stirling2nd function computes the Stirling numbers of the second kind.}
\description{This function returns the estimated Stirling numbers of the second kind \emph{i.e.}, the number of ways of partitioning a set of \emph{n} objects into \emph{k} nonempty groups.}
\usage{Stirling2nd(n,k)}
\arguments{
	\item{n}{number of objects.}
    \item{k}{number of groups (\emph{i.e.} classes).}    
}

\value{The Stirling number of the 2nd kind for \emph{n} elements and \emph{k} groups or \emph{NaN} (if the Stirling number for those \emph{n} and \emph{k} is greater than 1e300).}

\examples{
  Stirling2nd(n=3,k=2)
  # Expected value=3
  Stirling2nd(n=300,k=20)
  # Expected value=NaN  
}

%\keyword{ }

